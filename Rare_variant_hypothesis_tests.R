# =============================================================================
# DUAL gRNA SCRIPT — Rare pools eCROPseq
# =============================================================================

library(sceptre)
library(Matrix)
library(Seurat)
library(scuttle)
library(scDblFinder)
library(dplyr)
library(ggplot2)
library(tidyr)

setwd("path/to/directory")

POOL_NAME   <- "pool_name"
POOL_PATH_1 <- "path/to/pool"
POOL_PATH_2 <- NULL
POOL_PATH_3 <- NULL
POOL_PATH_4 <- NULL
POOL_PATH_5 <- NULL

# =============================================================================
# 1. LOAD DATA
# =============================================================================

fix_feature_names <- function(mat) {
  rownames(mat) <- gsub("'", "", rownames(mat))
  mat
}

load_pool <- function(path, rep_label) {
  Read10X(path) %>%
    fix_feature_names() %>%
    CreateSeuratObject(counts = ., min.cells = 6, min.features = 200) %>%
    AddMetaData(metadata = rep_label, col.name = "Sample")
}

pool_paths <- Filter(Negate(is.null),
                     list(POOL_PATH_1, POOL_PATH_2, POOL_PATH_3,
                          POOL_PATH_4, POOL_PATH_5))
rep_labels <- c("Rep1", "Rep2", "Rep3", "Rep4", "Rep5")
n_pools    <- length(pool_paths)

seurat_list <- mapply(load_pool, pool_paths, rep_labels[seq_len(n_pools)],
                      SIMPLIFY = FALSE)

if (n_pools == 1) {
  data <- seurat_list[[1]]
  message("Single pool loaded.")
} else {
  data <- merge(seurat_list[[1]], y = seurat_list[-1],
                add.cell.ids = rep_labels[seq_len(n_pools)])
  message(sprintf("Merged %d pools.", n_pools))
}
rm(seurat_list); gc()
data <- JoinLayers(data)
gc()

data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
message(sprintf("Total cells before QC: %d", ncol(data)))

# =============================================================================
# 2. DOUBLET DETECTION
# =============================================================================

data_for_doublets <- subset(data, subset = percent.mt < 35)
sce               <- as.SingleCellExperiment(data_for_doublets)
sce               <- scDblFinder(sce, samples = "Sample")

doublet_calls <- data.frame(
  cell       = colnames(sce),
  is_doublet = sce$scDblFinder.class == "doublet",
  stringsAsFactors = FALSE
)

data$is_doublet <- ifelse(
  is.na(doublet_calls$is_doublet[match(colnames(data), doublet_calls$cell)]),
  TRUE,
  doublet_calls$is_doublet[match(colnames(data), doublet_calls$cell)]
)

message(sprintf("Doublets detected: %d (%.1f%%)",
                sum(data$is_doublet), 100 * mean(data$is_doublet)))
data <- data[, !data$is_doublet]
message(sprintf("Cells after doublet removal: %d", ncol(data)))
rm(sce, data_for_doublets); gc()

# =============================================================================
# 3. nFEATURE QC
# =============================================================================

qc_per_sample <- function(seurat_obj, sample_id, nmads = 2.5) {
  cells <- colnames(seurat_obj)[seurat_obj$Sample == sample_id]
  obj   <- seurat_obj[, cells]
  qc_lo <- isOutlier(obj$nFeature_RNA, log = TRUE, type = "lower",  nmads = nmads)
  qc_hi <- isOutlier(obj$nFeature_RNA, log = TRUE, type = "higher", nmads = nmads)
  message(sprintf("Sample %s: nFeature > %.1f, < %.1f | MT: covariate only",
                  sample_id,
                  attr(qc_lo, "thresholds")["lower"],
                  attr(qc_hi, "thresholds")["higher"]))
  cells[obj$nFeature_RNA > attr(qc_lo, "thresholds")["lower"] &
          obj$nFeature_RNA < attr(qc_hi, "thresholds")["higher"]]
}

suppressWarnings(
  print(VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                group.by = "Sample", ncol = 3))
)

samples    <- unique(data$Sample)
keep_cells <- unlist(lapply(samples, function(s) qc_per_sample(data, s)))
data       <- data[, keep_cells]
message(sprintf("Cells after QC: %d", ncol(data)))

# =============================================================================
# 4. BUILD GENE MAP FROM SNP LIST
# =============================================================================

all_gene_features     <- grep("-gene$", rownames(data), value = TRUE)
variant_ids_in_matrix <- sub("-gene$", "", all_gene_features)
base_variant_ids      <- sub("-[0-9]+$", "", variant_ids_in_matrix)

snp_list <- read.csv('snp_list.csv') %>%  
  mutate(across(everything(), ~trimws(as.character(.)))) %>%
  separate_longer_delim(Gene, delim = "/") %>%
  mutate(Gene = trimws(Gene), Gene = toupper(Gene)) %>%
  mutate(Gene = case_when(
    Gene == "IL8RBP" ~ "CXCR2",
    Gene == "EBI2"   ~ "GPR183",
    Gene == "TNFSF8" ~ "CD153",
    Gene == "PSGM1"  ~ "PSMG1",
    TRUE             ~ Gene
  )) %>%
  distinct(rs.ID, Gene, .keep_all = TRUE)

snp_matched <- snp_list %>% filter(rs.ID %in% base_variant_ids)
message(sprintf("Matched: %d variants across %d genes",
                nrow(snp_matched), length(unique(snp_matched$Gene))))
print(table(snp_matched$Gene))

grna_to_gene <- setNames(
  snp_matched$Gene[match(base_variant_ids, snp_matched$rs.ID)],
  variant_ids_in_matrix
)
grna_to_gene <- grna_to_gene[!is.na(grna_to_gene)]

gene_groups_grna <- split(names(grna_to_gene), grna_to_gene)

# =============================================================================
# 5. AUTO-DETECT CONTROLS
# =============================================================================

all_features  <- rownames(data)
nt_grnas      <- grep("no-control|non-targeting|rs2251039-ctrl|rs35675666-ctrl",
                      all_features, value = TRUE, ignore.case = TRUE)
posctrl_grnas <- grep("positive-ctrl", all_features, value = TRUE, ignore.case = TRUE)

message(sprintf("NT gRNAs: %d", length(nt_grnas)));       print(nt_grnas)
message(sprintf("Pos ctrl gRNAs: %d", length(posctrl_grnas))); print(posctrl_grnas)

grna_gene_blocks <- lapply(gene_groups_grna, function(v) paste0(v, "-gene"))
all_grna_names   <- unique(c(unlist(grna_gene_blocks), nt_grnas, posctrl_grnas))
posctrl_names    <- posctrl_grnas
nonposctrl_names <- setdiff(all_grna_names, posctrl_names)

# =============================================================================
# 6. IDENTIFY gRNA-POSITIVE CELLS
# =============================================================================

counts_all      <- GetAssayData(data, assay = "RNA", layer = "counts")
present_nonpos  <- intersect(nonposctrl_names, rownames(counts_all))
present_posctrl <- intersect(posctrl_names,    rownames(counts_all))

missing_grnas <- setdiff(all_grna_names, rownames(counts_all))
if (length(missing_grnas) > 0)
  message("WARNING — gRNAs not found:\n", paste(missing_grnas, collapse = "\n"))

has_nonpos_grna <- colSums(counts_all[present_nonpos, , drop = FALSE] > 0) > 0
has_posctrl <- if (length(present_posctrl) > 0) {
  colSums(counts_all[present_posctrl, , drop = FALSE] > 0) > 0
} else {
  rep(FALSE, ncol(counts_all))
}
has_nonpos_grna <- has_nonpos_grna[colnames(data)]
has_posctrl     <- has_posctrl[colnames(data)]

data <- data[, has_nonpos_grna & !has_posctrl]
message(sprintf("gRNA-positive cells: %d", ncol(data)))

# =============================================================================
# 7. NORMALIZE + SCALE
# =============================================================================

data <- data %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(.))

# =============================================================================
# 8. SPLIT MATRICES
# =============================================================================

grna_features  <- grep("-gene$", rownames(data), value = TRUE)
counts         <- GetAssayData(data, assay = "RNA", layer = "counts")
grna_matrix    <- counts[grna_features, ]
gene_matrix    <- counts[!rownames(counts) %in% grna_features, ]
lognorm_matrix <- GetAssayData(data, assay = "RNA", layer = "data")
lognorm_matrix <- lognorm_matrix[!rownames(lognorm_matrix) %in% grna_features, ]

# =============================================================================
# 9. BUILD grna_target_df
# =============================================================================

grna_target_df <- data.frame(
  grna_id     = rownames(grna_matrix),
  grna_target = sub("-gene$", "", rownames(grna_matrix)),
  stringsAsFactors = FALSE
)

posctrl_ids    <- sub("-gene$", "", posctrl_names)
grna_target_df <- grna_target_df[!grna_target_df$grna_target %in% posctrl_ids, , drop = FALSE]
grna_matrix    <- grna_matrix[rownames(grna_matrix) %in% grna_target_df$grna_id, ]

is_nt <- grepl("non-targeting|no-control|rs2251039-ctrl|rs35675666-ctrl",
               grna_target_df$grna_id, ignore.case = TRUE)
grna_target_df$grna_target[is_nt] <- "non-targeting"
grna_target_df <- grna_target_df %>%
  filter(!grepl("^cas9$|^puro$|^GFP$", grna_target, ignore.case = TRUE))
grna_matrix <- grna_matrix[rownames(grna_matrix) %in% grna_target_df$grna_id, ]

missing_genes <- setdiff(unique(snp_matched$Gene), rownames(gene_matrix))
if (length(missing_genes) > 0) {
  message(sprintf("WARNING — %d genes missing:", length(missing_genes)))
  print(missing_genes)
}

# =============================================================================
# 10. COVARIATES
# =============================================================================

n_batches  <- length(unique(data$Sample))
covariates <- data.frame(
  log_nUMI   = log1p(data$nCount_RNA),
  percent_mt = data$percent.mt,
  row.names  = colnames(data)
)
if (n_batches > 1) {
  covariates$batch <- factor(data$Sample)
  message(sprintf("Batch covariate added: %d samples", n_batches))
}

# =============================================================================
# 11. IDENTIFY CELL GROUPS
# =============================================================================

run_ttest_sceptre <- function(focal_gene,
                              treatment_cells,
                              all_sibling_grna_ids,   # used to exclude from control
                              focal_grna_target,       # label for SCEPTRE (single) or "combined"
                              focal_grna_ids,          # actual gRNA feature IDs for SCEPTRE
                              grna_target_df,
                              grna_matrix,
                              gene_matrix,
                              lognorm_matrix,
                              covariates,
                              test_label) {
  
  sibling_cells <- if (length(all_sibling_grna_ids) > 0)
    colnames(grna_matrix)[colSums(grna_matrix[
      intersect(all_sibling_grna_ids, rownames(grna_matrix)), , drop = FALSE] > 0) > 0]
  else character(0)
  
  control_cells <- setdiff(colnames(grna_matrix), union(treatment_cells, sibling_cells))
  
  n_trt  <- length(treatment_cells)
  n_ctrl <- length(control_cells)
  message(sprintf("  [%s] Treatment: %d | Siblings excluded: %d | Control: %d",
                  test_label, n_trt, length(sibling_cells), n_ctrl))
  
  if (n_trt < 3 || n_ctrl < 3) {
    message(sprintf("  [%s] SKIPPING — insufficient cells", test_label))
    return(list(ttest_p = NA_real_, ttest_log2fc = NA_real_,
                sceptre_p = NA_real_, sceptre_log2fc = NA_real_,
                n_treatment = n_trt, n_control = n_ctrl))
  }
  
  # ------------------------------------------------------------------
  # T-TEST (Welch) 
  # ------------------------------------------------------------------
  trt_expr  <- as.numeric(lognorm_matrix[focal_gene, treatment_cells])
  ctrl_expr <- as.numeric(lognorm_matrix[focal_gene, control_cells])
  
  ttest_p <- tryCatch(
    t.test(ctrl_expr, trt_expr, alternative = "two.sided", var.equal = FALSE)$p.value,
    error = function(e) NA_real_
  )
  
  ttest_log2fc <- log2((mean(expm1(trt_expr)) + 1) / (mean(expm1(ctrl_expr)) + 1))
  
  # Expressors only
  trt_pos    <- trt_expr[trt_expr   > 0]
  ctrl_pos   <- ctrl_expr[ctrl_expr > 0]
  n_trt_exp  <- length(trt_pos)
  n_ctrl_exp <- length(ctrl_pos)
  
  # Welch t-test expressors
  ttest_p_exp <- tryCatch(
    if (n_trt_exp >= 3 && n_ctrl_exp >= 3)
      t.test(ctrl_pos, trt_pos, alternative = "two.sided", var.equal = FALSE)$p.value
    else NA_real_,
    error = function(e) NA_real_
  )
  ttest_log2fc_exp <- if (n_trt_exp >= 1 && n_ctrl_exp >= 1)
    log2((mean(expm1(trt_pos)) + 1) / (mean(expm1(ctrl_pos)) + 1))
  else NA_real_
  
  # KS test — all cells
  ks_p <- tryCatch(
    ks.test(ctrl_expr, trt_expr)$p.value,
    error = function(e) NA_real_
  )
  
  # KS test — expressors only
  ks_p_exp <- tryCatch(
    if (n_trt_exp >= 3 && n_ctrl_exp >= 3)
      ks.test(ctrl_pos, trt_pos)$p.value
    else NA_real_,
    error = function(e) NA_real_
  )
  
  # ------------------------------------------------------------------
  # SCEPTRE — log_nUMI + percent_mt + batch
  # ------------------------------------------------------------------
  combined_target_name <- paste0(focal_grna_target)   
  
  grna_target_df_run <- grna_target_df
  
  
  other_targets <- setdiff(
    grna_target_df_run$grna_target[grna_target_df_run$grna_target != "non-targeting"],
    focal_grna_ids 
  )
  grna_target_df_run$grna_target[
    grna_target_df_run$grna_target %in% other_targets
  ] <- "non-targeting"
  
  grna_target_df_run$grna_target[
    grna_target_df_run$grna_target %in% focal_grna_ids
  ] <- combined_target_name
  
  discovery_pairs_run <- data.frame(
    grna_target = combined_target_name,
    response_id = focal_gene,
    stringsAsFactors = FALSE
  )
  
  sceptre_p      <- NA_real_
  sceptre_log2fc <- NA_real_
  
  if (combined_target_name %in% grna_target_df_run$grna_target) {
    result <- tryCatch({
      import_data(
        response_matrix        = gene_matrix,
        grna_matrix            = grna_matrix,
        grna_target_data_frame = grna_target_df_run[, c("grna_id", "grna_target")],
        moi                    = "low",
        extra_covariates       = covariates
      ) %>%
        set_analysis_parameters(
          discovery_pairs           = discovery_pairs_run,
          side                      = "both",
          control_group             = "complement",
          grna_integration_strategy = "singleton"
        ) %>%
        assign_grnas(method = "thresholding", threshold = 1) %>%
        run_qc(
          n_nonzero_trt_thresh     = 1L,
          n_nonzero_cntrl_thresh   = 1L,
          response_n_umis_range    = c(0, 1),
          response_n_nonzero_range = c(0, 1),
          p_mito_threshold         = 100
        ) %>%
        run_discovery_analysis(parallel = TRUE, n_processors = 8) %>%
        get_result("run_discovery_analysis")
    }, error = function(e) {
      message(sprintf("  [%s] SCEPTRE error: %s", test_label, e$message)); NULL
    })
    
    if (!is.null(result) && nrow(result) > 0) {
      sceptre_p      <- result$p_value[1]
      sceptre_log2fc <- result$log_2_fold_change[1]
    }
  }
  
  list(ttest_p         = ttest_p,
       ttest_log2fc    = ttest_log2fc,
       ttest_p_exp     = ttest_p_exp,
       ttest_log2fc_exp= ttest_log2fc_exp,
       ks_p            = ks_p,
       ks_p_exp        = ks_p_exp,
       n_trt_exp       = n_trt_exp,
       n_ctrl_exp      = n_ctrl_exp,
       sceptre_p       = sceptre_p,
       sceptre_log2fc  = sceptre_log2fc,
       n_treatment     = n_trt,
       n_control       = n_ctrl)
}

# =============================================================================
# 12. All variants
# =============================================================================

run_all_tests_for_variant <- function(base_rsid,
                                      focal_gene,
                                      gene_groups_grna,
                                      grna_target_df,
                                      grna_matrix,
                                      gene_matrix,
                                      lognorm_matrix,
                                      covariates) {
  
  message(sprintf("\n====== %s -> %s ======", base_rsid, focal_gene))
  
  all_variant_targets <- names(grna_to_gene)[
    grna_to_gene == focal_gene &
      grepl(paste0("^", base_rsid, "-"), names(grna_to_gene))
  ]
  
  grna1_target <- paste0(base_rsid, "-1")
  grna2_target <- paste0(base_rsid, "-2")
  grna1_id     <- paste0(grna1_target, "-gene")
  grna2_id     <- paste0(grna2_target, "-gene")
  
  has_grna1 <- grna1_id %in% rownames(grna_matrix)
  has_grna2 <- grna2_id %in% rownames(grna_matrix)
  
  if (!has_grna1 && !has_grna2) {
    message("  No gRNAs found in matrix — skipping"); return(NULL)
  }
  if (!focal_gene %in% rownames(lognorm_matrix)) {
    message("  Gene not in matrix — skipping"); return(NULL)
  }
  
  all_gene_grna_ids <- paste0(gene_groups_grna[[focal_gene]], "-gene")
  
  get_cells <- function(ids) {
    ids <- intersect(ids, rownames(grna_matrix))
    if (length(ids) == 0) return(character(0))
    colnames(grna_matrix)[colSums(grna_matrix[ids, , drop = FALSE] > 0) > 0]
  }
  
  cells_grna1 <- if (has_grna1) get_cells(grna1_id) else character(0)
  cells_grna2 <- if (has_grna2) get_cells(grna2_id) else character(0)
  cells_union <- union(cells_grna1, cells_grna2)
  cells_inter <- intersect(cells_grna1, cells_grna2)
  
  results <- list()
  
  if (has_grna1 && length(cells_grna1) >= 3) {
    siblings_excl <- all_gene_grna_ids_excl
    results$grna1 <- run_ttest_sceptre(
      focal_gene       = focal_gene,
      treatment_cells  = cells_grna1,
      all_sibling_grna_ids = siblings_excl,
      focal_grna_target = grna1_target,
      focal_grna_ids   = grna1_target,
      grna_target_df   = grna_target_df,
      grna_matrix      = grna_matrix,
      gene_matrix      = gene_matrix,
      lognorm_matrix   = lognorm_matrix,
      covariates       = covariates,
      test_label       = "gRNA1"
    )
  } else {
    results$grna1 <- list(ttest_p=NA_real_, ttest_log2fc=NA_real_,
                          sceptre_p=NA_real_, sceptre_log2fc=NA_real_,
                          n_treatment=length(cells_grna1), n_control=NA_integer_)
  }
  
  if (has_grna2 && length(cells_grna2) >= 3) {
    siblings_excl <- all_gene_grna_ids_excl   # same clean control pool
    results$grna2 <- run_ttest_sceptre(
      focal_gene       = focal_gene,
      treatment_cells  = cells_grna2,
      all_sibling_grna_ids = siblings_excl,
      focal_grna_target = grna2_target,
      focal_grna_ids   = grna2_target,
      grna_target_df   = grna_target_df,
      grna_matrix      = grna_matrix,
      gene_matrix      = gene_matrix,
      lognorm_matrix   = lognorm_matrix,
      covariates       = covariates,
      test_label       = "gRNA2"
    )
  } else {
    results$grna2 <- list(ttest_p=NA_real_, ttest_log2fc=NA_real_,
                          sceptre_p=NA_real_, sceptre_log2fc=NA_real_,
                          n_treatment=length(cells_grna2), n_control=NA_integer_)
  }
  
  if (has_grna1 && has_grna2 && length(cells_union) >= 3) {
    combined_label <- paste0(base_rsid, "-combined_union")
    results$combined_union <- run_ttest_sceptre(
      focal_gene       = focal_gene,
      treatment_cells  = cells_union,
      all_sibling_grna_ids = all_gene_grna_ids_excl,  # exclude all gene gRNAs from control
      focal_grna_target = combined_label,
      focal_grna_ids   = c(grna1_target, grna2_target),
      grna_target_df   = grna_target_df,
      grna_matrix      = grna_matrix,
      gene_matrix      = gene_matrix,
      lognorm_matrix   = lognorm_matrix,
      covariates       = covariates,
      test_label       = "gRNA1+gRNA2 union"
    )
  } else {
    results$combined_union <- list(ttest_p=NA_real_, ttest_log2fc=NA_real_,
                                   sceptre_p=NA_real_, sceptre_log2fc=NA_real_,
                                   n_treatment=length(cells_union), n_control=NA_integer_)
  }
  
  if (has_grna1 && has_grna2 && length(cells_inter) >= 3) {
    combined_label_inter <- paste0(base_rsid, "-combined_inter")
    results$combined_inter <- run_ttest_sceptre(
      focal_gene       = focal_gene,
      treatment_cells  = cells_inter,
      all_sibling_grna_ids = all_gene_grna_ids_excl,  # exclude all gene gRNAs from control
      focal_grna_target = combined_label_inter,
      focal_grna_ids   = c(grna1_target, grna2_target),
      grna_target_df   = grna_target_df,
      grna_matrix      = grna_matrix,
      gene_matrix      = gene_matrix,
      lognorm_matrix   = lognorm_matrix,
      covariates       = covariates,
      test_label       = "gRNA1+gRNA2 intersection"
    )
  } else {
    results$combined_inter <- list(ttest_p=NA_real_, ttest_log2fc=NA_real_,
                                   sceptre_p=NA_real_, sceptre_log2fc=NA_real_,
                                   n_treatment=length(cells_inter), n_control=NA_integer_)
  }
  
  g <- function(res, field, default = NA_real_) {
    val <- res[[field]]
    if (is.null(val) || length(val) == 0) return(default)
    val[[1]]
  }
  gi <- function(res, field) g(res, field, NA_integer_)
  
  g1 <- results$grna1
  g2 <- results$grna2
  cu <- results$combined_union
  ci <- results$combined_inter
  
  data.frame(
    base_variant               = base_rsid,
    gene                       = focal_gene,
    n_grnas                    = has_grna1 + has_grna2,
    # gRNA1
    n_trt_grna1                = gi(g1, "n_treatment"),
    ttest_p_grna1              = g(g1, "ttest_p"),
    ttest_log2fc_grna1         = g(g1, "ttest_log2fc"),
    ttest_p_exp_grna1          = g(g1, "ttest_p_exp"),
    ttest_log2fc_exp_grna1     = g(g1, "ttest_log2fc_exp"),
    ks_p_grna1                 = g(g1, "ks_p"),
    ks_p_exp_grna1             = g(g1, "ks_p_exp"),
    n_trt_exp_grna1            = gi(g1, "n_trt_exp"),
    sceptre_p_grna1            = g(g1, "sceptre_p"),
    sceptre_log2fc_grna1       = g(g1, "sceptre_log2fc"),
    # gRNA2
    n_trt_grna2                = gi(g2, "n_treatment"),
    ttest_p_grna2              = g(g2, "ttest_p"),
    ttest_log2fc_grna2         = g(g2, "ttest_log2fc"),
    ttest_p_exp_grna2          = g(g2, "ttest_p_exp"),
    ttest_log2fc_exp_grna2     = g(g2, "ttest_log2fc_exp"),
    ks_p_grna2                 = g(g2, "ks_p"),
    ks_p_exp_grna2             = g(g2, "ks_p_exp"),
    n_trt_exp_grna2            = gi(g2, "n_trt_exp"),
    sceptre_p_grna2            = g(g2, "sceptre_p"),
    sceptre_log2fc_grna2       = g(g2, "sceptre_log2fc"),
    # Combined union
    n_trt_union                = gi(cu, "n_treatment"),
    ttest_p_union              = g(cu, "ttest_p"),
    ttest_log2fc_union         = g(cu, "ttest_log2fc"),
    ttest_p_exp_union          = g(cu, "ttest_p_exp"),
    ks_p_union                 = g(cu, "ks_p"),
    ks_p_exp_union             = g(cu, "ks_p_exp"),
    sceptre_p_union            = g(cu, "sceptre_p"),
    sceptre_log2fc_union       = g(cu, "sceptre_log2fc"),
    # Combined intersection
    n_trt_inter                = gi(ci, "n_treatment"),
    ttest_p_inter              = g(ci, "ttest_p"),
    ttest_log2fc_inter         = g(ci, "ttest_log2fc"),
    ttest_p_exp_inter          = g(ci, "ttest_p_exp"),
    ks_p_inter                 = g(ci, "ks_p"),
    ks_p_exp_inter             = g(ci, "ks_p_exp"),
    sceptre_p_inter            = g(ci, "sceptre_p"),
    sceptre_log2fc_inter       = g(ci, "sceptre_log2fc"),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# 13. BUILD RUN LIST 
# =============================================================================

variant_run_list <- grna_target_df %>%
  filter(grepl("-[0-9]+$", grna_target),
         grna_target != "non-targeting") %>%
  mutate(
    base_rsid = sub("-[0-9]+$", "", grna_target),
    gene      = grna_to_gene[grna_target]
  ) %>%
  filter(!is.na(gene), gene %in% rownames(gene_matrix)) %>%
  distinct(base_rsid, gene)   # one entry per variant+gene — the function handles both gRNAs

message(sprintf("Running %d variant+gene pairs across %d genes",
                nrow(variant_run_list), length(unique(variant_run_list$gene))))

# =============================================================================
# 14. RUN
# =============================================================================

all_results <- mapply(
  run_all_tests_for_variant,
  base_rsid  = variant_run_list$base_rsid,
  focal_gene = variant_run_list$gene,
  MoreArgs = list(
    gene_groups_grna = gene_groups_grna,
    grna_target_df   = grna_target_df,
    grna_matrix      = grna_matrix,
    gene_matrix      = gene_matrix,
    lognorm_matrix   = lognorm_matrix,
    covariates       = covariates
  ),
  SIMPLIFY = FALSE
)
all_results <- Filter(Negate(is.null), all_results)
message(sprintf("Completed %d / %d variant+gene pairs",
                length(all_results), nrow(variant_run_list)))

# =============================================================================
# 15. VARIANT-LEVEL RESULTS + BH CORRECTION
# =============================================================================

variant_results <- bind_rows(all_results) %>%
  mutate(
    # Nominal sig flags (p < 0.05) per test
    ttest_sig_grna1        = !is.na(ttest_p_grna1)   & ttest_p_grna1   < 0.05,
    ttest_sig_grna2        = !is.na(ttest_p_grna2)   & ttest_p_grna2   < 0.05,
    ttest_sig_union        = !is.na(ttest_p_union)   & ttest_p_union   < 0.05,
    ttest_sig_inter        = !is.na(ttest_p_inter)   & ttest_p_inter   < 0.05,
    sceptre_sig_grna1      = !is.na(sceptre_p_grna1) & sceptre_p_grna1 < 0.05,
    sceptre_sig_grna2      = !is.na(sceptre_p_grna2) & sceptre_p_grna2 < 0.05,
    sceptre_sig_union      = !is.na(sceptre_p_union) & sceptre_p_union < 0.05,
    sceptre_sig_inter      = !is.na(sceptre_p_inter) & sceptre_p_inter < 0.05,
    
    # Expressor t-test sig flags
    ttest_sig_exp_grna1    = !is.na(ttest_p_exp_grna1) & ttest_p_exp_grna1 < 0.05,
    ttest_sig_exp_grna2    = !is.na(ttest_p_exp_grna2) & ttest_p_exp_grna2 < 0.05,
    ttest_sig_exp_union    = !is.na(ttest_p_exp_union)  & ttest_p_exp_union  < 0.05,
    ttest_sig_exp_inter    = !is.na(ttest_p_exp_inter)  & ttest_p_exp_inter  < 0.05,
    # KS sig flags
    ks_sig_grna1           = !is.na(ks_p_grna1)  & ks_p_grna1  < 0.05,
    ks_sig_grna2           = !is.na(ks_p_grna2)  & ks_p_grna2  < 0.05,
    ks_sig_union           = !is.na(ks_p_union)   & ks_p_union   < 0.05,
    ks_sig_inter           = !is.na(ks_p_inter)   & ks_p_inter   < 0.05,
    ks_sig_exp_grna1       = !is.na(ks_p_exp_grna1) & ks_p_exp_grna1 < 0.05,
    ks_sig_exp_grna2       = !is.na(ks_p_exp_grna2) & ks_p_exp_grna2 < 0.05,
    
    # Any test sig
    ttest_sig_any          = ttest_sig_grna1 | ttest_sig_grna2 | ttest_sig_union | ttest_sig_inter,
    ttest_sig_exp_any      = ttest_sig_exp_grna1 | ttest_sig_exp_grna2 | ttest_sig_exp_union | ttest_sig_exp_inter,
    ks_sig_any             = ks_sig_grna1 | ks_sig_grna2 | ks_sig_union | ks_sig_inter,
    sceptre_sig_any        = sceptre_sig_grna1 | sceptre_sig_grna2 | sceptre_sig_union | sceptre_sig_inter,
    
    # Minimum p across all tests (for BH correction)
    ttest_p_min            = pmin(ttest_p_grna1, ttest_p_grna2, ttest_p_union, ttest_p_inter,
                                  na.rm = TRUE),
    ttest_p_exp_min        = pmin(ttest_p_exp_grna1, ttest_p_exp_grna2, ttest_p_exp_union, ttest_p_exp_inter,
                                  na.rm = TRUE),
    ks_p_min               = pmin(ks_p_grna1, ks_p_grna2, ks_p_union, ks_p_inter,
                                  na.rm = TRUE),
    ks_p_exp_min           = pmin(ks_p_exp_grna1, ks_p_exp_grna2, na.rm = TRUE),
    sceptre_p_min          = pmin(sceptre_p_grna1, sceptre_p_grna2, sceptre_p_union, sceptre_p_inter,
                                  na.rm = TRUE),
    
    # Mean log2FC across available individual tests (grna1 + grna2, not combined)
    mean_ttest_log2fc      = rowMeans(cbind(ttest_log2fc_grna1, ttest_log2fc_grna2), na.rm = TRUE),
    mean_sceptre_log2fc    = rowMeans(cbind(sceptre_log2fc_grna1, sceptre_log2fc_grna2), na.rm = TRUE),
    same_direction         = sign(ttest_log2fc_grna1) == sign(ttest_log2fc_grna2) |
      is.na(ttest_log2fc_grna1) | is.na(ttest_log2fc_grna2)
  ) %>%
  mutate(
    ttest_p_min_adj        = p.adjust(ttest_p_min,     method = "BH"),
    ttest_p_exp_min_adj    = p.adjust(ttest_p_exp_min, method = "BH"),
    ks_p_min_adj           = p.adjust(ks_p_min,        method = "BH"),
    ks_p_exp_min_adj       = p.adjust(ks_p_exp_min,    method = "BH"),
    sceptre_p_min_adj      = p.adjust(sceptre_p_min,   method = "BH"),
    ttest_sig_bh_any       = ttest_p_min_adj     < 0.05,
    ttest_sig_exp_bh_any   = ttest_p_exp_min_adj < 0.05,
    ks_sig_bh_any          = ks_p_min_adj        < 0.05,
    ks_sig_exp_bh_any      = ks_p_exp_min_adj    < 0.05,
    sceptre_sig_bh_any     = sceptre_p_min_adj   < 0.05
  ) %>%
  arrange(ttest_p_min)

# =============================================================================
# 16. SUMMARY
# =============================================================================

cat("\n=== VARIANT-LEVEL RESULTS — nominal p < 0.05 ===\n")
cat(sprintf("Total variants tested: %d\n\n", nrow(variant_results)))

for (test_label in c("gRNA1", "gRNA2", "Union", "Inter")) {
  ttest_col   <- paste0("ttest_sig_",   tolower(test_label))
  sceptre_col <- paste0("sceptre_sig_", tolower(test_label))
  cat(sprintf("%-20s T-test: %d sig | SCEPTRE: %d sig\n",
              test_label,
              sum(variant_results[[ttest_col]],   na.rm = TRUE),
              sum(variant_results[[sceptre_col]], na.rm = TRUE)))
}
cat(sprintf("\n%-20s T-test: %d sig | SCEPTRE: %d sig\n",
            "ANY test",
            sum(variant_results$ttest_sig_any,   na.rm = TRUE),
            sum(variant_results$sceptre_sig_any, na.rm = TRUE)))

# =============================================================================
# 17. SAVE
# =============================================================================

write.csv(variant_results, paste0(POOL_NAME, "_variant_results.csv"), row.names = FALSE)
cat(sprintf("\nSaved: %s_variant_results.csv\n", POOL_NAME))

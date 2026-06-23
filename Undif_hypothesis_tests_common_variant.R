# =============================================================================
# SCEPTRE + T-TEST — eCROPseq Pool (multi-replicate)
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

pool_paths  <- list(POOL_PATH_1, POOL_PATH_2, POOL_PATH_3, POOL_PATH_4, POOL_PATH_5)
rep_labels  <- c("Rep1", "Rep2", "Rep3", "Rep4", "Rep5")
pool_paths  <- Filter(Negate(is.null), pool_paths)
n_pools     <- length(pool_paths)

seurat_list <- mapply(load_pool, pool_paths, rep_labels[seq_len(n_pools)],
                      SIMPLIFY = FALSE)

if (n_pools == 1) {
  data <- seurat_list[[1]]
  message("Single pool loaded.")
} else {
  data <- merge(seurat_list[[1]],
                y            = seurat_list[-1],
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
# 3. QC 
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

snp_list <- read.csv('~/path/to/snp_list.csv') %>%
  mutate(across(everything(), ~trimws(as.character(.)))) %>%
  separate_longer_delim(Gene, delim = "/") %>%
  mutate(Gene = trimws(Gene)) %>%
  filter(!is.na(Gene) & nchar(Gene) > 0) %>%
  mutate(Gene = toupper(Gene)) %>%
  mutate(Gene = case_when(
    Gene == "IL8RBP" ~ "CXCR2",
    Gene == "EBI2"   ~ "GPR183",
    Gene == "TNFSF8" ~ "CD153",
    Gene == "PSGM1"  ~ "PSMG1",
    TRUE             ~ Gene
  )) %>%
  distinct(rs.ID, Gene, .keep_all = TRUE)

snp_matched <- snp_list %>% filter(rs.ID %in% variant_ids_in_matrix)

message(sprintf("Variants in matrix: %d", length(variant_ids_in_matrix)))
message(sprintf("Matched: %d variants across %d genes",
                nrow(snp_matched), length(unique(snp_matched$Gene))))
print(table(snp_matched$Gene))

variant_to_gene  <- setNames(as.character(snp_matched$Gene),
                             as.character(snp_matched$rs.ID))
gene_groups      <- split(as.character(snp_matched$rs.ID),
                          as.character(snp_matched$Gene))
grna_gene_blocks <- lapply(gene_groups, function(v) paste0(v, "-gene"))

# =============================================================================
# 5. AUTO-DETECT CONTROLS
# =============================================================================

all_features  <- rownames(data)
nt_grnas      <- grep(
  "no-control|non-targeting-control|non-targeting-ctrl|rs2251039-ctrl|rs35675666-ctrl",
  all_features, value = TRUE, ignore.case = TRUE)
posctrl_grnas <- grep(
  "positive-ctrl", all_features, value = TRUE, ignore.case = TRUE)

message(sprintf("NT gRNAs: %d", length(nt_grnas)));      print(nt_grnas)
message(sprintf("Pos ctrl gRNAs: %d", length(posctrl_grnas))); print(posctrl_grnas)

grna_names       <- c(grna_gene_blocks, list(NT = nt_grnas, POSCTRL = posctrl_grnas))
all_grna_names   <- unique(unlist(grna_names, use.names = FALSE))
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
grna_target_df <- grna_target_df[
  !grna_target_df$grna_target %in% posctrl_ids, , drop = FALSE]
grna_matrix    <- grna_matrix[rownames(grna_matrix) %in% grna_target_df$grna_id, ]

is_nt <- grepl(
  "non-targeting|no-control|rs2251039-ctrl|rs35675666-ctrl|non-targeting-ctrl|non-targeting-control",
  grna_target_df$grna_id, ignore.case = TRUE)
grna_target_df$grna_target[is_nt] <- "non-targeting"

grna_target_df <- grna_target_df %>%
  filter(!grepl("^cas9$|^puro$|^GFP$", grna_target, ignore.case = TRUE))
grna_matrix <- grna_matrix[rownames(grna_matrix) %in% grna_target_df$grna_id, ]

missing_check <- setdiff(names(variant_to_gene), grna_target_df$grna_target)
cat("Variants missing from grna_target_df:",
    if (length(missing_check) == 0) "none\n"
    else paste(missing_check, collapse = ", "), "\n")

missing_genes <- setdiff(unique(snp_matched$Gene), rownames(gene_matrix))
if (length(missing_genes) > 0) {
  message(sprintf("WARNING — %d genes missing:", length(missing_genes)))
  print(missing_genes)
}

variant_to_gene <- variant_to_gene[variant_to_gene %in% rownames(gene_matrix)]
gene_groups     <- gene_groups[names(gene_groups) %in% rownames(gene_matrix)]

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

get_cell_groups <- function(focal_variant, focal_gene,
                            gene_groups, grna_target_df, grna_matrix) {
  focal_grna_id <- paste0(focal_variant, "-gene")
  if (!focal_grna_id %in% rownames(grna_matrix)) return(NULL)
  
  treatment_cells  <- colnames(grna_matrix)[grna_matrix[focal_grna_id, ] > 0]
  sibling_variants <- setdiff(gene_groups[[focal_gene]], focal_variant)
  sibling_grna_ids <- intersect(paste0(sibling_variants, "-gene"), rownames(grna_matrix))
  
  sibling_cells <- if (length(sibling_grna_ids) > 0)
    colnames(grna_matrix)[colSums(grna_matrix[sibling_grna_ids, , drop = FALSE] > 0) > 0]
  else character(0)
  
  control_cells <- setdiff(colnames(grna_matrix), union(treatment_cells, sibling_cells))
  list(treatment = treatment_cells, control = control_cells, siblings = sibling_cells)
}

# =============================================================================
# 12. PER-VARIANT COMPARISON FUNCTION
# =============================================================================

run_comparison_for_variant <- function(focal_variant, focal_gene,
                                       gene_groups, grna_target_df,
                                       grna_matrix, gene_matrix,
                                       lognorm_matrix, covariates) {
  
  message(sprintf("\n=== %s (%s) ===", focal_variant, focal_gene))
  
  groups <- get_cell_groups(focal_variant, focal_gene,
                            gene_groups, grna_target_df, grna_matrix)
  if (is.null(groups)) { message("  SKIPPING — gRNA not found"); return(NULL) }
  
  n_trt  <- length(groups$treatment)
  n_ctrl <- length(groups$control)
  message(sprintf("  Treatment: %d | Control: %d", n_trt, n_ctrl))
  
  if (n_trt < 3 || n_ctrl < 3) { message("  SKIPPING — insufficient cells"); return(NULL) }
  if (!focal_gene %in% rownames(lognorm_matrix)) { message("  SKIPPING — gene not in matrix"); return(NULL) }
  
  # ------------------------------------------------------------------
  # T-TEST 
  # ------------------------------------------------------------------
  
  trt_expr  <- as.numeric(lognorm_matrix[focal_gene, groups$treatment])
  ctrl_expr <- as.numeric(lognorm_matrix[focal_gene, groups$control])
  
  ttest_p <- tryCatch(
    t.test(ctrl_expr, trt_expr, alternative = "two.sided")$p.value,
    error = function(e) NA_real_
  )
  
  trt_norm  <- expm1(trt_expr)
  ctrl_norm <- expm1(ctrl_expr)
  ttest_log2fc <- log2((mean(trt_norm) + 1) / (mean(ctrl_norm) + 1))
  
  # ------------------------------------------------------------------
  # SCEPTRE — percent_mt + log_nUMI + batch as covariates
  # ------------------------------------------------------------------
  
  other_gene_variants <- unlist(gene_groups[names(gene_groups) != focal_gene])
  grna_target_df_run  <- grna_target_df
  grna_target_df_run$grna_target[
    grna_target_df_run$grna_target %in% other_gene_variants
  ] <- "non-targeting"
  
  discovery_pairs_run <- data.frame(
    grna_target = focal_variant,
    response_id = focal_gene,
    stringsAsFactors = FALSE
  )
  
  sceptre_p      <- NA_real_
  sceptre_log2fc <- NA_real_
  
  if (focal_variant %in% grna_target_df_run$grna_target) {
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
      message(sprintf("  SCEPTRE error: %s", e$message))
      NULL
    })
    
    if (!is.null(result) && nrow(result) > 0) {
      sceptre_p      <- result$p_value[1]
      sceptre_log2fc <- result$log_2_fold_change[1]
    }
  } else {
    message("  SKIPPING SCEPTRE — variant not in grna_target_df")
  }
  
  data.frame(
    variant           = focal_variant,
    gene              = focal_gene,
    n_treatment_cells = n_trt,
    n_control_cells   = n_ctrl,
    n_expressing_trt  = sum(trt_expr  > 0),
    n_expressing_ctrl = sum(ctrl_expr > 0),
    ttest_p           = ttest_p,
    ttest_log2fc      = ttest_log2fc,
    ttest_sig_raw     = !is.na(ttest_p) & ttest_p < 0.05,
    sceptre_p         = sceptre_p,
    sceptre_log2fc    = sceptre_log2fc,
    sceptre_sig_raw   = !is.na(sceptre_p) & sceptre_p < 0.05,
    stringsAsFactors  = FALSE
  )
}

# =============================================================================
# 13. RUN ALL VARIANTS
# =============================================================================

variant_list <- data.frame(
  variant = names(variant_to_gene),
  gene    = as.character(variant_to_gene),
  stringsAsFactors = FALSE
)

message(sprintf("Running %d per-variant comparisons...", nrow(variant_list)))

all_results <- mapply(
  run_comparison_for_variant,
  focal_variant = variant_list$variant,
  focal_gene    = variant_list$gene,
  MoreArgs = list(
    gene_groups    = gene_groups,
    grna_target_df = grna_target_df,
    grna_matrix    = grna_matrix,
    gene_matrix    = gene_matrix,
    lognorm_matrix = lognorm_matrix,
    covariates     = covariates
  ),
  SIMPLIFY = FALSE
)

all_results <- Filter(Negate(is.null), all_results)
message(sprintf("Completed %d / %d", length(all_results), nrow(variant_list)))

# =============================================================================
# 14. COMBINE + BH CORRECTION
# =============================================================================

results_combined <- bind_rows(all_results) %>%
  mutate(
    ttest_p_adj    = p.adjust(ttest_p,   method = "BH"),
    ttest_sig_bh   = ttest_p_adj  < 0.05,
    sceptre_p_adj  = p.adjust(sceptre_p, method = "BH"),
    sceptre_sig_bh = sceptre_p_adj < 0.05
  ) %>%
  arrange(sceptre_p)

print(results_combined)
write.csv(results_combined,paste0(POOL_NAME, "_sceptre_vs_ttest.csv"),row.names = FALSE)

cat("\n=== DISCOVERY SUMMARY (BH < 0.05) ===\n")
cat(sprintf("SCEPTRE sig: %d\n", sum(results_combined$sceptre_sig_bh, na.rm = TRUE)))
cat(sprintf("T-test sig:  %d\n", sum(results_combined$ttest_sig_bh,   na.rm = TRUE)))
cat(sprintf("Both sig:    %d\n",
            sum(results_combined$sceptre_sig_bh &
                  results_combined$ttest_sig_bh,   na.rm = TRUE)))
cat(sprintf("\nSECPTRE NA rate: %.1f%% (%d / %d)\n",
            100 * mean(is.na(results_combined$sceptre_p)),
            sum(is.na(results_combined$sceptre_p)),
            nrow(results_combined)))

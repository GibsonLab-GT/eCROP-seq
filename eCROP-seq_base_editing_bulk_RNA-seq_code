library(ensembldb)
library(AnnotationHub)
library(DESeq2)
library(tximport)
library(biomaRt)
library(pheatmap)
library(clusterProfiler)
library(msigdbr)
library(ggplot2)
library(EnhancedVolcano)
library(dplyr)
library(tidyverse)
library(edgeR)
library(limma)

ah <- AnnotationHub()

human_datasets <- query(ah, c("EnsDb.Hsapiens"))

edb <- ah[["AH98047"]]

tx2gene <- transcripts(edb, return.type = "DataFrame")[, c("tx_id", "gene_id")]
tx2gene <- as.data.frame(tx2gene)

##################### INDIVIDUAL COMPARSION #######################
s1.1 <- '/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot/salmon_quant/output_quant_sample1-1/quant.sf'
s1.2 <- '/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot/salmon_quant/output_quant_sample1-2/quant.sf'
s2.1 <- '/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot/salmon_quant/output_quant_sample2-1/quant.sf'
s2.2 <- '/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot/salmon_quant/output_quant_sample2-2/quant.sf'
s3.1 <- '/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot/salmon_quant/output_quant_sample3-1/quant.sf'
s3.2 <- '/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot/salmon_quant/output_quant_sample3-2/quant.sf'
s4.1 <- '/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot/salmon_quant/output_quant_sample4-1/quant.sf'
s4.2 <- '/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot/salmon_quant/output_quant_sample4-2/quant.sf'
s5.1 <- '/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot/salmon_quant/output_quant_sample5-1/quant.sf'
s5.2 <- '/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot/salmon_quant/output_quant_sample5-2/quant.sf'
s6.1 <- '/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot/salmon_quant/output_quant_sample6-1/quant.sf'
s6.2 <- '/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot/salmon_quant/output_quant_sample6-2/quant.sf'
s7.1 <- '/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot/salmon_quant/output_quant_sample7-1/quant.sf'
s7.2 <- '/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot/salmon_quant/output_quant_sample7-2/quant.sf'
s8.1 <- '/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot/salmon_quant/output_quant_sample8-1/quant.sf'
s8.2 <- '/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot/salmon_quant/output_quant_sample8-2/quant.sf'

txi_ADCY3_N1 <- tximport(c(s4.1,s4.2,s3.1,s3.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_ADCY3_1 <- tximport(c(s5.1,s5.2,s3.1,s3.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_IFNAR2_1 <- tximport(c(s8.1,s8.2,s3.1,s3.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_IFNAR2_2 <- tximport(c(s6.1,s6.2,s3.1,s3.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_IFNAR2_N1 <- tximport(c(s7.1,s7.2,s3.1,s3.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_IFNAR2_N2 <- tximport(c(s1.1,s1.2,s3.1,s3.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_all <- tximport(c(s1.1,s1.2,s2.1, s2.2, s3.1,s3.2, s4.1, s4.2, s5.1, s5.2, s6.1, s6.2, s7.1, s7.2, s8.1, s8.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

sampleTable_ADCY3_N1 <- data.frame(
  samples = c("ADCY3_N1_1", "ADCY3_N1_2", "NC1_1", "NC1_2"),
  condition = c("ADCY3_N1", "ADCY3_N1", "NC1", "NC1")
)
rownames(sampleTable_ADCY3_N1) <- sampleTable_ADCY3_N1$samples

sampleTable_ADCY3_1 <- data.frame(
  samples = c("ADCY3_1_1", "ADCY3_1_2", "NC1_1", "NC1_2"),
  condition = c("ADCY3_1", "ADCY3_1", "NC1", "NC1")
)
rownames(sampleTable_ADCY3_1) <- sampleTable_ADCY3_1$samples

sampleTable_IFNAR2_1 <- data.frame(
  samples = c("IFNAR2_1_1", "IFNAR2_1_2", "NC1_1", "NC1_2"),
  condition = c("IFNAR2_1", "IFNAR2_1", "NC1", "NC1")
)
rownames(sampleTable_IFNAR2_1) <- sampleTable_IFNAR2_1$samples

sampleTable_IFNAR2_2 <- data.frame(
  samples = c("IFNAR2_2_1", "IFNAR2_2_2", "NC1_1", "NC1_2"),
  condition = c("IFNAR2_2", "IFNAR2_2", "NC1", "NC1")
)
rownames(sampleTable_IFNAR2_2) <- sampleTable_IFNAR2_2$samples

sampleTable_IFNAR2_N1 <- data.frame(
  samples = c("IFNAR2_N1_1", "IFNAR2_N1_2", "NC1_1", "NC1_2"),
  condition = c("IFNAR2_N1", "IFNAR2_N1", "NC1", "NC1")
)
rownames(sampleTable_IFNAR2_N1) <- sampleTable_IFNAR2_N1$samples

sampleTable_IFNAR2_N2 <- data.frame(
  samples = c("IFNAR2_N2_1", "IFNAR2_N2_2", "NC1_1", "NC1_2"),
  condition = c("IFNAR2_N2", "IFNAR2_N2", "NC1", "NC1")
)
rownames(sampleTable_IFNAR2_N2) <- sampleTable_IFNAR2_N2$samples

txi_names <- list(txi_ADCY3_N1,txi_ADCY3_1, txi_IFNAR2_1, txi_IFNAR2_2, txi_IFNAR2_N1, txi_IFNAR2_N2)
res_names <- c("res_ADCY3_N1", "res_ADCY3_1", "res_IFNAR2_1", "res_IFNAR2_2", "res_IFNAR2_N1", "res_IFNAR2_N2")
dds_names <- c("dds_ADCY3_N1", "dds_ADCY3_1", "dds_IFNAR2_1", "dds_IFNAR2_2", "dds_IFNAR2_N1", "dds_IFNAR2_N2")
res_df_names <- c("res_df_ADCY3_N1", "res_df_ADCY3_1", "res_df_IFNAR2_1", "res_df_IFNAR2_2", "res_df_IFNAR2_N1", "res_df_IFNAR2_N2")
sample_names <- list(
  sampleTable_ADCY3_N1,
  sampleTable_ADCY3_1,
  sampleTable_IFNAR2_1,
  sampleTable_IFNAR2_2,
  sampleTable_IFNAR2_N1,
  sampleTable_IFNAR2_N2
)

## DEG analysis
i = 1
for (x in txi_names){

  dds <- DESeqDataSetFromTximport(x, colData = sample_names[i], design = ~ condition)
  dds$condition <- relevel(dds$condition, ref = "NC1")
  dds <- DESeq(dds)
  res <- results(dds)

  res$padj <- p.adjust(res$pvalue, method="BH")

  assign(res_names[i], res)
  assign(dds_names[i], dds)

  i = i +1

}

res_vars <- list(
  res_ADCY3_N1,
  res_ADCY3_1,
  res_IFNAR2_1,
  res_IFNAR2_2,
  res_IFNAR2_N1,
  res_IFNAR2_N2
)

dds_vars <- list(
  dds_ADCY3_N1,
  dds_ADCY3_1,
  dds_IFNAR2_1,
  dds_IFNAR2_2,
  dds_IFNAR2_N1,
  dds_IFNAR2_N2
)

i = 1
for (x in txi_names){

  res_table_nam <- paste("/media/ming/Extra_SSD_4TB1/Linlin_bulkRNA/DESeq2_results", res_names[i], '.csv',sep = '_')
  path <- "/media/ming/Extra_SSD_4TB1/Linlin_bulkRNA/"
  write.csv(as.data.frame(res_vars[i]), file = res_table_nam)

  i = i +1
}

i = 1
for (x in res_vars){
  plotMA(x, main=res_names[i])
  i = i +1
}

i = 1
for (x in dds_vars){
  vsd <- vst(x)
  plotPCA(vsd, intgroup = "condition")
  i = i +1
}

## output tables with gene name
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

i = 1
for (x in res_vars) {
  deg_genes <- rownames(x)

  gene_conversion <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                           filters = "ensembl_gene_id",
                           values = deg_genes,
                           mart = ensembl)

  res_df <- as.data.frame(x)
  res_df$ensembl_gene_id <- rownames(res_df)
  res_df <- merge(res_df, gene_conversion, by = "ensembl_gene_id", all.x = TRUE)

  res_table_nam <- paste('/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot',
                         res_names[i], "with_gene.csv", sep = "_")
  write.csv(res_df, file = res_table_nam)

  assign(res_df_names[i], res_df)

  i = i + 1
}
subset(res_df_ADCY3_N1, hgnc_symbol == "ADCY3")
subset(res_df_ADCY3_1, hgnc_symbol == "ADCY3")
subset(res_df_IFNAR2_1, hgnc_symbol == "IFNAR2")
subset(res_df_IFNAR2_2, hgnc_symbol == "IFNAR2")
subset(res_df_IFNAR2_N1, hgnc_symbol == "IFNAR2")
subset(res_df_IFNAR2_N2, hgnc_symbol == "IFNAR2")

############################## vs NC1-1 only #############################
txi_ADCY3_N1 <- tximport(c(s4.1,s4.2,s3.1), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_ADCY3_1 <- tximport(c(s5.1,s5.2,s3.1), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_IFNAR2_1 <- tximport(c(s8.1,s8.2,s3.1), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_IFNAR2_2 <- tximport(c(s6.1,s6.2,s3.1), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_IFNAR2_N1 <- tximport(c(s7.1,s7.2,s3.1), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_IFNAR2_N2 <- tximport(c(s1.1,s1.2,s3.1), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

sampleTable_ADCY3_N1 <- data.frame(
  samples = c("ADCY3_N1_1", "ADCY3_N1_2", "NC1_1"),
  condition = c("ADCY3_N1", "ADCY3_N1", "NC1" )
)
rownames(sampleTable_ADCY3_N1) <- sampleTable_ADCY3_N1$samples

sampleTable_ADCY3_1 <- data.frame(
  samples = c("ADCY3_1_1", "ADCY3_1_2", "NC1_1"),
  condition = c("ADCY3_1", "ADCY3_1", "NC1" )
)
rownames(sampleTable_ADCY3_1) <- sampleTable_ADCY3_1$samples

sampleTable_IFNAR2_1 <- data.frame(
  samples = c("IFNAR2_1_1", "IFNAR2_1_2", "NC1_1"),
  condition = c("IFNAR2_1", "IFNAR2_1", "NC1" )
)
rownames(sampleTable_IFNAR2_1) <- sampleTable_IFNAR2_1$samples

sampleTable_IFNAR2_2 <- data.frame(
  samples = c("IFNAR2_2_1", "IFNAR2_2_2", "NC1_1"),
  condition = c("IFNAR2_2", "IFNAR2_2", "NC1" )
)
rownames(sampleTable_IFNAR2_2) <- sampleTable_IFNAR2_2$samples

sampleTable_IFNAR2_N1 <- data.frame(
  samples = c("IFNAR2_N1_1", "IFNAR2_N1_2", "NC1_1"),
  condition = c("IFNAR2_N1", "IFNAR2_N1", "NC1" )
)
rownames(sampleTable_IFNAR2_N1) <- sampleTable_IFNAR2_N1$samples

sampleTable_IFNAR2_N2 <- data.frame(
  samples = c("IFNAR2_N2_1", "IFNAR2_N2_2", "NC1_1"),
  condition = c("IFNAR2_N2", "IFNAR2_N2", "NC1" )
)
rownames(sampleTable_IFNAR2_N2) <- sampleTable_IFNAR2_N2$samples

txi_names <- list(txi_ADCY3_N1,txi_ADCY3_1, txi_IFNAR2_1, txi_IFNAR2_2, txi_IFNAR2_N1, txi_IFNAR2_N2)
res_names <- c("res_ADCY3_N1", "res_ADCY3_1", "res_IFNAR2_1", "res_IFNAR2_2", "res_IFNAR2_N1", "res_IFNAR2_N2")
dds_names <- c("dds_ADCY3_N1", "dds_ADCY3_1", "dds_IFNAR2_1", "dds_IFNAR2_2", "dds_IFNAR2_N1", "dds_IFNAR2_N2")
res_df_names <- c("res_df_ADCY3_N1", "res_df_ADCY3_1", "res_df_IFNAR2_1", "res_df_IFNAR2_2", "res_df_IFNAR2_N1", "res_df_IFNAR2_N2")
sample_names <- list(
  sampleTable_ADCY3_N1,
  sampleTable_ADCY3_1,
  sampleTable_IFNAR2_1,
  sampleTable_IFNAR2_2,
  sampleTable_IFNAR2_N1,
  sampleTable_IFNAR2_N2
)

## DEG analysis
i = 1
for (x in txi_names){

  dds <- DESeqDataSetFromTximport(x, colData = sample_names[i], design = ~ condition)
  dds$condition <- relevel(dds$condition, ref = "NC1")
  dds <- DESeq(dds)
  res <- results(dds)

  res$padj <- p.adjust(res$pvalue, method="BH")

  assign(res_names[i], res)
  assign(dds_names[i], dds)

  i = i +1

}

res_vars <- list(
  res_ADCY3_N1,
  res_ADCY3_1,
  res_IFNAR2_1,
  res_IFNAR2_2,
  res_IFNAR2_N1,
  res_IFNAR2_N2
)

dds_vars <- list(
  dds_ADCY3_N1,
  dds_ADCY3_1,
  dds_IFNAR2_1,
  dds_IFNAR2_2,
  dds_IFNAR2_N1,
  dds_IFNAR2_N2
)

i = 1
for (x in res_vars){
  plotMA(x, main=res_names[i])
  i = i +1
}

i = 1
for (x in dds_vars){
  vsd <- vst(x)
  plotPCA(vsd, intgroup = "condition")
  i = i +1
}

## output tables with gene name
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

i = 1
for (x in res_vars) {
  deg_genes <- rownames(x)

  gene_conversion <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                           filters = "ensembl_gene_id",
                           values = deg_genes,
                           mart = ensembl)

  res_df <- as.data.frame(x)
  res_df$ensembl_gene_id <- rownames(res_df)
  res_df <- merge(res_df, gene_conversion, by = "ensembl_gene_id", all.x = TRUE)

  assign(res_df_names[i], res_df)

  i = i + 1
}
subset(res_df_ADCY3_N1, hgnc_symbol == "ADCY3")
subset(res_df_ADCY3_1, hgnc_symbol == "ADCY3")
subset(res_df_IFNAR2_1, hgnc_symbol == "IFNAR2")
subset(res_df_IFNAR2_2, hgnc_symbol == "IFNAR2")
subset(res_df_IFNAR2_N1, hgnc_symbol == "IFNAR2")
subset(res_df_IFNAR2_N2, hgnc_symbol == "IFNAR2")
############################## vs NC1-2 only #############################
txi_ADCY3_N1 <- tximport(c(s4.1,s4.2,s3.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_ADCY3_1 <- tximport(c(s5.1,s5.2,s3.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_IFNAR2_1 <- tximport(c(s8.1,s8.2,s3.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_IFNAR2_2 <- tximport(c(s6.1,s6.2,s3.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_IFNAR2_N1 <- tximport(c(s7.1,s7.2,s3.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_IFNAR2_N2 <- tximport(c(s1.1,s1.2,s3.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

sampleTable_ADCY3_N1 <- data.frame(
  samples = c("ADCY3_N1_1", "ADCY3_N1_2",  "NC1_2"),
  condition = c("ADCY3_N1", "ADCY3_N1", "NC1" )
)
rownames(sampleTable_ADCY3_N1) <- sampleTable_ADCY3_N1$samples

sampleTable_ADCY3_1 <- data.frame(
  samples = c("ADCY3_1_1", "ADCY3_1_2",  "NC1_2"),
  condition = c("ADCY3_1", "ADCY3_1", "NC1" )
)
rownames(sampleTable_ADCY3_1) <- sampleTable_ADCY3_1$samples

sampleTable_IFNAR2_1 <- data.frame(
  samples = c("IFNAR2_1_1", "IFNAR2_1_2",  "NC1_2"),
  condition = c("IFNAR2_1", "IFNAR2_1", "NC1" )
)
rownames(sampleTable_IFNAR2_1) <- sampleTable_IFNAR2_1$samples

sampleTable_IFNAR2_2 <- data.frame(
  samples = c("IFNAR2_2_1", "IFNAR2_2_2",  "NC1_2"),
  condition = c("IFNAR2_2", "IFNAR2_2", "NC1" )
)
rownames(sampleTable_IFNAR2_2) <- sampleTable_IFNAR2_2$samples

sampleTable_IFNAR2_N1 <- data.frame(
  samples = c("IFNAR2_N1_1", "IFNAR2_N1_2",  "NC1_2"),
  condition = c("IFNAR2_N1", "IFNAR2_N1", "NC1" )
)
rownames(sampleTable_IFNAR2_N1) <- sampleTable_IFNAR2_N1$samples

sampleTable_IFNAR2_N2 <- data.frame(
  samples = c("IFNAR2_N2_1", "IFNAR2_N2_2",  "NC1_2"),
  condition = c("IFNAR2_N2", "IFNAR2_N2", "NC1" )
)
rownames(sampleTable_IFNAR2_N2) <- sampleTable_IFNAR2_N2$samples

txi_names <- list(txi_ADCY3_N1,txi_ADCY3_1, txi_IFNAR2_1, txi_IFNAR2_2, txi_IFNAR2_N1, txi_IFNAR2_N2)
res_names <- c("res_ADCY3_N1", "res_ADCY3_1", "res_IFNAR2_1", "res_IFNAR2_2", "res_IFNAR2_N1", "res_IFNAR2_N2")
dds_names <- c("dds_ADCY3_N1", "dds_ADCY3_1", "dds_IFNAR2_1", "dds_IFNAR2_2", "dds_IFNAR2_N1", "dds_IFNAR2_N2")
res_df_names <- c("res_df_ADCY3_N1", "res_df_ADCY3_1", "res_df_IFNAR2_1", "res_df_IFNAR2_2", "res_df_IFNAR2_N1", "res_df_IFNAR2_N2")
sample_names <- list(
  sampleTable_ADCY3_N1,
  sampleTable_ADCY3_1,
  sampleTable_IFNAR2_1,
  sampleTable_IFNAR2_2,
  sampleTable_IFNAR2_N1,
  sampleTable_IFNAR2_N2
)

## DEG analysis
i = 1
for (x in txi_names){

  dds <- DESeqDataSetFromTximport(x, colData = sample_names[i], design = ~ condition)
  dds$condition <- relevel(dds$condition, ref = "NC1")
  dds <- DESeq(dds)
  res <- results(dds)

  res$padj <- p.adjust(res$pvalue, method="BH")

  assign(res_names[i], res)
  assign(dds_names[i], dds)

  i = i +1

}

res_vars <- list(
  res_ADCY3_N1,
  res_ADCY3_1,
  res_IFNAR2_1,
  res_IFNAR2_2,
  res_IFNAR2_N1,
  res_IFNAR2_N2
)

dds_vars <- list(
  dds_ADCY3_N1,
  dds_ADCY3_1,
  dds_IFNAR2_1,
  dds_IFNAR2_2,
  dds_IFNAR2_N1,
  dds_IFNAR2_N2
)

i = 1
for (x in res_vars){
  plotMA(x, main=res_names[i])
  i = i +1
}

i = 1
for (x in dds_vars){
  vsd <- vst(x)
  plotPCA(vsd, intgroup = "condition")
  i = i +1
}

## output tables with gene name
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

i = 1
for (x in res_vars) {
  deg_genes <- rownames(x)

  gene_conversion <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                           filters = "ensembl_gene_id",
                           values = deg_genes,
                           mart = ensembl)

  res_df <- as.data.frame(x)
  res_df$ensembl_gene_id <- rownames(res_df)
  res_df <- merge(res_df, gene_conversion, by = "ensembl_gene_id", all.x = TRUE)

  assign(res_df_names[i], res_df)

  i = i + 1
}
subset(res_df_ADCY3_N1, hgnc_symbol == "ADCY3")
subset(res_df_ADCY3_1, hgnc_symbol == "ADCY3")
subset(res_df_IFNAR2_1, hgnc_symbol == "IFNAR2")
subset(res_df_IFNAR2_2, hgnc_symbol == "IFNAR2")
subset(res_df_IFNAR2_N1, hgnc_symbol == "IFNAR2")
subset(res_df_IFNAR2_N2, hgnc_symbol == "IFNAR2")

############################ vs others ####################
txi_ADCY3_N1 <- tximport(c(s1.1, s1.2, s3.1, s3.2, s4.1, s4.2, s6.1, s6.2, s7.1, s7.2, s8.1, s8.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_ADCY3_1 <- tximport(c(s1.1, s1.2, s3.1, s3.2, s5.1, s5.2, s6.1, s6.2, s7.1, s7.2, s8.1, s8.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_IFNAR2_1 <- tximport(c(s3.1, s3.2, s4.1, s4.2, s5.1, s5.2, s8.1, s8.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_IFNAR2_2 <- tximport(c(s3.1, s3.2, s4.1, s4.2, s5.1, s5.2, s6.1, s6.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_IFNAR2_N1 <- tximport(c(s3.1, s3.2, s4.1, s4.2, s5.1, s5.2, s7.1, s7.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_IFNAR2_N2 <- tximport(c(s1.1, s1.2, s3.1, s3.2, s4.1, s4.2, s5.1, s5.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

sampleTable_ADCY3_N1 <- data.frame(
  samples = c("1-1", "1-2", "3-1", "3-2", "4-1", "4-2", "6-1","6-2","7-1","7-2","8-1","8-2"),
  condition = c("Other", "Other", "Other","Other","ADCY3_N1","ADCY3_N1",
                "Other", "Other","Other","Other","Other","Other")
)
rownames(sampleTable_ADCY3_N1) <- sampleTable_ADCY3_N1$samples

sampleTable_ADCY3_1 <- data.frame(
  samples = c("1-1", "1-2", "3-1", "3-2", "5-1","5-2","6-1","6-2","7-1","7-2","8-1","8-2"),
  condition = c("Other", "Other", "Other","Other",
                "ADCY3_1", "ADCY3_1", "Other", "Other","Other","Other","Other","Other")
)
rownames(sampleTable_ADCY3_1) <- sampleTable_ADCY3_1$samples

sampleTable_IFNAR2_1 <- data.frame(
  samples = c("3-1", "3-2", "4-1", "4-2", "5-1","5-2","8-1","8-2"),
  condition = c("Other","Other","Other","Other",
                "Other", "Other", "IFNAR2_1","IFNAR2_1")
)
rownames(sampleTable_IFNAR2_1) <- sampleTable_IFNAR2_1$samples

sampleTable_IFNAR2_2 <- data.frame(
  samples = c("3-1", "3-2", "4-1", "4-2", "5-1","5-2","6-1","6-2"),
  condition = c("Other","Other","Other","Other",
                "Other", "Other", "IFNAR2_2", "IFNAR2_2")
)
rownames(sampleTable_IFNAR2_2) <- sampleTable_IFNAR2_2$samples

sampleTable_IFNAR2_N1 <- data.frame(
  samples = c("3-1", "3-2", "4-1", "4-2", "5-1","5-2","7-1","7-2"),
  condition = c("Other","Other","Other","Other",
                "Other", "Other", "IFNAR2_N1","IFNAR2_N1")
)
rownames(sampleTable_IFNAR2_N1) <- sampleTable_IFNAR2_N1$samples

sampleTable_IFNAR2_N2 <- data.frame(
  samples = c("1-1", "1-2","3-1", "3-2", "4-1", "4-2", "5-1","5-2"),
  condition = c("IFNAR2_N2", "IFNAR2_N2", "Other","Other","Other","Other",
                "Other", "Other")
)
rownames(sampleTable_IFNAR2_N2) <- sampleTable_IFNAR2_N2$samples

txi_names <- list(txi_ADCY3_N1,txi_ADCY3_1, txi_IFNAR2_1, txi_IFNAR2_2, txi_IFNAR2_N1, txi_IFNAR2_N2)
res_names <- c("res_ADCY3_N1", "res_ADCY3_1", "res_IFNAR2_1", "res_IFNAR2_2", "res_IFNAR2_N1", "res_IFNAR2_N2")
dds_names <- c("dds_ADCY3_N1", "dds_ADCY3_1", "dds_IFNAR2_1", "dds_IFNAR2_2", "dds_IFNAR2_N1", "dds_IFNAR2_N2")
res_df_names <- c("res_df_ADCY3_N1", "res_df_ADCY3_1", "res_df_IFNAR2_1", "res_df_IFNAR2_2", "res_df_IFNAR2_N1", "res_df_IFNAR2_N2")
sample_names <- list(
  sampleTable_ADCY3_N1,
  sampleTable_ADCY3_1,
  sampleTable_IFNAR2_1,
  sampleTable_IFNAR2_2,
  sampleTable_IFNAR2_N1,
  sampleTable_IFNAR2_N2
)

## DEG analysis
i = 1
for (x in txi_names){

  dds <- DESeqDataSetFromTximport(x, colData = sample_names[i], design = ~ condition)
  dds$condition <- relevel(dds$condition, ref = "Other")

  dds <- DESeq(dds, test = "LRT", reduced = ~1)
  res <- results(dds)

  res$padj <- p.adjust(res$pvalue, method="BH")

  assign(res_names[i], res)
  assign(dds_names[i], dds)

  i = i +1

}

## DEG analysis with limma-voom
i <- 1
for (x in txi_names) {

  txi <- x
  sample_info <- sample_names[[i]]

  y <- DGEList(counts = txi$counts)
  y <- calcNormFactors(y)

  sample_info$condition <- relevel(factor(sample_info$condition), ref = "Other")
  design <- model.matrix(~ condition, data = sample_info)

  v <- voom(y, design, plot = FALSE)

  fit <- lmFit(v, design)
  fit <- eBayes(fit)

  res <- topTable(fit, coef = 2, number = Inf, sort.by = "none")

  assign(res_names[i], res)
  assign(dds_names[i], v)

  i <- i + 1
}

## limma-voom optimized
library(edgeR)
library(limma)
library(statmod)

i <- 1
for (x in txi_names) {

  txi <- x
  sample_info <- sample_names[[i]]

  keep <- rowSums(txi$counts >= 10) >= 2
  counts_filtered <- txi$counts[keep, ]

  y <- DGEList(counts = counts_filtered)
  y <- calcNormFactors(y)

  sample_info$condition <- relevel(factor(sample_info$condition), ref = "IFNAR2_N2")
  design <- model.matrix(~ condition, data = sample_info)

  v <- voomWithQualityWeights(y, design, plot = FALSE)

  fit <- lmFit(v, design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)

  res <- topTable(fit, coef = 2, number = Inf, sort.by = "none")

  assign(res_names[i], res)
  assign(dds_names[i], v)

  i <- i + 1
}

## DESeq2 lfcShrink
library(DESeq2)
library(apeglm)

i <- 1
for (x in txi_names) {

  txi <- x
  sample_info <- sample_names[[i]]

  sample_info$condition <- relevel(factor(sample_info$condition), ref = "Other")

  dds <- DESeqDataSetFromTximport(txi, colData = sample_info, design = ~ condition)

  dds <- DESeq(dds)

  res <- lfcShrink(dds, coef = 2, type = "normal")

  assign(res_names[i], res)
  assign(dds_names[i], dds)

  i <- i + 1
}

## apeglm
i <- 1
for (x in txi_names) {

  txi <- x
  sample_info <- sample_names[[i]]
  sample_info$condition <- relevel(factor(sample_info$condition), ref = "Other")

  dds <- DESeqDataSetFromTximport(txi, colData = sample_info, design = ~ condition)

  keep <- rowSums(counts(dds) >= 10) >= 2
  dds <- dds[keep, ]

  dds <- DESeq(dds)

  coef_name <- resultsNames(dds)[2]

  res <- tryCatch({
    lfcShrink(dds, coef = coef_name, type = "apeglm")
  }, error = function(e) {
    message(paste("apeglm failed for", res_names[i], "- using normal shrinkage"))
    lfcShrink(dds, coef = coef_name, type = "normal")
  })

  assign(res_names[i], res)
  assign(dds_names[i], dds)

  i <- i + 1
}

## edgeR QLF
library(edgeR)

i <- 1
for (x in txi_names) {

  txi <- x
  sample_info <- sample_names[[i]]

  sample_info$condition <- relevel(factor(sample_info$condition), ref = "Other")

  y <- DGEList(counts = txi$counts)
  y <- calcNormFactors(y)

  design <- model.matrix(~ condition, data = sample_info)

  y <- estimateDisp(y, design)

  fit <- glmQLFit(y, design)

  qlf <- glmQLFTest(fit, coef = 2)

  res <- topTags(qlf, n = Inf)$table

  assign(res_names[i], res)
  assign(dds_names[i], fit)

  i <- i + 1
}

## edgeR QLF opt
library(edgeR)

summary_log <- list()

for (i in seq_along(txi_names)) {

  txi <- txi_names[[i]]
  sample_info <- sample_names[[i]]

  sample_info$condition <- relevel(factor(sample_info$condition), ref = "NC1")

  keep <- rowSums(txi$counts >= 10) >= 2
  counts_filtered <- txi$counts[keep, ]

  y <- DGEList(counts = counts_filtered)
  y <- calcNormFactors(y)

  design <- model.matrix(~ condition, data = sample_info)

  y <- estimateDisp(y, design)

  fit <- glmQLFit(y, design)

  coef_to_test <- if (ncol(design) >= 2) 2 else {
    warning(paste("Design matrix has <2 columns in iteration", i))
    next
  }

  qlf <- glmQLFTest(fit, coef = coef_to_test)
  res <- topTags(qlf, n = Inf)$table

  assign(res_names[i], res)
  assign(dds_names[i], fit)

  summary_log[[res_names[i]]] <- list(
    DEGs_padj05 = sum(res$FDR < 0.05, na.rm = TRUE),
    mean_logFC = mean(abs(res$logFC), na.rm = TRUE)
  )
}

## ShrinkBayes
library(ShrinkBayes)
library(BiocParallel)

i <- 1
for (x in txi_names) {

  txi <- x
  sample_info <- sample_names[[i]]
  counts <- as.matrix(txi$counts)

  group_factor <- relevel(factor(sample_info$condition), ref = "Other")

  dat <- list(y = counts, group = group_factor)

  sb <- ShrinkSeq(dat = dat,
                  formula = "y ~ group",
                  parameter = "group",
                  prior = "spike",
                  ngroups = 2,
                  model = "nbinom",
                  ncpus = 1)

  res <- SummarizePost(sb)

  assign(res_names[i], res)
  assign(dds_names[i], sb)

  i <- i + 1
}

res_vars <- list(
  res_ADCY3_N1,
  res_ADCY3_1,
  res_IFNAR2_1,
  res_IFNAR2_2,
  res_IFNAR2_N1,
  res_IFNAR2_N2
)

dds_vars <- list(
  dds_ADCY3_N1,
  dds_ADCY3_1,
  dds_IFNAR2_1,
  dds_IFNAR2_2,
  dds_IFNAR2_N1,
  dds_IFNAR2_N2
)

## output tables with gene name
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

i = 1
for (x in res_vars) {
  deg_genes <- rownames(x)

  gene_conversion <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                           filters = "ensembl_gene_id",
                           values = deg_genes,
                           mart = ensembl)

  res_df <- as.data.frame(x)
  res_df$ensembl_gene_id <- rownames(res_df)
  res_df <- merge(res_df, gene_conversion, by = "ensembl_gene_id", all.x = TRUE)

  assign(res_df_names[i], res_df)

  i = i + 1
}
subset(res_df_ADCY3_N1, hgnc_symbol == "ADCY3")
subset(res_df_ADCY3_1, hgnc_symbol == "ADCY3")
subset(res_df_IFNAR2_1, hgnc_symbol == "IFNAR2")
subset(res_df_IFNAR2_2, hgnc_symbol == "IFNAR2")
subset(res_df_IFNAR2_N1, hgnc_symbol == "IFNAR2")
subset(res_df_IFNAR2_N2, hgnc_symbol == "IFNAR2")

## exclude NC1-1
txi_ADCY3_N1 <- tximport(c(s1.1, s1.2, s3.2, s4.1, s4.2, s6.1, s6.2, s7.1, s7.2, s8.1, s8.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_ADCY3_1 <- tximport(c(s1.1, s1.2, s3.2, s5.1, s5.2, s6.1, s6.2, s7.1, s7.2, s8.1, s8.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_IFNAR2_1 <- tximport(c(s3.2, s4.1, s4.2, s5.1, s5.2, s8.1, s8.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_IFNAR2_2 <- tximport(c(s3.2, s4.1, s4.2, s5.1, s5.2, s6.1, s6.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_IFNAR2_N1 <- tximport(c(s3.2, s4.1, s4.2, s5.1, s5.2, s7.1, s7.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_IFNAR2_N2 <- tximport(c(s1.1, s1.2, s3.2, s4.1, s4.2, s5.1, s5.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

sampleTable_ADCY3_N1 <- data.frame(
  samples = c("1-1", "1-2", "3-2", "4-1", "4-2", "6-1","6-2","7-1","7-2","8-1","8-2"),
  condition = c("Other", "Other","Other","ADCY3_N1","ADCY3_N1",
                "Other", "Other","Other","Other","Other","Other")
)
rownames(sampleTable_ADCY3_N1) <- sampleTable_ADCY3_N1$samples

sampleTable_ADCY3_1 <- data.frame(
  samples = c("1-1", "1-2", "3-2", "5-1","5-2","6-1","6-2","7-1","7-2","8-1","8-2"),
  condition = c("Other", "Other","Other",
                "ADCY3_1", "ADCY3_1", "Other", "Other","Other","Other","Other","Other")
)
rownames(sampleTable_ADCY3_1) <- sampleTable_ADCY3_1$samples

sampleTable_IFNAR2_1 <- data.frame(
  samples = c("3-2", "4-1", "4-2", "5-1","5-2","8-1","8-2"),
  condition = c("Other","Other","Other",
                "Other", "Other", "IFNAR2_1","IFNAR2_1")
)
rownames(sampleTable_IFNAR2_1) <- sampleTable_IFNAR2_1$samples

sampleTable_IFNAR2_2 <- data.frame(
  samples = c("3-2", "4-1", "4-2", "5-1","5-2","6-1","6-2"),
  condition = c("Other","Other","Other",
                "Other", "Other", "IFNAR2_2", "IFNAR2_2")
)
rownames(sampleTable_IFNAR2_2) <- sampleTable_IFNAR2_2$samples

sampleTable_IFNAR2_N1 <- data.frame(
  samples = c("3-2", "4-1", "4-2", "5-1","5-2","7-1","7-2"),
  condition = c("Other","Other","Other",
                "Other", "Other", "IFNAR2_N1","IFNAR2_N1")
)
rownames(sampleTable_IFNAR2_N1) <- sampleTable_IFNAR2_N1$samples

sampleTable_IFNAR2_N2 <- data.frame(
  samples = c("1-1", "1-2","3-2", "4-1", "4-2", "5-1","5-2"),
  condition = c("IFNAR2_N2", "IFNAR2_N2", "Other","Other","Other",
                "Other", "Other")
)
rownames(sampleTable_IFNAR2_N2) <- sampleTable_IFNAR2_N2$samples

txi_names <- list(txi_ADCY3_N1,txi_ADCY3_1, txi_IFNAR2_1, txi_IFNAR2_2, txi_IFNAR2_N1, txi_IFNAR2_N2)
res_names <- c("res_ADCY3_N1", "res_ADCY3_1", "res_IFNAR2_1", "res_IFNAR2_2", "res_IFNAR2_N1", "res_IFNAR2_N2")
dds_names <- c("dds_ADCY3_N1", "dds_ADCY3_1", "dds_IFNAR2_1", "dds_IFNAR2_2", "dds_IFNAR2_N1", "dds_IFNAR2_N2")
res_df_names <- c("res_df_ADCY3_N1", "res_df_ADCY3_1", "res_df_IFNAR2_1", "res_df_IFNAR2_2", "res_df_IFNAR2_N1", "res_df_IFNAR2_N2")
sample_names <- list(
  sampleTable_ADCY3_N1,
  sampleTable_ADCY3_1,
  sampleTable_IFNAR2_1,
  sampleTable_IFNAR2_2,
  sampleTable_IFNAR2_N1,
  sampleTable_IFNAR2_N2
)

## exclude NC1-2
txi_ADCY3_N1 <- tximport(c(s1.1, s1.2, s3.1, s4.1, s4.2, s6.1, s6.2, s7.1, s7.2, s8.1, s8.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_ADCY3_1 <- tximport(c(s1.1, s1.2, s3.1, s5.1, s5.2, s6.1, s6.2, s7.1, s7.2, s8.1, s8.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_IFNAR2_1 <- tximport(c(s3.1, s4.1, s4.2, s5.1, s5.2, s8.1, s8.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_IFNAR2_2 <- tximport(c(s3.1, s4.1, s4.2, s5.1, s5.2, s6.1, s6.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_IFNAR2_N1 <- tximport(c(s3.1, s4.1, s4.2, s5.1, s5.2, s7.1, s7.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_IFNAR2_N2 <- tximport(c(s1.1, s1.2, s3.1, s4.1, s4.2, s5.1, s5.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_IFNAR2_com <- tximport(c(s1.1, s1.2, s8.1, s8.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

sampleTable_ADCY3_N1 <- data.frame(
  samples = c("1-1", "1-2", "3-1", "4-1", "4-2", "6-1","6-2","7-1","7-2","8-1","8-2"),
  condition = c("Other", "Other","Other","ADCY3_N1","ADCY3_N1",
                "Other", "Other","Other","Other","Other","Other")
)
rownames(sampleTable_ADCY3_N1) <- sampleTable_ADCY3_N1$samples

sampleTable_ADCY3_1 <- data.frame(
  samples = c("1-1", "1-2", "3-1", "5-1","5-2","6-1","6-2","7-1","7-2","8-1","8-2"),
  condition = c("Other", "Other","Other",
                "ADCY3_1", "ADCY3_1", "Other", "Other","Other","Other","Other","Other")
)
rownames(sampleTable_ADCY3_1) <- sampleTable_ADCY3_1$samples

sampleTable_IFNAR2_1 <- data.frame(
  samples = c("3-1", "4-1", "4-2", "5-1","5-2","8-1","8-2"),
  condition = c("Other","Other","Other",
                "Other", "Other", "IFNAR2_1","IFNAR2_1")
)
rownames(sampleTable_IFNAR2_1) <- sampleTable_IFNAR2_1$samples

sampleTable_IFNAR2_2 <- data.frame(
  samples = c("3-1", "4-1", "4-2", "5-1","5-2","6-1","6-2"),
  condition = c("Other","Other","Other",
                "Other", "Other", "IFNAR2_2", "IFNAR2_2")
)
rownames(sampleTable_IFNAR2_2) <- sampleTable_IFNAR2_2$samples

sampleTable_IFNAR2_N1 <- data.frame(
  samples = c("3-1", "4-1", "4-2", "5-1","5-2","7-1","7-2"),
  condition = c("Other","Other","Other",
                "Other", "Other", "IFNAR2_N1","IFNAR2_N1")
)
rownames(sampleTable_IFNAR2_N1) <- sampleTable_IFNAR2_N1$samples

sampleTable_IFNAR2_N2 <- data.frame(
  samples = c("1-1", "1-2","3-1", "4-1", "4-2", "5-1","5-2"),
  condition = c("IFNAR2_N2", "IFNAR2_N2", "Other","Other","Other",
                "Other", "Other")
)
rownames(sampleTable_IFNAR2_N2) <- sampleTable_IFNAR2_N2$samples

sampleTable_IFNAR2_com <- data.frame(
  samples = c("1-1", "1-2","8-1", "8-2"),
  condition = c("IFNAR2_N2", "IFNAR2_N2", "IFNAR2_1","IFNAR2_1")
)
rownames(sampleTable_IFNAR2_com) <- sampleTable_IFNAR2_com$samples

txi_names <- list(txi_ADCY3_N1,txi_ADCY3_1, txi_IFNAR2_1, txi_IFNAR2_2, txi_IFNAR2_N1, txi_IFNAR2_N2, txi_IFNAR2_com)
res_names <- c("res_ADCY3_N1", "res_ADCY3_1", "res_IFNAR2_1", "res_IFNAR2_2", "res_IFNAR2_N1", "res_IFNAR2_N2", "res_IFNAR2_com")
dds_names <- c("dds_ADCY3_N1", "dds_ADCY3_1", "dds_IFNAR2_1", "dds_IFNAR2_2", "dds_IFNAR2_N1", "dds_IFNAR2_N2", "dds_IFNAR2_com")
res_df_names <- c("res_df_ADCY3_N1", "res_df_ADCY3_1", "res_df_IFNAR2_1", "res_df_IFNAR2_2", "res_df_IFNAR2_N1", "res_df_IFNAR2_N2", "res_df_IFNAR2_com")
sample_names <- list(
  sampleTable_ADCY3_N1,
  sampleTable_ADCY3_1,
  sampleTable_IFNAR2_1,
  sampleTable_IFNAR2_2,
  sampleTable_IFNAR2_N1,
  sampleTable_IFNAR2_N2,
  sampleTable_IFNAR2_com
)

txi_names <- list(txi_IFNAR2_com)
res_names <- c("res_IFNAR2_com")
dds_names <- c("dds_IFNAR2_com")
res_df_names <- c("res_df_IFNAR2_com")
sample_names <- list(
  sampleTable_IFNAR2_com
)

res_vars <- list(
  res_IFNAR2_com
)

dds_vars <- list(
  dds_IFNAR2_com
)

## output tables with gene name
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

i = 1
for (x in res_vars) {
  deg_genes <- rownames(x)

  gene_conversion <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                           filters = "ensembl_gene_id",
                           values = deg_genes,
                           mart = ensembl)

  res_df <- as.data.frame(x)
  res_df$ensembl_gene_id <- rownames(res_df)
  res_df <- merge(res_df, gene_conversion, by = "ensembl_gene_id", all.x = TRUE)

  assign(res_df_names[i], res_df)

  i = i + 1
}

subset(res_df_IFNAR2_com, hgnc_symbol == "IFNAR2")

##################################### pilot batch ########################
### No UMI
s1.1 <- '/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot/salmon_quant_noUMI/output_quant_sample1-1/quant.sf'
s1.2 <- '/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot/salmon_quant_noUMI/output_quant_sample1-2/quant.sf'
s2.1 <- '/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot/salmon_quant_noUMI/output_quant_sample2-1/quant.sf'
s2.2 <- '/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot/salmon_quant_noUMI/output_quant_sample2-2/quant.sf'
s3.1 <- '/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot/salmon_quant_noUMI/output_quant_sample3-1/quant.sf'
s3.2 <- '/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot/salmon_quant_noUMI/output_quant_sample3-2/quant.sf'
s4.1 <- '/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot/salmon_quant_noUMI/output_quant_sample4-1/quant.sf'
s4.2 <- '/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot/salmon_quant_noUMI/output_quant_sample4-2/quant.sf'
s5.1 <- '/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot/salmon_quant_noUMI/output_quant_sample5-1/quant.sf'
s5.2 <- '/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot/salmon_quant_noUMI/output_quant_sample5-2/quant.sf'
s6.1 <- '/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot/salmon_quant_noUMI/output_quant_sample6-1/quant.sf'
s6.2 <- '/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot/salmon_quant_noUMI/output_quant_sample6-2/quant.sf'
s7.1 <- '/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot/salmon_quant_noUMI/output_quant_sample7-1/quant.sf'
s7.2 <- '/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot/salmon_quant_noUMI/output_quant_sample7-2/quant.sf'
s8.1 <- '/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot/salmon_quant_noUMI/output_quant_sample8-1/quant.sf'
s8.2 <- '/media/ming/Extra_SSD_4TB1/eCROPseq/20250725_bulk_cell_RNA_seq_pilot/salmon_quant_noUMI/output_quant_sample8-2/quant.sf'

txi_ADCY3_N1 <- tximport(c(s1.1, s1.2, s3.1, s3.2, s4.1, s4.2, s6.1, s6.2, s7.1, s7.2, s8.1, s8.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_ADCY3_1 <- tximport(c(s1.1, s1.2, s3.1, s3.2, s5.1, s5.2, s6.1, s6.2, s7.1, s7.2, s8.1, s8.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_IFNAR2_1 <- tximport(c(s3.1, s3.2, s4.1, s4.2, s5.1, s5.2, s8.1, s8.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_IFNAR2_2 <- tximport(c(s3.1, s3.2, s4.1, s4.2, s5.1, s5.2, s6.1, s6.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_IFNAR2_N1 <- tximport(c(s3.1, s3.2, s4.1, s4.2, s5.1, s5.2, s7.1, s7.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_IFNAR2_N2 <- tximport(c(s1.1, s1.2, s3.1, s3.2, s4.1, s4.2, s5.1, s5.2), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")

sampleTable_ADCY3_N1 <- data.frame(
  samples = c("1-1", "1-2", "3-1", "3-2", "4-1", "4-2", "6-1","6-2","7-1","7-2","8-1","8-2"),
  condition = c("Other", "Other", "Other","Other","ADCY3_N1","ADCY3_N1",
                "Other", "Other","Other","Other","Other","Other")
)
rownames(sampleTable_ADCY3_N1) <- sampleTable_ADCY3_N1$samples

sampleTable_ADCY3_1 <- data.frame(
  samples = c("1-1", "1-2", "3-1", "3-2", "5-1","5-2","6-1","6-2","7-1","7-2","8-1","8-2"),
  condition = c("Other", "Other", "Other","Other",
                "ADCY3_1", "ADCY3_1", "Other", "Other","Other","Other","Other","Other")
)
rownames(sampleTable_ADCY3_1) <- sampleTable_ADCY3_1$samples

sampleTable_IFNAR2_1 <- data.frame(
  samples = c("3-1", "3-2", "4-1", "4-2", "5-1","5-2","8-1","8-2"),
  condition = c("Other","Other","Other","Other",
                "Other", "Other", "IFNAR2_1","IFNAR2_1")
)
rownames(sampleTable_IFNAR2_1) <- sampleTable_IFNAR2_1$samples

sampleTable_IFNAR2_2 <- data.frame(
  samples = c("3-1", "3-2", "4-1", "4-2", "5-1","5-2","6-1","6-2"),
  condition = c("Other","Other","Other","Other",
                "Other", "Other", "IFNAR2_2", "IFNAR2_2")
)
rownames(sampleTable_IFNAR2_2) <- sampleTable_IFNAR2_2$samples

sampleTable_IFNAR2_N1 <- data.frame(
  samples = c("3-1", "3-2", "4-1", "4-2", "5-1","5-2","7-1","7-2"),
  condition = c("Other","Other","Other","Other",
                "Other", "Other", "IFNAR2_N1","IFNAR2_N1")
)
rownames(sampleTable_IFNAR2_N1) <- sampleTable_IFNAR2_N1$samples

sampleTable_IFNAR2_N2 <- data.frame(
  samples = c("1-1", "1-2","3-1", "3-2", "4-1", "4-2", "5-1","5-2"),
  condition = c("IFNAR2_N2", "IFNAR2_N2", "Other","Other","Other","Other",
                "Other", "Other")
)
rownames(sampleTable_IFNAR2_N2) <- sampleTable_IFNAR2_N2$samples

txi_names <- list(txi_ADCY3_N1,txi_ADCY3_1, txi_IFNAR2_1, txi_IFNAR2_2, txi_IFNAR2_N1, txi_IFNAR2_N2)
res_names <- c("res_ADCY3_N1", "res_ADCY3_1", "res_IFNAR2_1", "res_IFNAR2_2", "res_IFNAR2_N1", "res_IFNAR2_N2")
dds_names <- c("dds_ADCY3_N1", "dds_ADCY3_1", "dds_IFNAR2_1", "dds_IFNAR2_2", "dds_IFNAR2_N1", "dds_IFNAR2_N2")
res_df_names <- c("res_df_ADCY3_N1", "res_df_ADCY3_1", "res_df_IFNAR2_1", "res_df_IFNAR2_2", "res_df_IFNAR2_N1", "res_df_IFNAR2_N2")
sample_names <- list(
  sampleTable_ADCY3_N1,
  sampleTable_ADCY3_1,
  sampleTable_IFNAR2_1,
  sampleTable_IFNAR2_2,
  sampleTable_IFNAR2_N1,
  sampleTable_IFNAR2_N2
)

i <- 1
for (x in txi_names) {

  txi <- x
  sample_info <- sample_names[[i]]

  keep <- rowSums(txi$counts >= 10) >= 2
  counts_filtered <- txi$counts[keep, ]

  y <- DGEList(counts = counts_filtered)
  y <- calcNormFactors(y)

  sample_info$condition <- relevel(factor(sample_info$condition), ref = "Other")
  design <- model.matrix(~ condition, data = sample_info)

  v <- voomWithQualityWeights(y, design, plot = FALSE)

  fit <- lmFit(v, design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)

  res <- topTable(fit, coef = 2, number = Inf, sort.by = "none")

  assign(res_names[i], res)
  assign(dds_names[i], v)

  i <- i + 1
}

res_vars <- list(
  res_ADCY3_N1,
  res_ADCY3_1,
  res_IFNAR2_1,
  res_IFNAR2_2,
  res_IFNAR2_N1,
  res_IFNAR2_N2
)

dds_vars <- list(
  dds_ADCY3_N1,
  dds_ADCY3_1,
  dds_IFNAR2_1,
  dds_IFNAR2_2,
  dds_IFNAR2_N1,
  dds_IFNAR2_N2
)

## output tables with gene name
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

i = 1
for (x in res_vars) {
  deg_genes <- rownames(x)

  gene_conversion <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                           filters = "ensembl_gene_id",
                           values = deg_genes,
                           mart = ensembl)

  res_df <- as.data.frame(x)
  res_df$ensembl_gene_id <- rownames(res_df)
  res_df <- merge(res_df, gene_conversion, by = "ensembl_gene_id", all.x = TRUE)

  assign(res_df_names[i], res_df)

  i = i + 1
}
subset(res_df_ADCY3_N1, hgnc_symbol == "ADCY3")
subset(res_df_ADCY3_1, hgnc_symbol == "ADCY3")
subset(res_df_IFNAR2_1, hgnc_symbol == "IFNAR2")
subset(res_df_IFNAR2_2, hgnc_symbol == "IFNAR2")
subset(res_df_IFNAR2_N1, hgnc_symbol == "IFNAR2")
subset(res_df_IFNAR2_N2, hgnc_symbol == "IFNAR2")

#################################### GENE EXPRESSIONS #########################
selected_genes <- c("ADCY3", "AKAP11", "ATG16L1", "CEBPB", "CARD9", "CARS1", "DAP", "ERAP2", "HLA-C", "IFNAR2", "CEBPB_V1", "MTMR3", "PTPN22", "RIPK2", "SNAKPC4")

filtered_df <- res_df_ADCY3_N1[res_df_ADCY3_N1$hgnc_symbol %in% selected_genes, ]
print(filtered_df)

selected_symbols <- filtered_df$ensembl_gene_id

txi_filtered <- txi_all$counts[rownames(txi_all$counts) %in% selected_symbols, ]

print(txi_filtered)

##################################### 2nd batch ###################################
samples <- 1:36
paths <- sprintf("/media/ming/Extra_SSD_4TB1/eCROPseq/20250903_bulk_cell_RNA_seq_2nd/salmon_quant/output_quant_sample%d_S%d/quant.sf", samples, samples)

txi_DAP_N4 <- tximport(c(paths[1:3],paths[16:36]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_DAP_N4_1 <- tximport(c(paths[1],paths[16:36]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_DAP_N4_2 <- tximport(c(paths[2],paths[16:36]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_DAP_N4_3 <- tximport(c(paths[3],paths[16:36]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")

txi_DAP_V4 <- tximport(c(paths[4:6],paths[16:36]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_DAP_V4_1 <- tximport(c(paths[4],paths[16:36]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_DAP_V4_2 <- tximport(c(paths[5],paths[16:36]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_DAP_V4_3 <- tximport(c(paths[6],paths[16:36]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")

txi_DAP_N2 <- tximport(c(paths[7:9],paths[16:36]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_DAP_N2_1 <- tximport(c(paths[7],paths[16:36]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_DAP_N2_2 <- tximport(c(paths[8],paths[16:36]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_DAP_N2_3 <- tximport(c(paths[9],paths[16:36]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")

txi_DAP_V2 <- tximport(c(paths[10:12],paths[16:36]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_DAP_V2_1 <- tximport(c(paths[10],paths[16:36]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_DAP_V2_2 <- tximport(c(paths[11],paths[16:36]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_DAP_V2_3 <- tximport(c(paths[12],paths[16:36]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")

txi_DAP_V1 <- tximport(c(paths[13:15],paths[16:36]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_DAP_V1_1 <- tximport(c(paths[13],paths[16:36]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_DAP_V1_2 <- tximport(c(paths[14],paths[16:36]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_DAP_V1_3 <- tximport(c(paths[15],paths[16:36]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")

txi_RIPK2_V1 <- tximport(c(paths[16:18],paths[1:15],paths[22:36]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_RIPK2_V1_1 <- tximport(c(paths[16],paths[1:15],paths[22:36]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_RIPK2_V1_2 <- tximport(c(paths[17],paths[1:15],paths[22:36]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_RIPK2_V1_3 <- tximport(c(paths[18],paths[1:15],paths[22:36]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")

txi_RIPK2_N2 <- tximport(c(paths[19:21],paths[1:15],paths[22:36]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_RIPK2_N2_1 <- tximport(c(paths[19],paths[1:15],paths[22:36]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_RIPK2_N2_2 <- tximport(c(paths[20],paths[1:15],paths[22:36]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_RIPK2_N2_3 <- tximport(c(paths[21],paths[1:15],paths[22:36]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")

txi_AKAP11_N1 <- tximport(c(paths[31:33],paths[1:30]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_AKAP11_N1_1 <- tximport(c(paths[31],paths[1:30]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_AKAP11_N1_2 <- tximport(c(paths[32],paths[1:30]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_AKAP11_N1_3 <- tximport(c(paths[33],paths[1:30]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")

txi_AKAP11_V1 <- tximport(c(paths[34:36],paths[1:30]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_AKAP11_V1_1 <- tximport(c(paths[34],paths[1:30]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_AKAP11_V1_2 <- tximport(c(paths[35],paths[1:30]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_AKAP11_V1_3 <- tximport(c(paths[36],paths[1:30]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")

txi_NC1 <- tximport(c(paths[22:24],paths[1:21],paths[25:36]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_HLAC_N1 <- tximport(c(paths[25:27],paths[1:24],paths[31:36]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_HLAC_V3 <- tximport(c(paths[28:30],paths[1:24],paths[31:36]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
samples_order <- as.character(1:36)

samples_order <- as.character(1:36)
# ==================== DAP ====================

sampleTable_DAP_N4 <- data.frame(
  samples = c(samples_order[c(1:3,16:36)]),
  condition = c(rep("DAP_N4", 3), rep("Other", 21))
)
rownames(sampleTable_DAP_N4) <- sampleTable_DAP_N4$samples

sampleTable_DAP_N4_1 <- data.frame(
  samples = c(samples_order[c(1,16:36)]),
  condition = c("DAP_N4", rep("Other", 21))
)
rownames(sampleTable_DAP_N4_1) <- sampleTable_DAP_N4_1$samples

sampleTable_DAP_N4_2 <- data.frame(
  samples = c(samples_order[c(2,16:36)]),
  condition = c("DAP_N4", rep("Other", 21))
)
rownames(sampleTable_DAP_N4_2) <- sampleTable_DAP_N4_2$samples

sampleTable_DAP_N4_3 <- data.frame(
  samples = c(samples_order[c(3,16:36)]),
  condition = c("DAP_N4", rep("Other", 21))
)
rownames(sampleTable_DAP_N4_3) <- sampleTable_DAP_N4_3$samples

sampleTable_DAP_V4 <- data.frame(
  samples = c(samples_order[c(4:6,16:36)]),
  condition = c(rep("DAP_V4", 3), rep("Other", 21))
)
rownames(sampleTable_DAP_V4) <- sampleTable_DAP_V4$samples

sampleTable_DAP_V4_1 <- data.frame(
  samples = c(samples_order[c(4,16:36)]),
  condition = c("DAP_V4", rep("Other", 21))
)
rownames(sampleTable_DAP_V4_1) <- sampleTable_DAP_V4_1$samples

sampleTable_DAP_V4_2 <- data.frame(
  samples = c(samples_order[c(5,16:36)]),
  condition = c("DAP_V4", rep("Other", 21))
)
rownames(sampleTable_DAP_V4_2) <- sampleTable_DAP_V4_2$samples

sampleTable_DAP_V4_3 <- data.frame(
  samples = c(samples_order[c(6,16:36)]),
  condition = c("DAP_V4", rep("Other", 21))
)
rownames(sampleTable_DAP_V4_3) <- sampleTable_DAP_V4_3$samples

sampleTable_DAP_N2 <- data.frame(
  samples = c(samples_order[c(7:9,16:36)]),
  condition = c(rep("DAP_N2", 3), rep("Other", 21))
)
rownames(sampleTable_DAP_N2) <- sampleTable_DAP_N2$samples

sampleTable_DAP_N2_1 <- data.frame(
  samples = c(samples_order[c(7,16:36)]),
  condition = c("DAP_N2", rep("Other", 21))
)
rownames(sampleTable_DAP_N2_1) <- sampleTable_DAP_N2_1$samples

sampleTable_DAP_N2_2 <- data.frame(
  samples = c(samples_order[c(8,16:36)]),
  condition = c("DAP_N2", rep("Other", 21))
)
rownames(sampleTable_DAP_N2_2) <- sampleTable_DAP_N2_2$samples

sampleTable_DAP_N2_3 <- data.frame(
  samples = c(samples_order[c(9,16:36)]),
  condition = c("DAP_N2", rep("Other", 21))
)
rownames(sampleTable_DAP_N2_3) <- sampleTable_DAP_N2_3$samples

sampleTable_DAP_V2 <- data.frame(
  samples = c(samples_order[c(10:12,16:36)]),
  condition = c(rep("DAP_V2", 3), rep("Other", 21))
)
rownames(sampleTable_DAP_V2) <- sampleTable_DAP_V2$samples

sampleTable_DAP_V2_1 <- data.frame(
  samples = c(samples_order[c(10,16:36)]),
  condition = c("DAP_V2", rep("Other", 21))
)
rownames(sampleTable_DAP_V2_1) <- sampleTable_DAP_V2_1$samples

sampleTable_DAP_V2_2 <- data.frame(
  samples = c(samples_order[c(11,16:36)]),
  condition = c("DAP_V2", rep("Other", 21))
)
rownames(sampleTable_DAP_V2_2) <- sampleTable_DAP_V2_2$samples

sampleTable_DAP_V2_3 <- data.frame(
  samples = c(samples_order[c(12,16:36)]),
  condition = c("DAP_V2", rep("Other", 21))
)
rownames(sampleTable_DAP_V2_3) <- sampleTable_DAP_V2_3$samples

sampleTable_DAP_V1 <- data.frame(
  samples = c(samples_order[c(13:15,16:36)]),
  condition = c(rep("DAP_V1", 3), rep("Other", 21))
)
rownames(sampleTable_DAP_V1) <- sampleTable_DAP_V1$samples

sampleTable_DAP_V1_1 <- data.frame(
  samples = c(samples_order[c(13,16:36)]),
  condition = c("DAP_V1", rep("Other", 21))
)
rownames(sampleTable_DAP_V1_1) <- sampleTable_DAP_V1_1$samples

sampleTable_DAP_V1_2 <- data.frame(
  samples = c(samples_order[c(14,16:36)]),
  condition = c("DAP_V1", rep("Other", 21))
)
rownames(sampleTable_DAP_V1_2) <- sampleTable_DAP_V1_2$samples

sampleTable_DAP_V1_3 <- data.frame(
  samples = c(samples_order[c(15,16:36)]),
  condition = c("DAP_V1", rep("Other", 21))
)
rownames(sampleTable_DAP_V1_3) <- sampleTable_DAP_V1_3$samples

# ==================== RIPK2 ====================

sampleTable_RIPK2_V1 <- data.frame(
  samples = c(samples_order[c(16:18,1:15,22:36)]),
  condition = c(rep("RIPK2_V1", 3), rep("Other", 30))
)
rownames(sampleTable_RIPK2_V1) <- sampleTable_RIPK2_V1$samples

sampleTable_RIPK2_V1_1 <- data.frame(
  samples = c(samples_order[c(16,1:15,22:36)]),
  condition = c("RIPK2_V1", rep("Other", 30))
)
rownames(sampleTable_RIPK2_V1_1) <- sampleTable_RIPK2_V1_1$samples

sampleTable_RIPK2_V1_2 <- data.frame(
  samples = c(samples_order[c(17,1:15,22:36)]),
  condition = c("RIPK2_V1", rep("Other", 30))
)
rownames(sampleTable_RIPK2_V1_2) <- sampleTable_RIPK2_V1_2$samples

sampleTable_RIPK2_V1_3 <- data.frame(
  samples = c(samples_order[c(18,1:15,22:36)]),
  condition = c("RIPK2_V1", rep("Other", 30))
)
rownames(sampleTable_RIPK2_V1_3) <- sampleTable_RIPK2_V1_3$samples

sampleTable_RIPK2_N2 <- data.frame(
  samples = c(samples_order[c(19:21,1:15,22:36)]),
  condition = c(rep("RIPK2_N2", 3), rep("Other", 30))
)
rownames(sampleTable_RIPK2_N2) <- sampleTable_RIPK2_N2$samples

sampleTable_RIPK2_N2_1 <- data.frame(
  samples = c(samples_order[c(19,1:15,22:36)]),
  condition = c("RIPK2_N2", rep("Other", 30))
)
rownames(sampleTable_RIPK2_N2_1) <- sampleTable_RIPK2_N2_1$samples

sampleTable_RIPK2_N2_2 <- data.frame(
  samples = c(samples_order[c(20,1:15,22:36)]),
  condition = c("RIPK2_N2", rep("Other", 30))
)
rownames(sampleTable_RIPK2_N2_2) <- sampleTable_RIPK2_N2_2$samples

sampleTable_RIPK2_N2_3 <- data.frame(
  samples = c(samples_order[c(21,1:15,22:36)]),
  condition = c("RIPK2_N2", rep("Other", 30))
)
rownames(sampleTable_RIPK2_N2_3) <- sampleTable_RIPK2_N2_3$samples

# ==================== AKAP11 ====================

sampleTable_AKAP11_N1 <- data.frame(
  samples = c(samples_order[c(31:33,1:30)]),
  condition = c(rep("AKAP11_N1", 3), rep("Other", 30))
)
rownames(sampleTable_AKAP11_N1) <- sampleTable_AKAP11_N1$samples

sampleTable_AKAP11_N1_1 <- data.frame(
  samples = c(samples_order[c(31,1:30)]),
  condition = c("AKAP11_N1", rep("Other", 30))
)
rownames(sampleTable_AKAP11_N1_1) <- sampleTable_AKAP11_N1_1$samples

sampleTable_AKAP11_N1_2 <- data.frame(
  samples = c(samples_order[c(32,1:30)]),
  condition = c("AKAP11_N1", rep("Other", 30))
)
rownames(sampleTable_AKAP11_N1_2) <- sampleTable_AKAP11_N1_2$samples

sampleTable_AKAP11_N1_3 <- data.frame(
  samples = c(samples_order[c(33,1:30)]),
  condition = c("AKAP11_N1", rep("Other", 30))
)
rownames(sampleTable_AKAP11_N1_3) <- sampleTable_AKAP11_N1_3$samples

sampleTable_AKAP11_V1 <- data.frame(
  samples = c(samples_order[c(34:36,1:30)]),
  condition = c(rep("AKAP11_V1", 3), rep("Other", 30))
)
rownames(sampleTable_AKAP11_V1) <- sampleTable_AKAP11_V1$samples

sampleTable_AKAP11_V1_1 <- data.frame(
  samples = c(samples_order[c(34,1:30)]),
  condition = c("AKAP11_V1", rep("Other", 30))
)
rownames(sampleTable_AKAP11_V1_1) <- sampleTable_AKAP11_V1_1$samples

sampleTable_AKAP11_V1_2 <- data.frame(
  samples = c(samples_order[c(35,1:30)]),
  condition = c("AKAP11_V1", rep("Other", 30))
)
rownames(sampleTable_AKAP11_V1_2) <- sampleTable_AKAP11_V1_2$samples

sampleTable_AKAP11_V1_3 <- data.frame(
  samples = c(samples_order[c(36,1:30)]),
  condition = c("AKAP11_V1", rep("Other", 30))
)
rownames(sampleTable_AKAP11_V1_3) <- sampleTable_AKAP11_V1_3$samples

# ==================== HLAC & NC1 (No _1/_2/_3 variants) ====================

sampleTable_HLAC_N1 <- data.frame(
  samples = c(samples_order[c(25:27,1:24,31:36)]),
  condition = c(rep("HLAC_N1", 3), rep("Other", 30))
)
rownames(sampleTable_HLAC_N1) <- sampleTable_HLAC_N1$samples

sampleTable_HLAC_V3 <- data.frame(
  samples = c(samples_order[c(28:30,1:24,31:36)]),
  condition = c(rep("HLAC_V3", 3), rep("Other", 30))
)
rownames(sampleTable_HLAC_V3) <- sampleTable_HLAC_V3$samples

sampleTable_NC1 <- data.frame(
  samples = c(samples_order[c(22:24,1:21,25:36)]),
  condition = c(rep("NC1", 3), rep("Other", 33))
)
rownames(sampleTable_NC1) <- sampleTable_NC1$samples

txi_names <- list(
  txi_DAP_N4,
  txi_DAP_N4_1,
  txi_DAP_N4_2,
  txi_DAP_N4_3,
  txi_DAP_V4,
  txi_DAP_V4_1,
  txi_DAP_V4_2,
  txi_DAP_V4_3,
  txi_DAP_N2,
  txi_DAP_N2_1,
  txi_DAP_N2_2,
  txi_DAP_N2_3,
  txi_DAP_V2,
  txi_DAP_V2_1,
  txi_DAP_V2_2,
  txi_DAP_V2_3,
  txi_DAP_V1,
  txi_DAP_V1_1,
  txi_DAP_V1_2,
  txi_DAP_V1_3,
  txi_RIPK2_V1,
  txi_RIPK2_V1_1,
  txi_RIPK2_V1_2,
  txi_RIPK2_V1_3,
  txi_RIPK2_N2,
  txi_RIPK2_N2_1,
  txi_RIPK2_N2_2,
  txi_RIPK2_N2_3,
  txi_NC1,
  txi_HLAC_N1,
  txi_HLAC_V3,
  txi_AKAP11_N1,
  txi_AKAP11_N1_1,
  txi_AKAP11_N1_2,
  txi_AKAP11_N1_3,
  txi_AKAP11_V1,
  txi_AKAP11_V1_1,
  txi_AKAP11_V1_2,
  txi_AKAP11_V1_3
)

res_names <- c(
  "res_DAP_N4",
  "res_DAP_N4_1",
  "res_DAP_N4_2",
  "res_DAP_N4_3",
  "res_DAP_V4",
  "res_DAP_V4_1",
  "res_DAP_V4_2",
  "res_DAP_V4_3",
  "res_DAP_N2",
  "res_DAP_N2_1",
  "res_DAP_N2_2",
  "res_DAP_N2_3",
  "res_DAP_V2",
  "res_DAP_V2_1",
  "res_DAP_V2_2",
  "res_DAP_V2_3",
  "res_DAP_V1",
  "res_DAP_V1_1",
  "res_DAP_V1_2",
  "res_DAP_V1_3",
  "res_RIPK2_V1",
  "res_RIPK2_V1_1",
  "res_RIPK2_V1_2",
  "res_RIPK2_V1_3",
  "res_RIPK2_N2",
  "res_RIPK2_N2_1",
  "res_RIPK2_N2_2",
  "res_RIPK2_N2_3",
  "res_NC1",
  "res_HLAC_N1",
  "res_HLAC_V3",
  "res_AKAP11_N1",
  "res_AKAP11_N1_1",
  "res_AKAP11_N1_2",
  "res_AKAP11_N1_3",
  "res_AKAP11_V1",
  "res_AKAP11_V1_1",
  "res_AKAP11_V1_2",
  "res_AKAP11_V1_3"
)

dds_names <- c(
  "dds_DAP_N4",
  "dds_DAP_N4_1",
  "dds_DAP_N4_2",
  "dds_DAP_N4_3",
  "dds_DAP_V4",
  "dds_DAP_V4_1",
  "dds_DAP_V4_2",
  "dds_DAP_V4_3",
  "dds_DAP_N2",
  "dds_DAP_N2_1",
  "dds_DAP_N2_2",
  "dds_DAP_N2_3",
  "dds_DAP_V2",
  "dds_DAP_V2_1",
  "dds_DAP_V2_2",
  "dds_DAP_V2_3",
  "dds_DAP_V1",
  "dds_DAP_V1_1",
  "dds_DAP_V1_2",
  "dds_DAP_V1_3",
  "dds_RIPK2_V1",
  "dds_RIPK2_V1_1",
  "dds_RIPK2_V1_2",
  "dds_RIPK2_V1_3",
  "dds_RIPK2_N2",
  "dds_RIPK2_N2_1",
  "dds_RIPK2_N2_2",
  "dds_RIPK2_N2_3",
  "dds_NC1",
  "dds_HLAC_N1",
  "dds_HLAC_V3",
  "dds_AKAP11_N1",
  "dds_AKAP11_N1_1",
  "dds_AKAP11_N1_2",
  "dds_AKAP11_N1_3",
  "dds_AKAP11_V1",
  "dds_AKAP11_V1_1",
  "dds_AKAP11_V1_2",
  "dds_AKAP11_V1_3"
)

res_df_names <- c(
  "res_df_DAP_N4",
  "res_df_DAP_N4_1",
  "res_df_DAP_N4_2",
  "res_df_DAP_N4_3",
  "res_df_DAP_V4",
  "res_df_DAP_V4_1",
  "res_df_DAP_V4_2",
  "res_df_DAP_V4_3",
  "res_df_DAP_N2",
  "res_df_DAP_N2_1",
  "res_df_DAP_N2_2",
  "res_df_DAP_N2_3",
  "res_df_DAP_V2",
  "res_df_DAP_V2_1",
  "res_df_DAP_V2_2",
  "res_df_DAP_V2_3",
  "res_df_DAP_V1",
  "res_df_DAP_V1_1",
  "res_df_DAP_V1_2",
  "res_df_DAP_V1_3",
  "res_df_RIPK2_V1",
  "res_df_RIPK2_V1_1",
  "res_df_RIPK2_V1_2",
  "res_df_RIPK2_V1_3",
  "res_df_RIPK2_N2",
  "res_df_RIPK2_N2_1",
  "res_df_RIPK2_N2_2",
  "res_df_RIPK2_N2_3",
  "res_df_NC1",
  "res_df_HLAC_N1",
  "res_df_HLAC_V3",
  "res_df_AKAP11_N1",
  "res_df_AKAP11_N1_1",
  "res_df_AKAP11_N1_2",
  "res_df_AKAP11_N1_3",
  "res_df_AKAP11_V1",
  "res_df_AKAP11_V1_1",
  "res_df_AKAP11_V1_2",
  "res_df_AKAP11_V1_3"
)

sample_names <- list(
  sampleTable_DAP_N4,
  sampleTable_DAP_N4_1,
  sampleTable_DAP_N4_2,
  sampleTable_DAP_N4_3,
  sampleTable_DAP_V4,
  sampleTable_DAP_V4_1,
  sampleTable_DAP_V4_2,
  sampleTable_DAP_V4_3,
  sampleTable_DAP_N2,
  sampleTable_DAP_N2_1,
  sampleTable_DAP_N2_2,
  sampleTable_DAP_N2_3,
  sampleTable_DAP_V2,
  sampleTable_DAP_V2_1,
  sampleTable_DAP_V2_2,
  sampleTable_DAP_V2_3,
  sampleTable_DAP_V1,
  sampleTable_DAP_V1_1,
  sampleTable_DAP_V1_2,
  sampleTable_DAP_V1_3,
  sampleTable_RIPK2_V1,
  sampleTable_RIPK2_V1_1,
  sampleTable_RIPK2_V1_2,
  sampleTable_RIPK2_V1_3,
  sampleTable_RIPK2_N2,
  sampleTable_RIPK2_N2_1,
  sampleTable_RIPK2_N2_2,
  sampleTable_RIPK2_N2_3,
  sampleTable_NC1,
  sampleTable_HLAC_N1,
  sampleTable_HLAC_V3,
  sampleTable_AKAP11_N1,
  sampleTable_AKAP11_N1_1,
  sampleTable_AKAP11_N1_2,
  sampleTable_AKAP11_N1_3,
  sampleTable_AKAP11_V1,
  sampleTable_AKAP11_V1_1,
  sampleTable_AKAP11_V1_2,
  sampleTable_AKAP11_V1_3
)

i <- 1
for (x in txi_names) {

  txi <- x
  sample_info <- sample_names[[i]]

  keep <- rowSums(txi$counts >= 10) >= 2
  counts_filtered <- txi$counts[keep, ]

  y <- DGEList(counts = counts_filtered)
  y <- calcNormFactors(y)

  sample_info$condition <- relevel(factor(sample_info$condition), ref = "Other")
  design <- model.matrix(~ condition, data = sample_info)

  v <- voomWithQualityWeights(y, design, plot = FALSE)

  fit <- lmFit(v, design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)

  res <- topTable(fit, coef = 2, number = Inf, sort.by = "none")

  assign(res_names[i], res)
  assign(dds_names[i], v)

  i <- i + 1
  print(i)
}

res_vars <- list(
  res_DAP_N4,
  res_DAP_N4_1,
  res_DAP_N4_2,
  res_DAP_N4_3,
  res_DAP_V4,
  res_DAP_V4_1,
  res_DAP_V4_2,
  res_DAP_V4_3,
  res_DAP_N2,
  res_DAP_N2_1,
  res_DAP_N2_2,
  res_DAP_N2_3,
  res_DAP_V2,
  res_DAP_V2_1,
  res_DAP_V2_2,
  res_DAP_V2_3,
  res_DAP_V1,
  res_DAP_V1_1,
  res_DAP_V1_2,
  res_DAP_V1_3,
  res_RIPK2_V1,
  res_RIPK2_V1_1,
  res_RIPK2_V1_2,
  res_RIPK2_V1_3,
  res_RIPK2_N2,
  res_RIPK2_N2_1,
  res_RIPK2_N2_2,
  res_RIPK2_N2_3,
  res_NC1,
  res_HLAC_N1,
  res_HLAC_V3,
  res_AKAP11_N1,
  res_AKAP11_N1_1,
  res_AKAP11_N1_2,
  res_AKAP11_N1_3,
  res_AKAP11_V1,
  res_AKAP11_V1_1,
  res_AKAP11_V1_2,
  res_AKAP11_V1_3
)

dds_vars <- list(
  dds_DAP_N4,
  dds_DAP_N4_1,
  dds_DAP_N4_2,
  dds_DAP_N4_3,
  dds_DAP_V4,
  dds_DAP_V4_1,
  dds_DAP_V4_2,
  dds_DAP_V4_3,
  dds_DAP_N2,
  dds_DAP_N2_1,
  dds_DAP_N2_2,
  dds_DAP_N2_3,
  dds_DAP_V2,
  dds_DAP_V2_1,
  dds_DAP_V2_2,
  dds_DAP_V2_3,
  dds_DAP_V1,
  dds_DAP_V1_1,
  dds_DAP_V1_2,
  dds_DAP_V1_3,
  dds_RIPK2_V1,
  dds_RIPK2_V1_1,
  dds_RIPK2_V1_2,
  dds_RIPK2_V1_3,
  dds_RIPK2_N2,
  dds_RIPK2_N2_1,
  dds_RIPK2_N2_2,
  dds_RIPK2_N2_3,
  dds_NC1,
  dds_HLAC_N1,
  dds_HLAC_V3,
  dds_AKAP11_N1,
  dds_AKAP11_N1_1,
  dds_AKAP11_N1_2,
  dds_AKAP11_N1_3,
  dds_AKAP11_V1,
  dds_AKAP11_V1_1,
  dds_AKAP11_V1_2,
  dds_AKAP11_V1_3
)

## output tables with gene name
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

i = 1
for (x in res_vars) {
  deg_genes <- rownames(x)

  gene_conversion <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                           filters = "ensembl_gene_id",
                           values = deg_genes,
                           mart = ensembl)

  res_df <- as.data.frame(x)
  res_df$ensembl_gene_id <- rownames(res_df)
  res_df <- merge(res_df, gene_conversion, by = "ensembl_gene_id", all.x = TRUE)

  assign(res_df_names[i], res_df)

  i = i + 1
}

# ==================== DAP ====================
subset(res_df_DAP_N4, hgnc_symbol == "DAP")
subset(res_df_DAP_N4_1, hgnc_symbol == "DAP")
subset(res_df_DAP_N4_2, hgnc_symbol == "DAP")
subset(res_df_DAP_N4_3, hgnc_symbol == "DAP")

subset(res_df_DAP_V4, hgnc_symbol == "DAP")
subset(res_df_DAP_V4_1, hgnc_symbol == "DAP")
subset(res_df_DAP_V4_2, hgnc_symbol == "DAP")
subset(res_df_DAP_V4_3, hgnc_symbol == "DAP")

subset(res_df_DAP_N2, hgnc_symbol == "DAP")
subset(res_df_DAP_N2_1, hgnc_symbol == "DAP")
subset(res_df_DAP_N2_2, hgnc_symbol == "DAP")
subset(res_df_DAP_N2_3, hgnc_symbol == "DAP")

subset(res_df_DAP_V2, hgnc_symbol == "DAP")
subset(res_df_DAP_V2_1, hgnc_symbol == "DAP")
subset(res_df_DAP_V2_2, hgnc_symbol == "DAP")
subset(res_df_DAP_V2_3, hgnc_symbol == "DAP")

subset(res_df_DAP_V1, hgnc_symbol == "DAP")
subset(res_df_DAP_V1_1, hgnc_symbol == "DAP")
subset(res_df_DAP_V1_2, hgnc_symbol == "DAP")
subset(res_df_DAP_V1_3, hgnc_symbol == "DAP")

# ==================== RIPK2 ====================
subset(res_df_RIPK2_V1, hgnc_symbol == "RIPK2")
subset(res_df_RIPK2_V1_1, hgnc_symbol == "RIPK2")
subset(res_df_RIPK2_V1_2, hgnc_symbol == "RIPK2")
subset(res_df_RIPK2_V1_3, hgnc_symbol == "RIPK2")

subset(res_df_RIPK2_N2, hgnc_symbol == "RIPK2")
subset(res_df_RIPK2_N2_1, hgnc_symbol == "RIPK2")
subset(res_df_RIPK2_N2_2, hgnc_symbol == "RIPK2")
subset(res_df_RIPK2_N2_3, hgnc_symbol == "RIPK2")

# ==================== HLAC ====================
subset(res_df_HLAC_N1, hgnc_symbol == "HLA-C")
subset(res_df_HLAC_V3, hgnc_symbol == "HLA-C")

# ==================== AKAP11 ====================
subset(res_df_AKAP11_N1, hgnc_symbol == "AKAP11")
subset(res_df_AKAP11_N1_1, hgnc_symbol == "AKAP11")
subset(res_df_AKAP11_N1_2, hgnc_symbol == "AKAP11")
subset(res_df_AKAP11_N1_3, hgnc_symbol == "AKAP11")

subset(res_df_AKAP11_V1, hgnc_symbol == "AKAP11")
subset(res_df_AKAP11_V1_1, hgnc_symbol == "AKAP11")
subset(res_df_AKAP11_V1_2, hgnc_symbol == "AKAP11")
subset(res_df_AKAP11_V1_3, hgnc_symbol == "AKAP11")

##################################### 3rd batch ###################################
samples <- 1:33
paths <- sprintf("/media/ming/Extra_SSD_4TB1/eCROPseq/20250926_bulk_cell_RNA_seq_3rd/salmon_quant/output_quant_sample%d_S%d/quant.sf", samples, samples)

txi_CARS1_N1 <- tximport(c(paths[1:3],paths[19],paths[21:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_CARS1_N1_1 <- tximport(c(paths[1],paths[19],paths[21:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_CARS1_N1_2 <- tximport(c(paths[2],paths[19],paths[21:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_CARS1_N1_3 <- tximport(c(paths[3],paths[19],paths[21:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")

txi_CARS1_V1 <- tximport(c(paths[4:6],paths[19],paths[21:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_CARS1_V1_1 <- tximport(c(paths[4],paths[19],paths[21:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_CARS1_V1_2 <- tximport(c(paths[5],paths[19],paths[21:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_CARS1_V1_3 <- tximport(c(paths[6],paths[19],paths[21:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")

txi_CARS1_V3 <- tximport(c(paths[7:9],paths[19],paths[21:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_CARS1_V3_1 <- tximport(c(paths[7],paths[19],paths[21:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_CARS1_V3_2 <- tximport(c(paths[8],paths[19],paths[21:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_CARS1_V3_3 <- tximport(c(paths[9],paths[19],paths[21:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")

txi_CARS1_V5 <- tximport(c(paths[10:12],paths[19],paths[21:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_CARS1_V5_1 <- tximport(c(paths[10],paths[19],paths[21:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_CARS1_V5_2 <- tximport(c(paths[11],paths[19],paths[21:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_CARS1_V5_3 <- tximport(c(paths[12],paths[19],paths[21:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")

txi_CARS1_V6 <- tximport(c(paths[13:15],paths[19],paths[21:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_CARS1_V6_1 <- tximport(c(paths[13],paths[19],paths[21:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_CARS1_V6_2 <- tximport(c(paths[14],paths[19],paths[21:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_CARS1_V6_3 <- tximport(c(paths[15],paths[19],paths[21:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")

txi_CARS1_V7 <- tximport(c(paths[16:18],paths[19],paths[21:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_CARS1_V7_1 <- tximport(c(paths[16],paths[19],paths[21:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_CARS1_V7_2 <- tximport(c(paths[17],paths[19],paths[21:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_CARS1_V7_3 <- tximport(c(paths[18],paths[19],paths[21:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")

txi_PTPN_N1 <- tximport(c(paths[19],paths[21],paths[1:18],paths[28:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_PTPN_N1_1 <- tximport(c(paths[19],paths[1:18],paths[28:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_PTPN_N1_2 <- tximport(c(paths[21],paths[1:18],paths[28:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_PTPN_N1_3 <- tximport(c(paths[21],paths[1:18],paths[28:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")

txi_PTPN_N2 <- tximport(c(paths[22:24],paths[1:18],paths[28:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_PTPN_N2_1 <- tximport(c(paths[22],paths[1:18],paths[28:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_PTPN_N2_2 <- tximport(c(paths[23],paths[1:18],paths[28:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_PTPN_N2_3 <- tximport(c(paths[24],paths[1:18],paths[28:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")

txi_PTPN_V1 <- tximport(c(paths[25:27],paths[1:18],paths[28:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_PTPN_V1_1 <- tximport(c(paths[25],paths[1:18],paths[28:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_PTPN_V1_2 <- tximport(c(paths[26],paths[1:18],paths[28:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_PTPN_V1_3 <- tximport(c(paths[27],paths[1:18],paths[28:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")

txi_RIPK_V1 <- tximport(c(paths[28:30],paths[1:19],paths[21:27],paths[31:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_RIPK_V1_1 <- tximport(c(paths[28],paths[1:19],paths[21:27],paths[31:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_RIPK_V1_2 <- tximport(c(paths[29],paths[1:19],paths[21:27],paths[31:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_RIPK_V1_3 <- tximport(c(paths[30],paths[1:19],paths[21:27],paths[31:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")

txi_NC1 <- tximport(c(paths[31:33],paths[1:19],paths[21:30]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
samples_order <- as.character(1:33)

# ==================== CARS1 ====================

sampleTable_CARS1_N1 <- data.frame(
  samples = c(samples_order[c(1:3, 19, 21:33)]),
  condition = c(rep("CARS1_N1", 3), rep("Other", 14))
)
rownames(sampleTable_CARS1_N1) <- sampleTable_CARS1_N1$samples

sampleTable_CARS1_N1_1 <- data.frame(
  samples = c(samples_order[c(1, 19, 21:33)]),
  condition = c("CARS1_N1", rep("Other", 14))
)
rownames(sampleTable_CARS1_N1_1) <- sampleTable_CARS1_N1_1$samples

sampleTable_CARS1_N1_2 <- data.frame(
  samples = c(samples_order[c(2, 19, 21:33)]),
  condition = c("CARS1_N1", rep("Other", 14))
)
rownames(sampleTable_CARS1_N1_2) <- sampleTable_CARS1_N1_2$samples

sampleTable_CARS1_N1_3 <- data.frame(
  samples = c(samples_order[c(3, 19, 21:33)]),
  condition = c("CARS1_N1", rep("Other", 14))
)
rownames(sampleTable_CARS1_N1_3) <- sampleTable_CARS1_N1_3$samples

sampleTable_CARS1_V1 <- data.frame(
  samples = c(samples_order[c(4:6, 19, 21:33)]),
  condition = c(rep("CARS1_V1", 3), rep("Other", 14))
)
rownames(sampleTable_CARS1_V1) <- sampleTable_CARS1_V1$samples

sampleTable_CARS1_V1_1 <- data.frame(
  samples = c(samples_order[c(4, 19, 21:33)]),
  condition = c("CARS1_V1", rep("Other", 14))
)
rownames(sampleTable_CARS1_V1_1) <- sampleTable_CARS1_V1_1$samples

sampleTable_CARS1_V1_2 <- data.frame(
  samples = c(samples_order[c(5, 19, 21:33)]),
  condition = c("CARS1_V1", rep("Other", 14))
)
rownames(sampleTable_CARS1_V1_2) <- sampleTable_CARS1_V1_2$samples

sampleTable_CARS1_V1_3 <- data.frame(
  samples = c(samples_order[c(6, 19, 21:33)]),
  condition = c("CARS1_V1", rep("Other", 14))
)
rownames(sampleTable_CARS1_V1_3) <- sampleTable_CARS1_V1_3$samples

sampleTable_CARS1_V3 <- data.frame(
  samples = c(samples_order[c(7:9, 19, 21:33)]),
  condition = c(rep("CARS1_V3", 3), rep("Other", 14))
)
rownames(sampleTable_CARS1_V3) <- sampleTable_CARS1_V3$samples

sampleTable_CARS1_V3_1 <- data.frame(
  samples = c(samples_order[c(7, 19, 21:33)]),
  condition = c("CARS1_V3", rep("Other", 14))
)
rownames(sampleTable_CARS1_V3_1) <- sampleTable_CARS1_V3_1$samples

sampleTable_CARS1_V3_2 <- data.frame(
  samples = c(samples_order[c(8, 19, 21:33)]),
  condition = c("CARS1_V3", rep("Other", 14))
)
rownames(sampleTable_CARS1_V3_2) <- sampleTable_CARS1_V3_2$samples

sampleTable_CARS1_V3_3 <- data.frame(
  samples = c(samples_order[c(9, 19, 21:33)]),
  condition = c("CARS1_V3", rep("Other", 14))
)
rownames(sampleTable_CARS1_V3_3) <- sampleTable_CARS1_V3_3$samples

sampleTable_CARS1_V5 <- data.frame(
  samples = c(samples_order[c(10:12, 19, 21:33)]),
  condition = c(rep("CARS1_V5", 3), rep("Other", 14))
)
rownames(sampleTable_CARS1_V5) <- sampleTable_CARS1_V5$samples

sampleTable_CARS1_V5_1 <- data.frame(
  samples = c(samples_order[c(10, 19, 21:33)]),
  condition = c("CARS1_V5", rep("Other", 14))
)
rownames(sampleTable_CARS1_V5_1) <- sampleTable_CARS1_V5_1$samples

sampleTable_CARS1_V5_2 <- data.frame(
  samples = c(samples_order[c(11, 19, 21:33)]),
  condition = c("CARS1_V5", rep("Other", 14))
)
rownames(sampleTable_CARS1_V5_2) <- sampleTable_CARS1_V5_2$samples

sampleTable_CARS1_V5_3 <- data.frame(
  samples = c(samples_order[c(12, 19, 21:33)]),
  condition = c("CARS1_V5", rep("Other", 14))
)
rownames(sampleTable_CARS1_V5_3) <- sampleTable_CARS1_V5_3$samples

sampleTable_CARS1_V6 <- data.frame(
  samples = c(samples_order[c(13:15, 19, 21:33)]),
  condition = c(rep("CARS1_V6", 3), rep("Other", 14))
)
rownames(sampleTable_CARS1_V6) <- sampleTable_CARS1_V6$samples

sampleTable_CARS1_V6_1 <- data.frame(
  samples = c(samples_order[c(13, 19, 21:33)]),
  condition = c("CARS1_V6", rep("Other", 14))
)
rownames(sampleTable_CARS1_V6_1) <- sampleTable_CARS1_V6_1$samples

sampleTable_CARS1_V6_2 <- data.frame(
  samples = c(samples_order[c(14, 19, 21:33)]),
  condition = c("CARS1_V6", rep("Other", 14))
)
rownames(sampleTable_CARS1_V6_2) <- sampleTable_CARS1_V6_2$samples

sampleTable_CARS1_V6_3 <- data.frame(
  samples = c(samples_order[c(15, 19, 21:33)]),
  condition = c("CARS1_V6", rep("Other", 14))
)
rownames(sampleTable_CARS1_V6_3) <- sampleTable_CARS1_V6_3$samples

sampleTable_CARS1_V7 <- data.frame(
  samples = c(samples_order[c(16:18, 19, 21:33)]),
  condition = c(rep("CARS1_V7", 3), rep("Other", 14))
)
rownames(sampleTable_CARS1_V7) <- sampleTable_CARS1_V7$samples

sampleTable_CARS1_V7_1 <- data.frame(
  samples = c(samples_order[c(16, 19, 21:33)]),
  condition = c("CARS1_V7", rep("Other", 14))
)
rownames(sampleTable_CARS1_V7_1) <- sampleTable_CARS1_V7_1$samples

sampleTable_CARS1_V7_2 <- data.frame(
  samples = c(samples_order[c(17, 19, 21:33)]),
  condition = c("CARS1_V7", rep("Other", 14))
)
rownames(sampleTable_CARS1_V7_2) <- sampleTable_CARS1_V7_2$samples

sampleTable_CARS1_V7_3 <- data.frame(
  samples = c(samples_order[c(18, 19, 21:33)]),
  condition = c("CARS1_V7", rep("Other", 14))
)
rownames(sampleTable_CARS1_V7_3) <- sampleTable_CARS1_V7_3$samples

samples_order <- as.character(1:33)

# ==================== PTPN_N1 ====================

sampleTable_PTPN_N1 <- data.frame(
  samples = c(samples_order[c(19, 21, 1:18, 28:33)]),
  condition = c(rep("PTPN_N1", 2), rep("Other", 24))
)
rownames(sampleTable_PTPN_N1) <- sampleTable_PTPN_N1$samples

sampleTable_PTPN_N1_1 <- data.frame(
  samples = c(samples_order[c(19, 1:18, 28:33)]),
  condition = c("PTPN_N1", rep("Other", 24))
)
rownames(sampleTable_PTPN_N1_1) <- sampleTable_PTPN_N1_1$samples

sampleTable_PTPN_N1_2 <- data.frame(
  samples = c(samples_order[c(21, 1:18, 28:33)]),
  condition = c("PTPN_N1", rep("Other", 24))
)
rownames(sampleTable_PTPN_N1_2) <- sampleTable_PTPN_N1_2$samples

sampleTable_PTPN_N1_3 <- data.frame(
  samples = c(samples_order[c(21, 1:18, 28:33)]),
  condition = c("PTPN_N1", rep("Other", 24))
)
rownames(sampleTable_PTPN_N1_3) <- sampleTable_PTPN_N1_3$samples

# ==================== PTPN_N2 ====================

sampleTable_PTPN_N2 <- data.frame(
  samples = c(samples_order[c(22:24, 1:18, 28:33)]),
  condition = c(rep("PTPN_N2", 3), rep("Other", 24))
)
rownames(sampleTable_PTPN_N2) <- sampleTable_PTPN_N2$samples

sampleTable_PTPN_N2_1 <- data.frame(
  samples = c(samples_order[c(22, 1:18, 28:33)]),
  condition = c("PTPN_N2", rep("Other", 24))
)
rownames(sampleTable_PTPN_N2_1) <- sampleTable_PTPN_N2_1$samples

sampleTable_PTPN_N2_2 <- data.frame(
  samples = c(samples_order[c(23, 1:18, 28:33)]),
  condition = c("PTPN_N2", rep("Other", 24))
)
rownames(sampleTable_PTPN_N2_2) <- sampleTable_PTPN_N2_2$samples

sampleTable_PTPN_N2_3 <- data.frame(
  samples = c(samples_order[c(24, 1:18, 28:33)]),
  condition = c("PTPN_N2", rep("Other", 24))
)
rownames(sampleTable_PTPN_N2_3) <- sampleTable_PTPN_N2_3$samples

# ==================== PTPN_V1 ====================

sampleTable_PTPN_V1 <- data.frame(
  samples = c(samples_order[c(25:27, 1:18, 28:33)]),
  condition = c(rep("PTPN_V1", 3), rep("Other", 24))
)
rownames(sampleTable_PTPN_V1) <- sampleTable_PTPN_V1$samples

sampleTable_PTPN_V1_1 <- data.frame(
  samples = c(samples_order[c(25, 1:18, 28:33)]),
  condition = c("PTPN_V1", rep("Other", 24))
)
rownames(sampleTable_PTPN_V1_1) <- sampleTable_PTPN_V1_1$samples

sampleTable_PTPN_V1_2 <- data.frame(
  samples = c(samples_order[c(26, 1:18, 28:33)]),
  condition = c("PTPN_V1", rep("Other", 24))
)
rownames(sampleTable_PTPN_V1_2) <- sampleTable_PTPN_V1_2$samples

sampleTable_PTPN_V1_3 <- data.frame(
  samples = c(samples_order[c(27, 1:18, 28:33)]),
  condition = c("PTPN_V1", rep("Other", 24))
)
rownames(sampleTable_PTPN_V1_3) <- sampleTable_PTPN_V1_3$samples

# ==================== RIPK_V1 ====================

sampleTable_RIPK_V1 <- data.frame(
  samples = c(samples_order[c(28:30, 1:19, 21:27, 31:33)]),
  condition = c(rep("RIPK_V1", 3), rep("Other", 29))
)
rownames(sampleTable_RIPK_V1) <- sampleTable_RIPK_V1$samples

sampleTable_RIPK_V1_1 <- data.frame(
  samples = c(samples_order[c(28, 1:19, 21:27, 31:33)]),
  condition = c("RIPK_V1", rep("Other", 29))
)
rownames(sampleTable_RIPK_V1_1) <- sampleTable_RIPK_V1_1$samples

sampleTable_RIPK_V1_2 <- data.frame(
  samples = c(samples_order[c(29, 1:19, 21:27, 31:33)]),
  condition = c("RIPK_V1", rep("Other", 29))
)
rownames(sampleTable_RIPK_V1_2) <- sampleTable_RIPK_V1_2$samples

sampleTable_RIPK_V1_3 <- data.frame(
  samples = c(samples_order[c(30, 1:19, 21:27, 31:33)]),
  condition = c("RIPK_V1", rep("Other", 29))
)
rownames(sampleTable_RIPK_V1_3) <- sampleTable_RIPK_V1_3$samples

# ==================== txi_names ====================
txi_names <- list(
  txi_CARS1_N1,
  txi_CARS1_N1_1,
  txi_CARS1_N1_2,
  txi_CARS1_N1_3,
  txi_CARS1_V1,
  txi_CARS1_V1_1,
  txi_CARS1_V1_2,
  txi_CARS1_V1_3,
  txi_CARS1_V3,
  txi_CARS1_V3_1,
  txi_CARS1_V3_2,
  txi_CARS1_V3_3,
  txi_CARS1_V5,
  txi_CARS1_V5_1,
  txi_CARS1_V5_2,
  txi_CARS1_V5_3,
  txi_CARS1_V6,
  txi_CARS1_V6_1,
  txi_CARS1_V6_2,
  txi_CARS1_V6_3,
  txi_CARS1_V7,
  txi_CARS1_V7_1,
  txi_CARS1_V7_2,
  txi_CARS1_V7_3,
  txi_PTPN_N1,
  txi_PTPN_N1_1,
  txi_PTPN_N1_2,
  txi_PTPN_N1_3,
  txi_PTPN_N2,
  txi_PTPN_N2_1,
  txi_PTPN_N2_2,
  txi_PTPN_N2_3,
  txi_PTPN_V1,
  txi_PTPN_V1_1,
  txi_PTPN_V1_2,
  txi_PTPN_V1_3,
  txi_RIPK_V1,
  txi_RIPK_V1_1,
  txi_RIPK_V1_2,
  txi_RIPK_V1_3,
  txi_NC1
)

# ==================== res_names ====================
res_names <- c(
  "res_CARS1_N1",
  "res_CARS1_N1_1",
  "res_CARS1_N1_2",
  "res_CARS1_N1_3",
  "res_CARS1_V1",
  "res_CARS1_V1_1",
  "res_CARS1_V1_2",
  "res_CARS1_V1_3",
  "res_CARS1_V3",
  "res_CARS1_V3_1",
  "res_CARS1_V3_2",
  "res_CARS1_V3_3",
  "res_CARS1_V5",
  "res_CARS1_V5_1",
  "res_CARS1_V5_2",
  "res_CARS1_V5_3",
  "res_CARS1_V6",
  "res_CARS1_V6_1",
  "res_CARS1_V6_2",
  "res_CARS1_V6_3",
  "res_CARS1_V7",
  "res_CARS1_V7_1",
  "res_CARS1_V7_2",
  "res_CARS1_V7_3",
  "res_PTPN_N1",
  "res_PTPN_N1_1",
  "res_PTPN_N1_2",
  "res_PTPN_N1_3",
  "res_PTPN_N2",
  "res_PTPN_N2_1",
  "res_PTPN_N2_2",
  "res_PTPN_N2_3",
  "res_PTPN_V1",
  "res_PTPN_V1_1",
  "res_PTPN_V1_2",
  "res_PTPN_V1_3",
  "res_RIPK_V1",
  "res_RIPK_V1_1",
  "res_RIPK_V1_2",
  "res_RIPK_V1_3",
  "res_NC1"
)

# ==================== dds_names ====================
dds_names <- c(
  "dds_CARS1_N1",
  "dds_CARS1_N1_1",
  "dds_CARS1_N1_2",
  "dds_CARS1_N1_3",
  "dds_CARS1_V1",
  "dds_CARS1_V1_1",
  "dds_CARS1_V1_2",
  "dds_CARS1_V1_3",
  "dds_CARS1_V3",
  "dds_CARS1_V3_1",
  "dds_CARS1_V3_2",
  "dds_CARS1_V3_3",
  "dds_CARS1_V5",
  "dds_CARS1_V5_1",
  "dds_CARS1_V5_2",
  "dds_CARS1_V5_3",
  "dds_CARS1_V6",
  "dds_CARS1_V6_1",
  "dds_CARS1_V6_2",
  "dds_CARS1_V6_3",
  "dds_CARS1_V7",
  "dds_CARS1_V7_1",
  "dds_CARS1_V7_2",
  "dds_CARS1_V7_3",
  "dds_PTPN_N1",
  "dds_PTPN_N1_1",
  "dds_PTPN_N1_2",
  "dds_PTPN_N1_3",
  "dds_PTPN_N2",
  "dds_PTPN_N2_1",
  "dds_PTPN_N2_2",
  "dds_PTPN_N2_3",
  "dds_PTPN_V1",
  "dds_PTPN_V1_1",
  "dds_PTPN_V1_2",
  "dds_PTPN_V1_3",
  "dds_RIPK_V1",
  "dds_RIPK_V1_1",
  "dds_RIPK_V1_2",
  "dds_RIPK_V1_3",
  "dds_NC1"
)

# ==================== res_df_names ====================
res_df_names <- c(
  "res_df_CARS1_N1",
  "res_df_CARS1_N1_1",
  "res_df_CARS1_N1_2",
  "res_df_CARS1_N1_3",
  "res_df_CARS1_V1",
  "res_df_CARS1_V1_1",
  "res_df_CARS1_V1_2",
  "res_df_CARS1_V1_3",
  "res_df_CARS1_V3",
  "res_df_CARS1_V3_1",
  "res_df_CARS1_V3_2",
  "res_df_CARS1_V3_3",
  "res_df_CARS1_V5",
  "res_df_CARS1_V5_1",
  "res_df_CARS1_V5_2",
  "res_df_CARS1_V5_3",
  "res_df_CARS1_V6",
  "res_df_CARS1_V6_1",
  "res_df_CARS1_V6_2",
  "res_df_CARS1_V6_3",
  "res_df_CARS1_V7",
  "res_df_CARS1_V7_1",
  "res_df_CARS1_V7_2",
  "res_df_CARS1_V7_3",
  "res_df_PTPN_N1",
  "res_df_PTPN_N1_1",
  "res_df_PTPN_N1_2",
  "res_df_PTPN_N1_3",
  "res_df_PTPN_N2",
  "res_df_PTPN_N2_1",
  "res_df_PTPN_N2_2",
  "res_df_PTPN_N2_3",
  "res_df_PTPN_V1",
  "res_df_PTPN_V1_1",
  "res_df_PTPN_V1_2",
  "res_df_PTPN_V1_3",
  "res_df_RIPK_V1",
  "res_df_RIPK_V1_1",
  "res_df_RIPK_V1_2",
  "res_df_RIPK_V1_3",
  "res_df_NC1"
)

# ==================== sample_names ====================
sample_names <- list(
  sampleTable_CARS1_N1,
  sampleTable_CARS1_N1_1,
  sampleTable_CARS1_N1_2,
  sampleTable_CARS1_N1_3,
  sampleTable_CARS1_V1,
  sampleTable_CARS1_V1_1,
  sampleTable_CARS1_V1_2,
  sampleTable_CARS1_V1_3,
  sampleTable_CARS1_V3,
  sampleTable_CARS1_V3_1,
  sampleTable_CARS1_V3_2,
  sampleTable_CARS1_V3_3,
  sampleTable_CARS1_V5,
  sampleTable_CARS1_V5_1,
  sampleTable_CARS1_V5_2,
  sampleTable_CARS1_V5_3,
  sampleTable_CARS1_V6,
  sampleTable_CARS1_V6_1,
  sampleTable_CARS1_V6_2,
  sampleTable_CARS1_V6_3,
  sampleTable_CARS1_V7,
  sampleTable_CARS1_V7_1,
  sampleTable_CARS1_V7_2,
  sampleTable_CARS1_V7_3,
  sampleTable_PTPN_N1,
  sampleTable_PTPN_N1_1,
  sampleTable_PTPN_N1_2,
  sampleTable_PTPN_N1_3,
  sampleTable_PTPN_N2,
  sampleTable_PTPN_N2_1,
  sampleTable_PTPN_N2_2,
  sampleTable_PTPN_N2_3,
  sampleTable_PTPN_V1,
  sampleTable_PTPN_V1_1,
  sampleTable_PTPN_V1_2,
  sampleTable_PTPN_V1_3,
  sampleTable_RIPK_V1,
  sampleTable_RIPK_V1_1,
  sampleTable_RIPK_V1_2,
  sampleTable_RIPK_V1_3,
  sampleTable_NC1
)

i <- 1
for (x in txi_names) {

  txi <- x
  sample_info <- sample_names[[i]]

  keep <- rowSums(txi$counts >= 10) >= 2
  counts_filtered <- txi$counts[keep, ]

  y <- DGEList(counts = counts_filtered)
  y <- calcNormFactors(y)

  sample_info$condition <- relevel(factor(sample_info$condition), ref = "Other")
  design <- model.matrix(~ condition, data = sample_info)

  v <- voomWithQualityWeights(y, design, plot = FALSE)

  fit <- lmFit(v, design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)

  res <- topTable(fit, coef = 2, number = Inf, sort.by = "none")

  assign(res_names[i], res)
  assign(dds_names[i], v)

  i <- i + 1
  print(i)
}

res_vars <- list(
  res_CARS1_N1,
  res_CARS1_N1_1,
  res_CARS1_N1_2,
  res_CARS1_N1_3,
  res_CARS1_V1,
  res_CARS1_V1_1,
  res_CARS1_V1_2,
  res_CARS1_V1_3,
  res_CARS1_V3,
  res_CARS1_V3_1,
  res_CARS1_V3_2,
  res_CARS1_V3_3,
  res_CARS1_V5,
  res_CARS1_V5_1,
  res_CARS1_V5_2,
  res_CARS1_V5_3,
  res_CARS1_V6,
  res_CARS1_V6_1,
  res_CARS1_V6_2,
  res_CARS1_V6_3,
  res_CARS1_V7,
  res_CARS1_V7_1,
  res_CARS1_V7_2,
  res_CARS1_V7_3,
  res_PTPN_N1,
  res_PTPN_N1_1,
  res_PTPN_N1_2,
  res_PTPN_N1_3,
  res_PTPN_N2,
  res_PTPN_N2_1,
  res_PTPN_N2_2,
  res_PTPN_N2_3,
  res_PTPN_V1,
  res_PTPN_V1_1,
  res_PTPN_V1_2,
  res_PTPN_V1_3,
  res_RIPK_V1,
  res_RIPK_V1_1,
  res_RIPK_V1_2,
  res_RIPK_V1_3,
  res_NC1
)

dds_vars <- list(
  dds_CARS1_N1,
  dds_CARS1_N1_1,
  dds_CARS1_N1_2,
  dds_CARS1_N1_3,
  dds_CARS1_V1,
  dds_CARS1_V1_1,
  dds_CARS1_V1_2,
  dds_CARS1_V1_3,
  dds_CARS1_V3,
  dds_CARS1_V3_1,
  dds_CARS1_V3_2,
  dds_CARS1_V3_3,
  dds_CARS1_V5,
  dds_CARS1_V5_1,
  dds_CARS1_V5_2,
  dds_CARS1_V5_3,
  dds_CARS1_V6,
  dds_CARS1_V6_1,
  dds_CARS1_V6_2,
  dds_CARS1_V6_3,
  dds_CARS1_V7,
  dds_CARS1_V7_1,
  dds_CARS1_V7_2,
  dds_CARS1_V7_3,
  dds_PTPN_N1,
  dds_PTPN_N1_1,
  dds_PTPN_N1_2,
  dds_PTPN_N1_3,
  dds_PTPN_N2,
  dds_PTPN_N2_1,
  dds_PTPN_N2_2,
  dds_PTPN_N2_3,
  dds_PTPN_V1,
  dds_PTPN_V1_1,
  dds_PTPN_V1_2,
  dds_PTPN_V1_3,
  dds_RIPK_V1,
  dds_RIPK_V1_1,
  dds_RIPK_V1_2,
  dds_RIPK_V1_3,
  dds_NC1
)

## output tables with gene name
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

i = 1
for (x in res_vars) {
  deg_genes <- rownames(x)

  gene_conversion <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                           filters = "ensembl_gene_id",
                           values = deg_genes,
                           mart = ensembl)

  res_df <- as.data.frame(x)
  res_df$ensembl_gene_id <- rownames(res_df)
  res_df <- merge(res_df, gene_conversion, by = "ensembl_gene_id", all.x = TRUE)

  assign(res_df_names[i], res_df)

  i = i + 1
}

# ==================== CARS1 ====================

subset(res_df_CARS1_N1,   hgnc_symbol == "CARS1")
subset(res_df_CARS1_N1_1, hgnc_symbol == "CARS1")
subset(res_df_CARS1_N1_2, hgnc_symbol == "CARS1")
subset(res_df_CARS1_N1_3, hgnc_symbol == "CARS1")

subset(res_df_CARS1_V1,   hgnc_symbol == "CARS1")
subset(res_df_CARS1_V1_1, hgnc_symbol == "CARS1")
subset(res_df_CARS1_V1_2, hgnc_symbol == "CARS1")
subset(res_df_CARS1_V1_3, hgnc_symbol == "CARS1")

subset(res_df_CARS1_V3,   hgnc_symbol == "CARS1")
subset(res_df_CARS1_V3_1, hgnc_symbol == "CARS1")
subset(res_df_CARS1_V3_2, hgnc_symbol == "CARS1")
subset(res_df_CARS1_V3_3, hgnc_symbol == "CARS1")

subset(res_df_CARS1_V5,   hgnc_symbol == "CARS1")
subset(res_df_CARS1_V5_1, hgnc_symbol == "CARS1")
subset(res_df_CARS1_V5_2, hgnc_symbol == "CARS1")
subset(res_df_CARS1_V5_3, hgnc_symbol == "CARS1")

subset(res_df_CARS1_V6,   hgnc_symbol == "CARS1")
subset(res_df_CARS1_V6_1, hgnc_symbol == "CARS1")
subset(res_df_CARS1_V6_2, hgnc_symbol == "CARS1")
subset(res_df_CARS1_V6_3, hgnc_symbol == "CARS1")

subset(res_df_CARS1_V7,   hgnc_symbol == "CARS1")
subset(res_df_CARS1_V7_1, hgnc_symbol == "CARS1")
subset(res_df_CARS1_V7_2, hgnc_symbol == "CARS1")
subset(res_df_CARS1_V7_3, hgnc_symbol == "CARS1")

# ==================== PTPN22 ====================

subset(res_df_PTPN_N1,   hgnc_symbol == "PTPN22")
subset(res_df_PTPN_N1_1, hgnc_symbol == "PTPN22")
subset(res_df_PTPN_N1_2, hgnc_symbol == "PTPN22")
subset(res_df_PTPN_N1_3, hgnc_symbol == "PTPN22")

subset(res_df_PTPN_N2,   hgnc_symbol == "PTPN22")
subset(res_df_PTPN_N2_1, hgnc_symbol == "PTPN22")
subset(res_df_PTPN_N2_2, hgnc_symbol == "PTPN22")
subset(res_df_PTPN_N2_3, hgnc_symbol == "PTPN22")

subset(res_df_PTPN_V1,   hgnc_symbol == "PTPN22")
subset(res_df_PTPN_V1_1, hgnc_symbol == "PTPN22")
subset(res_df_PTPN_V1_2, hgnc_symbol == "PTPN22")
subset(res_df_PTPN_V1_3, hgnc_symbol == "PTPN22")

# ==================== RIPK2 ====================

subset(res_df_RIPK_V1,   hgnc_symbol == "RIPK2")
subset(res_df_RIPK_V1_1, hgnc_symbol == "RIPK2")
subset(res_df_RIPK_V1_2, hgnc_symbol == "RIPK2")
subset(res_df_RIPK_V1_3, hgnc_symbol == "RIPK2")

# ==================== NC1 (Control) ====================

subset(res_df_NC1, hgnc_symbol == "NC1")

##################################### 4th batch ###################################
samples <- 1:33
paths <- sprintf("/media/ming/Extra_SSD_4TB1/eCROPseq/20251018_bulk_cell_RNA_seq_4th/salmon_quant/output_quant_sample%d_S%d/quant.sf", samples, samples)

txi_AKAP11_N1 <- tximport(c(paths[1:3],paths[10:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_AKAP11_N1_1 <- tximport(c(paths[1],paths[10:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_AKAP11_N1_2 <- tximport(c(paths[2],paths[10:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_AKAP11_N1_3 <- tximport(c(paths[3],paths[10:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_AKAP11_N2 <- tximport(c(paths[4:6],paths[10:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_AKAP11_N2_1 <- tximport(c(paths[4],paths[10:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_AKAP11_N2_2 <- tximport(c(paths[5],paths[10:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_AKAP11_N2_3 <- tximport(c(paths[6],paths[10:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_AKAP11_V1 <- tximport(c(paths[7:9],paths[10:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_AKAP11_V1_1 <- tximport(c(paths[7],paths[10:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_AKAP11_V1_2 <- tximport(c(paths[8],paths[10:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_AKAP11_V1_3 <- tximport(c(paths[9],paths[10:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_IFNAR2_N1 <- tximport(c(paths[10:12],paths[1:9],paths[19:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_IFNAR2_N1_1 <- tximport(c(paths[10],paths[1:9],paths[19:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_IFNAR2_N1_2 <- tximport(c(paths[11],paths[1:9],paths[19:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_IFNAR2_N1_3 <- tximport(c(paths[12],paths[1:9],paths[19:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")

txi_IFNAR2_V1 <- tximport(c(paths[13:15],paths[1:9],paths[19:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_IFNAR2_V1_1 <- tximport(c(paths[13],paths[1:9],paths[19:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_IFNAR2_V1_2 <- tximport(c(paths[14],paths[1:9],paths[19:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_IFNAR2_V1_3 <- tximport(c(paths[15],paths[1:9],paths[19:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")

txi_IFNAR2_V2 <- tximport(c(paths[16:18],paths[1:9],paths[19:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_IFNAR2_V2_1 <- tximport(c(paths[16],paths[1:9],paths[19:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_IFNAR2_V2_2 <- tximport(c(paths[17],paths[1:9],paths[19:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_IFNAR2_V2_3 <- tximport(c(paths[18],paths[1:9],paths[19:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")

txi_CARS1_V9 <- tximport(c(paths[19:21],paths[1:18],paths[22:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_CARS1_V9_1 <- tximport(c(paths[19],paths[1:18],paths[22:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_CARS1_V9_2 <- tximport(c(paths[20],paths[1:18],paths[22:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_CARS1_V9_3 <- tximport(c(paths[21],paths[1:18],paths[22:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_MFHAS1_V2 <- tximport(c(paths[22:24],paths[1:21],paths[25:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_CEBPB_V1 <- tximport(c(paths[25:27],paths[1:24],paths[31:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_CEBPB_V1_1 <- tximport(c(paths[25],paths[1:24],paths[31:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_CEBPB_V1_2 <- tximport(c(paths[26],paths[1:24],paths[31:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_CEBPB_V1_3 <- tximport(c(paths[27],paths[1:24],paths[31:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_CEBPB_N1 <- tximport(c(paths[28:30],paths[1:24],paths[31:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_CEBPB_N1_1 <- tximport(c(paths[28],paths[1:24],paths[31:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_CEBPB_N1_2 <- tximport(c(paths[29],paths[1:24],paths[31:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_CEBPB_N1_3 <- tximport(c(paths[30],paths[1:24],paths[31:33]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
txi_NC1 <- tximport(c(paths[31:33],paths[1:30]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")

samples_order <- as.character(1:33)

sampleTable_IFNAR2_N1_1 <- data.frame(
  samples = c(samples_order[10], samples_order[1:9], samples_order[19:33]),
  condition = c("ABE_IFNAR2_N1-1", rep("Other", 24))
)
rownames(sampleTable_IFNAR2_N1_1) <- sampleTable_IFNAR2_N1_1$samples

sampleTable_IFNAR2_N1_2 <- data.frame(
  samples = c(samples_order[11], samples_order[1:9], samples_order[19:33]),
  condition = c("ABE_IFNAR2_N1-2", rep("Other", 24))
)
rownames(sampleTable_IFNAR2_N1_2) <- sampleTable_IFNAR2_N1_2$samples

sampleTable_IFNAR2_N1_3 <- data.frame(
  samples = c(samples_order[12], samples_order[1:9], samples_order[19:33]),
  condition = c("ABE_IFNAR2_N1-3", rep("Other", 24))
)
rownames(sampleTable_IFNAR2_N1_3) <- sampleTable_IFNAR2_N1_3$samples

sampleTable_IFNAR2_V1_1 <- data.frame(
  samples = c(samples_order[13], samples_order[1:9], samples_order[19:33]),
  condition = c("ABE_IFNAR2_V1-1", rep("Other", 24))
)
rownames(sampleTable_IFNAR2_V1_1) <- sampleTable_IFNAR2_V1_1$samples

sampleTable_IFNAR2_V1_2 <- data.frame(
  samples = c(samples_order[14], samples_order[1:9], samples_order[19:33]),
  condition = c("ABE_IFNAR2_V1-2", rep("Other", 24))
)
rownames(sampleTable_IFNAR2_V1_2) <- sampleTable_IFNAR2_V1_2$samples

sampleTable_IFNAR2_V1_3 <- data.frame(
  samples = c(samples_order[15], samples_order[1:9], samples_order[19:33]),
  condition = c("ABE_IFNAR2_V1-3", rep("Other", 24))
)
rownames(sampleTable_IFNAR2_V1_3) <- sampleTable_IFNAR2_V1_3$samples

sampleTable_IFNAR2_V2_1 <- data.frame(
  samples = c(samples_order[16], samples_order[1:9], samples_order[19:33]),
  condition = c("ABE_IFNAR2_V2-1", rep("Other", 24))
)
rownames(sampleTable_IFNAR2_V2_1) <- sampleTable_IFNAR2_V2_1$samples

sampleTable_IFNAR2_V2_2 <- data.frame(
  samples = c(samples_order[17], samples_order[1:9], samples_order[19:33]),
  condition = c("ABE_IFNAR2_V2-2", rep("Other", 24))
)
rownames(sampleTable_IFNAR2_V2_2) <- sampleTable_IFNAR2_V2_2$samples

sampleTable_IFNAR2_V2_3 <- data.frame(
  samples = c(samples_order[18], samples_order[1:9], samples_order[19:33]),
  condition = c("ABE_IFNAR2_V2-3", rep("Other", 24))
)
rownames(sampleTable_IFNAR2_V2_3) <- sampleTable_IFNAR2_V2_3$samples

sampleTable_AKAP11_N1 <- data.frame(
  samples = c(samples_order[1:3],samples_order[10:33]),
  condition = c(rep("AKAP11_N1", 3), rep("Other", 24))
)
rownames(sampleTable_AKAP11_N1) <- sampleTable_AKAP11_N1$samples

sampleTable_AKAP11_N2 <- data.frame(
  samples = c(samples_order[4:6],samples_order[10:33]),
  condition = c(rep("AKAP11_N2", 3), rep("Other", 24))
)
rownames(sampleTable_AKAP11_N2) <- sampleTable_AKAP11_N2$samples

sampleTable_AKAP11_V1 <- data.frame(
  samples = c(samples_order[7:9],samples_order[10:33]),
  condition = c(rep("AKAP11_V1", 3), rep("Other", 24))
)
rownames(sampleTable_AKAP11_V1) <- sampleTable_AKAP11_V1$samples

sampleTable_IFNAR2_N1 <- data.frame(
  samples = c(samples_order[10:12],samples_order[1:9],samples_order[19:33]),
  condition = c(rep("IFNAR2_N1", 3), rep("Other", 24))
)
rownames(sampleTable_IFNAR2_N1) <- sampleTable_IFNAR2_N1$samples

sampleTable_IFNAR2_V1 <- data.frame(
  samples = c(samples_order[13:15],samples_order[1:9],samples_order[19:33]),
  condition = c(rep("IFNAR2_V1", 3), rep("Other", 24))
)
rownames(sampleTable_IFNAR2_V1) <- sampleTable_IFNAR2_V1$samples

sampleTable_IFNAR2_V2 <- data.frame(
  samples = c(samples_order[16:18],samples_order[1:9],samples_order[19:33]),
  condition = c(rep("IFNAR2_V2", 3), rep("Other", 24))
)
rownames(sampleTable_IFNAR2_V2) <- sampleTable_IFNAR2_V2$samples

sampleTable_CARS1_V9 <- data.frame(
  samples = c(samples_order[19:21],samples_order[1:18],samples_order[22:33]),
  condition = c(rep("CARS1_V9", 3), rep("Other", 30))
)
rownames(sampleTable_CARS1_V9) <- sampleTable_CARS1_V9$samples

sampleTable_MFHAS1_V2 <- data.frame(
  samples = c(samples_order[22:24],samples_order[1:21],samples_order[25:33]),
  condition = c(rep("MFHAS1_V2", 3), rep("Other", 30))
)
rownames(sampleTable_MFHAS1_V2) <- sampleTable_MFHAS1_V2$samples

sampleTable_CEBPB_V1 <- data.frame(
  samples = c(samples_order[25:27],samples_order[1:24],samples_order[31:33]),
  condition = c(rep("CEBPB_V1", 3), rep("Other", 27))
)
rownames(sampleTable_CEBPB_V1) <- sampleTable_CEBPB_V1$samples

## demo testing
sampleTable_CEBPB_V1_1 <- data.frame(
  samples = c(samples_order[25],samples_order[1:24],samples_order[31:33]),
  condition = c("CEBPB_V1-1", rep("Other", 27))
)
rownames(sampleTable_CEBPB_V1_1) <- sampleTable_CEBPB_V1_1$samples

sampleTable_CEBPB_V1_2 <- data.frame(
  samples = c(samples_order[26],samples_order[1:24],samples_order[31:33]),
  condition = c("CEBPB_V1-2", rep("Other", 27))
)
rownames(sampleTable_CEBPB_V1_2) <- sampleTable_CEBPB_V1_2$samples

sampleTable_CEBPB_V1_3 <- data.frame(
  samples = c(samples_order[27],samples_order[1:24],samples_order[31:33]),
  condition = c("CEBPB_V1-3", rep("Other", 27))
)
rownames(sampleTable_CEBPB_V1_3) <- sampleTable_CEBPB_V1_3$samples

sampleTable_CEBPB_N1 <- data.frame(
  samples = c(samples_order[28:30],samples_order[1:24],samples_order[31:33]),
  condition = c(rep("CEBPB_N1", 3), rep("Other", 27))
)
rownames(sampleTable_CEBPB_N1) <- sampleTable_CEBPB_N1$samples

sampleTable_NC1 <- data.frame(
  samples = c(samples_order[31:33],samples_order[1:30]),
  condition = c(rep("NC1", 3), rep("Other", 30))
)
rownames(sampleTable_NC1) <- sampleTable_NC1$samples

sampleTable_AKAP11_N1_1 <- data.frame(
  samples = c(samples_order[1], samples_order[10:33]),
  condition = c("ABE_AKAP11_N1-1", rep("Other", 24))
)
rownames(sampleTable_AKAP11_N1_1) <- sampleTable_AKAP11_N1_1$samples

sampleTable_AKAP11_N1_2 <- data.frame(
  samples = c(samples_order[2], samples_order[10:33]),
  condition = c("ABE_AKAP11_N1-2", rep("Other", 24))
)
rownames(sampleTable_AKAP11_N1_2) <- sampleTable_AKAP11_N1_2$samples

sampleTable_AKAP11_N1_3 <- data.frame(
  samples = c(samples_order[3], samples_order[10:33]),
  condition = c("ABE_AKAP11_N1-3", rep("Other", 24))
)
rownames(sampleTable_AKAP11_N1_3) <- sampleTable_AKAP11_N1_3$samples

sampleTable_AKAP11_N2_1 <- data.frame(
  samples = c(samples_order[4], samples_order[10:33]),
  condition = c("ABE_AKAP11_N2-1", rep("Other", 24))
)
rownames(sampleTable_AKAP11_N2_1) <- sampleTable_AKAP11_N2_1$samples

sampleTable_AKAP11_N2_2 <- data.frame(
  samples = c(samples_order[5], samples_order[10:33]),
  condition = c("ABE_AKAP11_N2-2", rep("Other", 24))
)
rownames(sampleTable_AKAP11_N2_2) <- sampleTable_AKAP11_N2_2$samples

sampleTable_AKAP11_N2_3 <- data.frame(
  samples = c(samples_order[6], samples_order[10:33]),
  condition = c("ABE_AKAP11_N2-3", rep("Other", 24))
)
rownames(sampleTable_AKAP11_N2_3) <- sampleTable_AKAP11_N2_3$samples

sampleTable_AKAP11_V1_1 <- data.frame(
  samples = c(samples_order[7], samples_order[10:33]),
  condition = c("ABE_AKAP11_V1-1", rep("Other", 24))
)
rownames(sampleTable_AKAP11_V1_1) <- sampleTable_AKAP11_V1_1$samples

sampleTable_AKAP11_V1_2 <- data.frame(
  samples = c(samples_order[8], samples_order[10:33]),
  condition = c("ABE_AKAP11_V1-2", rep("Other", 24))
)
rownames(sampleTable_AKAP11_V1_2) <- sampleTable_AKAP11_V1_2$samples

sampleTable_AKAP11_V1_3 <- data.frame(
  samples = c(samples_order[9], samples_order[10:33]),
  condition = c("ABE_AKAP11_V1-3", rep("Other", 24))
)
rownames(sampleTable_AKAP11_V1_3) <- sampleTable_AKAP11_V1_3$samples

sampleTable_CARS1_V9_1 <- data.frame(
  samples = c(samples_order[19], samples_order[1:18], samples_order[22:33]),
  condition = c("CBE_CARS1_V9-1", rep("Other", 30))
)
rownames(sampleTable_CARS1_V9_1) <- sampleTable_CARS1_V9_1$samples

sampleTable_CARS1_V9_2 <- data.frame(
  samples = c(samples_order[20], samples_order[1:18], samples_order[22:33]),
  condition = c("CBE_CARS1_V9-2", rep("Other", 30))
)
rownames(sampleTable_CARS1_V9_2) <- sampleTable_CARS1_V9_2$samples

sampleTable_CARS1_V9_3 <- data.frame(
  samples = c(samples_order[21], samples_order[1:18], samples_order[22:33]),
  condition = c("CBE_CARS1_V9-3", rep("Other", 30))
)
rownames(sampleTable_CARS1_V9_3) <- sampleTable_CARS1_V9_3$samples

sampleTable_CEBPB_N1_1 <- data.frame(
  samples = c(samples_order[28], samples_order[1:24], samples_order[31:33]),
  condition = c("CBE_CEBPB_N1-1", rep("Other", 27))
)
rownames(sampleTable_CEBPB_N1_1) <- sampleTable_CEBPB_N1_1$samples

sampleTable_CEBPB_N1_2 <- data.frame(
  samples = c(samples_order[29], samples_order[1:24], samples_order[31:33]),
  condition = c("CBE_CEBPB_N1-2", rep("Other", 27))
)
rownames(sampleTable_CEBPB_N1_2) <- sampleTable_CEBPB_N1_2$samples

sampleTable_CEBPB_N1_3 <- data.frame(
  samples = c(samples_order[30], samples_order[1:24], samples_order[31:33]),
  condition = c("CBE_CEBPB_N1-3", rep("Other", 27))
)
rownames(sampleTable_CEBPB_N1_3) <- sampleTable_CEBPB_N1_3$samples

txi_names <- list(
  txi_AKAP11_N1,
  txi_AKAP11_N1_1,
  txi_AKAP11_N1_2,
  txi_AKAP11_N1_3,
  txi_AKAP11_N2,
  txi_AKAP11_N2_1,
  txi_AKAP11_N2_2,
  txi_AKAP11_N2_3,
  txi_AKAP11_V1,
  txi_AKAP11_V1_1,
  txi_AKAP11_V1_2,
  txi_AKAP11_V1_3,
  txi_IFNAR2_N1,
  txi_IFNAR2_N1_1,
  txi_IFNAR2_N1_2,
  txi_IFNAR2_N1_3,
  txi_IFNAR2_V1,
  txi_IFNAR2_V1_1,
  txi_IFNAR2_V1_2,
  txi_IFNAR2_V1_3,
  txi_IFNAR2_V2,
  txi_IFNAR2_V2_1,
  txi_IFNAR2_V2_2,
  txi_IFNAR2_V2_3,
  txi_CARS1_V9,
  txi_CARS1_V9_1,
  txi_CARS1_V9_2,
  txi_CARS1_V9_3,
  txi_MFHAS1_V2,
  txi_CEBPB_V1,
  txi_CEBPB_V1_1,
  txi_CEBPB_V1_2,
  txi_CEBPB_V1_3,
  txi_CEBPB_N1,
  txi_CEBPB_N1_1,
  txi_CEBPB_N1_2,
  txi_CEBPB_N1_3,
  txi_NC1
)
res_names <- c(
  "res_AKAP11_N1",
  "res_AKAP11_N1_1",
  "res_AKAP11_N1_2",
  "res_AKAP11_N1_3",
  "res_AKAP11_N2",
  "res_AKAP11_N2_1",
  "res_AKAP11_N2_2",
  "res_AKAP11_N2_3",
  "res_AKAP11_V1",
  "res_AKAP11_V1_1",
  "res_AKAP11_V1_2",
  "res_AKAP11_V1_3",
  "res_IFNAR2_N1",
  "res_IFNAR2_N1_1",
  "res_IFNAR2_N1_2",
  "res_IFNAR2_N1_3",
  "res_IFNAR2_V1",
  "res_IFNAR2_V1_1",
  "res_IFNAR2_V1_2",
  "res_IFNAR2_V1_3",
  "res_IFNAR2_V2",
  "res_IFNAR2_V2_1",
  "res_IFNAR2_V2_2",
  "res_IFNAR2_V2_3",
  "res_CARS1_V9",
  "res_CARS1_V9_1",
  "res_CARS1_V9_2",
  "res_CARS1_V9_3",
  "res_MFHAS1_V2",
  "res_CEBPB_V1",
  "res_CEBPB_V1_1",
  "res_CEBPB_V1_2",
  "res_CEBPB_V1_3",
  "res_CEBPB_N1",
  "res_CEBPB_N1_1",
  "res_CEBPB_N1_2",
  "res_CEBPB_N1_3",
  "res_NC1"
)
dds_names <- c(
  "dds_AKAP11_N1",
  "dds_AKAP11_N1_1",
  "dds_AKAP11_N1_2",
  "dds_AKAP11_N1_3",
  "dds_AKAP11_N2",
  "dds_AKAP11_N2_1",
  "dds_AKAP11_N2_2",
  "dds_AKAP11_N2_3",
  "dds_AKAP11_V1",
  "dds_AKAP11_V1_1",
  "dds_AKAP11_V1_2",
  "dds_AKAP11_V1_3",
  "dds_IFNAR2_N1",
  "dds_IFNAR2_N1_1",
  "dds_IFNAR2_N1_2",
  "dds_IFNAR2_N1_3",
  "dds_IFNAR2_V1",
  "dds_IFNAR2_V1_1",
  "dds_IFNAR2_V1_2",
  "dds_IFNAR2_V1_3",
  "dds_IFNAR2_V2",
  "dds_IFNAR2_V2_1",
  "dds_IFNAR2_V2_2",
  "dds_IFNAR2_V2_3",
  "dds_CARS1_V9",
  "dds_CARS1_V9_1",
  "dds_CARS1_V9_2",
  "dds_CARS1_V9_3",
  "dds_MFHAS1_V2",
  "dds_CEBPB_V1",
  "dds_CEBPB_V1_1",
  "dds_CEBPB_V1_2",
  "dds_CEBPB_V1_3",
  "dds_CEBPB_N1",
  "dds_CEBPB_N1_1",
  "dds_CEBPB_N1_2",
  "dds_CEBPB_N1_3",
  "dds_NC1"
)
res_df_names <- c(
  "res_df_AKAP11_N1",
  "res_df_AKAP11_N1_1",
  "res_df_AKAP11_N1_2",
  "res_df_AKAP11_N1_3",
  "res_df_AKAP11_N2",
  "res_df_AKAP11_N2_1",
  "res_df_AKAP11_N2_2",
  "res_df_AKAP11_N2_3",
  "res_df_AKAP11_V1",
  "res_df_AKAP11_V1_1",
  "res_df_AKAP11_V1_2",
  "res_df_AKAP11_V1_3",
  "res_df_IFNAR2_N1",
  "res_df_IFNAR2_N1_1",
  "res_df_IFNAR2_N1_2",
  "res_df_IFNAR2_N1_3",
  "res_df_IFNAR2_V1",
  "res_df_IFNAR2_V1_1",
  "res_df_IFNAR2_V1_2",
  "res_df_IFNAR2_V1_3",
  "res_df_IFNAR2_V2",
  "res_df_IFNAR2_V2_1",
  "res_df_IFNAR2_V2_2",
  "res_df_IFNAR2_V2_3",
  "res_df_CARS1_V9",
  "res_df_CARS1_V9_1",
  "res_df_CARS1_V9_2",
  "res_df_CARS1_V9_3",
  "res_df_MFHAS1_V2",
  "res_df_CEBPB_V1",
  "res_df_CEBPB_V1_1",
  "res_df_CEBPB_V1_2",
  "res_df_CEBPB_V1_3",
  "res_df_CEBPB_N1",
  "res_df_CEBPB_N1_1",
  "res_df_CEBPB_N1_2",
  "res_df_CEBPB_N1_3",
  "res_df_NC1"
)
sample_names <- list(
  sampleTable_AKAP11_N1,
  sampleTable_AKAP11_N1_1,
  sampleTable_AKAP11_N1_2,
  sampleTable_AKAP11_N1_3,
  sampleTable_AKAP11_N2,
  sampleTable_AKAP11_N2_1,
  sampleTable_AKAP11_N2_2,
  sampleTable_AKAP11_N2_3,
  sampleTable_AKAP11_V1,
  sampleTable_AKAP11_V1_1,
  sampleTable_AKAP11_V1_2,
  sampleTable_AKAP11_V1_3,
  sampleTable_IFNAR2_N1,
  sampleTable_IFNAR2_N1_1,
  sampleTable_IFNAR2_N1_2,
  sampleTable_IFNAR2_N1_3,
  sampleTable_IFNAR2_V1,
  sampleTable_IFNAR2_V1_1,
  sampleTable_IFNAR2_V1_2,
  sampleTable_IFNAR2_V1_3,
  sampleTable_IFNAR2_V2,
  sampleTable_IFNAR2_V2_1,
  sampleTable_IFNAR2_V2_2,
  sampleTable_IFNAR2_V2_3,
  sampleTable_CARS1_V9,
  sampleTable_CARS1_V9_1,
  sampleTable_CARS1_V9_2,
  sampleTable_CARS1_V9_3,
  sampleTable_MFHAS1_V2,
  sampleTable_CEBPB_V1,
  sampleTable_CEBPB_V1_1,
  sampleTable_CEBPB_V1_2,
  sampleTable_CEBPB_V1_3,
  sampleTable_CEBPB_N1,
  sampleTable_CEBPB_N1_1,
  sampleTable_CEBPB_N1_2,
  sampleTable_CEBPB_N1_3,
  sampleTable_NC1
)

i <- 1
for (x in txi_names) {

  txi <- x
  sample_info <- sample_names[[i]]

  keep <- rowSums(txi$counts >= 10) >= 2
  counts_filtered <- txi$counts[keep, ]

  y <- DGEList(counts = counts_filtered)
  y <- calcNormFactors(y)

  sample_info$condition <- relevel(factor(sample_info$condition), ref = "Other")
  design <- model.matrix(~ condition, data = sample_info)

  v <- voomWithQualityWeights(y, design, plot = FALSE)

  fit <- lmFit(v, design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)

  res <- topTable(fit, coef = 2, number = Inf, sort.by = "none")

  assign(res_names[i], res)
  assign(dds_names[i], v)

  i <- i + 1
  print(i)

}

res_vars <- list(
  res_AKAP11_N1,
  res_AKAP11_N1_1,
  res_AKAP11_N1_2,
  res_AKAP11_N1_3,
  res_AKAP11_N2,
  res_AKAP11_N2_1,
  res_AKAP11_N2_2,
  res_AKAP11_N2_3,
  res_AKAP11_V1,
  res_AKAP11_V1_1,
  res_AKAP11_V1_2,
  res_AKAP11_V1_3,
  res_IFNAR2_N1,
  res_IFNAR2_N1_1,
  res_IFNAR2_N1_2,
  res_IFNAR2_N1_3,
  res_IFNAR2_V1,
  res_IFNAR2_V1_1,
  res_IFNAR2_V1_2,
  res_IFNAR2_V1_3,
  res_IFNAR2_V2,
  res_IFNAR2_V2_1,
  res_IFNAR2_V2_2,
  res_IFNAR2_V2_3,
  res_CARS1_V9,
  res_CARS1_V9_1,
  res_CARS1_V9_2,
  res_CARS1_V9_3,
  res_MFHAS1_V2,
  res_CEBPB_V1,
  res_CEBPB_V1_1,
  res_CEBPB_V1_2,
  res_CEBPB_V1_3,
  res_CEBPB_N1,
  res_CEBPB_N1_1,
  res_CEBPB_N1_2,
  res_CEBPB_N1_3,
  res_NC1
)

dds_vars <- list(
  dds_AKAP11_N1,
  dds_AKAP11_N1_1,
  dds_AKAP11_N1_2,
  dds_AKAP11_N1_3,
  dds_AKAP11_N2,
  dds_AKAP11_N2_1,
  dds_AKAP11_N2_2,
  dds_AKAP11_N2_3,
  dds_AKAP11_V1,
  dds_AKAP11_V1_1,
  dds_AKAP11_V1_2,
  dds_AKAP11_V1_3,
  dds_IFNAR2_N1,
  dds_IFNAR2_N1_1,
  dds_IFNAR2_N1_2,
  dds_IFNAR2_N1_3,
  dds_IFNAR2_V1,
  dds_IFNAR2_V1_1,
  dds_IFNAR2_V1_2,
  dds_IFNAR2_V1_3,
  dds_IFNAR2_V2,
  dds_IFNAR2_V2_1,
  dds_IFNAR2_V2_2,
  dds_IFNAR2_V2_3,
  dds_CARS1_V9,
  dds_CARS1_V9_1,
  dds_CARS1_V9_2,
  dds_CARS1_V9_3,
  dds_MFHAS1_V2,
  dds_CEBPB_V1,
  dds_CEBPB_V1_1,
  dds_CEBPB_V1_2,
  dds_CEBPB_V1_3,
  dds_CEBPB_N1,
  dds_CEBPB_N1_1,
  dds_CEBPB_N1_2,
  dds_CEBPB_N1_3,
  dds_NC1
)

## output tables with gene name
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

i = 1
for (x in res_vars) {
  deg_genes <- rownames(x)

  gene_conversion <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                           filters = "ensembl_gene_id",
                           values = deg_genes,
                           mart = ensembl)

  res_df <- as.data.frame(x)
  res_df$ensembl_gene_id <- rownames(res_df)
  res_df <- merge(res_df, gene_conversion, by = "ensembl_gene_id", all.x = TRUE)

  assign(res_df_names[i], res_df)

  i = i + 1
}

# =================== AKAP11 ===================
subset(res_df_AKAP11_N1, hgnc_symbol == "AKAP11")
subset(res_df_AKAP11_N1_1, hgnc_symbol == "AKAP11")
subset(res_df_AKAP11_N1_2, hgnc_symbol == "AKAP11")
subset(res_df_AKAP11_N1_3, hgnc_symbol == "AKAP11")

subset(res_df_AKAP11_N2, hgnc_symbol == "AKAP11")
subset(res_df_AKAP11_N2_1, hgnc_symbol == "AKAP11")
subset(res_df_AKAP11_N2_2, hgnc_symbol == "AKAP11")
subset(res_df_AKAP11_N2_3, hgnc_symbol == "AKAP11")

subset(res_df_AKAP11_V1, hgnc_symbol == "AKAP11")
subset(res_df_AKAP11_V1_1, hgnc_symbol == "AKAP11")
subset(res_df_AKAP11_V1_2, hgnc_symbol == "AKAP11")
subset(res_df_AKAP11_V1_3, hgnc_symbol == "AKAP11")

# =================== IFNAR2 ===================
# ==================== IFNAR2 ====================

subset(res_df_IFNAR2_N1,   hgnc_symbol == "IFNAR2")
subset(res_df_IFNAR2_N1_1, hgnc_symbol == "IFNAR2")
subset(res_df_IFNAR2_N1_2, hgnc_symbol == "IFNAR2")
subset(res_df_IFNAR2_N1_3, hgnc_symbol == "IFNAR2")

subset(res_df_IFNAR2_V1,   hgnc_symbol == "IFNAR2")
subset(res_df_IFNAR2_V1_1, hgnc_symbol == "IFNAR2")
subset(res_df_IFNAR2_V1_2, hgnc_symbol == "IFNAR2")
subset(res_df_IFNAR2_V1_3, hgnc_symbol == "IFNAR2")

subset(res_df_IFNAR2_V2,   hgnc_symbol == "IFNAR2")
subset(res_df_IFNAR2_V2_1, hgnc_symbol == "IFNAR2")
subset(res_df_IFNAR2_V2_2, hgnc_symbol == "IFNAR2")
subset(res_df_IFNAR2_V2_3, hgnc_symbol == "IFNAR2")

# =================== CARS1 ===================
subset(res_df_CARS1_V9, hgnc_symbol == "CARS1")
subset(res_df_CARS1_V9_1, hgnc_symbol == "CARS1")
subset(res_df_CARS1_V9_2, hgnc_symbol == "CARS1")
subset(res_df_CARS1_V9_3, hgnc_symbol == "CARS1")

# =================== MFHAS1 ===================
subset(res_df_MFHAS1_V2, hgnc_symbol == "MFHAS1")

# =================== CEBPB ===================
subset(res_df_CEBPB_V1, hgnc_symbol == "CEBPB")
subset(res_df_CEBPB_V1_1, hgnc_symbol == "CEBPB")
subset(res_df_CEBPB_V1_2, hgnc_symbol == "CEBPB")
subset(res_df_CEBPB_V1_3, hgnc_symbol == "CEBPB")

subset(res_df_CEBPB_N1, hgnc_symbol == "CEBPB")
subset(res_df_CEBPB_N1_1, hgnc_symbol == "CEBPB")
subset(res_df_CEBPB_N1_2, hgnc_symbol == "CEBPB")
subset(res_df_CEBPB_N1_3, hgnc_symbol == "CEBPB")

#### Back-up code: manual log2fc and individual replicate

compare_CEBPB_V1 <- compare_limma_vs_manual(
  res_df_CEBPB_V1,
  stats_res_CEBPB_V1,
  contrast_name = "CEBPB_V1"
)

subset(compare_CEBPB_V1, hgnc_symbol == "CEBPB")

compare_limma_vs_manual <- function(res_df, stats_df, contrast_name = NULL) {
  library(dplyr)
  library(tidyr)

  merged <- res_df %>%
    select(ensembl_gene_id, hgnc_symbol, logFC, AveExpr, P.Value, adj.P.Val) %>%
    left_join(
      stats_df %>%
        select(gene, hgnc_symbol, condition, mean_log2, sd_log2),
      by = c("ensembl_gene_id" = "gene", "hgnc_symbol" = "hgnc_symbol")
    )

  merged_summary <- merged %>%
    pivot_wider(
      id_cols = c(ensembl_gene_id, hgnc_symbol, logFC, AveExpr, P.Value, adj.P.Val),
      names_from = condition,
      values_from = c(mean_log2, sd_log2)
    )

  cond_cols <- colnames(merged_summary)
  conds <- gsub("^mean_log2_", "", cond_cols[grepl("^mean_log2_", cond_cols)])
  conds <- conds[!conds %in% c("Other", "Control")]
  exp_cond <- ifelse(length(conds) > 0, conds[1], "Exp")

  merged_summary <- merged_summary %>%
    mutate(
      manual_log2FC = get(paste0("mean_log2_", exp_cond)) - mean_log2_Other,
      diff_vs_model = manual_log2FC - logFC
    )

  merged_summary <- merged_summary %>%
    relocate(
      hgnc_symbol, ensembl_gene_id, logFC, manual_log2FC, diff_vs_model,
      mean_log2_Other, sd_log2_Other,
      !!sym(paste0("mean_log2_", exp_cond)),
      !!sym(paste0("sd_log2_", exp_cond)),
      AveExpr, P.Value, adj.P.Val
    )

  if (!is.null(contrast_name)) {
    message("✅ Comparison summary for: ", contrast_name)
  }

  return(merged_summary)
}

i <- 1
for (x in txi_names) {

  txi <- x
  sample_info <- sample_names[[i]]

  keep <- rowSums(txi$counts >= 10) >= 2
  counts_filtered <- txi$counts[keep, ]

  y <- DGEList(counts = counts_filtered)
  y <- calcNormFactors(y)

  sample_info$condition <- relevel(factor(sample_info$condition), ref = "Other")
  design <- model.matrix(~ condition, data = sample_info)

  v <- voomWithQualityWeights(y, design, plot = FALSE)

  fit <- lmFit(v, design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)

  res <- topTable(fit, coef = 2, number = Inf, sort.by = "none")

  ############### Individual replicates ##############

  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

  expr_mat <- v$E
  expr_long <- melt(expr_mat)
  colnames(expr_long) <- c("gene", "sample", "log2expr")
  expr_long$sample <- gsub("Sample", "", expr_long$sample)
  expr_long <- merge(expr_long, sample_info, by.x = "sample", by.y = "samples")

  stats_by_cond <- expr_long %>%
    group_by(gene, condition) %>%
    summarise(
      mean_log2 = mean(log2expr, na.rm = TRUE),
      sd_log2   = sd(log2expr, na.rm = TRUE),
      .groups = "drop"
    )

  gene_conversion <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = unique(stats_by_cond$gene),
    mart = ensembl
  )

  stats_by_cond <- stats_by_cond %>%
    left_join(gene_conversion, by = c("gene" = "ensembl_gene_id")) %>%
    relocate(hgnc_symbol, .after = gene)

  assign(paste0("stats_", res_names[i]), stats_by_cond)

  assign(res_names[i], res)
  assign(dds_names[i], v)

  ############### Individual replicates ##############

  i <- i + 1
  print(i)

}

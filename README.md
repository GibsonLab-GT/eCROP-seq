# eCROP-seq

Expression CROPseq is a process that uses CRISPR/cas9 to edit cells from a pool of gRNA and measure the impact on gene expression from scRNAseq. Guide RNA, corresponding to a unique SNP knockout, are detected as transcripts enabling us to differentiate which cell got which perturbation, and thus, which perturbation causes an effect. Published pilot experiments ([DOI: 10.1093/biomethods/bpaa008](https://pubmed.ncbi.nlm.nih.gov/32665975/)) show this working well for two genes (*PARK7* and *CISD1*). Here, we use screens with the function of rare and common variants (manuscript in submission). Uploaded files to perform hypothesis testing of eCROP-seq experimental data in R and Phyton and are described below. Analyses report on work that was conducted by Emily Greenwood (Georgia Tech) and Mingming Cao (Rice University). A preprint of this work is available at https://doi.org/10.1101/2025.01.30.635675.

------------------
**Pre-Processing**
------------------

**Pre-processing Common Variant:**

Genes associated with Inflammatory Bowel Disease (IBD) were taken from https://www.opentargets.org/.
Expression QTL for these genes were isolated via all-but-one conditional analysis by Maggie Brown (Georgia Tech) and described in previous publication. Please see https://github.com/GibsonLab-GT/All-but-One conditional analysis and additional publication [DOI: 10.1093/genetics/iyac162](https://pubmed.ncbi.nlm.nih.gov/36321965/). 

Variants in eQTL were then selected based on significance, linkage disequilibrium with lead variants, and overlap with ATACseq peaks. Then they were brought forth to gRNA design and were targeted by one unique gRNA.

**Pre-processing Rare Variant:**

Rare variants that were predicted to be functional by Watershed (Ferraro et al. DOI: 10.1126/science.aaz5900), indicated by a posterior greater than 0.9, for genes expressed in HL-60 cell line, were selected for eCROPseq testing. Controls with a posterior less than 0.01 were also selected. Rare variants were targeted by two unique gRNA and controls by one.

**gRNA Design:**

This code has been modified from original work (https://github.com/pdpdpd/eCROP-kit) created by author's (Yidan Pan, and in collaboration with Ruoyu Tian and Ciaran Lee) in initial eCROP-seq analysis ([DOI: 10.1093/biomethods/bpaa008](https://pubmed.ncbi.nlm.nih.gov/32665975/)). Modifications were made to fit our study design. 

1) gRNA_assessment_hg19.py and gRNA_assessment_hg38.py

File takes in a data frame titled SNP.csv with columns ‘RSID Chr’, ‘RSID BP’, and ‘Focal_RSID’ and outputs gRNA with a cut site within 10bp and GC content between 20-80. Output file with all information,  is labeled ‘gRNA_to_process.txt’ and will contain the columns: Focal_RSID, gRNA sequence, SNP-gRNA distance, Strand. Another file label 'no_gRNA_SNP.csv' will contain SNP names without a gRNA fitting specified criteria. 
If positions in SNP.csv file are in GrCH37/Hg19 human reference genome format, use gRNA_assessment_hg19.py. If positions are in GrCH38/Hg38 human reference genome format, use gRNA_assessment_hg38.py. Run script with Python v 2.7.18 with the following command: 
  
  	python gRNA_assessment_hg19.py -i SNP.csv. 

**Supplementing reference genome**

Fasta and GTF so for each gRNA sequence as a chromosome were created following protocol laid out in Datlinger et al. [DOI: https://doi.org/10.1038/nmeth.4177](https://www.nature.com/articles/nmeth.4177). To create reference genomes, we followed https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_mr by first copying the Human reference genome's (https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz) fasta and GTF files into a new folder. Then concatenated the gRNA fasta and GTF file to the Hg38 fasta and GTF file, respectively. 
For samples processed with 10x genomics, we then ran cell ranger mkref to create the reference. For samples processed with Fluent Biosciences, we ran Star v.2.7.11a with the following command:
	
 	./STAR --runMode genomeGenerate --genomeDir ~/path/for/file --runThreadN 4 --genomeFastaFiles ~/path/to/fasta --sjdbGTFfile ~/path/to/gtf --sjdbOverhang 69 --runThreadN 7 --genomeSAsparseD 3.

**Aligning Reads to Supplemented Reference Genome**

With 10x genomics samples, reads were aligned to supplemented reference genome with the following command: 

	cellranger count --id=NAME --fastqs=/path/to/fastqs --transcriptome=/path/to/supplemented/genome --expect-cells=10000.

With Fluent Biosciences samples, reads were aligned to supplemented reference genome with the following command: 
	
 	PIPseeker full --chemistry v4 --fastq /path/to/fastq/. --star-index-path /path/to/gtf --output-path /output/path.

**Quality Control in Seurat**

Standard quality control was followed as described in https://satijalab.org/seurat/articles/pbmc3k_tutorial. Feature and mitochondrial counts cutoffsfor all pools were selected within 2.5 median absolute deviations (MADS) surrounding the median for feature counr and 35% mito count. This was down through the R package scuttle (https://rdrr.io/bioc/scuttle/man/isOutlier.html) with the commands below.

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

 
----------------------
**Hypothesis Testing**
----------------------

**Hypothesis Testing Common Variant:**

2) Initial_hypothesis_tests_common_variant.R
3) Undif_hypothesis_tests_common_variant.R
4) Macro_hypothesis_tests_common_variant.R
5) Neutro_hypothesis_tests_common_variant.R

These scripts take in a path to barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz for a scRNAseq aligned pool. This script also takes in a SNP list with two important columns: Gene, rs.ID 

**Hypothesis Testing Rare Variant:**

4) Rare_variant_hypothesis_tests.R
This script is slightly modified from above due to rare variants being targeted by two gRNA, named SNP_1 for cells with gRNA 1 and SNP_2 for cells with gRNA 2. This script takes in a SNP list with columns: rs.ID,grna_target,Gene,grna_number as well as a path to barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz for a scRNAseq aligned pool. Note this script performs 3 different groupings for hypothesis tests:

  A) cells with gRNA 1 to those without gRNA 1 or gRNA
  
  B) cells with gRNA 2 to those without gRNA 2 or gRNA2
  
  C) cells with either gRNA 1 or 2 to those without either 
  
Tests performed by this script include a Welch’s t-test, KS test, and SCEPTRE. As well as all of those comparisons for only cells expressing target gene and repeats. SCEPTRE and Welch results were the results reported in the manuscript. 






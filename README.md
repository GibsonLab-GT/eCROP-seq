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


**gRNA Naming**

The hypothesis-testing scripts identify gRNA features by name, so the supplemented reference must follow these conventions:

Common variant gRNA: <rs.ID>-gene
Rare variant gRNA 1: <rs.ID>-1-gene
Rare variant gRNA 2: <rs.ID>-2-gene
Non-targeting controls name matching: non-targeting, non-targeting-control, non-targeting-ctrl or no-control
Previously validated controls: rs2251039-ctrl (CISD1), rs35675666-ctrl (PARK7)
Survival positive controls name matchin:  positive-ctrl (RUNX1, TUBB)


**Quality Control in Seurat**

Applied identically by every hypothesis-testing script, in this order:

1. Object creation: CreateSeuratObject(min.cells = 6, min.features = 200). Pools split across multiple sequencing reactions are merged, one Sample label per reaction.
2. Doublet and high-mitochondrial removal: scDblFinder is run per Sample on cells with percent.mt < 35. Cells above that threshold are therefore absent from the doublet calls and are flagged is_doublet = TRUE, so the single step removes both predicted doublets and cells with > 35% mitochondrial transcripts.
3. nFeature outliers: per Sample, cells outside ±2.5 median absolute deviations of the median nFeature_RNA are removed, on the log scale, using scuttle::isOutlier. No further mitochondrial filter is applied; percent.mt is retained as a SCEPTRE covariate.
4. gRNA requirement: cells are retained only if they carry at least one non-positive-control gRNA. Cells carrying a survival positive-control gRNA (RUNX1, TUBB) are removed.
5. Normalization:LogNormalize, scale factor 10,000, followed by FindVariableFeatures (vst, 2,000 features) and ScaleData. All t-tests and KS tests are computed on the log-normalized layer; SCEPTRE receives the raw count matrix.
   
Example code:

		qc_per_sample <- function(seurat_obj, sample_id, nmads = 2.5) {
		  cells <- colnames(seurat_obj)[seurat_obj$Sample == sample_id]
		  obj   <- seurat_obj[, cells]
		
		  qc_lo <- isOutlier(obj$nFeature_RNA, log = TRUE, type = "lower", nmads = nmads)
		
		  qc_hi <- isOutlier(obj$nFeature_RNA, log = TRUE, type = "higher", nmads = nmads)
		
		  message(
		    sprintf(
		      "Sample %s: nFeature > %.1f, < %.1f",
		      sample_id,
		      attr(qc_lo, "thresholds")["lower"],
		      attr(qc_hi, "thresholds")["higher"]
		    )
		  )
		
		  cells[
		    obj$nFeature_RNA > attr(qc_lo, "thresholds")["lower"] &
		    obj$nFeature_RNA < attr(qc_hi, "thresholds")["higher"]
		  ]
		}


**Cell Type Selection**

After normalization, the macrophage and neutrophil scripts retain cells expressing either marker at or above the pool mean:

-Macrophage: CD14 or CD68
-Neutrophil: CSF3R or S100A8

The initial and undifferentiated scripts apply no such filter.

----------------------
**Hypothesis Testing**
----------------------

**Hypothesis Testing Common Variant:**

*Scripts:**
- Initial_hypothesis_tests_common_variant.R - initial screen 
- Undif_hypothesis_tests_common_variant.R - undifferentiated retest
- Macro_hypothesis_tests_common_variant.R - macrophage retest
- Neutro_hypothesis_tests_common_variant.R - neutrophil retest

*Inputs:**
- For each aligned scRNA-seq pool:
- barcodes.tsv.gz
- features.tsv.gz
- matrix.mtx.gz

Variant annotation file (snp_list.csv) containing:
- Gene
- rs.ID

*Comparison:**
For each variant-gene pair, cells carrying the focal gRNA are compared with a control group of all remaining cells except those carrying any other gRNA targeting the same gene. 


**Hypothesis Testing Rare Variant:**

4) Rare_variant_hypothesis_tests.R

Because rare variants were targeted by two independent gRNAs, cells are classified as:
- SNP_1: cells containing gRNA 1
- SNP_2: cells containing gRNA 2

*Inputs for each aligned scRNA-seq pool:**
- barcodes.tsv.gz
- features.tsv.gz
- matrix.mtx.gz

Variant annotation file (snp_list.csv) containing:
- rs.ID
- grna_target
- Gene
- grna_number

*Comparison:** 
Treatment cells: 
A. gRNA 1 (cells carrying gRNA 1)
B. gRNA 2 (cells carrying gRNA 2)
C. Union (cells carrying gRNA 1 or gRNA 2)
D. Intersection (cells carrying gRNA 1 and gRNA 2)

In all four, the control group is every remaining cell except those carrying any gRNA targeting the focal gene, so the same control pool is used throughout. Control variants carry a single gRNA, so only one of A and B is populated for them and C and D are empty.

**Statistical Tests:**
The following analyses are performed for each comparison:

- Welch's t-test
- Kolmogorov-Smirnov (KS) test
- SCEPTRE with log_nUMI, percent_mt and sequencing batch as covariates

The Welch's t-test and SCEPTRE results at nominal p < 0.05 are those reported in the manuscript. Benjamini-Hochberg adjusted values are included in the output for reference.


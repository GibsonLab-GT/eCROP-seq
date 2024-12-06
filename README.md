# eCROP-seq

Expression CROPseq is a process that uses CRISPR/cas9 to edit cells from a pool of gRNA and measure the impact on gene expression from scRNAseq. Guide RNA, corresponding to a unique SNP knockout, are detected as transcripts enabling us to differentiate which cell got which perturbation, and thus, which perturbation causes an effect. Published pilot experiments (DOI: 10.1093/biomethods/bpaa008) show this working well for two genes (PARK7 and CISD1). Here, we use screens with the function of rare and common variants (manuscript in submission). Uploaded files to perform hypothesis testing of eCROP-seq experimental data in R, Phyton, and Bash and are described below. Analyses report on work that was conducted by Emily Greenwood (Georgia Tech) and Mingming Cao (Rice University).

**Pre-processing Common Variant:**

Genes associated with Inflammatory Bowel Disease (IBD) were taken from https://www.opentargets.org/.
Expression QTL for these genes were isolated via all-but-one conditional analysis by Maggie Brown (Georgia Tech) and described in previous publication. Please see https://github.com/GibsonLab-GT/All-but-One conditional analysis and additional publication DOI: 10.1093/genetics/iyac162. 

Variants in eQTL were then selected based on significance, linkage disequilibrium with lead variants, and overlap with ATACseq peaks. Then they were brought forth to gRNA design and were targeted by one unique gRNA.

**Pre-processing Rare Variant:**

Rare variants that were predicted to be functional by Watershed (Ferraro et al. DOI: 10.1126/science.aaz5900), indicated by a posterior greater than 0.9, for genes expressed in HL-60 cell line, were selected for eCROPseq testing. Controls with a posterior less than 0.01 were also selected. Rare variants were targeted by two unique gRNA and controls by one.

**gRNA Design:**

This code has been modified from original work created by author's (Yidan Pan, Ruoyu Tian and Ciaran Lee) in initial eCROP-seq analysis (DOI: 10.1093/biomethods/bpaa008). Modifications were made to fit our study design. 

1) gRNA_assessment_hg19.py and gRNA_assessment_hg38.py

File takes in a data frame titled SNP.csv with columns ‘RSID Chr’, ‘RSID BP’, and ‘Focal_RSID’ and outputs gRNA with a cut site within 10bp and GC content between 20-80. Output file with all information,  is labeled ‘gRNA_to_process.txt’ and will contain the columns: Focal_RSID, gRNA sequence, SNP-gRNA distance, Strand. Another file label 'no_gRNA_SNP.csv' will contain SNP names without a gRNA fitting specified criteria. 
If positions in SNP.csv file are in GrCH37/Hg19 human reference genome format, use gRNA_assessment_hg19.py. If positions are in GrCH38/Hg38 human reference genome format, use gRNA_assessment_hg38.py. Run script with Python v 2.7.18 with the following command: 
  python gRNA_assessment_hg19.py -i SNP.csv. 

**Supplementing reference genome**

Fasta and GTF so for each gRNA sequence as a chromosome were created following protocol laid out in Datlinger et al. DOI: https://doi.org/10.1038/nmeth.4177. To create reference genomes, we followed https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_mr by first copying the Human reference genome's (https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz) fasta and GTF files into a new folder. Then concatenated the gRNA fasta and GTF file to the Hg38 fasta and GTF file, respectively. 
For samples processed with 10x genomics, we then ran cell ranger mkref to create the reference. For samples processed with Fluent Biosciences, we ran Star v.2.7.11a with the following command:
	./STAR --runMode genomeGenerate --genomeDir ~/path/for/file --runThreadN 4 --genomeFastaFiles ~/path/to/fasta --sjdbGTFfile ~/path/to/gtf --sjdbOverhang 69 --runThreadN 7 --genomeSAsparseD 3.

**Aligning Reads to Supplemented Reference Genome**

With 10x genomics samples, reads were aligned to supplemented reference genome with the following command: cellranger count --id=NAME --fastqs=/path/to/fastqs --transcriptome=/path/to/supplemented/genome --expect-cells=10000.

With Fluent Biosciences samples, reads were aligned to supplemented reference genome with the following command: PIPseeker full --chemistry v4 --fastq /path/to/fastq/. --star-index-path /path/to/gtf --output-path /output/path.

**Quality Control in Seurat**

Standard quality control was followed as described in https://satijalab.org/seurat/articles/pbmc3k_tutorial. The only difference being in how the feature and mitochondrial counts cutoffs were selected in rare variant and common validation analyses. 

**Hypothesis Testing Common Variant:**

2) Target_gene_hypothesis_tests_common_variant_initial.R
This script takes in a file that has a column containing the variant’s name in Seurat Object titled ‘RSID’, the status of each variant (i.e. negative control, to be tested, positive control) titled ‘STATUS’, and the gene targeted titled ‘GENE’. As well as a quality-controlled Seurat object titled ‘data’. Proximal and random gene analyses are conducted with this script, instead read in a file with RSID, status, and proximal or random, respectively, instead of target gene.
Additional information not essential to this script can be included by simply adding ‘New  info col’ when initializing the info data frame and later when adding a row to the info data frame. It is recommended to do so for proximal and random gene analysis, by adding an additional column for the original target gene, even though the expression change of the proximal/random gene is the one being tested. Unless otherwise specified, the generated csv file will have the following information: Variant name, status, target gene, number of cells with gRNA, number of cells with gRNA and expressing target, p-value for cells with gRNA, p-value when subsetting further for cells with gRNA and expressing target, hypothesis test. 

3) Target_gene_hypothesis_tests_common_variant_validation.R
This file takes the same information as 2) (i.e. QC Seurat object and file with RSID, target gene, and status information). Pre-processing before hypothesis testing involves filtering the above files for the gene of interest, specifying each gRNA in a pool, extracting the cells with those gRNA, removing cells without a gRNA, then subsetted into two groups 1) cells without and gRNA targeting gene of interest and 2) cells with gRNA targeting genes of interest. Note the script runs for a specific gene at a time, so for a different target gene, go back to ‘SUBSET’ and repeat. 

**Hypothesis Testing Rare Variant:**

4) Target_gene_hypothesis_tests_rare_variant.R
This script is slightly modified from above due to rare variants being targeted by two gRNA, named SNP_1 for cells with gRNA 1 and SNP_2 for cells with gRNA 2. Example of the file in csv format:
SNP_1, SNP_2, SNP,  STATUS, GENE
chr6:143501293-1, NA, chr6:143501293,  NegativeControl, FUCA2
chr6:143507791-1, chr6:143507791-2, chr6:143507791,  Rare, FUCA2

Thus, this script takes in a file that has a column containing the variant targeted by gRNA 1 in Seurat Object titled ‘SNP_1’, the variant targeted by gRNA 2 titled ‘SNP_2’, the general SNP name targeted with either guide titled ‘SNP’, the status of each variant (i.e. negative control, rare to be tested, positive control) titled ‘STATUS’, and the gene targeted titled ‘GENE’. As well as a quality-controlled Seurat object titled ‘data’. 
As above, proximal and random gene analyses are conducted with this script, instead read in a file with SNP_1, SNP_2, SNP, status, and proximal or random, respectively, instead of target gene.
Additional information not essential to this script can be included by simply adding ‘New  info col’ when initializing the info data frame and later when adding a row to the info data frame. It is recommended to do so for proximal and random gene analysis, by adding an additional column for the original target gene, even though the expression change of the proximal/random gene is the one being tested. Unless otherwise specified, the generated csv file will have the following information: SNP_1, SNP_2, SNP, status, target gene, number of cells with gRNA, number of cells with gRNA and expressing target, p-value for cells with gRNA, p-value when subsetting further for cells with gRNA and expressing target, hypothesis test. Note this script performs 5 different groupings for hypothesis tests:

  A) cells with gRNA 1 to those without gRNA 1
  
  B) cells with gRNA 2 to those without gRNA 2
  
  C) cells with either gRNA 1 or 2 to those without either 
  
  D) cells with gRNA 1 to those without gRNA 1 or 2
  
  E) cells with gRNA 2 to those without gRNA 1 or 2
  
And performs a student’s t-test and KS test. As well as all of those comparisons for only cells expressing target gene and repeats.

**After Output from Hypothesis Testing:**

5) Filter_for_sig_common_variants.R
This script takes a combined data frame of all hypothesis test files generated from 2) 3) or 4). If across multiple files and all in one folder, there is a portion of the script that will combine all files for samples set (e.g. all genes from HL60 samples), then sort based on a p-value less than equal to 0.05 in 5 or more cells. This was the criteria for validation of common variant and rare variant experiments. For initial common variant experiments, the criteria was 0.005 in 5 or more cells. This can be tailored to the precise needs of differing experimental designs. Code filters columns with names of ‘X.p.value.all. , X.p.value.test. , X.cells.w.guide.and.target. , and X.cells.w.guide.’ which will be the output in 2) 3) or 4). 








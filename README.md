# eCROP-seq
Uploaded files to perform hypothesis testing of eCROP-seq experimental data in R are described below. Analyses report on work that was conducted by Emily Greenwood (Georgia Tech) and Mingming Cao (Rice Univeristy).

**Pre-processing Common Variant:**
Genes associated with Inflammatory Bowel Disease (IBD) were taken from https://www.opentargets.org/.
Expression QTL for these genes were isolated via all-but-one conditional analysis by Maggie Brown (Georgia Tech) and described in previous publication. Please see https://github.com/GibsonLab-GT/All-but-One conditional analysis and additional publication DOI: 10.1093/genetics/iyac162. 

Variants in eQTL were then selected based on significance, linkage disequilibrium with lead variant, and overlap with ATACseq peaks. Then they were brought forth to gRNA design and were targeted by one unique gRNA.

**Pre-processing Rare Variant:**
Rare variants that were predicted to be functional by Watershed (DOI: 10.1126/science.aaz5900), indicated by a posterior greater than 0.9, for genes expressed in HL-60 cell line, were selected for eCROPseq testing. Controls with a postieror less than 0.01 were also selected. Rare variants were targeted by two unique gRNA and controls by one.

**gRNA Design:**
These files were originally created from author's (Yidan Pan, Ruoyu Tian and Ciaran Lee) in initial eCROP-seq analysis (DOI: 10.1093/biomethods/bpaa008). We have made a few modifications to fit out study design. 

1) File takes 

**Hypothesis Testing Common Variant:**

1) Target_gene_hypothesis_tests_common_variant_initial
This script takes in a file that has a column containing the variant’s name in Seurat Object titled ‘RSID’, the status of each variant (i.e. negative control, to be tested, positive control) titled ‘STATUS’, and the gene targeted titled ‘GENE’. As well as a quality-controlled Seurat object titled ‘data’. Proximal and random gene analyses are conducted with this script, instead read in a file with RSID, status, and proximal or random, respectively, instead of target gene.
Additional information not essential to this script can be included by simply adding ‘New  info col’ when initializing the info data frame and later when adding a row to the info data frame. It is recommended to do so for proximal and random gene analysis, by adding an additional column for the original target gene, even though the expression change of the proximal/random gene is the one being tested. Unless otherwise specified, the generated csv file will have the following information: Variant name, status, target gene, number of cells with gRNA, number of cells with gRNA and expressing target, p-value for cells with gRNA, p-value when subsetting further for cells with gRNA and expressing target, hypothesis test. 

2) Target_gene_hypothesis_tests_common_variant_validation
This file takes in the same information as 1) (i.e. QC Seurat object and file with RSID, target gene, and status information). Pre-processing before hypothesis testing involves filtering the above files for the gene of interest, specifying each gRNA in a pool, extracted the cells with those gRNA, removing cells without a gRNA, then subsetted into two groups 1) cells without and gRNA targeting gene of interest and 2) cells with gRNA targeting genes of interest. Note the script runs for a specific gene at a time, so for a different target gene, go back to ‘SUBSET’ and repeat. 

**Hypothesis Testing Rare Variant:**

3) Target_gene_hypothesis_tests_rare_variant
This script is slightly modified from above due to rare variants being targeted by two gRNA, named SNP_1 for cells with gRNA 1 and SNP_2 for cells with gRNA 2. Example of the file in csv format:
SNP_1, SNP_2, SNP,  STATUS, GENE
chr6:143501293-1, NA, chr6:143501293,  NegativeControl, FUCA2
chr6:143507791-1, chr6:143507791-2, chr6:143507791,  Rare, FUCA2

Thus, this script takes in a file that has a column containing the variant targeted by gRNA 1 in Seurat Object titled ‘SNP_1’, the variant targeted by gRNA 2 titled ‘SNP_2’, the general SNP name targeted with either guide titled ‘SNP’, the status of each variant (i.e. negative control, rare to be tested, positive control) titled ‘STATUS’, and the gene targeted titled ‘GENE’. As well as a quality-controlled Seurat object titled ‘data’. 
As above, proximal and random gene analyses are conducted with this script, instead read in a file with SNP_1, SNP_2, SNP, status, and proximal or random, respectively, instead of target gene.
Additional information not essential to this script can be included by simply adding ‘New  info col’ when initializing the info data frame and later when adding a row to the info data frame. It is recommended to do so for proximal and random gene analysis, by adding an additional column for the original target gene, even though the expression change of the proximal/random gene is the one being tested. Unless otherwise specified, the generated csv file will have the following information: SNP_1, SNP_2, SNP, status, target gene, number of cells with gRNA, number of cells with gRNA and expressing target, p-value for cells with gRNA, p-value when subsetting further for cells with gRNA and expressing target, hypothesis test. Note this script performs 5 different groupings for hypothesis tests:
1) cells with gRNA 1 to those without gRNA 1
2) cells with gRNA 2 to those without gRNA 1
3) cells with either gRNA 1 or 2 to those without either 
4) cells with gRNA 1 to those without gRNA 1 or 2
5) cells with gRNA 2 to those without gRNA 1 or 2
And performs a student’s t-test and KS test. As well as all of those comparisons for only cells expressing target gene and repeats.

**After Output from Hypothesis Testing:**
4) Filter_for_sig_common_variants.R
This script takes a combined data frame of all hypothesis test files generated from 1) 2) or 3). If across multiple files and all in one folder, there is a portion of the script that will combine all files for samples set (e.g. all genes from HL60 samples), then sort based on a p-value less than equal to 0.05 in 5 or more cells. This was the criteria for validation of common variant and rare variant experiments. For initial common variant experiments, the criteria was 0.005 in 5 or more cells. This can be tailored to the precise needs of differing experimental designs. Code filters columns with names of ‘X.p.value.all. , X.p.value.test. , X.cells.w.guide.and.target. , and X.cells.w.guide.’ which will be outputted in 1) 2) or 3). 






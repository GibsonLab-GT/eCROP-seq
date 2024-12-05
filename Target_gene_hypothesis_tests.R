library(dplyr)
library(Seurat)

#read QC data from eCROPseq experiments in as data

##### PRE-PROCESSING #####
info <- data.frame('Gene','gRNA','STATUS','CHR','TARGET','cells.w.guide', 'cells.w.guide.and.target', 'p.value.all', 'p.value.test','Test')

#filter file for target gene of interest
setwd("Path to file with gRNA and target gene information")
to_test<-read.csv('Name of File.csv')
to_test <- to_test %>% filter(GENE=='Target gene of interest')

#specify all gRNA in pool
GENE1= 'RSID_1'
GENE2= 'RSID_2'
GENE3= 'RSID_3'

#extract cells with gRNA in pool
EXPR1=GetAssayData(object=data, assay='RNA',slot='data')[GENE1,]
EXPR2=GetAssayData(object=data, assay='RNA',slot='data')[GENE2,]
EXPR3=GetAssayData(object=data, assay='RNA',slot='data')[GENE3,]

#remove cells that do not have gRNA targeting any target gene
test=data.frame(positive= EXPR1>0 | EXPR2>0 | EXPR3>0)
data <- SetIdent(data, value = test$positive)
data <- subset(x = data, idents = "TRUE")

#specify gRNA corresponding to target gene of interest 
test=data.frame(positive= EXPR1==0 & EXPR2==0)
data <- SetIdent(data, value = test$positive)

# and separate in group without any gRNA targeting target gene
data_sub_pos <- subset(x = data, idents = "TRUE")

# or with one of the gRNA
data_sub_neg <- subset(x = data, idents = "FALSE")


##### Perform Hypothesis Tests (Student's t-test and KS)
i=1
for (i in 1:dim(to_test)[1]) {{
  tryCatch({
    GENE=to_test$RSID[i]
    gRNA=to_test$gRNA[i]
    STATUS=to_test$STATUS[i]
    CHR=to_test$CHR[i]
    TARGET=to_test$GENE[i]
    empty='NA'
    
###### grna 1 to cells without gRNA 1 in subsetted data ######
## create dataframes for target gene expression and gRNA presence in cells with gRNA targeting target gene 
EXPR = GetAssayData(object=data_sub_neg, assay="RNA",slot="data")[GENE,]
EXPR_target = GetAssayData(object=data_sub_neg, assay="RNA",slot="data")[TARGET,]
EXPR_df2=data.frame(EXPR, EXPR_target)
rm(EXPR,EXPR_target)

## create dataframes for target gene expression in cells without any gRNA targeting target gene
EXPR_target_NULL = GetAssayData(object=data_sub_pos, assay="RNA",slot="data")[TARGET,]
EXPR_target_NULL=data.frame(EXPR_target_NULL)
EXPR_target_NULL <- EXPR_target_NULL %>% rename('EXPR_target'='EXPR_target_NULL')

## extract target gene expression in all cells (i.e. not further subsetted for expression of target gene)
df_nogRNA <- EXPR_target_NULL %>%select(EXPR_target)
df_gRNA <- EXPR_df2 %>%filter(EXPR != 0)%>%select(EXPR_target)

## Subset further for cells that express target gene and extract target gene expression
test <- df_gRNA %>%filter(EXPR_target > 0)
no_g_test <- df_nogRNA %>%filter(EXPR_target > 0)

##ttest
p1 <- t.test(df_nogRNA, df_gRNA, alterntive= c("two.sided"), paired=FALSE, conf.level = 0.95)$p.value
p2 <- t.test(no_g_test, test, alterntive= c("two.sided"), paired=FALSE, conf.level = 0.95)$p.value
a<- nrow(df_gRNA)
b<- nrow(test)
TEST='T'
info[nrow(info) + 1,] = c(GENE,gRNA,STATUS,CHR,TARGET,a,b,p1,p2,TEST)

##ks
TEST='KS'
p1 <- ks.test(df_nogRNA$EXPR_target, df_gRNA$EXPR_target, alterntive= c("two.sided"))$p.value
p2 <- ks.test(no_g_test$EXPR_target, test$EXPR_target, alterntive= c("two.sided"))$p.value
info[nrow(info) + 1,] = c(GENE,gRNA,STATUS,CHR,TARGET,a,b,p1,p2,TEST)

write.csv(info,'Change_name_of_file.csv')

}, error=function(e){})
}
i= i + 1
}

library(dplyr)
library(Seurat)

#read QC data in as data

##### PRE-PROCESSING #####
info <- data.frame('Gene','gRNA','STATUS','CHR','TARGET','cells.w.guide', 'cells.w.guide.and.target', 'p.value.all', 'p.value.test','Test')

#filter file for gene of interest
setwd("C:/Users/green/GaTech Dropbox/Emily Greenwood/Sequencing runs 4-2024/Files to run")
to_test<-read.csv('pool12_hyptest_file.csv')
to_test <- to_test %>% filter(GENE=='FADS1')

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

# and separate in group with one of those gRNA
data_sub_pos <- subset(x = data, idents = "TRUE")

# or without any of them
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
    
######grna 1 to cells without gRNA 1###########################
EXPR = GetAssayData(object=data_sub_neg, assay="RNA",slot="data")[GENE,]
EXPR_target = GetAssayData(object=data_sub_neg, assay="RNA",slot="data")[TARGET,]
EXPR_df2=data.frame(EXPR, EXPR_target)
rm(EXPR,EXPR_target)

EXPR_target_NULL = GetAssayData(object=data_sub_pos, assay="RNA",slot="data")[TARGET,]
EXPR_target_NULL=data.frame(EXPR_target_NULL)
EXPR_target_NULL <- EXPR_target_NULL %>% rename('EXPR_target'='EXPR_target_NULL')
df_nogRNA <- EXPR_target_NULL %>%select(EXPR_target)
df_gRNA <- EXPR_df2 %>%filter(EXPR != 0)%>%select(EXPR_target)
test <- df_gRNA %>%filter(EXPR_target > 0)
no_g_test <- df_nogRNA %>%filter(EXPR_target > 0)

##ttest
p1 <- t.test(df_nogRNA, df_gRNA, alterntive= c("two.sided"), paired=FALSE, conf.level = 0.95)$p.value
p2 <- t.test(no_g_test, test, alterntive= c("two.sided"), paired=FALSE, conf.level = 0.95)$p.value
a<- nrow(df_gRNA)
b<- nrow(test)
TEST='T'
info[nrow(info) + 1,] = c(GENE,gRNA,STATUS,CHR,TARGET,a,b,p1,p2,TEST)

##kns
TEST='KS'
p1 <- ks.test(df_nogRNA$EXPR_target, df_gRNA$EXPR_target, alterntive= c("two.sided"))$p.value
p2 <- ks.test(no_g_test$EXPR_target, test$EXPR_target, alterntive= c("two.sided"))$p.value
info[nrow(info) + 1,] = c(GENE,gRNA,STATUS,CHR,TARGET,a,b,p1,p2,TEST)

#write.csv(info,'Change_name_of_file.csv')

}, error=function(e){})
}
i= i + 1
}

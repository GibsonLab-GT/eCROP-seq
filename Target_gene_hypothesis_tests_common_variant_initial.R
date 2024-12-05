library(dplyr)
library(Seurat)

#read QC Seurat object in as data
#initialize data frame
info <- data.frame('SNP','STATUS','TARGET','cells.w.guide', 'cells.w.guide.and.target', 'p.value.all', 'p.value.test','Test')

#read in file with gRNA and target information
setwd("Path to file with gRNA and target gene information")
to_test<-read.csv('Name of File.csv')

i=1
for (i in 1:dim(sig)[1]) {{
  tryCatch({
    RSID=sig$RSID[i]
    STATUS=sig$STATUS[i]
    CHR=sig$CHR[i]
    TARGET=sig$GENE[i]
    
######grna 1 to cells without gRNA 1###########################
EXPR = GetAssayData(object=data, assay="RNA",slot="data")[RSID,]
EXPR_target = GetAssayData(object=data, assay="RNA",slot="data")[TARGET,]
EXPR_df2=data.frame(EXPR, EXPR_target)
rm(EXPR,EXPR_target)

df_nogRNA <- EXPR_df2 %>%filter(EXPR == 0)%>%select(EXPR_target)
df_gRNA <- EXPR_df2 %>%filter(EXPR != 0)%>%select(EXPR_target)
test <- df_gRNA %>%filter(EXPR_target > 0)
no_g_test <- df_nogRNA %>%filter(EXPR_target > 0)

##ttest
p1 <- t.test(df_nogRNA, df_gRNA, alterntive= c("two.sided"), paired=FALSE, conf.level = 0.95)$p.value
p2 <- t.test(no_g_test, test, alterntive= c("two.sided"), paired=FALSE, conf.level = 0.95)$p.value
a<- nrow(df_gRNA)
b<- nrow(test)
TEST='T'
info[nrow(info) + 1,] = c(RSID,STATUS,TARGET,a,b,p1,p2,TEST)

##ks
TEST='KS'
p1 <- ks.test(df_nogRNA$EXPR_target, df_gRNA$EXPR_target, alterntive= c("two.sided"))$p.value
p2 <- ks.test(no_g_test$EXPR_target, test$EXPR_target, alterntive= c("two.sided"))$p.value
info[nrow(info) + 1,] = c(RSID,STATUS,TARGET,a,b,p1,p2,TEST)

write.csv(info, 'Name of hypothesis test file.csv')

  }, error=function(e){})
}
  i= i + 1
}

  

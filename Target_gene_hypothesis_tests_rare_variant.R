library(dplyr)
library(Seurat)
#read QC seurat object in as data

#initialize data frame
info <- data.frame('SNP_1','SNP_2','SNP','STATUS','Target','cells.w.guide', 'cells.w.guide.and.target', 'p.value.all', 'p.value.test','test')

#read in file with target gene of interest with corresponding gRNA information
setwd("Path to file with gRNA and target gene information")
to_test<-read.csv('Name of File.csv')
to_test <- to_test %>% filter(GENE=='Target gene of interest')

##### Perform Hypothesis Tests (Student's t-test and KS)
i=1
for (i in 1:dim(sig)[1]) {{
  tryCatch({
    SNP_1=sig$SNP_1[i]
    SNP_2=sig$SNP_2[i]
    SNP=sig$SNP[i]
    STATUS=sig$STATUS[i]
    TARGET=sig$GENE[i]
    empty='NA'

######grna 1 to cells without gRNA 1###########################
EXPR = GetAssayData(object=data, assay="RNA",slot="data")[SNP_1,]
EXPR_target = GetAssayData(object=data, assay="RNA",slot="data")[TARGET,]
EXPR_df2=data.frame(EXPR, EXPR_target)
rm(EXPR,EXPR_target)
df_nogRNA <- EXPR_df2 %>%filter(EXPR == 0)%>%select(EXPR_target)
df_gRNA <- EXPR_df2 %>%filter(EXPR != 0)%>%select(EXPR_target)
test <- df_gRNA %>%filter(EXPR_target > 0)
no_g_test <- df_nogRNA %>%filter(EXPR_target > 0)

TEST='t-test'
p1 <- t.test(df_nogRNA, df_gRNA, alterntive= c("two.sided"), paired=FALSE, conf.level = 0.95)$p.value
p2 <- t.test(no_g_test, test, alterntive= c("two.sided"), paired=FALSE, conf.level = 0.95)$p.value
a<- nrow(df_gRNA)
b<- nrow(test)
info[nrow(info) + 1,] = c(SNP_1,empty,SNP,STATUS,TARGET,a,b,p1,p2,TEST)

TEST='ks'
p1 <- ks.test(df_nogRNA$EXPR_target, df_gRNA$EXPR_target, alterntive= c("two.sided"))$p.value
p2 <- ks.test(no_g_test$EXPR_target, test$EXPR_target, alterntive= c("two.sided"))$p.value
info[nrow(info) + 1,] = c(SNP_1,empty,SNP,STATUS,TARGET,a,b,p1,p2,TEST)

#############gRNA 2 to cells without gRNA2##############################################
EXPR = GetAssayData(object=data, assay="RNA",slot="data")[SNP_2,]
EXPR_target = GetAssayData(object=data, assay="RNA",slot="data")[TARGET,]
EXPR_df2=data.frame(EXPR, EXPR_target)
rm(EXPR,EXPR_target)
df_nogRNA <- EXPR_df2 %>%filter(EXPR == 0)%>%select(EXPR_target)
df_gRNA <- EXPR_df2 %>%filter(EXPR != 0)%>%select(EXPR_target)
test <- df_gRNA %>%filter(EXPR_target > 0)
no_g_test <- df_nogRNA %>%filter(EXPR_target > 0)

TEST='t-test'
p1 <- t.test(df_nogRNA, df_gRNA, alterntive= c("two.sided"), paired=FALSE, conf.level = 0.95)$p.value
p2 <- t.test(no_g_test, test, alterntive= c("two.sided"), paired=FALSE, conf.level = 0.95)$p.value
a<- nrow(df_gRNA)
b<- nrow(test)
info[nrow(info) + 1,] = c(empty,SNP_2,SNP,STATUS,TARGET,a,b,p1,p2, TEST)

TEST='ks'
p1 <- ks.test(df_nogRNA$EXPR_target, df_gRNA$EXPR_target, alterntive= c("two.sided"))$p.value
p2 <- ks.test(no_g_test$EXPR_target, test$EXPR_target, alterntive= c("two.sided"))$p.value
info[nrow(info) + 1,] = c(empty,SNP_2,SNP,STATUS,TARGET,a,b,p1,p2, TEST)

################gRNA 1 or 2 to cells without 1 or 2###############################
EXPR1 = GetAssayData(object=data, assay="RNA",slot="data")[SNP_1,]
EXPR2 = GetAssayData(object=data, assay="RNA",slot="data")[SNP_2,]
EXPR_target = GetAssayData(object=data, assay="RNA",slot="data")[TARGET,]
EXPR_df2=data.frame(EXPR1,EXPR2,EXPR_target)
rm(EXPR1,EXPR2,EXPR_target)
df_nogRNA <- EXPR_df2 %>%filter(EXPR1 == 0)%>% filter(EXPR2 ==0) %>%select(EXPR_target)
df_gRNA1 <- EXPR_df2 %>%filter(EXPR1 != 0)%>%select(EXPR_target)
df_gRNA2 <- EXPR_df2 %>%filter(EXPR2 != 0)%>%select(EXPR_target)
df_gRNA <- rbind(df_gRNA1,df_gRNA2)

test <- df_gRNA %>%filter(EXPR_target > 0)
no_g_test <- df_nogRNA %>%filter(EXPR_target > 0)

TEST='t-test'
p1 <- t.test(df_nogRNA, df_gRNA, alterntive= c("two.sided"), paired=FALSE, conf.level = 0.95)$p.value
p2 <- t.test(no_g_test, test, alterntive= c("two.sided"), paired=FALSE, conf.level = 0.95)$p.value
a<- nrow(df_gRNA)
b<- nrow(test)
info[nrow(info) + 1,] = c(SNP_1, SNP_2,SNP,STATUS,TARGET,a,b,p1,p2, TEST)

TEST='ks'
p1 <- ks.test(df_nogRNA$EXPR_target, df_gRNA$EXPR_target, alterntive= c("two.sided"))$p.value
p2 <- ks.test(no_g_test$EXPR_target, test$EXPR_target, alterntive= c("two.sided"))$p.value
info[nrow(info) + 1,] = c(SNP_1, SNP_2,SNP,STATUS,TARGET,a,b,p1,p2, TEST)


######gRNA1 to cells without either gRNA#########
EXPR1 = GetAssayData(object=data, assay="RNA",slot="data")[SNP_1,]
EXPR2 = GetAssayData(object=data, assay="RNA",slot="data")[SNP_2,]
EXPR_target = GetAssayData(object=data, assay="RNA",slot="data")[TARGET,]
EXPR_df2=data.frame(EXPR1,EXPR2,EXPR_target)
rm(EXPR1,EXPR2,EXPR_target)
df_nogRNA <- EXPR_df2 %>%filter(EXPR1 == 0)%>% filter(EXPR2 ==0) %>%select(EXPR_target)
df_gRNA <- EXPR_df2 %>%filter(EXPR1 != 0)%>%select(EXPR_target)

test <- df_gRNA %>%filter(EXPR_target > 0)
no_g_test <- df_nogRNA %>%filter(EXPR_target > 0)

TEST='t-test'
p1 <- t.test(df_nogRNA, df_gRNA, alterntive= c("two.sided"), paired=FALSE, conf.level = 0.95)$p.value
p2 <- t.test(no_g_test, test, alterntive= c("two.sided"), paired=FALSE, conf.level = 0.95)$p.value
a<- nrow(df_gRNA)
b<- nrow(test)
info[nrow(info) + 1,] = c(SNP_1, empty,SNP,STATUS,TARGET,a,b,p1,p2, TEST)

TEST='ks'
p1 <- ks.test(df_nogRNA$EXPR_target, df_gRNA$EXPR_target, alterntive= c("two.sided"))$p.value
p2 <- ks.test(no_g_test$EXPR_target, test$EXPR_target, alterntive= c("two.sided"))$p.value
info[nrow(info) + 1,] = c(SNP_1, empty,SNP,STATUS,TARGET,a,b,p1,p2, TEST)

######gRNA2 to cells without either gRNA##########
EXPR1 = GetAssayData(object=data, assay="RNA",slot="data")[SNP_1,]
EXPR2 = GetAssayData(object=data, assay="RNA",slot="data")[SNP_2,]
EXPR_target = GetAssayData(object=data, assay="RNA",slot="data")[TARGET,]
EXPR_df2=data.frame(EXPR1,EXPR2,EXPR_target)
rm(EXPR1,EXPR2,EXPR_target)

df_nogRNA <- EXPR_df2 %>%filter(EXPR1 == 0)%>% filter(EXPR2 ==0) %>%select(EXPR_target)
df_gRNA <- EXPR_df2 %>%filter(EXPR2 != 0)%>%select(EXPR_target)

test <- df_gRNA %>%filter(EXPR_target > 0)
no_g_test <- df_nogRNA %>%filter(EXPR_target > 0)

TEST='t-test'
p1 <- t.test(df_nogRNA, df_gRNA, alterntive= c("two.sided"), paired=FALSE, conf.level = 0.95)$p.value
p2 <- t.test(no_g_test, test, alterntive= c("two.sided"), paired=FALSE, conf.level = 0.95)$p.value
a<- nrow(df_gRNA)
b<- nrow(test)
info[nrow(info) + 1,] = c(empty,SNP_2,SNP,STATUS,TARGET,a,b,p1,p2,TEST)

TEST='ks'
p1 <- ks.test(df_nogRNA$EXPR_target, df_gRNA$EXPR_target, alterntive= c("two.sided"))$p.value
p2 <- ks.test(no_g_test$EXPR_target, test$EXPR_target, alterntive= c("two.sided"))$p.value
info[nrow(info) + 1,] = c(empty,SNP_2,SNP,STATUS,TARGET,a,b,p1,p2,TEST)

write.csv(info2, 'Name of hypothesis test file.csv')

  }, error=function(e){})
}
  i= i + 1
}

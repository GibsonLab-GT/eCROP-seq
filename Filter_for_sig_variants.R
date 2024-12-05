library(dplyr)
library(readr)
library(data.table)
library(stringr)
library(purrr)


# If multiple, read in all hypothesis test files that were outputted from "Target_gene_hypothesis_tests.R"
# Combine to one file
# Must be in same folder to start
# Note, this was modified from: https://www.gerkelab.com/blog/2018/09/import-directory-csv-purrr-readr/

path_alps="path_to_hypothesis_test_output_files"
DF <- list.files(path_alps, pattern = "*.csv", full.names = TRUE) %>%
  map_df(function(x) read.csv(x,) %>%
           mutate(X = as.character(x)) %>%
           mutate(filename=gsub(".csv","", basename(x))))


write.csv(DF, 'Combined file name.csv')


# Next filter for significant hits (criteria here is with a p-value less than or equal to 0.05 in 5 or more cells)
i<-info
info<-read.csv('C:/Users/green/GaTech Dropbox/Emily Greenwood/FILES FOR PAPER/BEST ANALYSIS Updated 10-28-2024 USED FOR PAPER/Updated_isig_idis/All_combined_undif_target.csv')
i$X.p.value.all. <- as.numeric(i$X.p.value.all.)
i$X.p.value.test. <- as.numeric(i$X.p.value.test.)
i$X.cells.w.guide.and.target. <- as.numeric(i$X.cells.w.guide.and.target.)
i$X.cells.w.guide. <- as.numeric(i$X.cells.w.guide.)

i_sig <- i %>%filter(X.cells.w.guide.and.target. >= 5 & X.p.value.test. <= 0.05 | X.cells.w.guide. >= 5 & X.p.value.all. <= 0.05)
i_sig<- i_sig %>% distinct(X.SNP., .keep_all = T) #to get accurate count of controls which were included in each pool, add "filename" here

i_dis <- i %>% distinct(X.SNP.,.keep_all = T) #to get accurate count of controls which were included in each pool, add "filename" here
table(i_dis$X.STATUS.)
table(i_sig$X.STATUS.)

write.csv(i_dis, 'All variants tested file.csv')
write.csv(i_sig, 'All variants significant file.csv')



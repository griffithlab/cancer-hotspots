library(plyr)
library(ggplot2)

setwd("~/Desktop/GBM")
readcount_files <-list.files("~/Desktop/GBM/download")
data <- read.table("chr6-g.152098787T-A_9273_classII_TUMOR.all_epitopes.tsv", header = TRUE, sep = '\t')
data2 <- read.table("chr17-g.7673802C-A_TUMOR.all_epitopes.tsv", header = TRUE, sep = '\t')
readcount_data <-lapply(readcount_files[3],read.table, header = TRUE, sep = '\t')
all_samples_correlation <-list()
read_count_vaf <-data.frame()
for(i in 1:length(readcount_data)){
  names <-colnames(readcount_data[[i]])
  sub_names <-names[c(12,10,11,13,14,15,19,33,34,35,25,26,27,28,31,32)] ## same order/headers as pvacseq_results_ranked ##
  gene_name <- readcount_data[[i]][[12]]
  mutation <-readcount_data[[i]][[10]]
  protein_pos <- readcount_data[[i]][[11]]
  HGVSc <- readcount_data[[i]][[13]]
  HGVSp <- readcount_data[[i]][[14]]
  hla_allele <-readcount_data[[i]][[15]]
  mt_epitope_seq <- readcount_data[[i]][[19]]
  mt_ic50 <- readcount_data[[i]][[33]]
  wt_ic50 <- readcount_data[[i]][[34]]
  fold_change <-readcount_data[[i]][[35]]
  DNA_depth <- readcount_data[[i]][[25]]
  DNA_VAF <- readcount_data[[i]][[26]]
  RNA_depth <- readcount_data[[i]][[27]]
  RNA_VAF <-readcount_data[[i]][[28]]
  gene_exp <- readcount_data[[i]][[31]]
  txpt_exp <- readcount_data[[i]][[32]]
  read_count_vaf <- as.data.frame(cbind(gene_name, mutation, protein_pos, HGVSc, HGVSp, hla_allele, mt_epitope_seq, mt_ic50, wt_ic50, fold_change, DNA_depth, DNA_VAF, RNA_depth, RNA_VAF, gene_exp, txpt_exp))
  colnames(read_count_vaf) <- sub_names
  all_samples_correlation[[i]] <-read_count_vaf
}
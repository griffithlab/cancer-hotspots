library(readxl)

hs <- read_tsv("~/Downloads/class1_hla_freq_dist.tsv", col_names = F)
hs2 <- read_tsv("~/Downloads/hla_freq_class2_dist.tsv", col_names = F)

hs3 <- rbind(hs, hs2)

uniq_genes <- length(unique(filter(hs3, X2 > 0)$X1))

ch <- read_excel("~/Desktop/Cancer_Hotspots/hotspots_summary_final.xlsx")
ch <- ch[c("Gene Name", "Neo ClassI Predictions", "Neo ClassII Predictions")]
ch$Total_neogs <- rowSums(ch[,-1])

cht <- ch %>% group_by(`Gene Name`) %>% summarize(Total_neogs = sum(Total_neogs))

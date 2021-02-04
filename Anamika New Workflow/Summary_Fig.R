library(ggplot2)
library(plotly)
library(gapminder)
library(stringr)

data1 <- read.csv("/gscmnt/gc2142/griffithlab/abasu/Summary_Fig_Filtered_Data.tsv", sep = "\t", header=FALSE)
colnames(data1) <- c("Gene", "HLA_Freq", "Neoantigens")

p <- 
  ggplot(data1, aes(Neoantigens, HLA_Freq, label = Gene)) +
  geom_point() +
  labs(title = "Maximum HLA Allele Frequency vs. Total Neoantigens per Gene",
       x = "Total Neoantigens per Gene",
       y = "HLA Allele Frequency(%)") 


options(scipen=999) # scientific notation
ggplotly(p) 
ggsave("Summary_Figure.png", width = 10, height = 6)
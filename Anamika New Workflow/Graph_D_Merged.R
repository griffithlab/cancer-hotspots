library(ggplot2)
library(gridExtra)
data <- read.csv("/gscmnt/gc2142/griffithlab/abasu/Graph_D_Class1_Data_Filtered.tsv", sep = "\t", header=FALSE)
colnames(data) <- c("Gene", "Count", "Group")

ggp <- ggplot(data, aes(reorder(Gene,Count), Count, fill = factor(Group, levels=c("50-500 nM","0-50 nM" )))) +  # Create stacked bar chart
  geom_bar(stat = "identity") +
  labs(title = "Graph D with Class I HLA Alleles", 
       y = "Neoantigen Count",
       fill = "Binding Affinity") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

data2 <- read.csv("/gscmnt/gc2142/griffithlab/abasu/Graph_D_Class2_Data_Filtered.tsv", sep = "\t", header=FALSE)
colnames(data2) <- c("Gene", "Count", "Group")

ggp2 <- ggplot(data2, aes(reorder(Gene,Count), Count, fill = factor(Group, levels=c("50-500 nM","0-50 nM" )))) +  # Create stacked bar chart
  geom_bar(stat = "identity") +
  labs(title = "Graph D with Class II HLA Alleles", 
       y = "Neoantigen Count",
       fill = "Binding Affinity") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())



ggsave("Graph_D_Merged.png", plot = grid.arrange(ggp, ggp2, nrow = 2), width = 10, height = 6)
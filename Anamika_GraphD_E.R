library(ggplot2)

# Graph E Class I
data <- read.csv("/Users/anamikabasu/Code/Heat_Map/Cancer_Hotspots_Graphs/Graph_E_Class1_Data.tsv", sep = "\t", header=FALSE)
colnames(data) <- c("Gene", "Count", "Group")

ggp <- ggplot(data, aes(reorder(Gene,Count), Count, fill = Group, label = Count)) +  # Create stacked bar chart
  geom_bar(stat = "identity") +
  labs(title = "Graph E with Class I HLA Alleles", 
       y = "Neoantigen Count") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave("/Users/anamikabasu/Code/Heat_Map/Cancer_Hotspots_Graphs/Graph_E_Class1_5_bins.png", plot = ggp,  width = 10, height = 6)

# Graph E Class II

data2 <- read.csv("/Users/anamikabasu/Code/Heat_Map/Cancer_Hotspots_Graphs/Graph_E_Class2_Data.tsv", sep = "\t", header=FALSE)
colnames(data2) <- c("Gene", "Count", "Group")

ggp2 <- ggplot(data2, aes(reorder(Gene,Count), Count, fill = Group, label = Count)) +  # Create stacked bar chart
  geom_bar(stat = "identity") +
  labs(title = "Graph E with Class II HLA Alleles", 
       y = "Neoantigen Count") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave("/Users/anamikabasu/Code/Heat_Map/Cancer_Hotspots_Graphs/Graph_E_Class2_5_bins.png", plot = ggp2, width = 10, height = 6)

# Graph D Class I

data3 <- read.csv("/Users/anamikabasu/Code/Heat_Map/Cancer_Hotspots_Graphs/Graph_D_Class1_Data.tsv", sep = "\t", header=FALSE)
colnames(data3) <- c("Gene", "Count", "Group")

ggp3 <- ggplot(data3, aes(reorder(Gene,Count), Count, fill = factor(Group, levels=c("50-500 nM","0-50 nM" )))) +  # Create stacked bar chart
  geom_bar(stat = "identity") +
  labs(title = "Graph D with Class I HLA Alleles", 
       y = "Neoantigen Count",
       fill = "Binding Affinity") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave("/Users/anamikabasu/Code/Heat_Map/Cancer_Hotspots_Graphs/Graph_D_Class1.png", plot = ggp3,  width = 10, height = 6)

# Graph D Class II

data4 <- read.csv("/Users/anamikabasu/Code/Heat_Map/Cancer_Hotspots_Graphs/Graph_D_Class2_Data.tsv", sep = "\t", header=FALSE)
colnames(data4) <- c("Gene", "Count", "Group")

ggp4 <- ggplot(data4, aes(reorder(Gene,Count), Count, fill = factor(Group, levels=c("50-500 nM","0-50 nM" )))) +  # Create stacked bar chart
  geom_bar(stat = "identity") +
  labs(title = "Graph D with Class II HLA Alleles", 
       y = "Neoantigen Count",
       fill = "Binding Affinity") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave("/Users/anamikabasu/Code/Heat_Map/Cancer_Hotspots_Graphs/Graph_D_Class2.png", plot = ggp4, width = 10, height = 6)

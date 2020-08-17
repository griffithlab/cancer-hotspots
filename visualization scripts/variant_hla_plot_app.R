# load packages
library(tidyverse)
library(hash)
library(ggsci)
library(wesanderson)
library(grid)
library(readxl)
library(ggpmisc)

# plot overview:
# data=results
# x=median binding affinity
# y=hla_frequency
# put in quadrant
# size of dot=number of peptides
# color of dot=hla class ? or make separate plot for class I/II

# functions
add_class <- function(allele) {
    if (grepl("HLA-A|HLA-B|HLA-C", allele)) {
        class <- "ClassI"
    } else { 
        class <- "ClassII"
    }
    return(class)
}



# setwd
setwd("~/Desktop/Cancer_Hotspots/")
tsv_files <- list.files(path="~/Desktop/Cancer_Hotspots", pattern = "*_TUMOR.all_epitopes.tsv", full.names = T)

# load in tsv_files for plotting
results <- read_tsv("chr1-g.11109318A-C_TUMOR.all_epitopes.tsv", col_types = cols(Reference = col_character(), Variant = col_character()))
# table with only neoantigens
neo_results <- filter(results, `Median MT Score` <= 500)
# create a table with neo hla alleles, frequency
hla_table <- as.data.frame(table(neo_results$`HLA Allele`), stringsAsFactors = F)
names(hla_table) <- c("HLA Allele", "Neoantigens per HLA Allele")
# add the class of each hla
hla_table$`HLA Class` <- sapply(hla_table$`HLA Allele`, add_class)
# add hla population frequency for each hla
hla_freq <- data.frame(read_excel("HLA_allele_frequencies/PVAC_HLA_Freq.xlsx"))
hla_table$`HLA Frequency` <- sapply(hla_table$`HLA Allele`, function(x) return((filter(hla_freq, HLA.Allele == x)$Frequency)*100))
# add median binding affinity for all neoantigens each hla
hla_table$`Median Binding Affinity` <- sapply(hla_table$`HLA Allele`, function(x) return(median(filter(neo_results, `HLA Allele` == x)$`Median MT Score`)))

# plot
# palette:
my_pal <- c(pal_npg("nrc", alpha = 1)(10), pal_npg("nrc", alpha = 1)(10))

# with and without hla frequency = 0
ggplot(data=hla_table) + geom_point(mapping=aes(x=`Median Binding Affinity`, y=`HLA Frequency`, color=`Neoantigens per HLA Allele`), size=3) + facet_wrap( ~ `HLA Class`) + theme_bw() + theme(legend.position = "bottom", panel.border = element_rect(color = "black", fill=NA, size=1), text = element_text(size=12)) + ylab("HLA Frequency (%)") #+ coord_cartesian(expand = FALSE)

       
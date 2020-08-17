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

# setwd
setwd("~/Desktop/Cancer_Hotspots/")


# load in results file
results <- read_tsv("chr1-g.114713911C-T_TUMOR.all_epitopes.tsv")
## special case for this file - change "TRUE" back to "T" ##
results$Variant <- "T"

# step 1: create a new neo df with columns: c("HLA Allele", "HLA Frequency", "HLA Class", Neoantigens per HLA Allele", "Median Binding Affinity")
# table with only neoantigens
neo_results <- filter(results, `Median MT Score` <= 500)
# create a table with neo HLA alleles, frequency
hla_table <- as.data.frame(plyr::count(neo_results$`HLA Allele`))
colnames(hla_table) <- c("HLA Allele", "HLA Frequency")
hla_table$`HLA Allele` <- as.character(hla_table$`HLA Allele`)
# get median binding affinity for each HLA
# get hla frequency for each HLA
hla_freq <- as.data.frame(read_excel("HLA_allele_frequencies/PVAC_HLA_Freq.xlsx"))
plot_df <- data.frame()
for ( row in 1:nrow(hla_table) ) {
    hla <- hla_table[row, `HLA Allele`]
    peptide_count <- hla_table[row, "HLA Frequency"]
    frequency <- filter(hla_freq, `HLA Allele` == hla)$`Frequency`
    hla_results <- filter(neo_results, `HLA Allele` == hla)
    if ( grepl("HLA-A|HLA-B|HLA-C", unique(hla_results$`HLA Allele`)) ) {
        class <- "Class I"
    } else {
        class <- "Class II"
    }
    median_binding <- median(hla_results$`Median MT Score`)
    df_row <- data.frame("HLA Allele" = hla, "HLA Frequency" = frequency, "HLA Class" = class, "Neoantigens per HLA Allele" = peptide_count, "Median Binding Affinity" = median_binding, stringsAsFactors = FALSE)
    colnames(df_row) <- gsub('\\.', " ", colnames(df_row))
    plot_df <- rbind(plot_df, df_row)
}
plot_df$`HLA Frequency` <- sapply(plot_df$`HLA Frequency`, function(x) return(x*100))

# step 2: plot!
#palettes:
my_pal <- c(pal_npg("nrc", alpha = 1)(10), pal_npg("nrc", alpha = 1)(10))

# with and without hla frequency = 0
ggplot(data=filter(plot_df, `HLA Frequency` > 0)) + geom_point(mapping=aes(x=`Median Binding Affinity`, y=`HLA Frequency`, color=`Neoantigens per HLA Allele`), size=3) + facet_wrap( ~ `HLA Class`) + theme_bw() + theme(legend.position = "bottom", panel.border = element_rect(color = "black", fill=NA, size=1), text = element_text(size=12)) + ylab("HLA Frequency (%)") #+ coord_cartesian(expand = FALSE)

# next plot:
# flip hla allele and neoantigen sequence


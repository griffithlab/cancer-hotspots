# load packages
library(tidyverse)
library(hash)
library(ggsci)
library(wesanderson)
library(grid)

# setwd
setwd("~/Desktop/Cancer_Hotspots/")

# amino acid abbreviations
amino_acids <- hash(keys=c('Ala', 'Arg', 'Asn', 'Asp', 'Asx', 'Cys', 'Glu', 'Gln', 'Glx', 'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val'), values=c('A', 'R', 'N', 'D', 'B', 'C', 'E', 'Q', 'Z', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'))

#### functions ####
replace_abbreviations <- function(item) {
    mutation <- unlist(strsplit(item, split="\\."))[3]
    replace_items <- unlist(str_extract_all(mutation, paste(keys(amino_acids), collapse = '|')))
    for ( index in 1:length(replace_items) ) {
        if ( index == 1 ) {
            final_name <- gsub(replace_items[index], values(amino_acids[replace_items[index]], USE.NAMES=F), mutation)
        } else if ( index > 1 ) {
            final_name <- gsub(replace_items[index], values(amino_acids[replace_items[index]], USE.NAMES=F), final_name)
        }
    }
    return(final_name)
}

create_hotspots_id <- function(df) {
    new_id_df <- c()
    for ( index in 1:nrow(df) ) {
        gene_name <- df[index, "Gene.Name"]
        mutation <- replace_abbreviations(df[index, "HGVSp"])
        new_id <- paste(gene_name, ":", mutation, sep = "")
        new_id_df <- c(new_id_df, new_id)
    }
    return(new_id_df)
}

#### initial data set-up ####
# load in results file
results <- read.csv("1-5_test.csv")
# change NA values to 0
results[is.na(results)] = 0
# normalize all Chromosome values
results$Chromosome <- sapply(results$Chromosome, function(x) return(gsub("chr", "", x)))
# convert factor columns to character
results %>% mutate_if(is.factor, as.character) -> results
# create unique hotspot ID for each variant
results$Hotspot_ID <- create_hotspots_id(results)
# check if any IDs are duplicated 
dups <- results$Hotspot_ID[duplicated(results$Hotspot_ID)]
if ( length(dups) != 0 ) {
    dup_df <- filter(results, Hotspot_ID %in% dups)
}

#format column names correctly
# remove dots from column names
colnames(results) <- gsub('\\.', " ", colnames(results))
# first, format column names ("<" and "-" were replaced with ".")
colnames(results)[grep("50nM", colnames(results))] <- str_replace_all(colnames(results)[grep("50nM", colnames(results))], "  ", " <")
colnames(results)[grep("50 500nM", colnames(results))] <- str_replace_all(colnames(results)[grep("50 500nM", colnames(results))], "50 500", "50-500")
colnames(results)[grep("500 1000nM", colnames(results))] <- str_replace_all(colnames(results)[grep("500 1000nM", colnames(results))], "500 1000", "500-1000")

# order samples by total neos
results$`Total Score` <- rowSums(results[, c("Total Score < 50nM", "Total Score 50-500nM", "Total Score 500-1000nM")])
variant_order <- results[order(results$`Total Score`), ]$Hotspot_ID
# get total/classI/classII all formatted
results.gather <- gather(results, key=Score, value=Count, colnames(results)[grep("nM", colnames(results))])
# order on hotspots ID
results.gather$Hotspot_ID <- factor(results.gather$Hotspot_ID, levels = variant_order)
results.gather.order1 <- results.gather[order(results.gather$Hotspot_ID), ]
# create score type vector for facet grid
results.gather.order1$`Score Type` <- sapply(results.gather.order1$Score, function(x) return(unlist(strsplit(x, " "))[1]))
results.gather.order1$Score <- sapply(results.gather.order1$Score, function(x) return(str_replace(x, "Total |ClassII |ClassI |Score ", "")))
# order on scores
results.gather.order1$Score <- factor(results.gather.order1$Score, levels = unique(sort(results.gather.order1$Score, decreasing = T)))
results.gather.order12 <- results.gather.order1[order(results.gather.order1$Score), ]
# change factor levels for Score Type
results.gather.order12$`Score Type` <- factor(results.gather.order12$`Score Type`, levels = c("Total", "ClassI", "ClassII"))
colnames(results.gather.order12)[47] <- "Binding Affinity"

#### visualizations ####
# palettes
my_pal <- c("lightskyblue1",pal_npg("nrc", alpha = 1)(10)[c(4,1)])
pal = wes_palette("Zissou1", 100, type = "continuous")

# stacked bars
# y=Count, fill=Score
ggplot(data=results.gather.order12, aes(x=Hotspot_ID, y=Count, fill=`Binding Affinity`)) + geom_bar(stat="identity") + facet_wrap(~ `Score Type`, nrow=3, scales = "free") + scale_fill_manual(values = my_pal) + theme_bw() + theme(legend.position = "bottom", axis.title.x=element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.border = element_rect(color = "black", fill=NA, size=1), text = element_text(size=10), axis.ticks.y = element_line(), panel.grid.minor = element_blank()) + coord_cartesian(expand = FALSE)
 
# y=Score, fill=Chromosome
ggplot(data=results_scores.ordered2, aes(x=Hotspot_ID, y=Median.MT.Score.Min, fill=Chromosome)) + geom_bar(stat="identity") + theme(legend.position = "bottom", legend.title = element_blank(), axis.title.x=element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.border = element_rect(color = "black", fill=NA, size=1), text = element_text(size=12)) + coord_cartesian(expand = FALSE)



# visualization - length enrichment in neoantigens
# need to figure out enrichment first - how to represent this data?
results_length <- results[, c(colnames(results)[str_detect(colnames(results), "^Seq|^Neo")], "Hotspot_ID")]
results_length.gather <- gather(results_length, key="Lengths", value="Value", c(colnames(results)[str_detect(colnames(results), "^Seq|^Neo")]))
# define an enrichment score - just neo_length/seq_length ?
ggplot(results_length.gather, aes(x=Lengths, y=Value, fill=Lengths)) + geom_violin(scale="width", trim=F) + scale_fill_manual(values = my_pal) + theme_bw() + theme(legend.position = "bottom", legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_rect(color = "black", fill=NA, size=1), text = element_text(size=12))


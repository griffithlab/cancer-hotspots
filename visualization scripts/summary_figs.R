library(tidyverse)
library(RColorBrewer)
library(gridExtra)

setwd("~/Desktop/Cancer_Hotspots/")

### process data into df ###

# read in df
df <- read_excel("hotspots_summary.xlsx", col_types = c("text", rep("guess", 66))) 
# reorder Chromosome since I specified it as character
df$Chromosome <- factor(df$Chromosome, levels = c("X", rev(seq(1,22,1))))
df <- df[order(df$Chromosome), ]
# remove counts columns
df <- df[,c(!grepl("Length Count", names(df)))]
# remove any columns with all NA values
df <- df[1:24]

### summarize data ###

# hs per chr - variant type 
sum_1 <- df %>% group_by(Chromosome, `Variant Type`) %>% tally() %>% filter(`Variant Type` != "FS")
ggplot(data = sum_1, aes(x = Chromosome, y = n, fill = `Variant Type`)) + geom_col() + coord_flip() + theme_classic() + ylab("Hotspot Count") + theme(legend.position = "none")
# neo / total pred per chr - variant type
pred_type <- "Total Neoantigens"
sum_1.5 <- df[,c("Chromosome", pred_type, "Variant Type")] %>% filter(`Variant Type` != "FS")
sum_1.5.test <- aggregate(sum_1.5[,pred_type], by=list(Variant.Type = sum_1.5$`Variant Type`, Chromosome = sum_1.5$Chromosome), FUN=sum)
ggplot(data = sum_1.5.test, aes(x = Chromosome, y = `Total Predictions`, fill = Variant.Type)) + geom_col() + coord_flip() + theme_classic() + ylab(paste(pred_type, "Count", sep = ' ')) + theme(axis.title.y = element_blank(), legend.position = "none") 
# what percent of predictions are neoantigens?
percent_df <- df[,c("Chromosome", "Total Predictions", "Total Neoantigens")]
percent_df <- aggregate(list(Total.Predictions = percent_df$`Total Predictions`, Total.Neoantigens = percent_df$`Total Neoantigens`), by=list(Chromosome = percent_df$Chromosome), FUN=sum)
percent_df$Percent <- percent_df$Total.Neoantigens / percent_df$Total.Predictions * 100
percent_df$Percent <- sapply(percent_df$Percent, function(x) round(x, digits = 2))
#merged_df <- merge(sum_1.5.test, percent_df, by="Chromosome")
ggplot(data = percent_df, aes(x = Chromosome, y = Percent)) + geom_col() + coord_flip() + theme_classic() + ylim(0,10) + ylab("% Neoantigens from Total Predictions")


#+ theme_classic() + ylab(paste(pred_type, "Count", sep = ' ')) + theme(axis.title.y = element_blank(), legend.position = "none") 
# hotspots per gene
sum_2 <- df %>% group_by(`Gene Name`) %>% tally()
sum_2 <- sum_2[order(sum_2$n),]
sum_2$`Gene Name` <- factor(sum_2$`Gene Name`, levels = sum_2$`Gene Name`)
ggplot(data = sum_2, aes(x = `Gene Name`, y = n)) + geom_col() + theme_classic() + ylab("Hotspot Count") + xlab("Gene") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + coord_cartesian(expand=F)

# total predictions
sum_3 <- sum(df$`Total Predictions`)
sum_4 <- sum(df$`Total Neoantigens`)

### allele freqs ###
allele_df <- as.data.frame(read_excel("HLA_allele_frequencies/PVAC_HLA_Freq.xlsx"))
# how many zero freq?
zero_count <- nrow(allele_df %>% filter(Frequency >= 0.01))
# df with no zero values
allele_count <- allele_df %>% filter(Frequency > 0)
class_est <- function(x) {
    if (grepl("HLA", x)) {
        p <- "ClassI"
    } else {
        p <- "ClassII"
    }
    return(p)
}
allele_count$`HLA Class` <- sapply(allele_count$`HLA Allele`, class_est) 
allele_count <- allele_count[order(allele_count$Frequency),]
allele_count$`HLA Allele` <- factor(allele_count$`HLA Allele`, levels=c(allele_count$`HLA Allele`))
ggplot(data = allele_count, aes(x = `HLA Allele`, y = Frequency, fill = `HLA Class`)) + geom_col() + theme_classic() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + coord_cartesian(expand=F) + geom_hline(yintercept = c(0.01, 0.001), linetype="dashed")

### best neoantigens plot ###
sum_df <- na.omit(df[,c("Chromosome", "Gene Name", "HGVSp", "HLA Allele Min", "Median MT Score Min")])
sum_df$`HLA Frequency` <- sapply(sum_df$`HLA Allele Min`, function(x) filter(allele_df, `HLA Allele` == x)$Frequency)
sum_df$`HLA Frequency` <- as.numeric(sum_df$`HLA Frequency`)
sum_df$`HLA Class` <- sapply(sum_df$`HLA Allele Min`, class_est)
ggplot(data = sum_df, aes(x = `Median MT Score Min`, y = `HLA Frequency`, col = `HLA Class`)) + geom_point() + geom_hline(yintercept = c(0.001, 0.01, 0.05, 0.10), linetype = "dashed") + facet_wrap(~ `HLA Class`) + theme_classic()
ggplot(data = sum_df, aes(x = `HLA Frequency`, fill = `HLA Class`)) + geom_histogram(alpha = 0.5)
# subset
sum_df_ag <- sum_df %>% group_by(Chromosome, `HLA Class`) %>% tally()
ggplot(data = sum_df_ag, aes(x = Chromosome, y = n, fill = `HLA Class`)) + geom_col() + coord_flip() + theme_classic() + ylab("Median MT Score Min")

### updated summary csv for sep class I / II info ###
class_df <- read_excel("variant_summary_new.xlsx", col_types = c("text", "text", rep("guess", 51)))
# reorder Chromosome since I specified it as character
class_df$Chromosome <- factor(class_df$Chromosome, levels = c("X", rev(seq(1,22,1))))
class_df <- class_df[order(class_df$Chromosome), ]
# remove counts columns
class_df <- class_df[,c(!grepl("Length Count", names(class_df)))]
neo_I <- sum(class_df$`Neo ClassI Predictions`)
neo_II <- sum(class_df$`Neo ClassII Predictions`)

class_neo <- class_df[,c("Chromosome", "ClassI Predictions", "ClassII Predictions")]
class_neo_ag <- aggregate(list(ClassI = class_df$`ClassI Predictions`, ClassII = class_df$`ClassII Predictions`), by=list(Chromosome = class_neo$Chromosome), FUN=sum)
class_neo_ag <- class_neo_ag %>% gather(key = "HLA Class", value = "Total Count", -Chromosome)
ggplot(data = class_neo_ag, aes(x = Chromosome, y = `Total Count`, fill = `HLA Class`)) + geom_bar(stat = "identity") + coord_flip() + theme_classic() + ylab("Total Predictions Count") + theme(legend.position = "none")



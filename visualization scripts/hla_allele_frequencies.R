# 11.1.19
# parse HLA frequencies data from 1000 genome project
# this table is from the paper: HLA diversity in the 1000 genomes dataset (Gourrand et al., PLOS ONE, 2014)

# load packages
library(tidyverse)
library(readxl)
library(openxlsx)

# set wd
setwd("~/Desktop/Cancer_Hotspots/")

# provide path to excel file
hla_dataset <- "1000_genome_project_hla_frequencies.xlsx"

# create blank excel workbook
out <- createWorkbook()

# read in each sheet of original excel file and write results to "out"
for ( sheet in excel_sheets(hla_dataset)[1:5] ) {
  percent_df <- data.frame()
  print(sheet)
  initial_df <- as.data.frame(read_excel(hla_dataset, sheet = sheet))
  for ( row in 1:nrow(initial_df) ) {
    new_line <- c()
    allele <- initial_df[row, 1]
    allele <- gsub("g", "", allele)
    for ( item in initial_df[row, 2:ncol(initial_df)] ) {
      if ( grepl("NA", item) == F ) {
        percent <- str_extract_all(item, "\\([^()]+\\)")[[1]]
        new_item <- substring(percent, 2, nchar(percent)-2)
        new_line <- c(new_line, new_item)
      } else {
        new_item <- 0
        new_line <- c(new_line, new_item)
      }
    }
    new_line <- as.numeric(new_line)
    avg_fraction <- mean(new_line)
    line_df <- data.frame(allele, avg_fraction)
    colnames(line_df) <- c("Allele", "1000 Genomes Frequency (%)")
    percent_df <- rbind(percent_df, line_df)
  }
  ordered_df <- percent_df[order(percent_df$`1000 Genomes Frequency (%)`, decreasing = T), ]
  assign(sheet, ordered_df)
  addWorksheet(out, sheetName = sheet)
  writeData(out, sheet = sheet, x = ordered_df)
}

# save workbook "out" to excel file
saveWorkbook(out, "average_hla_frequencies.xlsx")

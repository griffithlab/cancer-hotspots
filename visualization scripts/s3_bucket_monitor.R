library(ggplot2)
library(dplyr)
library(tidyr)
setwd("~/Desktop/Cancer_Hotspots/R_scripts")

# read in file
tsv_file <- read.csv('../s3_status_hourly.txt', header=FALSE)

# create vectors from data in file
dates <- grep(" 2019", tsv_file$V1, value=TRUE)
total_objects <- grep("Objects:", tsv_file$V1, value=TRUE)
total_size <- grep("Size:", tsv_file$V1, value=TRUE)

#modify objects in vector
total_objects <- as.integer(sub("Total Objects: ", "",total_objects))
total_size <- as.integer(sub("Total Size: ", "",total_size))
#create df
plot_df <- data.frame(dates, total_objects, total_size)
#subset df
first_run_df <- plot_df[1:57,]
second_run_df <- plot_df[73:nrow(plot_df),]
# final df for plot
run_df <- rbind(first_run_df, second_run_df)
run_df$hours <- c(1:nrow(run_df))
run_df$Instance_Count <- c(rep("40", 57), rep("68", nrow(run_df)-57))
run_df$size_MB <- run_df$total_size/1000000

# objects plot
ggplot(run_df, aes(x=hours, y=total_objects, color=Instance_Count)) + geom_line(size=2) + ggtitle("AWS Hotspots Analysis s3 Bucket") + theme(plot.title = element_text(hjust = 0.5), text = element_text(size=15)) + ylab("Total Objects") + xlab("Hours")

# size plot
ggplot(run_df, aes(x=hours, y=size_MB, color=Instance_Count)) + geom_line(size=2) + ggtitle("AWS Hotspots Analysis s3 Bucket") + theme(plot.title = element_text(hjust = 0.5), text = element_text(size=15)) + ylab("Total Size (MB)") + xlab("Hours")


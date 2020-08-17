library(ggplot2)
library(dplyr)
library(tidyr)
setwd("~/Desktop/Cancer_Hotspots/R_scripts")

# read in file
tsv_file <- read.csv('../classI_monitor_plot_3.21.tsv', header=FALSE)

# create vectors from data in file
dates <- grep(" 2019", tsv_file$V1, value=TRUE)
total_objects <- grep("Total", tsv_file$V1, value=TRUE)


#modify objects in vector
total_objects <- as.integer(sub("Total class I jobs submitted: ", "",total_objects))
#create df
plot_df <- data.frame(dates, total_objects)

for ( row in 1:nrow(plot_df) ) {
  date <- plot_df[row, "dates"]
  print(date)
}
  



# objects plot
ggplot(run_df, aes(x=hours, y=total_objects, color=Instance_Count)) + geom_line(size=2) + ggtitle("AWS Hotspots Analysis s3 Bucket") + theme(plot.title = element_text(hjust = 0.5), text = element_text(size=15)) + ylab("Total Objects") + xlab("Hours")
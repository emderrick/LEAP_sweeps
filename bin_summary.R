library(tidyverse)

bin_files <- list.files(pattern = ".*bins_summary.txt", full.names = T)

all_bins <- data.frame()
for(i in 1:length(bin_files)){
  bin_pond <- read.table(bin_files[i], sep="\t", header = T)
  bin_pond <- bin_pond[, c(1:7)]
  pond <- bin_files[i] %>% substr(3,4)
  bin_pond <- cbind(bin_pond, pond = rep(pond, nrow(bin_pond)))
  all_bins <- rbind(all_bins, bin_pond)
}

all_bins$total <- 1
bins_sum <- all_bins %>% group_by(pond) %>% summarize(total_bins = sum(total))
write.csv(bins_sum, "bin_summary.csv", row.names = F)

final_bin_files <- list.files(pattern = ".*final_summary.txt", full.names = T)

final_bins <- data.frame()
for(i in 1:length(final_bin_files)){
  bin_pond <- read.table(final_bin_files[i], sep="\t", header = T)
  bin_pond <- bin_pond[, c(1:7)]
  pond <- final_bin_files[i] %>% substr(3,4)
  bin_pond <- cbind(bin_pond, pond = rep(pond, nrow(bin_pond)))
  final_bins <- rbind(final_bins, bin_pond)
}

final_bins$total <- 1
final_bins_sum <- final_bins %>% group_by(pond) %>% summarize(total_bins = sum(total))
write.csv(final_bins_sum, "final_bin_summary.csv", row.names = F)

bins_70 <- subset(final_bins, percent_completion >= 70 & percent_redundancy <= 10)
bins_70$total <- 1
bins_70_sum <- bins_70 %>% group_by(pond) %>% summarize(total_bins = sum(total))
test <- subset(bins_70, pond == "I4")

library(tidyverse)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

SNV_files <- list.files("data files/", recursive = T, pattern = "T1_filt_SNVs", full.names = T)

SNV_count <- data.frame()
for(i in 1:length(SNV_files)){
  SNV_sample <- read_csv(SNV_files[i])
  SNV_sample$Name_Time <- paste(SNV_sample$Name, SNV_sample$time, sep = " ")
  SNV_sample <- SNV_sample[, c("Name_Time", "mag", "mag_coverage", "med_cov", "mag_breadth", "mag_length", "class", "mutation_type", "group")]
  SNV_sample$class <- SNV_sample$class %>% str_remove("con_")
  SNV_sample$class <- SNV_sample$class %>% str_remove("pop_")
  SNV_count <- rbind(SNV_count, SNV_sample)
}

write.csv(SNV_count, "data files/T1_limited_SNV_info.csv", row.names = F)

mutation_counts <- SNV_count %>% group_by(mag, Name_Time, mag_coverage, med_cov, mag_breadth, mag_length, mutation_type)  %>% count(class)
mutation_counts <- pivot_wider(mutation_counts, values_from = n, names_from = mutation_type)
mutation_counts[is.na(mutation_counts)] <- 0
mutation_counts$NS_ratio <- mutation_counts$N / mutation_counts$S
mutation_counts$total <- rowSums(mutation_counts[, c("I","M","N","S","NA")], na.rm = T)
mutation_counts <- mutation_counts[, c(1:6,7,10,14,13)]
mutation_counts <- pivot_wider(mutation_counts, names_from = class, values_from = c(total, NS_ratio, N))
mutation_counts$SNSs_Mbp <- (mutation_counts$total_SNS/ mutation_counts$mag_length) * 10^6
mutation_counts$SNVs_Mbp <- (mutation_counts$total_SNV / mutation_counts$mag_length) * 10^6
mutation_counts$N_SNSs_Mbp <- (mutation_counts$N_SNS/ mutation_counts$mag_length) * 10^6
mutation_counts$N_SNVs_Mbp <- (mutation_counts$N_SNV / mutation_counts$mag_length) * 10^6
mag_info <- read_csv("data files/T1_mag_info.csv")
mag_info$Name_Time <- paste(mag_info$Name, mag_info$time, sep = " ")
mag_info <- mag_info[, c(1,2,9,10)]
mutation_counts <- full_join(mutation_counts, mag_info)
mutation_counts$Treatment <- ifelse(grepl("CTRL", mutation_counts$Name_Time), "Control", "GBH")
mutation_counts$Name <- mutation_counts$Name_Time %>% str_sub(end = -2)
mutation_counts$Time <- ifelse(grepl("1", mutation_counts$Name_Time), "Day 0", "Day 28")
mutation_counts$Treatment_Time <- paste(mutation_counts$Treatment, mutation_counts$Time, sep = " ")
write.csv(mutation_counts, "data files/T1_SNV_summary_MAG.csv", row.names = F)

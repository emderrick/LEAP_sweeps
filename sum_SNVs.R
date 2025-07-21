library(tidyverse)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

SNV_files <- list.files("data files/", recursive = T, pattern = "T1_filtered_SNVs", full.names = T)

SNV_count <- data.frame()
for(i in 1:length(SNV_files)){
  SNV_sample <- read_csv(SNV_files[i])
  SNV_sample$Name_Time <- paste(SNV_sample$Name, SNV_sample$time, sep = " ")
  SNV_sample <- SNV_sample[, c("Name_Time", "mag", "mag_coverage", "med_cov", "mag_breadth", "mag_length", "class","mutation_type", "group")]
  SNV_sample$class <- SNV_sample$class %>% str_remove("con_")
  SNV_sample$class <- SNV_sample$class %>% str_remove("pop_")
  SNV_count <- rbind(SNV_count, SNV_sample)
}

write.csv(SNV_count, "data files/T1_limited_SNV_info.csv", row.names = F)

SNV_count_MAG <- SNV_count %>% group_by(mag, Name_Time, mag_coverage, med_cov, mag_breadth, mag_length) %>% count(class)
SNV_count_MAG <- subset(SNV_count_MAG, class == "SNV")
SNV_count_MAG$SNVs_Mbp <- (SNV_count_MAG$n / SNV_count_MAG$mag_length) * 10^6
SNV_count_MAG$Treatment <- ifelse(grepl("CTRL", SNV_count_MAG$Name_Time), "Control", "GBH")
SNV_count_MAG$Name <- SNV_count_MAG$Name_Time %>% str_sub(end = -2)
SNV_count_MAG$Time <- ifelse(grepl("1", SNV_count_MAG$Name_Time), "Day 0", "Day 28")
SNV_count_MAG$Treatment_Time <- paste(SNV_count_MAG$Treatment, SNV_count_MAG$Time, sep = " ")
write.csv(SNV_count_MAG, "data files/T1_SNV_summary_MAG.csv", row.names = F)

SNS_count_MAG <- SNV_count %>% group_by(mag, Name_Time, mag_coverage, med_cov, mag_breadth, mag_length) %>% count(class)
SNS_count_MAG <- subset(SNS_count_MAG, class == "SNS")
SNS_count_MAG$SNSs_Mbp <- (SNS_count_MAG$n / SNS_count_MAG$mag_length) * 10^6
SNS_count_MAG$Treatment <- ifelse(grepl("CTRL", SNS_count_MAG$Name_Time), "Control", "GBH")
SNS_count_MAG$Name <- SNS_count_MAG$Name_Time %>% str_sub(end = -2)
SNS_count_MAG$Time <- ifelse(grepl("1", SNS_count_MAG$Name_Time), "Day 0", "Day 28")
SNS_count_MAG$Treatment_Time <- paste(SNS_count_MAG$Treatment, SNS_count_MAG$Time, sep = " ")
write.csv(SNS_count_MAG, "data files/T1_SNS_summary_MAG.csv", row.names = F)

NS_count <- SNV_count %>% group_by(mag, Name_Time, mag_coverage, med_cov, mag_breadth, mag_length, mutation_type)  %>% count(class)
NS_count_wide <- pivot_wider(NS_count, values_from = n, names_from = mutation_type)
NS_count_wide$NS_ratio <- NS_count_wide$N / NS_count_wide$S
NS_count_wide$Treatment <- ifelse(grepl("CTRL", NS_count_wide$Name_Time), "Control", "GBH")
NS_count_wide$Name <- NS_count_wide$Name_Time %>% str_sub(end = -2)
NS_count_wide$Time <- ifelse(grepl("1", NS_count_wide$Name_Time), "Day 0", "Day 28")
NS_count_wide$Treatment_Time <- paste(NS_count_wide$Treatment, NS_count_wide$Time, sep = " ")
write.csv(NS_count_wide, "data files/T1_NS_summary_MAG.csv", row.names = F)

  



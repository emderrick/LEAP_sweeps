library(tidyverse)
setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

sample_names <- read_csv("data files/chapter_1_sample_names.csv")
genome_files <- list.files("data files/T1_subsamp_inStrain", recursive = T, pattern = ".*genome_info.tsv", full.names = T)

all_mags <- data.frame()
for(i in 1:length(genome_files)){
  pond_time_mags <- read.table(genome_files[i], sep="\t", header = T)
  pond_time_mags$Sample <- genome_files[i] %>% substr(32,43)
  all_mags <- rbind(all_mags, pond_time_mags)
}

all_mags <- left_join(all_mags, sample_names)
all_mags$time <- all_mags$Pond_Time %>% substr(4,4)
all_mags <- all_mags %>% rename("mag_length" = "length")
all_mags <- all_mags %>% rename("mag" = "genome")
all_mags <- all_mags %>% rename("mag_coverage" = "coverage")
all_mags <- all_mags %>% rename("mag_breadth" = "breadth")
all_mags$mag <- all_mags$mag %>% str_sub(end = -4)
mapping_info <- all_mags[, c(1:3,5,26,27,30,31)]
write.csv(mapping_info, "data files/mag_subsamp_mapping_info.csv", row.names = F)

scaffold_files <- list.files("data files/T1_subsamp_inStrain", recursive = T, pattern = ".*scaffold_info.tsv", full.names = T)

all_scaffolds <- data.frame()
for(i in 1:length(scaffold_files)){
  pond_time_scaffolds <- read.table(scaffold_files[i], sep = "\t", header = T)
  pond_time_scaffolds$Sample <- scaffold_files[i] %>% substr(32,43)
  all_scaffolds <- rbind(all_scaffolds, pond_time_scaffolds)
} 

mag_info <- all_mags[, c(1:3,5,8,31,33,34)]

mag_scaffold_info <- read_tsv("data files/T1_refined.stb", col_names = F)
colnames(mag_scaffold_info) <-  c("scaffold", "mag")
all_scaffolds <- left_join(all_scaffolds, mag_scaffold_info)
all_scaffolds <- all_scaffolds[, c(1:3,22,23)]
all_scaffolds$mag <- all_scaffolds$mag %>% str_sub(end = -4)
mag_med_cov <- all_scaffolds %>% group_by(mag, Sample) %>% summarise(med_cov = median(coverage))
all_scaffolds <- left_join(all_scaffolds, mag_med_cov)
scaffold_mag <- left_join(all_scaffolds, mag_info)

mag_info <- left_join(mag_info, mag_med_cov)
write.csv(mag_info, "data files/T1_subsamp_mag_info.csv", row.names = F)

SNV_files <- list.files("data files/T1_subsamp_inStrain", recursive = T, pattern = ".*SNVs.tsv", full.names = T)

for(i in 1:length(SNV_files)){
  pond_time_SNV <- read.table(SNV_files[i], sep = "\t", header = T)
  Sample <- SNV_files[i] %>% substr(32,43)
  pond_time_SNV$Sample <- Sample
  SNV_mag <- left_join(pond_time_SNV, scaffold_mag)
  SNV_mag$pos_from_end <- SNV_mag$length - SNV_mag$position
  SNV_mag <- subset(SNV_mag, position > 100)
  SNV_mag <- subset(SNV_mag, pos_from_end > 100)
  SNV_mag <- SNV_mag %>% subset(!(position_coverage > mag_coverage * 3)) 
  SNV_mag <- SNV_mag %>% subset(!(position_coverage < mag_coverage / 3))
  SNV_mag$position <- str_pad(SNV_mag$position, 7, pad = "0")
  SNV_mag$group <- paste(SNV_mag$mag, SNV_mag$scaffold, SNV_mag$position, sep = "_")
  write.csv(SNV_mag, paste("data files/T1_subsamp_filt_SNVs_", Sample, ".csv", sep = ""), row.names = F)
}  

###

filt_SNV_files <- list.files("data files/", recursive = T, pattern = "subsamp_filt_SNVs", full.names = T)

SNV_count <- data.frame()
for(i in 1:length(filt_SNV_files)){
  SNV_sample <- read_csv(filt_SNV_files[i])
  SNV_sample$Name_Time <- paste(SNV_sample$Name, SNV_sample$time, sep = " ")
  SNV_sample <- SNV_sample[, c("Name_Time", "mag", "mag_coverage", "med_cov", "mag_breadth", "mag_length", "class", "mutation_type", "group")]
  SNV_sample$class <- SNV_sample$class %>% str_remove("con_")
  SNV_sample$class <- SNV_sample$class %>% str_remove("pop_")
  SNV_count <- rbind(SNV_count, SNV_sample)
}

write.csv(SNV_count, "data files/T1_subsamp_limited_SNV_info.csv", row.names = F)

mutation_counts <- SNV_count %>% group_by(mag, Name_Time, mag_coverage, med_cov, mag_breadth, mag_length, mutation_type)  %>% count(class)
mutation_counts <- pivot_wider(mutation_counts, values_from = n, names_from = mutation_type)
mutation_counts[is.na(mutation_counts)] <- 0
mutation_counts$NS_ratio <- mutation_counts$N / mutation_counts$S
mutation_counts$total <- rowSums(mutation_counts[, c("I","M","N","S","NA")], na.rm = T)
mutation_counts <- mutation_counts[, c(1:6,7,10,14,13)]
mutation_counts <- pivot_wider(mutation_counts, names_from = class, values_from = c(total, NS_ratio, N))
mutation_counts$total_SNS[is.na(mutation_counts$total_SNS)] <- 0
mutation_counts$N_SNS[is.na(mutation_counts$N_SNS)] <- 0
mutation_counts$SNSs_Mbp <- (mutation_counts$total_SNS/ mutation_counts$mag_length) * 10^6
mutation_counts$SNVs_Mbp <- (mutation_counts$total_SNV / mutation_counts$mag_length) * 10^6
mutation_counts$N_SNSs_Mbp <- (mutation_counts$N_SNS/ mutation_counts$mag_length) * 10^6
mutation_counts$N_SNVs_Mbp <- (mutation_counts$N_SNV / mutation_counts$mag_length) * 10^6
mutation_counts$Treatment <- ifelse(grepl("CTRL", mutation_counts$Name_Time), "Control", "GBH")
mutation_counts$Name <- mutation_counts$Name_Time %>% str_sub(end = -2)
mutation_counts$Time <- ifelse(grepl("1", mutation_counts$Name_Time), "Day 0", "Day 28")
mutation_counts$Treatment_Time <- paste(mutation_counts$Treatment, mutation_counts$Time, sep = " ")
write.csv(mutation_counts, "data files/T1_subsamp_SNV_summary_MAG.csv", row.names = F)


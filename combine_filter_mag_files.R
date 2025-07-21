library(tidyverse)
setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")
options(scipen=999)

sample_names <- read_csv("data files/chapter_1_sample_names.csv")

genome_files <- list.files("data files/inStrain_T1_50_output/", recursive = T, pattern = ".*genome_info.tsv", full.names = T)

all_mags <- data.frame()
for(i in 1:length(genome_files)){
  pond_time_mags <- read.table(genome_files[i], sep="\t", header = T)
  Sample <- genome_files[i] %>% substr(35,46)
  pond_time_mags$Sample <- Sample
  all_mags <- rbind(all_mags, pond_time_mags)
}

all_mags <- left_join(all_mags, sample_names)
all_mags$time <- all_mags$Pond_Time %>% substr(4,4)
all_mags <- all_mags %>% rename("mag_length" = "length")
all_mags <- all_mags %>% rename("mag" = "genome")
all_mags <- all_mags %>% rename("mag_coverage" = "coverage")
all_mags <- all_mags %>% rename("mag_breadth" = "breadth")
all_mags$mag <- all_mags$mag %>% str_sub(end = -4)

mag_info <- all_mags[, c(1:3,5,8,31,33,34)]
write.csv(all_mags, "data files/T1_mag_info.csv", row.names = F)

scaffold_files <- list.files("data files/inStrain_T1_50_output/", recursive = T, pattern = ".*scaffold_info.tsv", full.names = T)

all_scaffolds <- data.frame()
for(i in 1:length(scaffold_files)){
  pond_time_scaffolds <- read.table(scaffold_files[i], sep = "\t", header = T)
  Sample <- scaffold_files[i] %>% substr(35,46)
  pond_time_scaffolds$Sample <- Sample
  all_scaffolds <- rbind(all_scaffolds, pond_time_scaffolds)
} 

mag_scaffold_info <- read_tsv("data files/T1_50_MAGs.stb", col_names = F)
colnames(mag_scaffold_info) <-  c("scaffold", "mag")
all_scaffolds <- left_join(all_scaffolds, mag_scaffold_info)
all_scaffolds <- all_scaffolds[, c(1:3,22,23)]
all_scaffolds$mag <- all_scaffolds$mag %>% str_sub(end = -4)
mag_med_cov <- all_scaffolds %>% group_by(mag, Sample) %>% summarise(med_cov = median(coverage))
all_scaffolds <- left_join(all_scaffolds, mag_med_cov)
scaffold_mag <- left_join(all_scaffolds, mag_info)

SNV_files <- list.files("data files/inStrain_T1_50_output/", recursive = T, pattern = ".*SNVs.tsv", full.names = T)

for(i in 1:length(SNV_files)){
  pond_time_SNV <- read.table(SNV_files[i], sep = "\t", header = T)
  Sample <- SNV_files[i] %>% substr(35,46)
  pond_time_SNV$Sample <- Sample
  SNV_mag <- left_join(pond_time_SNV, scaffold_mag)
  SNV_mag$pos_from_end <- SNV_mag$length - SNV_mag$position
  SNV_mag <- subset(SNV_mag, position > 100)
  SNV_mag <- subset(SNV_mag, pos_from_end > 100)
  SNV_mag <- SNV_mag %>% subset(!(position_coverage > med_cov*3)) 
  SNV_mag <- SNV_mag %>% subset(!(position_coverage < med_cov/3))
  SNV_mag$position <- str_pad(SNV_mag$position, 7, pad = "0")
  SNV_mag$group <- paste(SNV_mag$mag, SNV_mag$scaffold, SNV_mag$position, sep = "_")
  write.csv(SNV_mag, paste("data files/T1_filtered_SNVs_", Sample, ".csv", sep = ""), row.names = F)
}  


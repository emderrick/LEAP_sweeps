library(tidyverse)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

mag_checkm <- read_tsv("data files/T1_refined_checkm.tsv")
mag_checkm <- subset(mag_checkm, Completeness >= 50 & Contamination < 10)
colnames(mag_checkm)[1] <- "mag"
mag_checkm$mag <- mag_checkm$mag %>% str_remove(".fa")

mag_tax <- read_tsv("data files/gtdbtk.bac120.summary.tsv")
colnames(mag_tax)[1] <- "mag"

mag_info <- full_join(mag_checkm, mag_tax)
mag_info <- mag_info[, c(1,12,13,15)]

mag_stats <- read.table("data files/refined_MAG_stats.tsv", header = T)
colnames(mag_stats)[1] <- "mag"
mag_stats$mag <- mag_stats$mag %>% str_remove(".fa")

mag_info <- full_join(mag_info, mag_stats[, c(1,4,5,13,18)])
mag_info$sum_len  <- gsub(",", "", mag_info$sum_len) %>% as.numeric()
mag_info$N50  <- gsub(",", "", mag_info$N50) %>% as.numeric()
colnames(mag_info) <- c("MAG", "Completeness", "Contamination", "GTDB Classification", "# Contigs", "Genome Length", "N50", "GC %")

avg_comp <- mean(mag_info$Completeness)
avg_contam <- mean(mag_info$Contamination)
avg_len <- mean(mag_info$`Genome Length`) / 1000000

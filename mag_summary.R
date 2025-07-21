library(tidyverse)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

v_mag_info <- read_tsv("data files/drep_vamb_T1_bins/data/checkM/checkM_outdir/results.tsv")
v_mag_stats <- read_csv("data files/vamb_T1_MAGs_stats.csv")
colnames(v_mag_info)[1] <- "mag"
colnames(v_mag_stats)[1] <- "mag"

v_mag_info <- right_join(v_mag_info, v_mag_stats)
v_mag_info <- v_mag_info[, c("mag", "Marker lineage", "Completeness", "Contamination", "num_seqs", "sum_len", "N50", "GC(%)")]

m_mag_info <- read_tsv("data files/checkM_T1_50_bins/data/checkM/checkM_outdir/results.tsv")
m_mag_stats <- read_csv("data files/T1_MAGs_50_stats.csv")
colnames(m_mag_info)[1] <- "mag"
colnames(m_mag_stats)[1] <- "mag"

m_mag_info <- right_join(m_mag_info, m_mag_stats)
m_mag_info <- m_mag_info[, c("mag", "Marker lineage", "Completeness", "Contamination", "num_seqs", "sum_len", "N50", "GC(%)")]

all_mag_info <- rbind(m_mag_info, v_mag_info)
all_mag_info$binner <- ifelse(grepl(".fna", all_mag_info$mag), "vamb", "metabat2")

ggplot(all_mag_info, aes(x = Completeness, y = Contamination, colour = binner))+
  geom_point()


mag_info_summary <- as.data.frame(sapply(mag_info[, c(2:6)], mean))
colnames(mag_info_summary)[1] <- "Mean"
mag_info_summary$category <- rownames(mag_info_summary)
write.csv(mag_info_summary, "data files/mag_info_summary.csv", row.names = F)

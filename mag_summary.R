library(tidyverse)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

v_mag_info <- read_tsv("data files/drep_vamb_T1_bins/data/checkM/checkM_outdir/results.tsv")
v_mag_stats <- read_csv("data files/vamb_T1_MAGs_stats.csv")
colnames(v_mag_info)[1] <- "mag"
colnames(v_mag_stats)[1] <- "mag"

v_mag_info <- right_join(v_mag_info, v_mag_stats)
v_mag_info <- v_mag_info[, c("mag", "Marker lineage", "Completeness", "Contamination", "num_seqs", "sum_len", "N50", "GC(%)")]

m_bins <- read_tsv("data files/all_T1_MAGs.stb", col_names = F)
m_bins$X2 <- m_bins$X2 %>% str_replace("\\.", "_")
m_bins$X2 <- m_bins$X2 %>% str_remove(".fa")
write_tsv(m_bins, "data files/metabat2_bins.txt")


m_mag_info <- read_tsv("data files/checkM_T1_50_bins/data/checkM/checkM_outdir/results.tsv")
m_mag_info$`Bin Id` <- m_mag_info$`Bin Id` %>% str_replace("\\.", "_")
m_mag_info$`Bin Id` <- m_mag_info$`Bin Id` %>% str_remove(".fa")

m_mag_stats <- read_csv("data files/T1_MAGs_50_stats.csv")
colnames(m_mag_info)[1] <- "mag"
colnames(m_mag_stats)[1] <- "mag"

anvi <- read_csv("data files/manually_refining_info.csv")
anvi_checkM <- full_join(anvi, m_mag_info, by = c("old_bin_name" = "Bin Id"))
anvi_checkM <- anvi_checkM[, c(1:6, 8, 18, 19)]

m_mag_info <- left_join(m_mag_info, m_mag_stats)
m_mag_info <- m_mag_info[, c("mag", "Marker lineage", "Completeness", "Contamination", "num_seqs", "sum_len", "N50", "GC(%)")]
m_mag_info <- subset(m_mag_info, Completeness >= 50)

all_mag_info <- rbind(m_mag_info, v_mag_info)
all_mag_info$binner <- ifelse(grepl(".fna", all_mag_info$mag), "vamb", "metabat2")

ggplot(all_mag_info, aes(x = Completeness, y = Contamination, colour = binner))+
  geom_point()


mag_info_summary <- as.data.frame(sapply(mag_info[, c(2:6)], mean))
colnames(mag_info_summary)[1] <- "Mean"
mag_info_summary$category <- rownames(mag_info_summary)
write.csv(mag_info_summary, "data files/mag_info_summary.csv", row.names = F)




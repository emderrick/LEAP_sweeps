library(tidyverse)
library(dplyr)
library(ggplot2)

naila_gene_files<-list.files("for_Emma/",recursive = T, pattern=".*gene_info.*",full.names = T)
naila_gene <- data.frame()
for(i in 1:length(naila_gene_files)){
  pond_time_mags <- read.table(naila_gene_files[i], sep="\t", header=T)
  timepoint <- naila_gene_files[i] %>% substr(32,40)
  pond_time_mags <- cbind(pond_time_mags, timepoint = rep(timepoint, nrow(pond_time_mags)))
  naila_gene<-rbind(naila_gene, pond_time_mags)
}

naila_gene$pond <- naila_gene$timepoint %>% substr(1,2)
naila_gene$time <- naila_gene$timepoint %>% substr(9,9)
genes_D8 <- subset(naila_gene, pond == "D8")
genes_D8$mag <- genes_D8$scaffold %>% substr(6,17)
genes_D8_C4_MAG_00010 <- subset(genes_D8, mag == "C4_MAG_00010")
T1_gene <- subset(genes_D8_C4_MAG_00010, time == 1)
T1_gene_cov <- median(T1_gene$coverage)
T2_gene <- subset(genes_D8_C4_MAG_00010, time == 2)
T2_gene_cov <- median(T2_gene$coverage)
T3_gene <- subset(genes_D8_C4_MAG_00010, time == 3)
T3_gene_cov <- median(T3_gene$coverage)
genes_D8_C4_MAG_00010$gene_cov <- with(genes_D8_C4_MAG_00010, ifelse(time == 1, T1_gene_cov, ifelse(time == 2, T2_gene_cov, T3_gene_cov)))
genes_D8_C4_MAG_00010$mag_cov <- with(genes_D8_C4_MAG_00010, ifelse(time == 1, 45, ifelse(time == 2, 103, 6)))
genes_D8_C4_MAG_00010$gene_cov_norm <- genes_D8_C4_MAG_00010$coverage / genes_D8_C4_MAG_00010$gene_cov
D8_C4_MAG_00010 <- genes_D8_C4_MAG_00010[, c("gene", "gene_cov_norm", "time")]
D8_C4_MAG_00010_wide <- pivot_wider(D8_C4_MAG_00010, values_from = gene_cov_norm, names_from = time)
D8_C4_MAG_00010_wide <- mutate_all(D8_C4_MAG_00010_wide, ~replace_na(., 0))
D8_C4_MAG_00010_wide$cov_dif <- abs(D8_C4_MAG_00010_wide$'3' - D8_C4_MAG_00010_wide$'1')
cov_dif_mean <- mean(D8_C4_MAG_00010_wide$cov_dif)
cov_dif_SD <- sd(D8_C4_MAG_00010_wide$cov_dif)
cutoff <- cov_dif_mean + 2*cov_dif_SD

D8_C4_MAG_00010_long <- pivot_longer(D8_C4_MAG_00010_wide, cols = c(2:4), values_to = "gene_cov_norm", names_to = "time" )
ggplot(subset(D8_C4_MAG_00010_long, cov_dif >= 0.7 & gene_cov_norm < 4), aes(x = time, y = gene_cov_norm, group = gene))+
  geom_point()+
  geom_line()

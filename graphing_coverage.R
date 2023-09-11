library(tidyverse)
library(ggplot2)
library(dplyr)
library(viridis)
library(stringr)
library(cowplot)

mag_list <- list("I4_MAG_00006", "I4_MAG_00065", "L2_MAG_00052", "L3_MAG_00058", "L4_MAG_00099",
                 "L7_MAG_00020", "L7_MAG_00028", "L7_MAG_00043", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")
mag_labs <- c(I4_MAG_00006 = "SJAQ100 sp016735685 assembled from Control B", I4_MAG_00065 = "Roseomonas sp. assembled from Control B", L2_MAG_00052 = "Erythrobacter sp. assembled from GBH A", 
              L3_MAG_00058 = "Prosthecobacter sp. assembled from Control C", L4_MAG_00099 = "Bosea sp001713455 assembled from Control D", L7_MAG_00020 = "Sphingorhabdus_B sp. assembled from GBH C",
              L7_MAG_00028 = "SYFN01 sp. assembled from GBH C", L7_MAG_00043 = "Luteolibacter sp. assembled from GBH C", L8_MAG_00011 = "UBA953 sp. assembled from GBH D", 
              L8_MAG_00019 = "UA16 family assembled from GBH D", L8_MAG_00042 = "UBA4660 sp. assembled from GBH D")

#load mag snv info
all_snv <- read_csv("filtered_ANI_95_mag_SNVs.csv")

mag_scaf_cov <- all_snv %>% group_by(scaffold, new_name, length, coverage, mag) %>% summarize(SNV_SNS_tot = sum(number_divergent))


all_MAG_scaf_cov <-  ggplot(mag_scaf_cov, aes(x = coverage, y = log10((SNV_SNS_tot/length)*10^6), colour = new_name)) + 
  geom_point()+
  scale_colour_viridis(discrete = T)+
  labs(y ="log10 SNVs / MBp", x="Coverage (x)", colour= "Pond") +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  facet_wrap(~mag, nrow = 4, ncol = 4, scales = "free_x", labeller = labeller(mag = mag_labs))

save_plot("MAG_scaf_cov_SNV.jpeg", all_MAG_scaf_cov, ncol = 4, nrow = 4, dpi = 300)


mag_cov <- all_snv %>% group_by(full_group, mag, new_name, mag_coverage, breadth.y, mag_length) %>% summarize(SNV_SNS_tot = sum(number_divergent))
mag_cov <- mag_cov %>% rename(mag_breadth = breadth.y)
mag_cov$snvs_per_mbp <- (mag_cov$SNV_SNS_tot / mag_cov$mag_length) *10^6
mag_cov$percent_id_to_ref <- (1 - mag_cov$SNV_SNS_tot / mag_cov$mag_length) *100
write.csv(mag_cov, "mag_cov_summary_info.csv")

all_MAG_cov <-  ggplot(mag_cov, aes(x = mag_coverage, y = log10((SNV_SNS_tot/mag_length)*10^6), colour = new_name)) + 
  geom_point()+
  scale_colour_viridis(discrete = T)+
  labs(y = "log10 SNVs / MBp", x = "Coverage (x)", colour = "Pond") +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  facet_wrap(~mag, nrow = 4, ncol = 4, scales = "free_x", labeller = labeller(mag = mag_labs))

save_plot("MAG_cov_SNV.jpeg", all_MAG_cov, ncol = 4, nrow = 4, dpi = 300)

library(tidyverse)
library(ggplot2)
library(cowplot)
library(reporter)

mag_list <- list("I4_MAG_00006", "I4_MAG_00065", "L2_MAG_00052", "L3_MAG_00058", "L4_MAG_00099",
                 "L7_MAG_00020", "L7_MAG_00028", "L7_MAG_00043", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")

mag_labs <- c(I4_MAG_00006 = "SJAQ100 sp016735685", I4_MAG_00065 = "Roseomonas sp.", L2_MAG_00052 = "Erythrobacter sp.", 
              L3_MAG_00058 = "Prosthecobacter sp.", L4_MAG_00099 = "Bosea sp001713455", L7_MAG_00020 = "Sphingorhabdus_B sp.",
              L7_MAG_00028 = "SYFN01 sp.", L7_MAG_00043 = "Luteolibacter sp.", L8_MAG_00011 = "UBA953 sp.", 
              L8_MAG_00019 = "UA16", L8_MAG_00042 = "UBA4660 sp.")

all_snv <- read_csv("filtered_ANI_95_mag_SNVs.csv") ###CHANGE WHAT FILE I"M USING####
all_snv <- subset(all_snv, new_time == 2)
mag_scaf_cov <- all_snv %>% group_by(scaffold, name, length, coverage, mag) %>% summarize(SNV_SNS_tot = sum(number_divergent))
mag_scaf_cov$mag_order = factor(mag_scaf_cov$mag, levels=c('I4_MAG_00006', 'I4_MAG_00065', 'L3_MAG_00058', 'L7_MAG_00020', 'L8_MAG_00011', 'L8_MAG_00019',
                                                           'L2_MAG_00052', 'L4_MAG_00099',  'L7_MAG_00028', 'L7_MAG_00043', 'L8_MAG_00042'))

all_MAG_scaf_cov <- ggplot(mag_scaf_cov, aes(x = coverage, y = log10((SNV_SNS_tot/length)*10^6), colour = name)) + 
  geom_point()+
  scale_colour_manual(values = c("#002C3D", "#005F73", "#0A9396", "#94D2BD", "#EE9B00", "#CA6702", "#BB3E03", "#AE2012", "#9B2226"))+
  labs(y = paste("log", {subsc('10')}, " SNVs / Mbp"), x="Coverage (X)", colour= "Pond") +
  theme_classic()+
  theme(text = element_text(size = 30, colour = 'black'),
        axis.text = element_text(colour = "black"),
        axis.title = element_text(size = 35), 
        axis.ticks.length = unit(.4, "cm"),
        axis.ticks = element_line(linewidth = 1.5, colour = "black"),
        axis.text.x = element_text(vjust = -1),
        axis.title.x = element_text(vjust = -1, margin = margin(t = 20, r = 0, b = 20, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.line = element_line(linewidth = 1),
        strip.text.x.top = element_text(size = 30, margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), face = "bold"),
        strip.background = element_rect(linewidth = 2),
        legend.title = element_text(size = 35),
        legend.position = "bottom",
        legend.text = element_text(size = 30),
        legend.background = element_blank(),
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.spacing.y = unit(1.5, "lines"))+
  scale_y_continuous(limits=c(0,5))+
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 10)))+
  facet_wrap(~mag_order, nrow = 2, ncol = 6, scales = "free", labeller = labeller(mag_order = mag_labs))

save_plot("MAG_scaf_cov_SNV.jpeg", all_MAG_scaf_cov, base_height = 7, base_width = 6, ncol = 6, nrow = 2, dpi = 300, limitsize = F)

#summary table of MAG info
mag_cov <- all_snv %>% group_by(mag, name, mag_length, mag_coverage, breadth.y, detected_scaffolds, conANI_reference.y) %>% summarize(SNV_SNS_tot = sum(number_divergent))
mag_cov <- mag_cov %>% rename(mag_breadth = breadth.y)
mag_cov$snv_sns__per_mbp <- (mag_cov$SNV_SNS_tot / mag_cov$mag_length) *10^6
mag_cov$mag_name <- with(mag_cov, ifelse(mag_cov$mag == "I4_MAG_00006", "SJAQ100 sp016735685", NA)) 
mag_cov$mag_name <- with(mag_cov, ifelse(mag_cov$mag == "I4_MAG_00065", "Roseomonas sp.", mag_name)) 
mag_cov$mag_name <- with(mag_cov, ifelse(mag_cov$mag == "L2_MAG_00052", "Erythrobacter sp.", mag_name)) 
mag_cov$mag_name <- with(mag_cov, ifelse(mag_cov$mag == "L3_MAG_00058", "Prosthecobacter sp.", mag_name)) 
mag_cov$mag_name <- with(mag_cov, ifelse(mag_cov$mag == "L4_MAG_00099", "Bosea sp001713455", mag_name)) 
mag_cov$mag_name <- with(mag_cov, ifelse(mag_cov$mag == "L7_MAG_00020", "Sphingorhabdus_B sp.", mag_name)) 
mag_cov$mag_name <- with(mag_cov, ifelse(mag_cov$mag == "L7_MAG_00028", "SYFN01 sp.", mag_name)) 
mag_cov$mag_name <- with(mag_cov, ifelse(mag_cov$mag == "L7_MAG_00043", "Luteolibacter sp.", mag_name)) 
mag_cov$mag_name <- with(mag_cov, ifelse(mag_cov$mag == "L8_MAG_00011", "UBA953 sp.", mag_name)) 
mag_cov$mag_name <- with(mag_cov, ifelse(mag_cov$mag == "L8_MAG_00019", "UA16", mag_name)) 
mag_cov$mag_name <- with(mag_cov, ifelse(mag_cov$mag == "L8_MAG_00042", "UBA4660 sp.", mag_name)) 
write.csv(mag_cov, "mag_cov_summary_info.csv", row.names = F)

library(tidyverse)
library(ggplot2)
library(cowplot)

mag_labs <- c(I4_MAG_00006 = "Burkholderiaceae 1", I4_MAG_00065 = "Roseomonas_A", L2_MAG_00052 = "Erythrobacter", 
              L3_MAG_00058 = "Prosthecobacter", L4_MAG_00099 = "Bosea sp001713455", L7_MAG_00020 = "Sphingorhabdus_B",
              L7_MAG_00028 = "Burkholderiaceae 2", L7_MAG_00043 = "Luteolibacter", L8_MAG_00011 = "Verrucomicrobiae", 
              L8_MAG_00019 = "Flavobacteriales 1", L8_MAG_00042 = "Flavobacteriales 2")

all_snv <- read_csv("all_MAG_SNVs_med_Apr9.csv")
all_snv <- subset(all_snv, new_time == 2)
all_snv$test <- 1
mag_scaf_sum <- all_snv %>% group_by(mag, scaffold, new_name, length, coverage.y) %>% summarize(SNV_SNS_tot = sum(test))

mag_scaf_cov_sum$mag_order = factor(mag_scaf_cov_sum$mag, levels=c('I4_MAG_00006','L7_MAG_00028','L8_MAG_00011', 'L8_MAG_00019', 'L8_MAG_00042',
                                                           'I4_MAG_00065', 'L3_MAG_00058', 'L7_MAG_00020',  'L7_MAG_00043',
                                                           'L2_MAG_00052', 'L4_MAG_00099'))

all_MAG_scaf_cov <- ggplot(mag_scaf_cov_sum, aes(x = coverage.y, y = log10((SNV_SNS_tot/length)*10^6), colour = new_name)) + 
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

save_plot("MAG_scaf_cov_SNV_sub.jpeg", all_MAG_scaf_cov, base_height = 7, base_width = 6, ncol = 6, nrow = 2, dpi = 300, limitsize = F)
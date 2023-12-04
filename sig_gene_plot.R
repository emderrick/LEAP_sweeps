library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

sig_summary <- read_csv("significant_genes_summary.csv")
sweep_mags <- list("I4_MAG_00006", "I4_MAG_00065", "L3_MAG_00058", "L7_MAG_00020", "L8_MAG_00011", "L8_MAG_00019")

sig_summary$sweep <- with(sig_summary, ifelse(mag %in% sweep_mags, "Potential Sweep", "No sweep"))

sig_summary_long <- pivot_longer(sig_summary, cols = 2:6, names_to = "method", values_to = "num_genes")
sig_summary_long <- subset(sig_summary_long, (method == "strict_genes" | method == "allele_frequency_genes" | method == "strict_allele_frequency"))
sig_summary_long$method <- factor(sig_summary_long$method, levels = c("allele_frequency_genes", "strict_genes", "strict_allele_frequency"))

sig_summary_long$mag_name <- with(sig_summary_long, ifelse(sig_summary_long$mag == "I4_MAG_00006", "SJAQ100 sp016735685", NA)) 
sig_summary_long$mag_name <- with(sig_summary_long, ifelse(sig_summary_long$mag == "I4_MAG_00065", "Roseomonas sp.", mag_name)) 
sig_summary_long$mag_name <- with(sig_summary_long, ifelse(sig_summary_long$mag == "L2_MAG_00052", "Erythrobacter sp.", mag_name)) 
sig_summary_long$mag_name <- with(sig_summary_long, ifelse(sig_summary_long$mag == "L3_MAG_00058", "Prosthecobacter sp.", mag_name)) 
sig_summary_long$mag_name <- with(sig_summary_long, ifelse(sig_summary_long$mag == "L4_MAG_00099", "Bosea sp001713455", mag_name)) 
sig_summary_long$mag_name <- with(sig_summary_long, ifelse(sig_summary_long$mag == "L7_MAG_00020", "Sphingorhabdus_B sp.", mag_name)) 
sig_summary_long$mag_name <- with(sig_summary_long, ifelse(sig_summary_long$mag == "L7_MAG_00028", "SYFN01 sp.", mag_name)) 
sig_summary_long$mag_name <- with(sig_summary_long, ifelse(sig_summary_long$mag == "L7_MAG_00043", "Luteolibacter sp.", mag_name)) 
sig_summary_long$mag_name <- with(sig_summary_long, ifelse(sig_summary_long$mag == "L8_MAG_00011", "UBA953 sp.", mag_name)) 
sig_summary_long$mag_name <- with(sig_summary_long, ifelse(sig_summary_long$mag == "L8_MAG_00019", "UA16", mag_name)) 
sig_summary_long$mag_name <- with(sig_summary_long, ifelse(sig_summary_long$mag == "L8_MAG_00042", "UBA4660 sp.", mag_name)) 

sig_gene <- ggplot(sig_summary_long, aes(x = method, y = num_genes))+
  geom_boxplot(fill = "grey90", outlier.shape = NA)+
  geom_jitter(aes(colour = sweep), width = 0.2, size = 2)+
  scale_colour_manual(values = c("#9B2226", "#005F73"))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  labs(y = "Number of Genes", x = NULL, colour = NULL)+
  theme_classic()+
  theme(text = element_text(size = 20, colour = 'black'),
        axis.text = element_text(colour = "black"))

save_plot("sig_gene_summary.jpeg", sig_gene, base_height = 8, base_width = 10, dpi = 300, limitsize = F)

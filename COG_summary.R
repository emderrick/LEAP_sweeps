library(tidyverse)
library(dplyr)
library(ggplot2)
library(cowplot)

COG_threshold <- read_csv("COG_enrichment_eggnog_threshold_significant_genes.csv")
COG_loose <- read_csv("COG_enrichment_eggnog_significant_genes_loose.csv")
COG_strict <- read_csv("COG_enrichment_eggnog_significant_genes_strict.csv")
COG_increase <- read_csv("COG_enrichment_eggnog_gene_cov_sig_increase.csv")
COG_decrease <- read_csv("COG_enrichment_eggnog_gene_cov_sig_decrease.csv")

# COG_threshold <- read_csv("COG_enrichment_eggnog_sig_genes_threshold_C1.csv")
# COG_loose <- read_csv("COG_enrichment_eggnog_sig_genes_loose_C1.csv")
# COG_strict <- read_csv("COG_enrichment_eggnog_sig_genes_strict_C1.csv")
# COG_increase <- read_csv("COG_enrichment_eggnog_sig_genes_increase_C1.csv")
# COG_decrease <- read_csv("COG_enrichment_eggnog_sig_genes_decrease_C1.csv")
# 
# COG_threshold <- read_csv("COG_enrichment_eggnog_sig_genes_threshold_C2.csv")
# COG_loose <- read_csv("COG_enrichment_eggnog_sig_genes_loose_C2.csv")
# COG_strict <- read_csv("COG_enrichment_eggnog_sig_genes_strict_C2.csv")
# COG_increase <- read_csv("COG_enrichment_eggnog_sig_genes_increase_C2.csv")
# COG_decrease <- read_csv("COG_enrichment_eggnog_sig_genes_decrease_C2.csv")
# 
# COG_threshold <- read_csv("COG_enrichment_eggnog_sig_genes_threshold_sweep.csv")
# COG_loose <- read_csv("COG_enrichment_eggnog_sig_genes_loose_sweep.csv")
# COG_strict <- read_csv("COG_enrichment_eggnog_sig_genes_strict_sweep.csv")
# COG_increase <- read_csv("COG_enrichment_eggnog_sig_genes_increase_sweep.csv")
# COG_decrease <- read_csv("COG_enrichment_eggnog_sig_genes_decrease_sweep.csv")
# 
# COG_threshold <- read_csv("COG_enrichment_eggnog_sig_genes_threshold_no_sweep.csv")
# COG_loose <- read_csv("COG_enrichment_eggnog_sig_genes_loose_no_sweep.csv")
# COG_strict <- read_csv("COG_enrichment_eggnog_sig_genes_strict_no_sweep.csv")
# COG_increase <- read_csv("COG_enrichment_eggnog_sig_genes_increase_no_sweep.csv")
# COG_decrease <- read_csv("COG_enrichment_eggnog_sig_genes_decrease_no_sweep.csv")

COG_threshold$method <- "Allele Frequency Change"
COG_loose$method <- "SNV Count Loose"
COG_strict$method <- "SNV Count Strict"
COG_increase$method <- "Gene Copy Increase"
COG_decrease$method <- "Gene Copy Decrease"

all_COG <- rbind(COG_threshold, COG_loose, COG_strict, COG_increase, COG_decrease)
all_COG <- all_COG[, c(1, 6:9)]
all_COG$method <- factor(all_COG$method, levels = c("Allele Frequency Change", "SNV Count Loose", "SNV Count Strict", "Gene Copy Increase", "Gene Copy Decrease"))
all_COG$OR_sig <- with(all_COG, ifelse(OR < 1 & p < 0.05, "p < 0.05 & depleted", "p > 0.05"))
all_COG$OR_sig <- with(all_COG, ifelse(OR > 1 & p < 0.05, "p < 0.05 & enriched", OR_sig))
all_COG$category <-  factor(all_COG$category, levels = c("J", "K", "L", "D", "V", "T", "M", "N", "W", "U", "O", "C", "G", "E", "F", "H", "I", "P", "Q", "X", "R", "S"))
all_COG_sig <- subset(all_COG, OR_sig != "p > 0.05" & fdr < 0.1)
all_COG_sig$fdr_sig <- "FDR < 0.1"


COG_summary <- ggplot(all_COG, aes(y = method, x = category, fill = OR_sig))+
  geom_tile(height=0.95, width=0.85)+
  geom_tile(data = all_COG_sig, aes(colour = fdr_sig), height=0.9, width=0.85, lwd = 2)+
  scale_fill_manual(values = c("#0A9396", "#EE9B00", "grey92"))+
  scale_colour_manual(values = c("black"))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 12, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 12, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.margin =  margin(t = 0, r = 0, b = 0, l = 10),
        legend.text = element_text(size = 16),
        legend.position = "right",
        legend.box = "veritcal",
        legend.key.size = unit(1, "cm"),
        plot.margin = unit(c(1, 1, 1, 1), "cm"))+
    labs(fill = "P-value significant")+
    scale_y_discrete(expand = expansion(mult = c(0.025, 0.05)))+
    scale_x_discrete(expand = c(0, 0.5))+
    guides(colour = guide_legend(order = 1))

save_plot("COG_summary_plot.jpeg", base_height = 7, base_width = 16, COG_summary, dpi = 200)


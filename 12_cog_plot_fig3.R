library(tidyverse)
library(ggpubr)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

cog_files <- list.files("data files", pattern = "enrichment.csv", full.names = T)

for(i in 1: length(cog_files)){
  file_name <- cog_files[i] %>% str_remove("data files/")
  file_name <- file_name %>% str_remove("_enrichment.csv")
  
  cog_results <- read_csv(cog_files[i])
  cog_results$category <-  factor(cog_results$category, levels = c("J", "K", "L", "D", "V", "T", "M", "N", "W", "U", "O", "C", "G", "E", "F", "H", "I", "P", "Q", "X", "R", "S"))
  cog_results$OR_sig <- with(cog_results, ifelse(OR < 1 & fdr < 0.05, "significantly depleted", "FDR > 0.05"))
  cog_results$OR_sig <- with(cog_results, ifelse(OR > 1 & fdr < 0.05, "significantly enriched", OR_sig))
  cog_results$logp <- -log10(cog_results$p)

  COG_summary_plot <- ggplot(cog_results, aes(x = category, y = -log10(p), fill = OR_sig))+
    geom_col()+
    scale_fill_manual(values = c("FDR > 0.05" = "grey92", "significantly depleted" = "#0A9396", "significantly enriched" ="#EE9B00"))+
    theme_bw()+
    theme(text = element_text(size = 18), axis.text = element_text(colour = "black"))+
    labs(x = "", y = "-log10(p-value)", fill = "")

  ggsave(paste("figures/", file_name, "_summary_plot.pdf", sep = ""), COG_summary_plot, height = 6, width = 14)

}

snv_freq <- read_csv("data files/nonsyn_allele_shifts_07_COG_enrichment.csv")
snv_freq$category <-  factor(snv_freq$category, levels = c("J", "K", "L", "D", "V", "T", "M", "N", "W", "U", "O", "C", "G", "E", "F", "H", "I", "P", "Q", "X", "R", "S"))
snv_freq$OR_sig <- with(snv_freq, ifelse(OR < 1 & fdr < 0.05, "significantly depleted", "FDR > 0.05"))
snv_freq$OR_sig <- with(snv_freq, ifelse(OR > 1 & fdr < 0.05, "significantly enriched", OR_sig))
snv_freq$logp <- -log10(snv_freq$p)

gbh_loss <- read_csv("data files/GBH_gene_loss_COG_enrichment.csv")
gbh_loss$category <-  factor(gbh_loss$category, levels = c("J", "K", "L", "D", "V", "T", "M", "N", "W", "U", "O", "C", "G", "E", "F", "H", "I", "P", "Q", "X", "R", "S"))
gbh_loss$OR_sig <- with(gbh_loss, ifelse(OR < 1 & fdr < 0.05, "significantly depleted", "FDR > 0.05"))
gbh_loss$OR_sig <- with(gbh_loss, ifelse(OR > 1 & fdr < 0.05, "significantly enriched", OR_sig))
gbh_loss$logp <- -log10(gbh_loss$p)

snv_plot <- ggplot(snv_freq, aes(x = category, y = -log10(p), fill = OR_sig))+
  geom_col()+
  scale_fill_manual(values = c("FDR > 0.05" = "grey92", "significantly depleted" = "#0A9396", "significantly enriched" ="#EE9B00"))+
  theme_classic()+
  theme(title = element_text(size = 10),
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 12, colour = "black"))+
  ylim(0, 10)+
  labs(x = "", y = "-log10(p-value)", fill = "", title = "SNV frequency shift")

gbh_loss_plot <- ggplot(gbh_loss, aes(x = category, y = -log10(p), fill = OR_sig))+
  geom_col()+
  scale_fill_manual(values = c("FDR > 0.05" = "grey92", "significantly depleted" = "#0A9396", "significantly enriched" ="#EE9B00"))+
  theme_classic()+
  theme(title = element_text(size = 10),
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 12, colour = "black"))+
  ylim(0, 10)+
  labs(x = "", y = "-log10(p-value)", fill = "", title = "Gene loss")

cog_plot <- ggarrange(snv_plot, gbh_loss_plot, common.legend = T, legend = "top", nrow = 2, labels = c("A", "B"))
ggsave("figures/snv_freq_gbh_loss_fig3.pdf", cog_plot, units = "cm", width = 15, height = 14)




library(tidyverse)
library(ggplot2)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

cog_files <- list.files("data files", pattern = "enrichment.csv", full.names = T)

for(i in 1: length(cog_files)){
  file_name <- cog_files[i] %>% str_remove("data files/")
  file_name <- file_name %>% str_remove("_enrichment.csv")
  
  cog_results <- read_csv(cog_files[i])
  cog_results$category <-  factor(cog_results$category, levels = c("J", "K", "L", "D", "V", "T", "M", "N", "W", "U", "O", "C", "G", "E", "F", "H", "I", "P", "Q", "X", "R", "S"))
  cog_results$OR_sig <- with(cog_results, ifelse(OR < 1 & fdr < 0.05, "significantly depleted", "p > 0.05"))
  cog_results$OR_sig <- with(cog_results, ifelse(OR > 1 & fdr < 0.05, "significantly enriched", OR_sig))
  cog_results$logp <- -log10(cog_results$p)

  COG_summary_plot <- ggplot(cog_results, aes(x = category, y = -log10(p), fill = OR_sig))+
    geom_col()+
    scale_fill_manual(values = c("grey92", "#0A9396", "#EE9B00"))+
    theme_bw()+
    theme(text = element_text(size = 18), axis.text = element_text(colour = "black"))+
    labs(x = "", y = "-log10(p-value)", fill = "")

  ggsave(paste("figures/", file_name, "_summary_plot.pdf", sep = ""), COG_summary_plot, height = 6, width = 14)

}
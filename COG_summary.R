library(tidyverse)
library(dplyr)
library(ggplot2)
library(cowplot)

test_groups <- list(all = list("COG_enrich_SNV_decrease_genes.csv", "COG_enrich_parallel_decrease_genes.csv", "COG_enrich_threshold_significant_genes_all_subsamp.csv", "COG_enrich_gene_cov_sig_increase_all_subsamp.csv", "COG_enrich_gene_cov_sig_decrease_all_subsamp.csv"),
                    class_1 = list("COG_enrich_sig_genes_loose_C1_subsamp.csv", "COG_enrich_sig_genes_strict_C1_subsamp.csv", "COG_enrich_sig_genes_threshold_C1_subsamp.csv", "COG_enrich_sig_genes_increase_C1_subsamp.csv", "COG_enrich_sig_genes_decrease_C1_subsamp.csv"),
                     class_2 = list( "COG_enrich_sig_genes_loose_C2_subsamp.csv", "COG_enrich_sig_genes_strict_C2_subsamp.csv", "COG_enrich_sig_genes_threshold_C2_subsamp.csv", "COG_enrich_sig_genes_increase_C2_subsamp.csv", "COG_enrich_sig_genes_decrease_C2_subsamp.csv"))

output_names <- names(test_groups)

for(name in output_names){
  COG_loose <- read_csv(test_groups[[name]][[1]])
  COG_strict <- read_csv(test_groups[[name]][[2]])
  COG_threshold <-  read_csv(test_groups[[name]][[3]])
  COG_increase <-  read_csv(test_groups[[name]][[4]])
  COG_decrease <-  read_csv(test_groups[[name]][[5]])
  
  COG_loose$method <- "SNV Number Decrease"
  COG_strict$method <- "Parallel SNV Number Decrease"
  COG_threshold$method <- "Allele Frequency Change"
  COG_increase$method <- "Gene Copy Number Increase"
  COG_decrease$method <- "Gene Copy Number Decrease"
  
  all_COG <- rbind(COG_loose, COG_strict, COG_threshold, COG_increase, COG_decrease)
  all_COG$OR_sig <- with(all_COG, ifelse(OR < 1 & p < 0.05, "p < 0.05 & depleted", "p > 0.05"))
  all_COG$OR_sig <- with(all_COG, ifelse(OR > 1 & p < 0.05, "p < 0.05 & enriched", OR_sig))
  write.csv(all_COG, paste(name, "COG.csv", sep ="_"), row.names = F)

}

all <- read_csv("all_COG.csv")
all$method <- factor(all$method, levels = c("Allele Frequency Change", "Parallel SNV Number Decrease", "SNV Number Decrease", "Gene Copy Number Increase", "Gene Copy Number Decrease"))
all$category <-  factor(all$category, levels = c("J", "K", "L", "D", "V", "T", "M", "N", "W", "U", "O", "C", "G", "E", "F", "H", "I", "P", "Q", "X", "R", "S"))
all_sig <- subset(all, OR_sig != "p > 0.05" & fdr < 0.1)
all_sig$fdr_sig <- "FDR < 0.1"

COG_summary_all <- ggplot(all, aes(y = method, x = category, fill = OR_sig))+
    geom_tile(height=0.95, width=0.85)+
    geom_tile(data = all_sig, aes(colour = fdr_sig), height=0.9, width=0.85, lwd = 2)+
    scale_fill_manual(values = c("#0A9396", "#EE9B00", "grey92"))+
    scale_colour_manual(values = c("black"))+
    theme_classic()+
    theme(axis.text.x = element_text(size = 16, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.text.y = element_text(size = 16, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          legend.title = element_blank(),
          legend.margin =  margin(t = 0, r = 0, b = 0, l = 10),
          legend.text = element_text(size = 16),
          legend.position = "right",
          legend.box = "veritcal",
          legend.key.size = unit(1, "cm"),
          plot.margin = unit(c(1, 1, 1, 1), "cm"))+
    scale_y_discrete(expand = expansion(mult = c(0.025, 0.05)))+
    scale_x_discrete(expand = c(0, 0.5))+
    guides(colour = guide_legend(order = 1))

save_plot("all_COG_summary_plot.jpeg", COG_summary_all, base_height = 6.5, base_width = 16, dpi = 200)


class_1 <- read_csv("class_1_COG.csv")
class_1$class <- "Class I - Sensitive"
class_2 <- read_csv("class_2_COG.csv")
class_2$class <- "Class II - Resistant"
all_class <- rbind(class_1, class_2)
all_class$method <- factor(all_class$method, levels = c("Allele Frequency Change", "Parallel SNV Number Decrease", "SNV Number Decrease", "Gene Copy Number Increase", "Gene Copy Number Decrease"))
all_class$category <-  factor(all_class$category, levels = c("J", "K", "L", "D", "V", "T", "M", "N", "W", "U", "O", "C", "G", "E", "F", "H", "I", "P", "Q", "X", "R", "S"))
all_class_sig <- subset(all_class, OR_sig != "p > 0.05" & fdr < 0.1)
all_class_sig$fdr_sig <- "FDR < 0.1"

class_COG <- ggplot(all_class, aes(y = method, x = category, fill = OR_sig))+
  geom_tile(height=0.95, width=0.85)+
  geom_tile(data = all_class_sig, aes(colour = fdr_sig), height=0.9, width=0.85, lwd = 2)+
  scale_fill_manual(values = c("#0A9396", "#EE9B00", "grey92"))+
  scale_colour_manual(values = c("black"))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 16, colour = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 16, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        strip.text.x.top = element_text(size = 20, margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"), face = "bold"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.margin =  margin(t = 0, r = 0, b = 0, l = 10),
        legend.text = element_text(size = 16),
        panel.spacing.y = unit(1, "lines"),
        legend.position = "right",
        legend.box = "veritcal",
        legend.key.size = unit(1, "cm"),
        plot.margin = unit(c(1, 1, 1, 1), "cm"))+
  scale_y_discrete(expand = expansion(mult = c(0.025, 0.05)))+
  scale_x_discrete(expand = c(0, 0.5))+
  guides(colour = guide_legend(order = 1))+
  facet_wrap(~class, ncol = 1, scales = "free")

save_plot("class_grouped_COG.jpeg", base_height = 14, base_width = 16, class_COG, dpi = 400)

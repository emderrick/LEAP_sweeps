library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(cowplot)

all_MAG_snvs <- read_csv("all_MAG_SNVs_med_Aug15.csv")
all_MAG_snvs <- subset(all_MAG_snvs, str_detect(new_name, "T2"))

MAG_snvs <- all_MAG_snvs[, c('mag', 'gene', 'mutation_type', 'new_name')] %>% subset(is.na(mutation_type) == F)
MAG_snvs$mutation <- 1
MAG_NS <- MAG_snvs_all %>% group_by(mag, new_name) %>% summarize(N = sum(mutation_type=="N"), S = sum(mutation_type=="S"), I = sum(mutation_type=="I"), M = sum(mutation_type=="M"))
MAG_NS$total <- MAG_NS$N + MAG_NS$S + MAG_NS$I + MAG_NS$M
MAG_NS$graph_name <- str_sub(MAG_NS$new_name, end = -6) 
MAG_NS$graph_name <- gsub('Control', 'Ctrl', MAG_NS$graph_name)
MAG_NS_long <- pivot_longer(MAG_NS, cols = c("N", "S", "I", "M"), names_to = "mutation_type", values_to = "count")

NS_plot <- ggplot(MAG_NS_long, aes(x = graph_name, y = count, fill = mutation_type))+
  geom_bar(stat = "identity")+
  theme_classic()+
  labs(y = "Fraction of SNVs", x = "pond", legend = "Mutation Type")+
  facet_wrap(~mag, nrow = 4, ncol = 4, scales="free")

save_plot("NS_MAG_plot.jpeg", NS_plot, ncol = 4, nrow = 4, dpi = 300)


# look at N:S ratio by gene
MAG_dnds_by_gene <- MAG_snvs %>% group_by(gene, new_name) %>% summarize(N = sum(mutation_type=="N"), S = sum(mutation_type=="S")) %>% 
  pivot_longer(names_to = "mutation_type", values_to = "count")
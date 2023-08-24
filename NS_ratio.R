library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(cowplot)

mag_labs <- c(I4_MAG_00006 = "SJAQ100 sp016735685 assembled from Control B", I4_MAG_00065 = "Roseomonas sp. assembled from Control B", L2_MAG_00052 = "Erythrobacter sp. assembled from GBH A", 
                 L3_MAG_00058 = "Prosthecobacter sp. assembled from Control C", L4_MAG_00099 = "Bosea sp001713455 assembled from Control D", L7_MAG_00020 = "Sphingorhabdus_B sp. assembled from GBH C",
                 L7_MAG_00028 = "SYFN01 sp. assembled from GBH C", L7_MAG_00043 = "Luteolibacter sp. assembled from GBH C", L8_MAG_00011 = "UBA953 sp. assembled from GBH D", 
                 L8_MAG_00019 = "UA16 family assembled from GBH D", L8_MAG_00042 = "UBA4660 sp. assembled from GBH D")

all_MAG_snvs <- read_csv("all_MAG_SNVs_med_Aug15.csv")
all_MAG_snvs <- subset(all_MAG_snvs, str_detect(new_name, "T2"))
MAG_snvs <- all_MAG_snvs[, c('mag', 'gene', 'mutation_type', 'new_name')] %>% subset(is.na(mutation_type) == F)
MAG_snvs$mutation <- 1

MAG_NS <- MAG_snvs %>% group_by(mag, new_name) %>% summarize(N = sum(mutation_type=="N"), S = sum(mutation_type=="S"), I = sum(mutation_type=="I"), M = sum(mutation_type=="M"))
MAG_NS$total <- MAG_NS$N + MAG_NS$S + MAG_NS$I + MAG_NS$M
MAG_NS$graph_name <- str_sub(MAG_NS$new_name, end = -6) 
MAG_NS$graph_name <- gsub('Control', 'Ctrl', MAG_NS$graph_name)
MAG_NS_long <- pivot_longer(MAG_NS, cols = c("N", "S", "I", "M"), names_to = "mutation_type", values_to = "count")

NS_plot <- ggplot(MAG_NS_long, aes(x = graph_name, y = count/total, fill = mutation_type))+
  geom_bar(stat = "identity")+
  theme_classic()+
  labs(y = "Fraction of SNVs", x = "pond", legend = "Mutation Type")+
  facet_wrap(~mag, nrow = 4, ncol = 4, scales="free", labeller = labeller(mag = mag_labs))

save_plot("NS_MAG_plot.jpeg", NS_plot, ncol = 4, nrow = 4, dpi = 300)


# look at N:S ratio by gene
MAG_dnds_by_gene <- MAG_snvs %>% group_by(gene, new_name, mag) %>% summarize(N = sum(mutation_type=="N"), S = sum(mutation_type=="S"), I = sum(mutation_type=="I"), M = sum(mutation_type=="M"))
MAG_dnds_by_gene$NS_ratio <- MAG_dnds_by_gene$N / MAG_dnds_by_gene$S
full_dnds_by_gene <- MAG_dnds_by_gene %>% group_by(mag) %>% complete(new_name, gene)


NS_gene_plot <- ggplot(full_dnds_by_gene, aes(x = new_name, y = NS_ratio, colour=new_name))+
  geom_boxplot(outlier.color = "black")+
  theme_classic()+
  labs(y = "N/S ratio", x = "Pond")+
  facet_wrap(~mag, nrow = 4, ncol = 4, scales="free", labeller = labeller(mag = mag_labs))

save_plot("NS_gene_plot.jpeg", NS_gene_plot, ncol = 4, nrow = 4, dpi = 300)

NS_by_gene_plot <- ggplot(full_dnds_by_gene, aes(x = gene, y = NS_ratio, colour=new_name))+
  geom_point()+
  theme_classic()+
  labs(y = "N/S ratio", x = "Pond")+
  facet_wrap(~mag, nrow = 4, ncol = 4, scales="free", labeller = labeller(mag = mag_labs))

save_plot("NS_by_gene_plot.jpeg", NS_by_gene_plot, ncol = 4, nrow = 4, dpi = 300)

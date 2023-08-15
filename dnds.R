library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(patchwork)

all_MAG_snvs <- read_csv("all_MAG_SNVs_med_Aug8.csv")
all_MAG_snvs <- subset(all_MAG_snvs, str_detect(new_name, "T2"))

MAG_snvs <- all_MAG_snvs[, c('mag', 'gene', 'mutation_type', 'new_name')] %>% subset(is.na(mutation_type) == F)
MAG_snvs$mutation <- 1
MAG_snvs_sum <- MAG_snvs %>% group_by(mag, new_name) %>% summarise(total = sum(is.na(mutation_type)==F))
MAG_snvs_all <- left_join(MAG_snvs, MAG_snvs_sum, by = c('mag', 'new_name'))
MAG_dnds_by_gene <- MAG_snvs %>% group_by(gene, new_name) %>% summarize(N = sum(mutation_type=="N"), S = (mutation_type=="S"))

dnds_plot <- ggplot(MAG_snvs_all, aes(x=new_name, y=mutation/total, fill=mutation_type))+
  geom_bar(stat="identity")+
  theme_classic()+
  labs(y= "Fraction of SNVs", x = "pond", legend = "Mutation Type")+
  facet_wrap(~mag, nrow=4, ncol=4, scales="free")
ggsave("dnds_MAG_plot.png", limitsize = F, dpi = 200)

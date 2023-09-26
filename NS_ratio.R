library(tidyverse)
library(dplyr)
library(ggplot2)
library(viridis)
library(cowplot)

mag_labs <- c(I4_MAG_00006 = "SJAQ100 sp016735685 assembled from Control B", I4_MAG_00065 = "Roseomonas sp. assembled from Control B", L2_MAG_00052 = "Erythrobacter sp. assembled from GBH A", 
                 L3_MAG_00058 = "Prosthecobacter sp. assembled from Control C", L4_MAG_00099 = "Bosea sp001713455 assembled from Control D", L7_MAG_00020 = "Sphingorhabdus_B sp. assembled from GBH C",
                 L7_MAG_00028 = "SYFN01 sp. assembled from GBH C", L7_MAG_00043 = "Luteolibacter sp. assembled from GBH C", L8_MAG_00011 = "UBA953 sp. assembled from GBH D", 
                 L8_MAG_00019 = "UA16 family assembled from GBH D", L8_MAG_00042 = "UBA4660 sp. assembled from GBH D")

select_mags <- list("I4_MAG_00006", "I4_MAG_00065", "L3_MAG_00058", "L7_MAG_00020", "L8_MAG_00011", "L8_MAG_00019")
EPSPS_class_1 <- list("I4_MAG_00006", "L7_MAG_00028", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")
EPSPS_unclass <- list("L2_MAG_00052", "L4_MAG_00099", "L7_MAG_00020")
all_MAG_snvs <- read_csv("all_MAG_SNVs_med_Aug15.csv")
all_MAG_snvs <- subset(all_MAG_snvs, str_detect(new_name, "T2"))
MAG_snvs <- all_MAG_snvs[, c('mag', 'gene', 'mutation_type', 'new_name', 'mag_length')] %>% subset(is.na(mutation_type) == F)

MAG_NS <- MAG_snvs %>% group_by(mag, mag_length, new_name) %>% summarize(N = sum(mutation_type=="N"), S = sum(mutation_type=="S"), I = sum(mutation_type=="I"), M = sum(mutation_type=="M"))
MAG_NS$total <- MAG_NS$N + MAG_NS$S + MAG_NS$I + MAG_NS$M
MAG_NS$graph_name <- str_sub(MAG_NS$new_name, end = -6) 
MAG_NS$graph_name <- gsub('Control', 'Ctrl', MAG_NS$graph_name)
MAG_NS_long <- pivot_longer(MAG_NS, cols = c("N", "S", "I", "M"), names_to = "mutation_type", values_to = "count")

MAG_NS$NS_ratio <-  MAG_NS$N / MAG_NS$S
MAG_NS$SNVs_MBp <- (MAG_NS$total / MAG_NS$mag_length) * 1000000
MAG_NS$treatment <- str_sub(MAG_NS$new_name, end = -8)
MAG_NS$diversity <- with(MAG_NS, ifelse(SNVs_MBp < 2500, 'Low (< 2500 SNV/MBp)', 'High (> 2500 SNV/MBp)'))
MAG_NS_wide <- MAG_NS[, c('mag', 'new_name', 'SNVs_MBp')] %>% pivot_wider(names_from = new_name, values_from = SNVs_MBp)
MAG_NS_wide <- MAG_NS_wide %>% rowwise() %>% mutate(max_SNVs_mag = max(c_across(contains("T")), na.rm = T))
MAG_NS_wide <- MAG_NS_wide %>% rowwise() %>% mutate(max_SNVs_mag_ctl = max(c_across(contains("Control")), na.rm = T))
MAG_NS_plotting <- pivot_longer(MAG_NS_wide, cols = contains("T2"), names_to = 'new_name', values_to = 'SNVs_MBp') %>% subset(is.na(SNVs_MBp) == F) 
MAG_NS_plotting <- right_join(MAG_NS, MAG_NS_plotting, by = c('mag', 'new_name', 'SNVs_MBp'))
MAG_NS_plotting$sweep <- with(MAG_NS_plotting, ifelse(mag %in% select_mags, "Yes", "No"))
MAG_NS_ctl_max <- subset(MAG_NS_plotting, SNVs_MBp == max_SNVs_mag_ctl)
MAG_NS_all_max <- subset(MAG_NS_plotting, SNVs_MBp == max_SNVs_mag)
MAG_NS_all_max$EPSPS_class <- with(MAG_NS_all_max, ifelse(mag %in% EPSPS_class_1, "Class I", "Class II"))
MAG_NS_all_max$EPSPS_class <- with(MAG_NS_all_max, ifelse(mag %in% EPSPS_unclass, "Unclassified", EPSPS_class))

NS_SNV_max_plot <- ggplot(MAG_NS_all_max, aes(x = NS_ratio, y = SNVs_MBp, colour = diversity))+
  geom_point(aes(shape=sweep), size=2.5)+
  scale_color_manual(values = c("#4c117a", "#cc3f71"))+
  labs(y = "SNVs / MBp", x = "N:S ratio", colour = "Diversity", shape = "Potential Sweep?")+
  theme_classic()+
  theme(text = element_text(family = "helvetica"), axis.title = element_text(face = "bold"), axis.text = element_text(face = "bold"), 
        legend.title = element_text(size=8, family="helvetica", face="bold"), legend.text = element_text(size=7, family="helvetica"))
save_plot("NS_SNV_max_plot.jpeg", NS_SNV_max_plot)

NS_SNV_max_plot_EPSPS <- ggplot(MAG_NS_all_max, aes(x = NS_ratio, y = SNVs_MBp, colour = EPSPS_class))+
  geom_point(aes(shape=sweep), size=2.5)+
  scale_color_manual(values = c("#4c117a", "#cc3f71", "darkorange"))+
  labs(y = "SNVs / MBp", x = "N:S ratio", colour = "EPSPS Class", shape = "Potential Sweep?")+
  theme_classic()+
  theme(text = element_text(family = "helvetica"), axis.title = element_text(face = "bold"), axis.text = element_text(face = "bold"), 
        legend.title = element_text(size=8, family="helvetica", face="bold"), legend.text = element_text(size=7, family="helvetica"))
save_plot("NS_SNV_max_EPSPS_plot.jpeg", NS_SNV_max_plot_EPSPS)

NS_SNV_max_log10_plot <- ggplot(MAG_NS_all_max, aes(x = NS_ratio, y = log10(SNVs_MBp), colour = diversity))+
  geom_point(aes(shape=sweep), size=2.5)+
  scale_color_manual(values = c("#4c117a", "#cc3f71"))+
  labs(y = "SNVs / MBp", x = "N:S ratio", colour = "Diversity", shape = "Potential Sweep?")+
  theme_classic()+
  theme(text = element_text(family = "helvetica"), axis.title = element_text(face = "bold"), axis.text = element_text(face = "bold"), 
        legend.title = element_text(size=8, family="helvetica", face="bold"), legend.text = element_text(size=7, family="helvetica"))
save_plot("NS_SNV_max_log10_plot.jpeg", NS_SNV_max_log10_plot)

NS_SNV_ctl_plot <- ggplot(MAG_NS_ctl_max, aes(x = NS_ratio, y = SNVs_MBp, colour = diversity))+
  geom_point(aes(shape=sweep), size=2.5)+
  scale_color_manual(values = c("#4c117a", "#cc3f71"))+
  labs(y = "SNVs / MBp", x = "N:S ratio", colour = "Diversity", shape = "Potential Sweep?")+
  theme_classic()+
  theme(text = element_text(family = "helvetica"), axis.title = element_text(face = "bold"), axis.text = element_text(face = "bold"), 
        legend.title = element_text(size=8, family="helvetica", face="bold"), legend.text = element_text(size=7, family="helvetica"))
save_plot("NS_SNV_ctl_plot.jpeg", NS_SNV_ctl_plot)

NS_SNV_ctl_log10_plot <- ggplot(MAG_NS_ctl_max, aes(x = NS_ratio, y = log10(SNVs_MBp), colour = diversity))+
  geom_point(aes(shape=sweep), size=2.5)+
  scale_color_manual(values = c("#4c117a", "#cc3f71"))+
  labs(y = "SNVs / MBp", x = "N:S ratio", colour = "Diversity", shape = "Potential Sweep?")+
  theme_classic()+
  theme(text = element_text(family = "helvetica"), axis.title = element_text(face = "bold"), axis.text = element_text(face = "bold"), 
        legend.title = element_text(size=8, family="helvetica", face="bold"), legend.text = element_text(size=7, family="helvetica"))
save_plot("NS_SNV_ctl_log10_plot.jpeg", NS_SNV_ctl_log10_plot)

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

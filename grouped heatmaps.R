library(tidyverse)
library(dplyr)
library(ggplot2)
library(viridis)
library(cowplot)
library(patchwork)

sens_mags <- list("I4_MAG_00006", "L7_MAG_00028", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")
res_mags <- list("I4_MAG_00065", "L3_MAG_00058", "L7_MAG_00020", "L7_MAG_00043")
unclass_mags <- list("L2_MAG_00052", "L4_MAG_00099")

sens_mag_labs <- c(I4_MAG_00006 = "Burkholderiaceae 1", L7_MAG_00028 = "Burkholderiaceae 2", L8_MAG_00011 = "Verrucomicrobiae", L8_MAG_00019 = "Flavobacteriales 1", L8_MAG_00042 = "Flavobacteriales 2")
res_mag_labs <-  c(I4_MAG_00065 = "Roseomonas_A",  L3_MAG_00058 = "Prosthecobacter", L7_MAG_00020 = "Sphingorhabdus_B", L7_MAG_00043 = "Luteolibacter")
unclass_mag_labs <- c(L2_MAG_00052 = "Erythrobacter", L4_MAG_00099 = "Bosea sp001713455")

all_snv <- read_csv("all_MAG_SNVs_med_Apr9.csv")
all_snv$graph_name <- gsub("Control", "CL", all_snv$new_name) %>% str_sub(end = -6)

sens_snv <- subset(all_snv, mag %in% sens_mags)
res_snv <- subset(all_snv, mag %in% res_mags & str_detect(all_snv$new_name, "T2"))
unclass_snv <- subset(all_snv, mag %in% unclass_mags)
all_L7 <- subset(all_snv, mag == "L7_MAG_00020")
all_L7 <- all_L7 %>% add_row(new_name = "GBH C at T2", final_ref_freq = NA, treatment = "glyphosate", name = "GBH C")
all_L7$graph_name <- all_L7$new_name %>% str_sub(start = -2)

all_sum <- read_csv("all_SNV_sum_subsamp.csv")
all_sum$graph_name <- gsub("Control", "CL", all_sum$new_name) %>% str_sub(end = -6)
all_sum$treatment <- with(all_sum, ifelse(str_detect(new_name, "Control"), "Control", "GBH"))
all_sum$treatment <- with(all_sum, ifelse(str_detect(new_name, "Control E"), "Phosphorus", treatment))

sens_sum <- subset(all_sum, mag %in% sens_mags)
res_sum <- subset(all_sum, mag %in% res_mags & str_detect(all_sum$new_name, "T2"))
unclass_sum <- subset(all_sum, mag %in% unclass_mags)
all_sum_L7 <- subset(all_sum, mag == "L7_MAG_00020")
all_sum_L7 <- all_sum_L7 %>% add_row(new_name = "GBH C at T2", SNV_Mbp = 0, treatment = "GBH")
all_sum_L7$graph_name <- all_sum_L7$new_name %>% str_sub(start = -2)
all_sum_L7$name <- all_sum_L7$new_name %>% str_sub(end = -7)

#FIGURE 2

sens_snv_heat <- ggplot(sens_snv, aes(x = graph_name, y = reorder(groups, all_mean), colour = final_ref_freq)) +
  geom_tile()+
  scale_colour_viridis(direction = -1, na.value = "white") +
  theme_classic() +
  theme(text = element_text(size = 6.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.title.y = element_text(size = 8), 
        axis.text.x = element_text(colour = "black", vjust = -1),
        axis.title.x = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.justification = "left",
        strip.text.x.top = element_text(size = 7, face = "bold"))+
  labs(colour = NULL, y = "Nucleotide Position")+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag, ncol = 5, scales = "free", labeller = labeller(mag = sens_mag_labs))

sens_snv_sum <- ggplot(sens_sum, aes(x = graph_name, y = SNV_Mbp, fill = treatment))+
  geom_bar(stat = "identity", colour = "black")+ 
  theme_classic()+
  scale_fill_manual(breaks = c('Control', 'Phosphorus', 'GBH'), values= c("white", "grey70", "grey30"))+
  theme(text = element_text(size = 6.5),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(vjust = -1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"),
        legend.justification = "left",
        legend.key.spacing.y = unit(0.1, "cm"),
        panel.spacing = unit(0.05, "cm"),
        strip.text.x = element_blank())+
  labs(y = "SNVs / Mbp", fill = NULL)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag, ncol = 5, scales = "free")

sens_all <- sens_snv_heat / sens_snv_sum
save_plot("sens_all.jpeg", sens_all, base_height = 3.25, base_width = 8.24, dpi = 600)

res_snv_heat <- ggplot(res_snv, aes(x = graph_name, y = reorder(groups, all_mean), colour = final_ref_freq)) +
  geom_tile()+
  scale_colour_viridis(direction = -1, na.value = "white") +
  theme_classic() +
  theme(text = element_text(size = 6.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.title.y = element_text(size = 8), 
        axis.text.x = element_text(colour = "black", vjust = -1),
        axis.title.x = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.justification = "left",
        strip.text.x.top = element_text(size = 7, face = "bold"))+
  labs(colour = NULL, y = "Nucleotide Position")+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag, ncol = 5, scales = "free", labeller = labeller(mag = res_mag_labs))

res_snv_sum <- ggplot(res_sum, aes(x = graph_name, y = SNV_Mbp, fill = treatment))+
  geom_bar(stat = "identity", colour = "black")+ 
  theme_classic()+
  scale_fill_manual(breaks = c('Control', 'Phosphorus', 'GBH'), values= c("white", "grey70", "grey30"))+
  theme(text = element_text(size = 6.5),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(vjust = -1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"),
        legend.justification = "left",
        legend.key.spacing.y = unit(0.1, "cm"),
        panel.spacing = unit(0.05, "cm"),
        strip.text.x = element_blank())+
  labs(y = "SNVs / Mbp", fill = NULL)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag, ncol = 5, scales = "free")

res_all <- res_snv_heat / res_snv_sum
save_plot("res_all.jpeg", res_all, base_height = 3.25, base_width = 6.9, dpi = 600)


#FIGURE 3
L7_snv_heat <- ggplot(all_L7, aes(x = graph_name, y = reorder(groups, all_mean), colour = final_ref_freq)) +
  geom_tile()+
  scale_colour_viridis(direction = -1, na.value = "white") +
  theme_classic() +
  theme(text = element_text(size = 6.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.title.y = element_text(size = 8), 
        axis.text.x = element_text(colour = "black", vjust = -1),
        axis.title.x = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.justification = "left",
        strip.text.x.top = element_text(size = 7, face = "bold"))+
  labs(colour = NULL, y = "Nucleotide Position")+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~name, ncol = 5, scales = "free")

L7_snv_sum <- ggplot(all_sum_L7, aes(x = graph_name, y = SNV_Mbp, fill = treatment))+
  geom_bar(stat = "identity", colour = "black")+ 
  theme_classic()+
  scale_fill_manual(breaks = c('Control', 'Phosphorus', 'GBH'), values= c("white", "grey70", "grey30"))+
  theme(text = element_text(size = 6.5),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(vjust = -1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"),
        legend.justification = "left",
        legend.key.spacing.y = unit(0.1, "cm"),
        panel.spacing = unit(0.05, "cm"),
        strip.text.x = element_blank())+
  labs(y = "SNVs / Mbp", fill = NULL)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 3000))+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~name, ncol = 5, scales = "free")

L7_all <- L7_snv_heat / L7_snv_sum
save_plot("L7_all.jpeg", L7_all, base_height = 3.25, base_width = 8.1, dpi = 600)


#FIGURE S2
unclass_snv_heat <- ggplot(unclass_snv, aes(x = graph_name, y = reorder(groups, all_mean), colour = final_ref_freq)) +
  geom_tile()+
  scale_colour_viridis(direction = -1, na.value = "white") +
  theme_classic() +
  theme(text = element_text(size = 6.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.title.y = element_text(size = 8), 
        axis.text.x = element_text(colour = "black", vjust = -1),
        axis.title.x = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.justification = "left",
        strip.text.x.top = element_text(size = 7, face = "bold"))+
  labs(colour = NULL, y = "Nucleotide Position")+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag, ncol = 2, scales = "free", labeller = labeller(mag = unclass_mag_labs))

unclass_snv_sum <- ggplot(unclass_sum, aes(x = graph_name, y = SNV_Mbp, fill = treatment))+
  geom_bar(stat = "identity", colour = "black")+ 
  theme_classic()+
  scale_fill_manual(breaks = c('Control', 'Phosphorus', 'GBH'), values= c("white", "grey70", "grey30"))+
  theme(text = element_text(size = 6.5),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(vjust = -1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"),
        legend.justification = "left",
        legend.key.spacing.y = unit(0.1, "cm"),
        panel.spacing = unit(0.05, "cm"),
        strip.text.x = element_blank())+
  labs(y = "SNVs / Mbp", fill = NULL)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag, ncol = 2, scales = "free")

unclass_all <- unclass_snv_heat / unclass_snv_sum
save_plot("unclass_all.jpeg", unclass_all, base_height = 3.25, base_width = 4.2, dpi = 600)

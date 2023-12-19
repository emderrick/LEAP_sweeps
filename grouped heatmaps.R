library(tidyverse)
library(dplyr)
library(ggplot2)
library(viridis)
library(cowplot)
library(patchwork)

sens_mags <- list("I4_MAG_00006", "L7_MAG_00028", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")
res_mags <- list("I4_MAG_00065", "L3_MAG_00058", "L7_MAG_00020", "L7_MAG_00043")
unclass_mags <- list("L2_MAG_00052", "L4_MAG_00099")

sens_mag_labs <- c(I4_MAG_00006 = "SJAQ100 sp016735685", L7_MAG_00028 = "SYFN01 sp.", L8_MAG_00011 = "UBA953 sp.", L8_MAG_00019 = "UA16", L8_MAG_00042 = "UBA4660 sp.")
res_mag_labs <-  c(I4_MAG_00065  = "Roseomonas sp.",  L3_MAG_00058 = "Prosthecobacter sp.", L7_MAG_00020 = "Sphingorhabdus_B sp.", L7_MAG_00043 = "Luteolibacter sp.")
unclass_mag_labs <- c(L2_MAG_00052 = "Erythrobacter sp.", L4_MAG_00099 = "Bosea sp001713455")

all_snv <- read_csv("all_MAG_SNVs_med_Dec7.csv")
all_snv$graph_name <- gsub("Control", "CTL", all_snv$new_name) %>% str_sub(end = -6)

sens_snv <- subset(all_snv, mag %in% sens_mags)
res_snv <- subset(all_snv, mag %in% res_mags & str_detect(all_snv$new_name, "T2"))
unclass_snv <- subset(all_snv, mag %in% unclass_mags)
all_L7 <- subset(all_snv, mag == "L7_MAG_00020")
all_L7$graph_name <- gsub("Control", "CTL", all_L7$new_name)
all_L7$graph_name <- gsub("at ","", all_L7$graph_name)
all_L7$pond <- all_L7$graph_name %>% str_sub(end = -3)

all_sum <- read_csv("all_SNV_sum_subsamp.csv")
all_sum$graph_name <- gsub("Control", "CTL", all_sum$new_name) %>% str_sub(end = -6)
all_sum$treatment <- with(all_sum, ifelse(str_detect(new_name, "Control"), "Control", "GBH"))
all_sum$treatment <- with(all_sum, ifelse(str_detect(new_name, "Control E"), "Phosphorus", treatment))

sens_sum <- subset(all_sum, mag %in% sens_mags)
res_sum <- subset(all_sum, mag %in% res_mags & str_detect(all_sum$new_name, "T2"))
unclass_sum <- subset(all_sum, mag %in% unclass_mags)
all_sum_L7 <- subset(all_sum, mag == "L7_MAG_00020")
all_sum_L7$graph_name <- gsub("Control", "CTL", all_sum_L7$new_name)
all_sum_L7$graph_name <- gsub("at ","", all_sum_L7$graph_name)
all_sum_L7$pond <- all_sum_L7$graph_name %>% str_sub(end = -3)

#sensitive
sens_snv_heat <- ggplot(sens_snv, aes(x = graph_name, y = reorder(groups, all_mean), colour = final_ref_freq)) +
  geom_tile()+
  scale_colour_viridis(direction = -1, na.value = "white") +
  theme_classic() +
  theme(text = element_text(size = 16),
        axis.line.x.bottom = element_line(linewidth = 0.85),
        axis.line.y.left = element_line(linewidth = 0.85),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.title.y = element_text(size = 16), 
        axis.text.x = element_text(colour = "black", vjust = -1),
        axis.title.x = element_blank(),
        strip.text.x.top = element_text(size = 16, margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"), face = "bold"),
        strip.background = element_rect(linewidth = 2),
        legend.key.size = unit(1, "cm"),
        legend.margin =  margin(t = 0, r = 0, b = 0, l = 40),
        legend.text = element_text(size = 12),
        plot.margin = unit(c(0,0,1,1), "cm"), 
        panel.spacing = unit(0.5, "cm"))+
  labs(colour = NULL, y = "Genome Position")+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag, ncol = 5, scales = "free", labeller = labeller(mag = sens_mag_labs))

sens_snv_sum <- ggplot(sens_sum, aes(x = graph_name, y = SNV_Mbp, fill = treatment))+
  geom_bar(stat = "identity", colour = "black")+ 
  theme_classic()+
  scale_fill_manual(breaks = c('Control', 'Phosphorus', 'GBH'), values= c("white", "grey70", "grey30"))+
  theme(text = element_text(size = 16),
        axis.line = element_line(linewidth = 0.85),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(vjust = -1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        strip.text.x = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0,0,1,1), "cm"), 
        panel.spacing = unit(0.05, "cm"))+
  labs(y = "SNVs / Mbp")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag, ncol = 5, scales = "free")

sens_snv_frac <- ggplot(sens_sum, aes(x = graph_name, y = SNS/(SNS+SNV), fill = treatment))+
  geom_bar(stat="identity", colour = "black")+ 
  theme_classic()+
  scale_fill_manual(breaks = c('Control', 'Phosphorus', 'GBH'), values= c("white", "grey70", "grey30"))+
  theme(text = element_text(size = 16),
        axis.line = element_line(linewidth = 0.85),
        axis.text = element_text(color = "black"),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_text(vjust = -1),
        strip.text.x = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0.2, "cm"))+
  labs(y = "SNSs / SNSs + SNVs")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0,1))+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag, ncol = 5, scales = "free")

sens_all <- sens_snv_heat / sens_snv_sum / sens_snv_frac
save_plot("sens_all.jpeg", sens_all, base_height = 10, base_width = 20, dpi = 400, limitsize = F)

#resistant
res_snv_heat <- ggplot(res_snv, aes(x = graph_name, y = reorder(groups, all_mean), colour = final_ref_freq)) +
  geom_tile()+
  scale_colour_viridis(direction = -1, na.value = "white") +
  theme_classic() +
  theme(text = element_text(size = 16),
        axis.line.x.bottom = element_line(linewidth = 0.85),
        axis.line.y.left = element_line(linewidth = 0.85),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.title.y = element_text(size = 16), 
        axis.text.x = element_text(colour = "black", vjust = -1),
        axis.title.x = element_blank(),
        strip.text.x.top = element_text(size = 16, margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"), face = "bold"),
        strip.background = element_rect(linewidth = 2),
        legend.key.size = unit(1, "cm"),
        legend.margin =  margin(t = 0, r = 0, b = 0, l = 40),
        legend.text = element_text(size = 12),
        plot.margin = unit(c(0,0,1,1), "cm"), 
        panel.spacing = unit(0.5, "cm"))+
  labs(colour = NULL, y = "Genome Position")+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag, ncol = 4, scales = "free", labeller = labeller(mag = res_mag_labs))

res_snv_sum <- ggplot(res_sum, aes(x = graph_name, y = SNV_Mbp, fill = treatment))+
  geom_bar(stat = "identity", colour = "black")+ 
  theme_classic()+
  scale_fill_manual(breaks = c('Control', 'Phosphorus', 'GBH'), values= c("white", "grey70", "grey30"))+
  theme(text = element_text(size = 16),
        axis.line = element_line(linewidth = 0.85),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(vjust = -1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        strip.text.x = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0,0,1,1), "cm"), 
        panel.spacing = unit(0.05, "cm"))+
  labs(y = "SNVs / Mbp")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag, ncol = 4, scales = "free")

res_snv_frac <- ggplot(res_sum, aes(x = graph_name, y = SNS/(SNS+SNV), fill = treatment))+
  geom_bar(stat="identity", colour = "black")+ 
  theme_classic()+
  scale_fill_manual(breaks = c('Control', 'Phosphorus', 'GBH'), values= c("white", "grey70", "grey30"))+
  theme(text = element_text(size = 16),
        axis.line = element_line(linewidth = 0.85),
        axis.text = element_text(color = "black"),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_text(vjust = -1),
        strip.text.x = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0.2, "cm"))+
  labs(y = "SNSs / SNSs + SNVs")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0,1))+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag, ncol = 4, scales = "free")

res_all <- res_snv_heat / res_snv_sum / res_snv_frac
save_plot("res_all.jpeg", res_all, base_height = 10, base_width = 17, dpi = 400, limitsize = F)

#unclassified
unclass_snv_heat <- ggplot(unclass_snv, aes(x = graph_name, y = reorder(groups, all_mean), colour = final_ref_freq)) +
  geom_tile()+
  scale_colour_viridis(direction = -1, na.value = "white") +
  theme_classic() +
  theme(text = element_text(size = 16),
        axis.line.x.bottom = element_line(linewidth = 0.85),
        axis.line.y.left = element_line(linewidth = 0.85),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.title.y = element_text(size = 16), 
        axis.text.x = element_text(colour = "black", vjust = -1),
        axis.title.x = element_blank(),
        strip.text.x.top = element_text(size = 16, margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"), face = "bold"),
        strip.background = element_rect(linewidth = 2),
        legend.key.size = unit(1, "cm"),
        legend.margin =  margin(t = 0, r = 0, b = 0, l = 40),
        legend.text = element_text(size = 12),
        plot.margin = unit(c(0,0,1,1), "cm"), 
        panel.spacing = unit(0.5, "cm"))+
  labs(colour = NULL, y = "Genome Position")+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag, ncol = 2, scales = "free", labeller = labeller(mag = unclass_mag_labs))

unclass_snv_sum <- ggplot(unclass_sum, aes(x = graph_name, y = SNV_Mbp, fill = treatment))+
  geom_bar(stat = "identity", colour = "black")+ 
  theme_classic()+
  scale_fill_manual(breaks = c('Control', 'Phosphorus', 'GBH'), values= c("white", "grey70", "grey30"))+
  theme(text = element_text(size = 16),
        axis.line = element_line(linewidth = 0.85),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(vjust = -1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        strip.text.x = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0,0,1,1), "cm"), 
        panel.spacing = unit(0.05, "cm"))+
  labs(y = "SNVs / Mbp")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag, ncol = 2, scales = "free")

unclass_snv_frac <- ggplot(unclass_sum, aes(x = graph_name, y = SNS/(SNS+SNV), fill = treatment))+
  geom_bar(stat="identity", colour = "black")+ 
  theme_classic()+
  scale_fill_manual(breaks = c('Control', 'Phosphorus', 'GBH'), values= c("white", "grey70", "grey30"))+
  theme(text = element_text(size = 16),
        axis.line = element_line(linewidth = 0.85),
        axis.text = element_text(color = "black"),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_text(vjust = -1),
        strip.text.x = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0.2, "cm"))+
  labs(y = "SNSs / SNSs + SNVs")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0,1))+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag, ncol = 2, scales = "free")

unclass_all <- unclass_snv_heat / unclass_snv_sum / unclass_snv_frac
save_plot("unclass_all.jpeg", unclass_all, base_height = 10, base_width = 10, dpi = 400, limitsize = F)

#L7_MAG_00020 heatmaps
L7_snv_heat <- ggplot(all_L7, aes(x = graph_name, y = reorder(groups, all_mean), colour = final_ref_freq)) +
  geom_tile()+
  scale_colour_viridis(direction = -1, na.value = "white") +
  theme_classic() +
  theme(text = element_text(size = 16),
        axis.line.x.bottom = element_line(linewidth = 0.85),
        axis.line.y.left = element_line(linewidth = 0.85),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.title.y = element_text(size = 16), 
        axis.text.x = element_text(colour = "black", vjust = -1),
        axis.title.x = element_blank(),
        strip.text.x.top = element_text(size = 16, margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"), face = "bold"),
        strip.background = element_rect(linewidth = 2),
        legend.key.size = unit(1, "cm"),
        legend.margin =  margin(t = 0, r = 0, b = 0, l = 40),
        legend.text = element_text(size = 12),
        plot.margin = unit(c(0,0,1,1), "cm"), 
        panel.spacing = unit(0.5, "cm"))+
  labs(colour = NULL, y = "Genome Position")+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~pond, ncol = 4, scales = "free")

L7_snv_sum <- ggplot(all_sum_L7, aes(x = graph_name, y = SNV_Mbp, fill = treatment))+
  geom_bar(stat = "identity", colour = "black")+ 
  theme_classic()+
  scale_fill_manual(breaks = c('Control', 'Phosphorus', 'GBH'), values= c("white", "grey70", "grey30"))+
  theme(text = element_text(size = 16),
        axis.line = element_line(linewidth = 0.85),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(vjust = -1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        strip.text.x = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0,0,1,1), "cm"), 
        panel.spacing = unit(0.05, "cm"))+
  labs(y = "SNVs / Mbp")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 2600))+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~pond, ncol = 4, scales = "free")

L7_snv_frac <- ggplot(all_sum_L7, aes(x = graph_name, y = SNS/(SNS+SNV), fill = treatment))+
  geom_bar(stat="identity", colour = "black")+ 
  theme_classic()+
  scale_fill_manual(breaks = c('Control', 'Phosphorus', 'GBH'), values= c("white", "grey70", "grey30"))+
  theme(text = element_text(size = 16),
        axis.line = element_line(linewidth = 0.85),
        axis.text = element_text(color = "black"),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_text(vjust = -1),
        strip.text.x = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0.2, "cm"))+
  labs(y = "SNSs / SNSs + SNVs")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0,1))+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~pond, ncol = 4, scales = "free")

L7_all <- L7_snv_heat / L7_snv_sum / L7_snv_frac
save_plot("L7_all.jpeg", L7_all, base_height = 10, base_width = 17, dpi = 400, limitsize = F)

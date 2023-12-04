library(tidyverse)
library(dplyr)
library(ggplot2)
library(viridis)
library(stringr)
library(cowplot)
library(patchwork)

select_mags <- list("I4_MAG_00006", "I4_MAG_00065", "L3_MAG_00058", "L7_MAG_00020", "L8_MAG_00011", "L8_MAG_00019")
other_mags <- list ("L2_MAG_00052", "L4_MAG_00099", "L7_MAG_00028", "L7_MAG_00043", "L8_MAG_00042")

all_snv <- read_csv("all_MAG_SNVs_med_Aug15.csv")
all_snv <- subset(all_snv, new_time == 2)
all_snv$graph_name <- gsub("Control", "CTL", all_snv$name)

select_snv <- subset(all_snv, mag %in% select_mags)
other_snv <- subset(all_snv, mag %in% other_mags)

all_sum <- read_csv("all_SNV_sum.csv")
all_sum <- subset(all_sum, str_detect(new_name, "2"))
all_sum$graph_name <- gsub("Control", "CTL", all_sum$new_name) %>% str_sub(end = -6)

all_sum$treatment <- with(all_sum, ifelse(str_detect(new_name, "Control"), "Control", "GBH"))
all_sum$treatment <- with(all_sum, ifelse(str_detect(new_name, "Control E"), "Phosphorus", treatment))

select_sum <- subset(all_sum, mag %in% select_mags)
other_sum <- subset(all_sum, mag %in% other_mags)

all_sum_long <- pivot_longer(all_sum, cols = contains("S"), names_to = "class", values_to = "divergent_sites", values_drop_na = F)
select_sum_long <- subset(all_sum_long, mag %in% select_mags)
other_sum_long <- subset(all_sum_long, mag %in% other_mags)

select_mag_labs <- (c(I4_MAG_00006 = "SJAQ100 sp016735685", I4_MAG_00065 = "Roseomonas sp.", L3_MAG_00058 = "Prosthecobacter sp.",
                      L7_MAG_00020 = "Sphingorhabdus_B sp. ", L8_MAG_00011 = "UBA953 sp.", L8_MAG_00019 = "UA16"))

other_mag_labs <- (c(L2_MAG_00052 = "Erythrobacter sp.", L4_MAG_00099 = "Bosea sp001713455", L7_MAG_00028 = "SYFN01 sp.",
                     L7_MAG_00043 = "Luteolibacter sp.", L8_MAG_00042 = "UBA4660 sp."))

#maybe sweep
select_snv_heat <- ggplot(select_snv, aes(x = graph_name, y = reorder(groups, all_mean), colour = final_ref_freq)) +
  geom_tile()+
  scale_colour_viridis(direction = -1, na.value = "white") +
  theme_classic() +
  theme(text = element_text(size = 16),
        axis.line.x.bottom = element_line(linewidth = 0.85),
        axis.line.y.left = element_line(linewidth = 0.85),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.title.y = element_text(size = 18), 
        axis.text.x = element_text(colour = "black", vjust = -1),
        axis.title.x = element_blank(),
        strip.text.x.top = element_text(size = 20, margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), face = "bold"),
        strip.background = element_rect(linewidth = 2),
        legend.title = element_text(size = 18),
        legend.key.size = unit(1, "cm"),
        legend.margin =  margin(t = 0, r = 0, b = 0, l = 40),
        legend.text = element_text(size = 12),
        plot.margin = unit(c(0,0,1,1), "cm"), 
        panel.spacing = unit(0.5, "cm"))+
  labs(colour = NULL, y = "Genome Position")+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag, nrow = 1, ncol = 6, scales = "free", labeller = labeller(mag = select_mag_labs))

select_snv_sum <- ggplot(subset(select_sum_long, class == "SNVs"), aes(x = graph_name, y=((divergent_sites/mag_length)*10^6), fill = treatment))+
  geom_bar(stat = "identity", colour = "black")+ 
  theme_classic()+
  scale_fill_manual(breaks = c('Control', 'Phosphorus', 'GBH'), values= c("white", "grey70", "grey30"))+
  theme(text = element_text(size = 16),
        axis.line = element_line(linewidth = 0.85),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(vjust = -1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        strip.text.x = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0,0,1,1), "cm"), 
        panel.spacing = unit(0.05, "cm"))+
  labs(y = "SNVs / Mbp")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag, nrow = 1, scales = "free")

select_snv_frac <- ggplot(select_sum, aes(x = graph_name, y = SNSs/(SNSs+SNVs), fill = treatment))+
  geom_bar(stat="identity", colour = "black")+ 
  theme_classic()+
  scale_fill_manual(breaks = c('Control', 'Phosphorus', 'GBH'), values= c("white", "grey70", "grey30"))+
  theme(text = element_text(size = 16),
        axis.line = element_line(linewidth = 0.85),
        axis.text = element_text(color = "black"),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.text.x = element_text(vjust = -1),
        strip.text.x = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0.2, "cm"))+
  labs(y = "SNSs / SNSs + SNVs")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0,1))+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag, nrow = 1, scales = "free")

select_all <- select_snv_heat / select_snv_sum / select_snv_frac
save_plot("select_all.jpeg", select_all, base_height = 2, base_width = 9.6, ncol = 3, nrow = 6, dpi = 400, limitsize = F)

#no sweep
other_snv_heat <- ggplot(other_snv, aes(x = graph_name, y = reorder(groups, all_mean), colour = final_ref_freq)) +
  geom_tile()+
  scale_colour_viridis(direction = -1, na.value = "white") +
  theme_classic() +
  theme(text = element_text(size = 16),
        axis.line.x.bottom = element_line(linewidth = 0.85),
        axis.line.y.left = element_line(linewidth = 0.85),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.title.y = element_text(size = 18), 
        axis.text.x = element_text(colour = "black", vjust = -1),
        axis.title.x = element_blank(),
        strip.text.x.top = element_text(size = 20, margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), face = "bold"),
        strip.background = element_rect(linewidth = 2),
        legend.title = element_text(size = 18),
        legend.key.size = unit(1, "cm"),
        legend.margin =  margin(t = 0, r = 0, b = 0, l = 40),
        legend.text = element_text(size = 12),
        plot.margin = unit(c(0,0,1,1), "cm"), 
        panel.spacing = unit(1, "cm"))+
  labs(colour = NULL, y = "Genome Position")+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag, nrow = 1, ncol = 6, scales = "free", labeller = labeller(mag = other_mag_labs))

other_snv_sum <- ggplot(subset(other_sum_long, class == "SNVs"), aes(x = graph_name, y=((divergent_sites/mag_length)*10^6), fill = treatment))+
  geom_bar(stat = "identity", colour = "black")+ 
  theme_classic()+
  scale_fill_manual(breaks = c('Control', 'Phosphorus', 'GBH'), values= c("white", "grey70", "grey30"))+
  theme(text = element_text(size = 16),
        axis.line = element_line(linewidth = 0.85),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(vjust = -1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        strip.text.x = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0,0,1,1), "cm"), 
        panel.spacing = unit(0.2, "cm"))+
  labs(y = "SNVs / Mbp")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag, nrow = 1, scales = "free")

other_snv_frac <- ggplot(other_sum, aes(x = graph_name, y = SNSs/(SNSs+SNVs), fill = treatment))+
  geom_bar(stat="identity", colour = "black")+ 
  theme_classic()+
  scale_fill_manual(breaks = c('Control', 'Phosphorus', 'GBH'), values= c("white", "grey70", "grey30"))+
  theme(text = element_text(size = 16),
        axis.line = element_line(linewidth = 0.85),
        axis.text = element_text(color = "black"),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.text.x = element_text(vjust = -1),
        strip.text.x = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0.2, "cm"))+
  labs(y = "SNSs / SNSs + SNVs")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0,1))+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag, nrow = 1, scales = "free")

other_all <- other_snv_heat / other_snv_sum / other_snv_frac
save_plot("other_all.jpeg", other_all, base_height = 2, base_width = 8, ncol = 3, nrow = 6, dpi = 400, limitsize = F)


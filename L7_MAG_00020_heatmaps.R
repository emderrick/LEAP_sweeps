library(tidyverse)
library(dplyr)
library(ggplot2)
library(viridis)
library(stringr)
library(cowplot)
library(patchwork)

all_snv <- read_csv("all_MAG_SNVs_med_Aug15.csv")

all_snv <- subset(all_snv, mag == "L7_MAG_00020")
all_snv$graph_name <- gsub("Control", "CTL", all_snv$new_name) %>% str_replace(" at", "") 

all_sum <- all_snv %>% group_by(mag, mag_length, new_name) %>% summarize(SNVs = sum(class %in% "SNV"), SNSs = sum(class %in% "SNS"))
all_sum$graph_name <- gsub("Control", "CTL", all_sum$new_name) %>% str_replace(" at", "")
all_sum$treatment <- with(all_sum, ifelse(str_detect(new_name, "Control"), "Control", "GBH"))
all_sum$treatment <- with(all_sum, ifelse(str_detect(new_name, "Control E"), "Phosphorus", treatment))
all_sum_long <- pivot_longer(all_sum, cols = contains("S"), names_to = "class", values_to = "divergent_sites", values_drop_na = F)

snv_heat <- ggplot(all_snv, aes(x = graph_name, y = reorder(groups, all_mean), colour = final_ref_freq)) +
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
        legend.title = element_text(size = 18),
        legend.key.size = unit(1, "cm"),
        legend.margin =  margin(t = 0, r = 0, b = 0, l = 40),
        legend.text = element_text(size = 12),
        plot.margin = unit(c(0,0,1,1), "cm"), 
        panel.spacing = unit(0.5, "cm"))+
  labs(colour = NULL, y = "Genome Position")+
  scale_x_discrete(expand = c(0, 0))

snv_sum <- ggplot(subset(all_sum_long, class == "SNVs"), aes(x = graph_name, y=((divergent_sites/mag_length)*10^6), fill = treatment))+
  geom_bar(stat = "identity", colour = "black")+ 
  theme_classic()+
  scale_fill_manual(breaks = c('Control', 'GBH'), values= c("white", "grey30"))+
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
  scale_x_discrete(expand = c(0, 0))

snv_frac <- ggplot(all_sum, aes(x = graph_name, y = SNSs/(SNSs+SNVs), fill = treatment))+
  geom_bar(stat="identity", colour = "black")+ 
  theme_classic()+
  scale_fill_manual(breaks = c('Control', 'GBH'), values= c("white", "grey30"))+
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
  scale_x_discrete(expand = c(0, 0))

L7_MAG_00020_all <- snv_heat / snv_sum / snv_frac
save_plot("L7_MAG_00020_all.jpeg", L7_MAG_00020_all, base_height = 16, base_width = 12, dpi = 400, limitsize = F)

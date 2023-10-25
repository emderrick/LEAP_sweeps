library(tidyverse)
library(dplyr)
library(ggplot2)
library(viridis)
library(stringr)
library(cowplot)
library(patchwork)

select_mags <- list("I4_MAG_00006", "I4_MAG_00065", "L3_MAG_00058", "L7_MAG_00020", "L8_MAG_00011", "L8_MAG_00019")
other_mags <- list ("L2_MAG_00052", "L4_MAG_00099", "L7_MAG_00028", "L7_MAG_00043", "L8_MAG_00042")
best_examples <- list("I4_MAG_00006", "L3_MAG_00058")

all_snv <- read_csv("all_MAG_SNVs_med_Aug15.csv")
all_snv$graph_name <- with(all_snv, ifelse(mag == "L7_MAG_00020", new_name, str_sub(all_snv$new_name, end = -6)))
all_snv$graph_name <- with(all_snv, ifelse(mag == "L7_MAG_00020", gsub('Control', 'Ctrl', graph_name), graph_name))
all_snv$graph_name <- with(all_snv, ifelse(mag == "L7_MAG_00020", gsub('at ', '', graph_name), graph_name))
select_snv <- subset(all_snv, mag %in% select_mags)
other_snv <- subset(all_snv, mag %in% other_mags)
best_snv <- subset(all_snv, mag %in% best_examples)

all_sum <- read_csv("all_snv_sum.csv")
all_sum$graph_name <- with(all_sum, ifelse(mag == "L7_MAG_00020", new_name, str_sub(all_sum$new_name, end = -6)))
all_sum$graph_name <- with(all_sum, ifelse(mag == "L7_MAG_00020", gsub('Control', 'Ctrl', graph_name), graph_name))
all_sum$graph_name <- with(all_sum, ifelse(mag == "L7_MAG_00020", gsub('at ', '', graph_name), graph_name))
all_sum$treatment <- with(all_sum, ifelse(graph_name == "GBH A T1" | str_detect(new_name, "Control"), "Control", "GBH"))
all_sum$treatment <- with(all_sum, ifelse(graph_name == "Control E ", "Phosphorus", treatment))

select_sum <- subset(all_sum, mag %in% select_mags)
other_sum <- subset(all_sum, mag %in% other_mags)
best_sum <- subset(all_sum, mag %in% best_examples)

all_sum_long <- pivot_longer(all_sum, cols = contains("S"), names_to = "class", values_to = "divergent_sites", values_drop_na = F)
select_sum_long <- subset(all_sum_long, mag %in% select_mags)
other_sum_long <- subset(all_sum_long, mag %in% other_mags)
best_sum_long <- subset(all_sum_long, mag %in% best_examples)

select_mag_labs <- (c(I4_MAG_00006 = "SJAQ100 sp016735685 assembled from Control B", I4_MAG_00065 = "Roseomonas sp. assembled from Control B", L3_MAG_00058 = "Prosthecobacter sp. assembled from Control C",
                      L7_MAG_00020 = "Sphingorhabdus_B sp. assembled from GBH C", L8_MAG_00011 = "UBA953 sp. assembled from GBH D", L8_MAG_00019 = "UA16 family assembled from GBH D"))

other_mag_labs <- (c(L2_MAG_00052 = "Erythrobacter sp. assembled from GBH A", L4_MAG_00099 = "Bosea sp001713455 assembled from Control D", L7_MAG_00028 = "SYFN01 sp. assembled from GBH C",
                     L7_MAG_00043 = "Luteolibacter sp. assembled from GBH C", L8_MAG_00042 = "UBA4660 sp. assembled from GBH D"))

best_mag_labs <- (c(I4_MAG_00006 = "SJAQ100 sp016735685 assembled from Control B", I4_MAG_00065 = "Roseomonas sp. assembled from Control B", L3_MAG_00058 = "Prosthecobacter sp. assembled from Control C",
                    L7_MAG_00020 = "Sphingorhabdus_B sp. assembled from GBH C", L8_MAG_00011 = "UBA953 sp. assembled from GBH D", L8_MAG_00019 = "UA16 family assembled from GBH D"))
#maybe sweep
select_snv_heat <- ggplot(select_snv, aes(x = graph_name, y = reorder(groups, all_mean), fill = final_ref_freq)) +
  geom_tile()+
  scale_fill_viridis(direction = -1, na.value = "white") +
  theme_classic() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_blank(), text = element_text(size = 30), axis.title.y = element_text(face = "bold"),
        strip.text.x.top = element_text(size = 25, face = "bold"), panel.spacing = unit(1, "cm"), axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"), legend.title = element_text(face = "bold"))+
  labs(fill = "Ref. Freq.", y = "Genome Position")+
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag, nrow = 1, ncol = 6, scales = "free", labeller = labeller(mag = select_mag_labs))
ggsave("select_snv_heat.png", limitsize = F, dpi = 100, width = 48, height = 10)

select_snv_sum <- ggplot(subset(select_sum_long, class == "SNVs"), aes(x = graph_name, y=((divergent_sites/mag_length)*10^6), fill = treatment))+
  geom_bar(stat = "identity")+ 
  theme_classic()+
  scale_fill_manual(values= c("grey", "darkred", "#FFB500"))+
  theme(axis.title.x = element_blank(), strip.text.x = element_blank(), axis.text = element_text(color = "black"), axis.title.y = element_text(face = "bold"),
        text = element_text(size = 30), panel.spacing = unit(0.05, "cm"), axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_text(face = "bold"))+
  labs(y = "SNVs / Mbp", fill = "Treatment")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag, nrow = 1, scales = "free")
ggsave("select_snv_sum.png", limitsize = F, width = 48, height = 10)

select_snv_frac <- ggplot(select_sum, aes(x = graph_name, y=SNSs/(SNSs+SNVs)))+
  geom_bar(stat="identity", fill= "#424242")+ 
  theme_classic()+
  theme(strip.text.x = element_blank(), axis.text = element_text(color = "black"), axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face = "bold"),
        text = element_text(size = 30), panel.spacing = unit(0.4, "cm"), axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank())+
  labs(y = "Fraction of SNVs dominated by a single allele", x = "Pond")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0,1))+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag, nrow = 1, scales = "free")
ggsave("select_snv_frac.png", limitsize = F, width = 48, height = 10)

select_all <- select_snv_heat / select_snv_sum / select_snv_frac
save_plot("select_all.jpeg", select_all, base_height = 5.5, base_width = 20, ncol = 3, nrow = 6, dpi = 300, limitsize = F)

#no sweep
other_snv_heat <- ggplot(other_snv, aes(x = graph_name, y = reorder(groups, all_mean), fill = final_ref_freq)) +
  geom_tile()+
  scale_fill_viridis(direction = -1, na.value = "white") +
  theme_classic() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_blank(), text = element_text(size = 25), axis.title.y = element_text(face = "bold"),
        strip.text.x.top = element_text(size = 25, face = "bold"), panel.spacing = unit(1.5, "cm"), axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"), legend.title = element_text(face = "bold"))+
  labs(fill = "Ref. Freq.", y = "Genome Position")+
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag, nrow = 1, ncol = 6, scales = "free", labeller = labeller(mag = other_mag_labs))
ggsave("other_snv_heat.png", limitsize = F, dpi = 100, width = 48, height = 10)

other_snv_sum <- ggplot(subset(other_sum_long, class == "SNVs"), aes(x = graph_name, y=((divergent_sites/mag_length)*10^6), fill = treatment))+
  geom_bar(stat = "identity")+ 
  theme_classic()+
  scale_fill_manual(values= c("grey", "darkred", "#FFB500"))+
  theme(axis.title.x = element_blank(), strip.text.x = element_blank(), axis.text = element_text(color = "black"), axis.title.y = element_text(face = "bold"),
        text = element_text(size = 25), panel.spacing = unit(0.5, "cm"), axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_text(face = "bold"))+
  labs(y = "SNVs / Mbp", fill = "Treatment")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag, nrow = 1, scales = "free")
ggsave("other_snv_sum.png", limitsize = F, width = 48, height = 10)

other_snv_frac <- ggplot(other_sum, aes(x = graph_name, y=SNSs/(SNSs+SNVs)))+
  geom_bar(stat="identity", fill= "#424242")+ 
  theme_classic()+
  theme(strip.text.x = element_blank(), axis.text = element_text(color = "black"), axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face = "bold"),
        text = element_text(size = 25), panel.spacing = unit(0.4, "cm"), axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank())+
  labs(y = "Fraction of SNVs dominated by a single allele", x = "Pond")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0,1))+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag, nrow = 1, scales = "free")
ggsave("other_snv_frac.png", limitsize = F, width = 48, height = 10)

other_all <- other_snv_heat / other_snv_sum / other_snv_frac
save_plot("other_all.jpeg", other_all, base_height = 5, base_width = 15, ncol = 3, nrow = 5, dpi = 300)


#best example
# best_heat <- ggplot(best_snv, aes(x = graph_name, y = reorder(groups, all_mean), fill = final_ref_freq)) +
#   geom_tile()+
#   scale_fill_viridis(direction = -1, na.value = "white") +
#   theme_classic() +
#   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), 
#         text = element_text(size = 17, face = "bold"), strip.text.x.top = element_text(size = 20), panel.spacing = unit(1, "cm"))+
#   labs(legend = "Reference Frequency")+
#   guides(fill = guide_legend(reverse = TRUE))+
#   scale_x_discrete(expand = c(0, 0))+
#   facet_wrap(~mag, nrow = 1, ncol = 6, scales = "free", labeller = labeller(mag = best_mag_labs))
# ggsave("best_snv_heat.png", limitsize = F, dpi = 100, width = 48, height = 8)
# 
# best_snv_sum <- ggplot(subset(best_sum_long, class == "SNVs"), aes(x = graph_name, y=((divergent_sites/mag_length)*10^6), fill = treatment))+
#   geom_bar(stat = "identity")+ 
#   theme_classic()+
#   scale_fill_manual(values= c("grey", "tomato3", "palegreen4"))+
#   theme(legend.title = element_blank(), axis.title.x = element_blank(), strip.background = element_blank(),
#         strip.text.x = element_blank(), text = element_text(size = 17, face = "bold"), panel.spacing = unit(0.05, "cm"))+
#   labs(y = "SNVs / Mbp")+
#   scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
#   scale_x_discrete(expand = c(0, 0))+
#   facet_wrap(~mag, nrow = 1, scales = "free")
# ggsave("best_snv_sum.png", limitsize = F, width = 48, height = 8)
# 
# best_snv_frac <- ggplot(best_sum, aes(x = graph_name, y=SNSs/(SNSs+SNVs)))+
#   geom_bar(stat="identity", fill= "#424242")+ 
#   theme_classic()+
#   theme(legend.title = element_blank(), strip.background = element_blank(), strip.text.x.top = element_blank(),
#         text = element_text(size = 17, face = "bold"), panel.spacing = unit(0.4, "cm"))+
#   labs(y = "fraction of SNVs dominated by a single allele", x = "Pond")+
#   scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0,1))+
#   scale_x_discrete(expand = c(0, 0))+
#   facet_wrap(~mag, nrow = 1, scales = "free")
# ggsave("best_snv_frac.png", limitsize = F, width = 48, height = 8)
# 
# best_all <- best_heat / best_snv_sum / best_snv_frac
# save_plot("best_all.jpeg", best_all, base_height = 25, base_width = 10, ncol = 2, nrow = 1, dpi = 300, limitsize = F)


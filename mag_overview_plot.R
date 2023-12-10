library(tidyverse)
library(dplyr)
library(ggplot2)
library(cowplot)


mag_list <- list("I4_MAG_00006", "I4_MAG_00065", "L2_MAG_00052", "L3_MAG_00058", "L4_MAG_00099",
                 "L7_MAG_00020", "L7_MAG_00028", "L7_MAG_00043", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")

mags <- read_csv("all_mags_subsamp.csv")
mags <- subset(mags, new_time == 2)
mags <- mags %>% mutate(mag_name = case_when(mag == "I4_MAG_00006" ~"SJAQ100 sp016735685", mag == "I4_MAG_00006" ~ "SJAQ100 sp016735685", mag == "I4_MAG_00065" ~ "Roseomonas sp.", mag == "L2_MAG_00052" ~ "Erythrobacter sp.",
                                             mag == "L3_MAG_00058" ~ "Prosthecobacter sp.", mag == "L4_MAG_00099" ~ "Bosea sp001713455", mag == "L7_MAG_00020" ~ "Sphingorhabdus_B sp.", mag == "L7_MAG_00028" ~ "SYFN01 sp.",
                                             mag == "L7_MAG_00043" ~ "Luteolibacter sp.", mag == "L8_MAG_00011" ~ "UBA953 sp.", mag == "L8_MAG_00019" ~ "UA16", mag == "L8_MAG_00042" ~ "UBA4660 sp."))

mags$treatment <- factor(mags$treatment, levels = c("control", "phosphorus", "glyphosate"))
mag_plot <- ggplot(mags, aes(x = mag_name, y = name))+
  geom_point(size = 5.5, colour = "black")+  
  geom_point(aes(colour = treatment), size = 4.5)+
  scale_colour_manual(labels = c('Control', 'Phosphorus', 'GBH'), values= c("white", "grey70", "grey30"))+
  theme_classic()+
  labs(y = "Pond")+
  scale_y_discrete(limits=rev)+
  theme(axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

save_plot("mag_overview_subsamp.jpeg", mag_plot, dpi = 500, base_height = 4, base_width = 8)

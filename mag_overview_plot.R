library(tidyverse)
library(dplyr)
library(ggplot2)
library(cowplot)

mags <- read_csv("all_mags_subsamp.csv")
mags <- subset(mags, new_time == 2)
mags <- mags %>% mutate(mag_name = case_when(mag == "I4_MAG_00006" ~ "Burkholderiaceae 1", mag == "I4_MAG_00065" ~ "Roseomonas_A", mag == "L2_MAG_00052" ~ "Erythrobacter",
                                             mag == "L3_MAG_00058" ~ "Prosthecobacter", mag == "L4_MAG_00099" ~ "Bosea sp001713455", mag == "L7_MAG_00020" ~ "Sphingorhabdus_B",
                                             mag == "L7_MAG_00028" ~ "Burkholderiaceae 2", mag == "L7_MAG_00043" ~ "Luteolibacter", mag == "L8_MAG_00011" ~ "Verrucomicrobiae",
                                             mag == "L8_MAG_00019" ~ "Flavobacteriales 1", mag == "L8_MAG_00042" ~ "Flavobacteriales 2"))

mags$mag_order = factor(mags$mag_name, levels=c('Burkholderiaceae 1', 'Burkholderiaceae 2', 'Verrucomicrobiae', 'Flavobacteriales 1', 'Flavobacteriales 2',
                                                'Roseomonas_A', 'Prosthecobacter', 'Sphingorhabdus_B',  'Luteolibacter',
                                                'Erythrobacter', 'Bosea sp001713455'))

mags$treatment <- factor(mags$treatment, levels = c("control", "phosphorus", "glyphosate"))
mag_plot <- ggplot(mags, aes(x = mag_order, y = name))+
  geom_point(size = 5.5, colour = "black")+  
  geom_point(aes(colour = treatment), size = 4.5)+
  scale_colour_manual(labels = c('Control', 'Phosphorus', 'GBH'), values= c("white", "grey70", "grey30"))+
  theme_classic()+
  labs(y = "Pond", x = "MAG")+
  scale_y_discrete(limits=rev)+
  theme(axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

save_plot("mag_overview_subsamp.jpeg", mag_plot, dpi = 500, base_height = 4, base_width = 8)

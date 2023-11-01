library(tidyverse)
library(dplyr)
library(ggplot2)
library(cowplot)


mag_list <- list("I4_MAG_00006", "I4_MAG_00065", "L2_MAG_00052", "L3_MAG_00058", "L4_MAG_00099",
                 "L7_MAG_00020", "L7_MAG_00028", "L7_MAG_00043", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")

mags <- read_csv("ANI_95_all_mags.csv")
mags <- subset(mags, mag %in% mag_list & new_time == 2)

mags$mag_name <- with(mags, ifelse(mags$mag == "I4_MAG_00006", "SJAQ100 sp016735685", NA)) 
mags$mag_name <- with(mags, ifelse(mags$mag == "I4_MAG_00065", "Roseomonas sp.", mag_name)) 
mags$mag_name <- with(mags, ifelse(mags$mag == "L2_MAG_00052", "Erythrobacter sp.", mag_name)) 
mags$mag_name <- with(mags, ifelse(mags$mag == "L3_MAG_00058", "Prosthecobacter sp.", mag_name)) 
mags$mag_name <- with(mags, ifelse(mags$mag == "L4_MAG_00099", "Bosea sp001713455", mag_name)) 
mags$mag_name <- with(mags, ifelse(mags$mag == "L7_MAG_00020", "Sphingorhabdus_B sp.", mag_name)) 
mags$mag_name <- with(mags, ifelse(mags$mag == "L7_MAG_00028", "SYFN01 sp.", mag_name)) 
mags$mag_name <- with(mags, ifelse(mags$mag == "L7_MAG_00043", "Luteolibacter sp.", mag_name)) 
mags$mag_name <- with(mags, ifelse(mags$mag == "L8_MAG_00011", "UBA953 sp.", mag_name)) 
mags$mag_name <- with(mags, ifelse(mags$mag == "L8_MAG_00019", "UA16 family", mag_name)) 
mags$mag_name <- with(mags, ifelse(mags$mag == "L8_MAG_00042", "UBA4660 sp.", mag_name)) 


mag_plot <- ggplot(mags, aes(x = mag_name, y = name, colour = treatment))+
  geom_point(size = 4)+  
  theme_classic()+
  labs(y = "Pond", x = "MAG", colour = "Treatment")+
  scale_colour_manual(labels = c("Control", "GBH (15mg/L)", "Phosphorus (320ug/L)"),
                      values = c("grey", "darkred", "#FFB500"))+
  scale_y_discrete(limits=rev)+
  theme(axis.text = element_text(colour = "black"), axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")


save_plot("mag_overview.jpeg", mag_plot, dpi = 300, base_asp = 3)


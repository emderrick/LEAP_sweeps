library(tidyverse)


options(scipen=999)
setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")
mag_info <- read_csv("data files/all_mag_info.csv")
mag_info <- subset(mag_info, mag_coverage >= 5 & mag_breadth >= 0.5)
mag_info$Treatment <- ifelse(grepl("CTRL", mag_info$Name), "Control", "GBH")


ggplot(mag_info, aes(x = popANI_reference))+
  geom_density()+
  theme_classic()+
  labs(x = "population ANI", y = "Number of MAGs")

ggplot(mag_info, aes(x = popANI_reference))+
  geom_density()+
  theme_classic()+
  labs(x = "population ANI", y = "Number of MAGs")+
  facet_grid(rows = vars(time), cols = vars(Treatment))

ggplot(mag_info, aes(x = conANI_reference))+
  geom_histogram()+
  theme_classic()+
  labs(x = "population ANI", y = "Number of MAGs")+
  facet_grid(rows = vars(time), cols = vars(Treatment))

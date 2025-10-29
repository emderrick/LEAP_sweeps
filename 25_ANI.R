library(tidyverse)
options(scipen=999)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

mag_list <- c("MAG_00097_1", "MAG_00110_1", "MAG_00179_1", "MAG_00194_1", "MAG_00197_1", "MAG_00201_1", "MAG_00674_1")

mag_info <- read_csv("data files/all_mag_info.csv")
mag_info <- subset(mag_info, mag_coverage >= 5 & mag_breadth >= 0.7)
mag_info$Treatment <- ifelse(grepl("CTRL", mag_info$Name), "Control", "GBH")

ggplot(mag_info, aes(x = popANI_reference))+
  geom_histogram()+
  theme_classic()+
  labs(x = "ANI", y = "Number of MAGs")

ggplot(mag_info, aes(x = popANI_reference))+
  geom_histogram()+
  theme_classic()+
  labs(x = "ANI", y = "Number of MAGs")+
  facet_grid(rows = vars(time), cols = vars(Treatment))

mag_info <- subset(mag_info, mag %in% mag_list)
ani <- mag_info %>% group_by(mag) %>% summarise(range(popANI_reference))

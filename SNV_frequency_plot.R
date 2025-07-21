library(tidyverse)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

mag_list <- c("T1_bin.305", "T1_bin.609", "T1_bin.676") 

snv_freq_wide <- read_csv("data files/mag_snvs_avg_frequency_bin.609.csv")
snv_freq_wide$GBH_change <- snv_freq_wide$GBH_28 - snv_freq_wide$GBH_0
snv_freq_wide$GBH_change_abs <- abs(snv_freq_wide$GBH_change)
snv_freq_wide$CTRL_change <- snv_freq_wide$CTRL_28 - snv_freq_wide$CTRL_0
snv_freq_wide$CTRL_change_abs <- abs(snv_freq_wide$CTRL_change)
snv_freq <- pivot_longer(snv_freq_wide, cols = c("CTRL_0", "GBH_0", "CTRL_28", "GBH_28"), names_to = "Pond_Treatment", values_to = "Reference_Frequency")
snv_freq$Treatment <- with (snv_freq, ifelse(grepl("CTRL", snv_freq$Pond_Treatment), "Control", "GBH"))

test <- subset(snv_freq, Treatment == "GBH")

ggplot(test, aes(x = Pond_Treatment, y = Reference_Frequency, colour = GBH_change, group = reorder(group, GBH_change_abs)))+
  geom_point()+
  geom_line()+
  scale_colour_gradientn(colours = brewer.pal(11, "RdBu"), limits = c(-1,1))

test_CTRL <- subset(snv_freq, Treatment == "Control")

ggplot(test_CTRL, aes(x = Pond_Treatment, y = Reference_Frequency, colour = CTRL_change, group = reorder(group, CTRL_change_abs)))+
  geom_point()+
  geom_line()+
  scale_colour_gradientn(colours = brewer.pal(11, "RdBu"), limits = c(-1,1))

library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(patchwork)

sweep <- read_csv("small_L4_MAG_00099_sweep.csv")
threshold_snvs <- read_csv("threshold_snvs.csv")

sweep_hor <- spread(sweep, key=new_name, value=final_ref_freq)
sweep_hor$control_mean <-rowMeans(sweep_hor[,c(8)], na.rm=T)
sweep_hor$GBH_mean <-rowMeans(sweep_hor[,c(9:12)], na.rm=T)
sweep_hor$abs_val<- abs(sweep_hor$control_mean - sweep_hor$GBH_mean)

L4_threshold <- subset(threshold_snvs, mag=="L4_MAG_00099")
sweep_test <- full_join(sweep_hor, L4_threshold, by = c('groups'))
sweep_test$pass <- with(sweep_test, ifelse(is.na(pass), 'no', pass))

L4_MAG_00099_sweep <- ggplot(sweep_hor, aes(x=groups, y=abs_val))+
  geom_point()+
  theme_classic()

ggsave("L4_MAG_00099_sweep.png", limitsize=F, width=64, height=12)

filtered_SNVs <- read_csv("filtered_ANI_95_mag_SNVs.csv")


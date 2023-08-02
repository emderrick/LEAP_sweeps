library(tidyverse)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggpubr)
library(stringr)
library(patchwork)

all_snv <- read_csv("all_MAG_SNVs_med_July25.csv")
small_snv <- all_snv[c(1:7,11)] %>% distinct()

L3_MAG_00058_all <- read_csv("merged_L3_MAG_00058_SNVs.csv")
L3_MAG_00058_sweep <- ggplot(subset(small_snv, mag=="L3_MAG_00058"), aes(x=groups, y=abs_val, group=1))+
  geom_line()+
  theme_classic()+
  scale_x_discrete(labels=all_snv$scaffold.x)

ggsave("L3_MAG_00058_sweep.png", limitsize=F, width=64, height=12)


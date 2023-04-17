library(tidyverse)
library(ggplot2)
library(dplyr)
library(viridis)
library(stringr)

setwd("/Users/Emma/OneDrive - McGill University/Documents/phd docs/Chapter 1/Aim 1A")

#create list of MAGs I'm using
mag_list<-list("L3_MAG_00058", "I8_MAG_00005", "L4_MAG_00099", "L8_MAG_00019", "L8_MAG_00011",
               "L7_MAG_00043", "L7_MAG_00028")

#PLOTTING SCAFFOLD COVERAGE 
#load scaffold summary
mag_scaffolds <- read_csv("mag_scaffolds.csv")
mag_scaf <- filter(mag_scaffolds, mag %in% mag_list)
mag_scaf<-filter(mag_scaf, coverage >=4)
mag_scaf$group<- paste(mag_scaf$mag, "in pond", mag_scaf$pond)

#manually filtering out timepoints I'm not using for certain MAGs
mag_scaf <- mag_scaf %>% filter(!(mag=="I8_MAG_00005" & time=="0"))
mag_scaf <- mag_scaf %>% filter(!(mag=="I8_MAG_00005" & time=="1"))
mag_scaf <- mag_scaf %>% filter(!(mag=="L3_MAG_00058" & time=="0"))
mag_scaf <- mag_scaf %>% filter(!(mag=="L3_MAG_00058" & time=="2"))
mag_scaf <- mag_scaf %>% filter(!(mag=="L7_MAG_00028" & time=="2"))
mag_scaf <- mag_scaf %>% filter(!(mag=="L7_MAG_00043" & time=="2"))
mag_scaf <- mag_scaf %>% filter(!(mag=="L8_MAG_00011" & time=="2"))

#scaffold coverage for all MAGs together
ggplot(mag_scaf, aes(x = coverage, y=((SNV_count/length)*10^6), colour = group)) +
  geom_point()+
  labs(y="SNVs / MBp", x="Coverage (x)") +
  theme_classic()
ggsave("no_L8_11_scaffold_coverage.png", limitsize = FALSE, width=20, height=4, dpi=500) 

#individual coverage graphs for each MAG
for (MAG in mag_list){
  ggplot(subset(mag_scaf, mag == MAG), aes(x = coverage, y=((SNV_count/length)*10^6), colour = group)) +
    geom_point()+
    labs(y="SNVs / MBp", x="Coverage (x)") +
    theme_classic()
  ggsave(filename=paste("scaffold_coverage_",MAG,".png",sep=""), limitsize = FALSE)
}
 

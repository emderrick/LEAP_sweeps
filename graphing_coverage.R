library(tidyverse)
library(ggplot2)
library(dplyr)
library(viridis)
library(stringr)
library(ggpubr)

setwd("/Users/Emma/OneDrive - McGill University/Documents/phd docs/Chapter 1/Aim 1A")

mag_list<-list("L3_MAG_00058", "I8_MAG_00005", "L4_MAG_00099", "L8_MAG_00019", "L8_MAG_00011",
               "L7_MAG_00043", "L7_MAG_00028", "I4_MAG_00006", "I4_MAG_00065", "L2_MAG_00052",
               "L7_MAG_00020", "L8_MAG_00042", "L2_MAG_00048")
#load mag snv info
mags <- read_csv("filtered_ANI_95_mag_SNVs.csv")

#sum up number of SNVs by scaffold
mag_cov <- aggregate(mags$number_SNVs, by=list(mag=mags$mag, scaffold=mags$scaffold, length=mags$length, group=mags$group, coverage=mags$coverage), FUN=sum) 
mag_cov <- mag_cov %>% rename("SNVs"="x")
mag_cov$snvbylength <- mag_cov$SNVs/mag_cov$length

#scaffold coverage for all together
ggplot(mag_cov, aes(x = coverage, y=snvbylength*10^6, group=scaffold)) + 
  geom_point(aes(colour=group))+
  scale_colour_viridis(discrete = T, option= "magma")+
  labs(y="SNVs / MBp", x="Coverage (x)")+
  theme_classic() +
  theme(legend.title=element_blank())
ggsave("all_scaffold_coverage.png", limitsize = FALSE, dpi=500) 

#individual coverage graphs
L7_MAG_00028_cov <- ggplot(subset(mag_cov, mag == "L7_MAG_00028"), aes(x = coverage, y=snvbylength*10^6)) + 
  geom_point(aes(colour=group))+
  scale_colour_viridis(discrete = T, option= "magma")+
  labs(y="SNVs / MBp", x="Coverage (x)") +
  ggtitle("L7_MAG_00028")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()+
  theme(legend.title=element_blank())
ggsave(filename="L7_MAG_00028_cov.png", limitsize = FALSE)

L8_MAG_00011_cov <- ggplot(subset(mag_cov, mag == "L8_MAG_00011"), aes(x = coverage, y=snvbylength*10^6)) + 
  geom_point(aes(colour=group))+
  scale_colour_viridis(discrete = T, option= "magma")+
  labs(y="SNVs / MBp", x="Coverage (x)") +
  ggtitle("L8_MAG_00011")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()+
  theme(legend.title=element_blank())
ggsave(filename="L8_MAG_00011_cov.png", limitsize = FALSE)

L8_MAG_00019_cov <- ggplot(subset(mag_cov, mag == "L8_MAG_00019"), aes(x = coverage, y=snvbylength*10^6)) + 
  geom_point(aes(colour=group))+
  scale_colour_viridis(discrete = T, option= "magma")+
  labs(y="SNVs / MBp", x="Coverage (x)") +
  ggtitle("L8_MAG_00019")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()+
  theme(legend.title=element_blank())
ggsave(filename="L8_MAG_00019_cov.png", limitsize = FALSE)

L7_MAG_00043_cov <-  ggplot(subset(mag_cov, mag == "L7_MAG_00043"), aes(x = coverage, y=snvbylength*10^6)) + 
  geom_point(aes(colour=group))+
  scale_colour_viridis(discrete = T, option= "magma")+
  labs(y="SNVs / MBp", x="Coverage (x)") +
  ggtitle("L7_MAG_00043")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()+
  theme(legend.title=element_blank())
ggsave(filename="L7_MAG_00043_cov.png", limitsize = FALSE)

L3_MAG_00058_cov <-  ggplot(subset(mag_cov, mag == "L3_MAG_00058"), aes(x = coverage, y=snvbylength*10^6)) + 
  geom_point(aes(colour=group))+
  scale_colour_viridis(discrete = T, option= "magma")+
  labs(y="SNVs / MBp", x="Coverage (x)") +
  ggtitle("L3_MAG_00058")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()+
  theme(legend.title=element_blank())
ggsave(filename="L3_MAG_00058_cov.png", limitsize = FALSE)

I8_MAG_00005_cov <-  ggplot(subset(mag_cov, mag == "I8_MAG_00005"), aes(x = coverage, y=snvbylength*10^6)) + 
  geom_point(aes(colour=group))+
  scale_colour_viridis(discrete = T, option= "magma")+
  labs(y="SNVs / MBp", x="Coverage (x)") +
  ggtitle("I8_MAG_00005")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()+
  theme(legend.title=element_blank())
ggsave(filename="I8_MAG_00005_cov.png", limitsize = FALSE)

L4_MAG_00099_cov <- ggplot(subset(mag_cov, mag == "L4_MAG_00099"), aes(x = coverage, y=snvbylength*10^6)) + 
  geom_point(aes(colour=group))+
  scale_colour_viridis(discrete = T, option= "magma")+
  labs(y="SNVs / MBp", x="Coverage (x)") +
  ggtitle("L4_MAG_00099")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()+
  theme(legend.title=element_blank())
ggsave(filename="L4_MAG_00099_cov.png", limitsize = FALSE)

I4_MAG_00006_cov <-  ggplot(subset(mag_cov, mag == "I4_MAG_00006"), aes(x = coverage, y=snvbylength*10^6)) + 
  geom_point(aes(colour=group))+
  scale_colour_viridis(discrete = T, option= "magma")+
  labs(y="SNVs / MBp", x="Coverage (x)") +
  ggtitle("I4_MAG_00006")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()+
  theme(legend.title=element_blank())
ggsave(filename="I4_MAG_00006_cov.png", limitsize = FALSE)

I4_MAG_00065_cov <-  ggplot(subset(mag_cov, mag == "I4_MAG_00065"), aes(x = coverage, y=snvbylength*10^6)) + 
  geom_point(aes(colour=group))+
  scale_colour_viridis(discrete = T, option= "magma")+
  labs(y="SNVs / MBp", x="Coverage (x)") +
  ggtitle("I4_MAG_00065")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()+
  theme(legend.title=element_blank())
ggsave(filename="I4_MAG_00065_cov.png", limitsize = FALSE)

L2_MAG_00052_cov <-  ggplot(subset(mag_cov, mag == "L2_MAG_00052"), aes(x = coverage, y=snvbylength*10^6)) + 
  geom_point(aes(colour=group))+
  scale_colour_viridis(discrete = T, option= "magma")+
  labs(y="SNVs / MBp", x="Coverage (x)") +
  ggtitle("L2_MAG_00052")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()+
  theme(legend.title=element_blank())
ggsave(filename="L2_MAG_00052_cov.png", limitsize = FALSE)

L7_MAG_00020_cov <-  ggplot(subset(mag_cov, mag == "L7_MAG_00020"), aes(x = coverage, y=snvbylength*10^6)) + 
  geom_point(aes(colour=group))+
  scale_colour_viridis(discrete = T, option= "magma")+
  labs(y="SNVs / MBp", x="Coverage (x)") +
  ggtitle("L7_MAG_00020")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()+
  theme(legend.title=element_blank())
ggsave(filename="L7_MAG_00020_cov.png", limitsize = FALSE)

L8_MAG_00042_cov <-  ggplot(subset(mag_cov, mag == "L8_MAG_00042"), aes(x = coverage, y=snvbylength*10^6)) + 
  geom_point(aes(colour=group))+
  scale_colour_viridis(discrete = T, option= "magma")+
  labs(y="SNVs / MBp", x="Coverage (x)") +
  ggtitle("L8_MAG_00042")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()+
  theme(legend.title=element_blank())
ggsave(filename="L8_MAG_00042_cov.png", limitsize = FALSE)

L2_MAG_00048_cov <-  ggplot(subset(mag_cov, mag == "L2_MAG_00048"), aes(x = coverage, y=snvbylength*10^6)) + 
  geom_point(aes(colour=group))+
  scale_colour_viridis(discrete = T, option= "magma")+
  labs(y="SNVs / MBp", x="Coverage (x)") +
  ggtitle("L2_MAG_00048")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()+
  theme(legend.title=element_blank())
ggsave(filename="L2_MAG_00048_cov.png", limitsize = FALSE)


ggarrange(L3_MAG_00058_cov, I8_MAG_00005_cov, L4_MAG_00099_cov, L8_MAG_00019_cov, L8_MAG_00011_cov,
          L7_MAG_00043_cov, L7_MAG_00028_cov, I4_MAG_00006_cov, I4_MAG_00065_cov, 
          L7_MAG_00020_cov, L8_MAG_00042_cov, L2_MAG_00048_cov, L2_MAG_00052_cov,
          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"),
          ncol = 3, nrow = 5)
ggsave(filename="95_individual_cov.png", limitsize = FALSE, width=16, height=10)


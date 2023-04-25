library(tidyverse)
library(ggplot2)
library(dplyr)
library(viridis)
library(ggpubr)
setwd("/Users/Emma/OneDrive - McGill University/Documents/phd docs/Chapter 1/Aim 1A")

#load mag info
summary_mags <- read_csv("ANI_95_mags.csv")
completed_mags <- complete(summary_mags, mag, pond)

#individual by MAG
L3_MAG_00058_pres <- ggplot(subset(completed_mags, mag=="L3_MAG_00058"), aes(x = new_time, y=pond, colour=treatment))+
  geom_point(size=4)+
  scale_colour_manual(values = c("control"="#F56B5CFF", "phosphorus"="#CC3F71FF", "glyphosate"="#3F0F72FF"))+
  theme_classic(base_size=13)+
  scale_x_discrete(limits=c("1", "2", "3"))+
  ggtitle("L3_MAG_00058")+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))

ggsave("L3_MAG_00058_time_pond.png", limitsize = FALSE)

I8_MAG_00005_pres <- ggplot(subset(completed_mags, mag=="I8_MAG_00005"), aes(x = new_time, y=pond, colour=treatment))+
  geom_point(size=4)+
  scale_colour_manual(values = c("control"="#F56B5CFF", "phosphorus"="#CC3F71FF", "glyphosate"="#3F0F72FF"))+
  theme_classic(base_size=13)+
  scale_x_discrete(limits=c("1", "2", "3"))+
  ggtitle("I8_MAG_00005")+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))

ggsave("I8_MAG_00005_time_pond.png", limitsize = FALSE)

L4_MAG_00099_pres <- ggplot(subset(completed_mags, mag=="L4_MAG_00099"), aes(x = new_time, y=pond, colour=treatment))+
  geom_point(size=4)+
  scale_colour_manual(values = c("control"="#F56B5CFF", "phosphorus"="#CC3F71FF", "glyphosate"="#3F0F72FF"))+
  theme_classic(base_size=13)+
  scale_x_discrete(limits=c("1", "2", "3"))+
  ggtitle("L4_MAG_00099")+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))

ggsave("L4_MAG_00099_time_pond.png", limitsize = FALSE)

L8_MAG_00019_pres <- ggplot(subset(completed_mags, mag=="L8_MAG_00019"), aes(x = new_time, y=pond, colour=treatment))+
  geom_point(size=4)+
  scale_colour_manual(values = c("control"="#F56B5CFF", "phosphorus"="#CC3F71FF", "glyphosate"="#3F0F72FF"))+
  theme_classic(base_size=13)+
  scale_x_discrete(limits=c("1", "2", "3"))+
  ggtitle("L8_MAG_00019")+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))

ggsave("L8_MAG_00019_time_pond.png", limitsize = FALSE)

L8_MAG_00011_pres <- ggplot(subset(completed_mags, mag=="L8_MAG_00011"), aes(x = new_time, y=pond, colour=treatment))+
  geom_point(size=4)+
  scale_colour_manual(values = c("control"="#F56B5CFF", "phosphorus"="#CC3F71FF", "glyphosate"="#3F0F72FF"))+
  theme_classic(base_size=13)+
  scale_x_discrete(limits=c("1", "2", "3"))+
  ggtitle("L8_MAG_00011")+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))

ggsave("L8_MAG_00011_time_pond.png", limitsize = FALSE)

L7_MAG_00043_pres <- ggplot(subset(completed_mags, mag=="L7_MAG_00043"), aes(x = new_time, y=pond, colour=treatment))+
  geom_point(size=4)+
  scale_colour_manual(values = c("control"="#F56B5CFF", "phosphorus"="#CC3F71FF", "glyphosate"="#3F0F72FF"))+
  theme_classic(base_size=13)+
  scale_x_discrete(limits=c("1", "2", "3"))+
  ggtitle("L7_MAG_00043")+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))

ggsave("L7_MAG_00043_time_pond.png", limitsize = FALSE)

L7_MAG_00028_pres <- ggplot(subset(completed_mags, mag=="L7_MAG_00028"), aes(x = new_time, y=pond, colour=treatment))+
  geom_point(size=4)+
  scale_colour_manual(values = c("control"="#F56B5CFF", "phosphorus"="#CC3F71FF", "glyphosate"="#3F0F72FF"))+
  theme_classic(base_size=13)+
  scale_x_discrete(limits=c("1", "2", "3"))+
  ggtitle("L7_MAG_00028")+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))

ggsave("L7_MAG_00028_time_pond.png", limitsize = FALSE)

I4_MAG_00006_pres <- ggplot(subset(completed_mags, mag=="I4_MAG_00006"), aes(x = new_time, y=pond, colour=treatment))+
  geom_point(size=4)+
  scale_colour_manual(values = c("control"="#F56B5CFF", "phosphorus"="#CC3F71FF", "glyphosate"="#3F0F72FF"))+
  theme_classic(base_size=13)+
  scale_x_discrete(limits=c("1", "2", "3"))+
  ggtitle("I4_MAG_00006")+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))

ggsave("I4_MAG_00006_time_pond.png", limitsize = FALSE)

I4_MAG_00065_pres <- ggplot(subset(completed_mags, mag=="I4_MAG_00065"), aes(x = new_time, y=pond, colour=treatment))+
  geom_point(size=4)+
  scale_colour_manual(values = c("control"="#F56B5CFF", "phosphorus"="#CC3F71FF", "glyphosate"="#3F0F72FF"))+
  theme_classic(base_size=13)+
  scale_x_discrete(limits=c("1", "2", "3"))+
  ggtitle("I4_MAG_00065")+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))

ggsave("I4_MAG_00065_time_pond.png", limitsize = FALSE)

L2_MAG_00052_pres <- ggplot(subset(completed_mags, mag=="L2_MAG_00052"), aes(x = new_time, y=pond, colour=treatment))+
  geom_point(size=4)+
  scale_colour_manual(values = c("control"="#F56B5CFF", "phosphorus"="#CC3F71FF", "glyphosate"="#3F0F72FF"))+
  theme_classic(base_size=13)+
  scale_x_discrete(limits=c("1", "2", "3"))+
  ggtitle("L2_MAG_00052")+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))

ggsave("L2_MAG_00052_time_pond.png", limitsize = FALSE)

L7_MAG_00020_pres <- ggplot(subset(completed_mags, mag=="L7_MAG_00020"), aes(x = new_time, y=pond, colour=treatment))+
  geom_point(size=4)+
  scale_colour_manual(values = c("control"="#F56B5CFF", "phosphorus"="#CC3F71FF", "glyphosate"="#3F0F72FF"))+
  theme_classic(base_size=13)+
  scale_x_discrete(limits=c("1", "2", "3"))+
  ggtitle("L7_MAG_00020")+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))

ggsave("L7_MAG_00020_time_pond.png", limitsize = FALSE)

L8_MAG_00042_pres < -ggplot(subset(completed_mags, mag=="L8_MAG_00042"), aes(x = new_time, y=pond, colour=treatment))+
  geom_point(size=4)+
  scale_colour_manual(values = c("control"="#F56B5CFF", "phosphorus"="#CC3F71FF", "glyphosate"="#3F0F72FF"))+
  theme_classic(base_size=13)+
  scale_x_discrete(limits=c("1", "2", "3"))+
  ggtitle("L8_MAG_00042")+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))

ggsave("L8_MAG_00042_time_pond.png", limitsize = FALSE)

L2_MAG_00048 <- ggplot(subset(completed_mags, mag=="L2_MAG_00048"), aes(x = new_time, y=pond, colour=treatment))+
  geom_point(size=4)+
  scale_colour_manual(values = c("control"="#F56B5CFF", "phosphorus"="#CC3F71FF", "glyphosate"="#3F0F72FF"))+
  theme_classic(base_size=13)+
  scale_x_discrete(limits=c("1", "2", "3"))+
  ggtitle("L2_MAG_00048")+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))

ggsave("L2_MAG_00048_time_pond.png", limitsize = FALSE)

ggarrange(L3_MAG_00058_pres, I8_MAG_00005_pres, L4_MAG_00099_pres, L8_MAG_00019_pres, L8_MAG_00011_pres,
          L7_MAG_00043_pres, L7_MAG_00028_pres, I4_MAG_00006_pres, I4_MAG_00065_pres,
          L7_MAG_00020_pres, L8_MAG_00042_pres, L2_MAG_00048_pres, L2_MAG_00052_pres,
          ncol = 5, nrow = 3)
ggsave(filename="individual_mag_pres.png", limitsize = FALSE, width=12, height=8)



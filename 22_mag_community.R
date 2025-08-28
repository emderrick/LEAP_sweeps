library(tidyverse)
library(ggplot2)
library(rstatix)
library(lme4)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

mags <- read_csv("data files/MAG_rel_abun.csv")
mags$Name_Time <- paste(mags$Name, mags$Time, sep = " ")
mags$Time <- ifelse(grepl("1", mags$Name_Time), "Day 0", "Day 28")

pond_community <- subset(mags, Relative_Abundance >= 0.01)
pond_community <- pond_community %>% group_by(Name, Time, Treatment) %>% count() %>% ungroup()
pond_community$Treatment_Time <- paste(pond_community$Treatment, pond_community$Time, sep = " ")

mag_community <- ggplot(pond_community, aes(x = Treatment, group = Treatment_Time, colour = Treatment_Time, fill = Treatment_Time, y = n))+
  geom_boxplot()+
  geom_point()+
  scale_fill_manual(values = c("#98BF64", "darkgreen", "#BE93D4", "darkmagenta"))+
  scale_colour_manual(values = c("darkgreen", "#0B4F02", "darkmagenta", "#5E0069"))+
  theme_classic()+
  theme(axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        text = element_text(size = 12, colour = "black"),
        legend.position = "none")+
  labs(y = "Number of Detected MAGs", x = "", fill = "", colour = "")+
  ylim(0,200)+
  facet_wrap(~Time)
ggsave("figures/mag_community.pdf", mag_community, limitsize = F, width = 7, height = 3)


mags_glm <- glm(n ~ Treatment * Time, family = poisson, data = pond_community)
summary(mags_glm)









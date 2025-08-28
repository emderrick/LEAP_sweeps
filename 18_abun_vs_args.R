library(tidyverse)
library(gghalves)
library(rstatix)
library(mgcv)
library(patchwork)

rel_abun <- read_csv("data files/MAG_avg_abun.csv")
rel_abun$avg_change <- rel_abun$Day_28_avg - rel_abun$Day_0_avg

all_efflux_abun <- ggplot(rel_abun, aes(x = as.factor(`antibiotic efflux`), y = log10(Day_28_avg), fill = Treatment))+
  geom_half_violin(side = "r")+
  geom_half_boxplot(side = "l")+
  scale_fill_manual(values = c("darkgreen", "darkmagenta"))+
  labs(x = " # Efflux ARGs", y = "log10(AVG MAG Abundance Day 28)")+
  theme_classic()+
  theme(axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        text = element_text(size = 12, colour = "black"),
        legend.position = "none")
ggsave("figures/all_efflux_rel_abun.pdf", all_efflux_abun, limitsize = F, width = 5, height = 4)

all_epsps_abun <- ggplot(rel_abun, aes(x = EPSPS_allele, y = log10(Day_28_avg)))+
  geom_half_violin(side = "r")+
  geom_half_boxplot(side = "l")+
  scale_fill_manual(values = c("darkgreen", "darkmagenta"))+
  labs(x = "", y = "")+
  theme_classic()+
  theme(axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        text = element_text(size = 12, colour = "black"))

ggsave("figures/all_epsps_rel_abun.pdf", all_epsps_abun, limitsize = F, width = 3.5, height = 4)

all_efflux_epsps <- all_efflux_abun + all_epsps_abun
ggsave("figures/all_efflux_epsps_abun.pdf", all_efflux_epsps, limitsize = F, width = 8.5, height = 4)

epsps_efflux_abun_stats <- glm(Day_28_avg ~ Treatment + EPSPS_allele + `antibiotic efflux`, data = rel_abun)
summary(epsps_efflux_abun_stats)

### change in abundance

all_efflux_change <- ggplot(rel_abun, aes(x = as.factor(`antibiotic efflux`), y = avg_change))+
  geom_half_violin(side = "r")+
  geom_half_boxplot(side = "l")+
  scale_fill_manual(values = c("darkgreen", "darkmagenta"))+
  labs(x = " # Efflux ARGs", y = "MAG abundance change")+
  theme_classic()+
  theme(axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        text = element_text(size = 12, colour = "black"),
        legend.position = "none")
ggsave("figures/all_efflux_change.pdf", all_efflux_change, limitsize = F, width = 5, height = 4)

all_epsps_change <- ggplot(rel_abun, aes(x = EPSPS_allele, y = avg_change))+
  geom_half_violin(side = "r")+
  geom_half_boxplot(side = "l")+
  scale_fill_manual(values = c("darkgreen", "darkmagenta"))+
  labs(x = "", y = "")+
  theme_classic()+
  theme(axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        text = element_text(size = 12, colour = "black"))
ggsave("figures/all_epsps_change.pdf", all_epsps_change, limitsize = F, width = 3.5, height = 4)

all_efflux_epsps_change <- all_efflux_change + all_epsps_change
ggsave("figures/all_efflux_epsps_change.pdf", all_efflux_epsps_change, limitsize = F, width = 8.5, height = 4)

epsps_efflux_change_stats <- glm(avg_change ~ Treatment + EPSPS_allele + `antibiotic efflux`, data = rel_abun)
summary(epsps_efflux_change_stats)
#plot(epsps_efflux_change_stats)




# with only GBH

rel_abun_GBH <- subset(rel_abun, Treatment == "GBH")

efflux_abun <- ggplot(rel_abun_GBH, aes(x = `antibiotic efflux`, y = log10(Day_28_avg)))+
  geom_violin(aes(group = `antibiotic efflux`), alpha = 0)+
  geom_jitter(aes(colour = EPSPS_allele), width = 0.2)+
  scale_colour_manual(values = c("darkmagenta", "darkgreen", "grey80"))+
  labs(x = " # Efflux ARGs", y = "", colour = "Classification")+
  theme_classic()+
  theme(axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        text = element_text(size = 12, colour = "black"))
ggsave("figures/efflux_rel_abun.pdf", efflux_abun, limitsize = F, width = 5, height = 4)

epsps_abun <- ggplot(rel_abun_GBH, aes(x = EPSPS_allele, y = log10(Day_28_avg)))+
  geom_violin(aes(group = EPSPS_allele), alpha = 0)+
  geom_jitter(aes(colour = EPSPS_allele), width = 0.2)+
  scale_colour_manual(values = c("darkmagenta", "darkgreen", "grey80"))+
  labs(x = "", y = "log10(AVG MAG Abundance Day 28)")+
  theme_classic()+
  theme(axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        text = element_text(size = 12, colour = "black"),
        legend.position = "none")
ggsave("figures/epsps_rel_abun.pdf", epsps_abun, limitsize = F, width = 3.5, height = 4)

efflux_epsps <- epsps_abun + efflux_abun
ggsave("figures/efflux_epsps_abun.pdf", efflux_epsps, limitsize = F, width = 8.5, height = 4)

#epsps_efflux_stats <- glm(Day_28_avg ~ EPSPS_allele + `antibiotic efflux`, data = rel_abun_GBH)
#summary(epsps_efflux_stats)

efflux_change <- ggplot(rel_abun_GBH, aes(x = `antibiotic efflux`, y = avg_change))+
  geom_violin(aes(group = `antibiotic efflux`), alpha = 0)+
  geom_jitter(aes(colour = EPSPS_allele), width = 0.2)+
  scale_colour_manual(values = c("darkmagenta", "darkgreen", "grey80"))+
  labs(x = " # Efflux ARGs", y = "", colour = "Classification")+
  theme_classic()+
  theme(axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        text = element_text(size = 12, colour = "black"))
ggsave("figures/efflux_change.pdf", efflux_change, limitsize = F, width = 5, height = 4)

epsps_change <- ggplot(rel_abun_GBH, aes(x = EPSPS_allele, y = avg_change))+
  geom_violin(aes(group = EPSPS_allele), alpha = 0)+
  geom_jitter(aes(colour = EPSPS_allele), width = 0.2)+
  scale_colour_manual(values = c("darkmagenta", "darkgreen", "grey80"))+
  labs(x = "", y = "MAG abundance change")+
  theme_classic()+
  theme(axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        text = element_text(size = 12, colour = "black"),
        legend.position = "none")
ggsave("figures/epsps_change.pdf", epsps_change, limitsize = F, width = 3.5, height = 4)

efflux_epsps_change <- epsps_change + efflux_change
ggsave("figures/efflux_epsps_change.pdf", efflux_epsps_change, limitsize = F, width = 8.5, height = 4)

#epsps_efflux_change_stats <- glm(avg_change ~ EPSPS_allele + `antibiotic efflux`, data = rel_abun_GBH)
#summary(epsps_efflux_change_stats)


library(tidyverse)
library(rstatix)
library(lme4)
library(patchwork)

rel_abun <- read_csv("data files/MAG_avg_abun.csv")
rel_abun$avg_change <- rel_abun$Day_28_avg - rel_abun$Day_0_avg
rel_abun$EPSPS_releveled <- factor(rel_abun$EPSPS_allele, levels = c("Sensitive","Unclassified","Resistant"))

epsps_abun_stats <- glm(Day_28_avg ~ Treatment * EPSPS_releveled, data = rel_abun)
summary(epsps_abun_stats)

epsps_change_stats <- glm(avg_change ~ Treatment * EPSPS_releveled, data = rel_abun)
summary(epsps_change_stats)

rel_abun_GBH <- subset(rel_abun, Treatment == "GBH")

epsps_abun_stats <- glm(Day_28_avg ~ EPSPS_releveled, data = rel_abun_GBH)
summary(epsps_abun_stats)

epsps_change_stats <- glm(avg_change ~ EPSPS_releveled, data = rel_abun_GBH)
summary(epsps_change_stats)


all_efflux_abun <- ggplot(rel_abun, aes(x = as.factor(`antibiotic efflux`), y = log10(Day_28_avg)))+
  geom_boxplot(aes(colour = Treatment), outliers = F)+
  geom_point(aes(colour = Treatment), position = position_jitterdodge(jitter.width = 0.2), size = 0.5)+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  labs(x = " # Efflux ARGs", y = "log10(AVG MAG Abundance Day 28)")+
  theme_classic()+
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 10, colour = "black"),
        legend.position = "none")
ggsave("figures/all_efflux_rel_abun.pdf", all_efflux_abun, units = "cm", width = 8, height = 7)

all_epsps_abun <- ggplot(rel_abun, aes(x = EPSPS_allele, y = log10(Day_28_avg)))+
  geom_boxplot(aes(colour = Treatment), outliers = T)+
  #geom_point(aes(colour = Treatment), position = position_jitterdodge(jitter.width = 0.2), size = 0.5)+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  labs(x = "", y = "")+
  theme_classic()+
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 10, colour = "black"))
ggsave("figures/all_epsps_rel_abun.pdf", all_epsps_abun, units = "cm", width = 9, height = 7)

all_efflux_epsps <- all_efflux_abun + all_epsps_abun
ggsave("figures/all_efflux_epsps_abun.pdf", all_efflux_epsps, units = "cm", width = 17, height = 7)

### change in abundance

all_efflux_change <- ggplot(rel_abun, aes(x = as.factor(`antibiotic efflux`), y = avg_change))+
  geom_boxplot(aes(colour = Treatment), outliers = F)+
  geom_point(aes(colour = Treatment), position = position_jitterdodge(jitter.width = 0.2), size = 0.5)+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  labs(x = " # Efflux ARGs", y = "MAG abundance change")+
  theme_classic()+
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 10, colour = "black"),
        legend.position = "none")
ggsave("figures/all_efflux_change.pdf", all_efflux_change, units = "cm", width = 17, height = 7)

all_epsps_change <- ggplot(rel_abun, aes(x = EPSPS_allele, y = avg_change))+
  geom_boxplot(aes(colour = Treatment), outliers = F)+
  geom_point(aes(colour = Treatment), position = position_jitterdodge(jitter.width = 0.2), size = 0.5)+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  labs(x = "", y = "")+
  theme_classic()+
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 10, colour = "black"))
ggsave("figures/all_epsps_change.pdf", all_epsps_change, units = "cm", width = 17, height = 7)

all_efflux_epsps_change <- all_efflux_change + all_epsps_change
ggsave("figures/all_efflux_epsps_change.pdf", all_efflux_epsps_change, units = "cm", width = 17, height = 7)


# with only GBH

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


library(tidyverse)
library(ggplot2)
library(rstatix)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

rel_abun <- read_csv("data files/MAG_rel_abun_change.csv")

mag_snvs <- read_csv("subsamp data files/T1_subsamp_SNV_summary_MAG.csv")
mag_snv_change <- mag_snvs[, c(1,14,17:19)]
mag_snv_change <- pivot_wider(mag_snv_change, names_from = Time, values_from = SNVs_Mbp)
mag_snv_change <- mag_snv_change[complete.cases(mag_snv_change ), ]
mag_snv_change$snv_change <- mag_snv_change$`Day 28` - mag_snv_change$`Day 0`

mag_sns_change <- mag_snvs[, c(1,13,17:19)]
mag_sns_change <- pivot_wider(mag_sns_change, names_from = Time, values_from = SNSs_Mbp)
mag_sns_change <- mag_sns_change[complete.cases(mag_sns_change ), ]
mag_sns_change$sns_change <- mag_sns_change$`Day 28` - mag_sns_change$`Day 0`

mag_n_snv_change <- mag_snvs[, c(1,16,17:19)]
mag_n_snv_change <- pivot_wider(mag_n_snv_change, names_from = Time, values_from = N_SNVs_Mbp)
mag_n_snv_change <- mag_n_snv_change[complete.cases(mag_n_snv_change ), ]
mag_n_snv_change$N_snv_change <- mag_n_snv_change$`Day 28` - mag_n_snv_change$`Day 0`

mag_n_sns_change <- mag_snvs[, c(1,15,17:19)]
mag_n_sns_change <- pivot_wider(mag_n_sns_change, names_from = Time, values_from = N_SNSs_Mbp)
mag_n_sns_change <- mag_n_sns_change[complete.cases(mag_n_sns_change ), ]
mag_n_sns_change$N_sns_change <- mag_n_sns_change$`Day 28` - mag_n_sns_change$`Day 0`

abun_snv_all <- left_join(rel_abun[, c(1:3,6:12)], mag_snv_change[,c(1,3,6)])
abun_snv_all <- left_join(abun_snv_all, mag_sns_change[, c(1,3,6)])
abun_snv_all <- left_join(abun_snv_all, mag_n_snv_change[, c(1,3,6)])
abun_snv_all <- left_join(abun_snv_all, mag_n_sns_change[, c(1,3,6)])
abun_snv_all <- abun_snv_all[complete.cases(abun_snv_all), ]
abun_snv_all$arg <- with(abun_snv_all, ifelse(total_hits >= 1, "Yes", "No"))

abun_snv_GBH <- subset(abun_snv_all, Treatment == "GBH")

gbh_abun_snv <- ggplot(abun_snv_GBH, aes(x = abun_change, y = snv_change, colour = EPSPS_allele))+
  geom_point()+
  geom_smooth(aes(group = Treatment), method = "glm", colour = "red", size = 0.5)+
  scale_colour_manual(values = c("darkmagenta", "darkgreen", "grey80"))+
  labs(x = "MAG abundance change", y = "SNV abundance change")+
  theme_classic()+
  theme(axis.text = element_text(size = 8, colour = "black"), legend.position = "none")

gbh_abun_sns <- ggplot(abun_snv_GBH, aes(x = abun_change, y = sns_change, colour = EPSPS_allele))+
  geom_point()+
  geom_smooth(aes(group = Treatment), method = "glm", colour = "red", size = 0.5)+
  scale_colour_manual(values = c("darkmagenta", "darkgreen", "grey80"))+
  labs(x = "MAG abundance change", y = "Fixed substitution abundance change", colour = "GBH resistance")+
  theme_classic()+
  theme(axis.text = element_text(size = 8, colour = "black"))

gbh_abun_snv_sns <- gbh_abun_snv + gbh_abun_sns
ggsave("figures/gbh_abun_snv_change.pdf", gbh_abun_snv_sns, limitsize = F, width = 8.5, height = 4)

gbh_abun_N_snv <- ggplot(abun_snv_GBH, aes(x = abun_change, y = N_snv_change, colour = EPSPS_allele))+
  geom_point()+
  geom_smooth(aes(group = Treatment), method = "glm", colour = "red", size = 0.5)+
  scale_colour_manual(values = c("darkmagenta", "darkgreen", "grey80"))+
  labs(x = "MAG abundance change", y = "Nonsyn SNV abundance change")+
  theme_classic()+
  theme(axis.text = element_text(size = 8, colour = "black"), legend.position = "none")

gbh_abun_N_sns <- ggplot(abun_snv_GBH, aes(x = abun_change, y = N_sns_change, colour = EPSPS_allele))+
  geom_point()+
  geom_smooth(aes(group = Treatment), method = "glm", colour = "red", size = 0.5)+
  scale_colour_manual(values = c("darkmagenta", "darkgreen", "grey80"))+
  labs(x = "MAG abundance change", y = "Nonsyn fixed substitution abundance change", colour = "GBH resistance")+
  theme_classic()+
  theme(axis.text = element_text(size = 8, colour = "black"))

gbh_abun_N_snv_sns <- gbh_abun_N_snv + gbh_abun_N_sns
ggsave("figures/gbh_nonsyn_abun_snv_change.pdf", gbh_abun_N_snv_sns, limitsize = F, width = 8.5, height = 4)


abun_snv_CTRL <- subset(abun_snv_all, Treatment == "Control")

ctrl_abun_snv <- ggplot(abun_snv_CTRL, aes(x = abun_change, y = snv_change, colour = EPSPS_allele))+
  geom_point()+
  geom_smooth(aes(group = Treatment), method = "glm", colour = "red", size = 0.5)+
  scale_colour_manual(values = c("darkmagenta", "darkgreen", "grey80"))+
  labs(x = "MAG abundance change", y = "SNV abundance change")+
  theme_classic()+
  theme(axis.text = element_text(size = 8, colour = "black"), legend.position = "none")

ctrl_abun_sns <- ggplot(abun_snv_CTRL, aes(x = abun_change, y = sns_change, colour = EPSPS_allele))+
  geom_point()+
  geom_smooth(aes(group = Treatment), method = "glm", colour = "red", size = 0.5)+
  scale_colour_manual(values = c("darkmagenta", "darkgreen", "grey80"))+
  labs(x = "MAG abundance change", y = "Fixed substitution abundance change", colour = "GBH resistance")+
  theme_classic()+
  theme(axis.text = element_text(size = 8, colour = "black"))

ctrl_abun_snv_sns <- ctrl_abun_snv + ctrl_abun_sns
ggsave("figures/ctrl_abun_snv_change.pdf", ctrl_abun_snv_sns, limitsize = F, width = 8.5, height = 4)

ctrl_abun_N_snv <- ggplot(abun_snv_CTRL, aes(x = abun_change, y = N_snv_change, colour = EPSPS_allele))+
  geom_point()+
  geom_smooth(aes(group = Treatment), method = "glm", colour = "red", size = 0.5)+
  scale_colour_manual(values = c("darkmagenta", "darkgreen", "grey80"))+
  labs(x = "MAG abundance change", y = "Nonsyn SNV abundance change")+
  theme_classic()+
  theme(axis.text = element_text(size = 8, colour = "black"), legend.position = "none")

ctrl_abun_N_sns <- ggplot(abun_snv_CTRL, aes(x = abun_change, y = N_sns_change, colour = EPSPS_allele))+
  geom_point()+
  geom_smooth(aes(group = Treatment), method = "glm", colour = "red", size = 0.5)+
  scale_colour_manual(values = c("darkmagenta", "darkgreen", "grey80"))+
  labs(x = "MAG abundance change", y = "Nonsyn fixed substitution abundance change", colour = "GBH resistance")+
  theme_classic()+
  theme(axis.text = element_text(size = 8, colour = "black"))

ctrl_abun_N_snv_sns <- ctrl_abun_N_snv + ctrl_abun_N_sns
ggsave("figures/ctrl_nonsyn_abun_snv_change.pdf", ctrl_abun_N_snv_sns, limitsize = F, width = 8.5, height = 4)



gbh_sns_abun <- glm(N_sns_change ~ abun_change, data = abun_snv_GBH)
summary(gbh_sns_abun)

ctrl_sns_abun <- glm(N_sns_change ~ abun_change, data = abun_snv_CTRL)
summary(ctrl_sns_abun)

sns_abun <- glm(N_sns_change ~ abun_change + Treatment, data = abun_snv_all)
summary(sns_abun)


gbh_snv_abun <- glm(N_snv_change ~ abun_change, data = abun_snv_GBH)
summary(gbh_snv_abun)

ctrl_snv_abun <- glm(N_snv_change ~ abun_change, data = abun_snv_CTRL)
summary(ctrl_snv_abun)

snv_abun <- glm(N_snv_change ~ abun_change + Treatment, data = abun_snv_all)
summary(snv_abun)




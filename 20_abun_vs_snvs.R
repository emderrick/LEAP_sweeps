library(tidyverse)
library(rstatix)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

rel_abun <- read_csv("data files/MAG_rel_abun_change.csv")

mag_snvs <- read_csv("data files/T1_subsamp_SNV_summary_MAG.csv")
mag_snvs$total_variants <- mag_snvs$SNSs_Mbp + mag_snvs$SNVs_Mbp
mag_snvs$N_total_variants <- mag_snvs$N_SNSs_Mbp + mag_snvs$N_SNVs_Mbp

mag_snv_change <- mag_snvs[, c(1,21,17:19)]
mag_snv_change <- pivot_wider(mag_snv_change, names_from = Time, values_from = total_variants)
mag_snv_change <- mag_snv_change[complete.cases(mag_snv_change ), ]
mag_snv_change$snv_change <- mag_snv_change$`Day 28` - mag_snv_change$`Day 0`

mag_n_snv_change <- mag_snvs[, c(1,22,17:19)]
mag_n_snv_change <- pivot_wider(mag_n_snv_change, names_from = Time, values_from = N_total_variants)
mag_n_snv_change <- mag_n_snv_change[complete.cases(mag_n_snv_change ), ]
mag_n_snv_change$N_snv_change <- mag_n_snv_change$`Day 28` - mag_n_snv_change$`Day 0`

abun_snv_all <- left_join(rel_abun[, c(1:3,6:12)], mag_snv_change[,c(1,3,6)])
abun_snv_all <- left_join(abun_snv_all, mag_n_snv_change[, c(1,3,6)])
abun_snv_all <- abun_snv_all[complete.cases(abun_snv_all), ]
abun_snv_all$arg <- with(abun_snv_all, ifelse(total_hits >= 1, "Yes", "No"))

snv_abun <- glm(N_snv_change ~ abun_change, data = abun_snv_all)
summary(snv_abun)

abun_snv_GBH <- subset(abun_snv_all, Treatment == "GBH")
abun_snv_CTRL <- subset(abun_snv_all, Treatment == "Control")

gbh_snv_abun <- glm(N_snv_change ~ abun_change, data = abun_snv_GBH)
summary(gbh_snv_abun)

ctrl_snv_abun <- glm(N_snv_change ~ abun_change, data = abun_snv_CTRL)
summary(ctrl_snv_abun)

ggplot(abun_snv_all, aes(x = abun_change, y = N_snv_change, colour = EPSPS_allele))+
  geom_point()+
  geom_smooth(aes(group = Treatment), method = "glm", colour = "red", size = 0.5)+
  scale_colour_manual(values = c("darkmagenta", "darkgreen", "grey80"))+
  labs(x = "MAG abundance change", y = "SNV abundance change")+
  theme_classic()+
  theme(axis.text = element_text(size = 8, colour = "black"), legend.position = "none")+
  facet_wrap(~Treatment, scales = "free")



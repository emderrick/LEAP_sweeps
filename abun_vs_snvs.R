library(tidyverse)
library(ggplot2)
library(rstatix)

rel_abun <- read_csv("refined data files/MAG_rel_abun_change.csv")

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
mag_n_sns_change$N_sns_change <- mag_n_sns_change$`Day 28` - mag_n_snv_change$`Day 0`

abun_snv <- left_join(rel_abun[, c(1:3,6:12)], mag_snv_change[,c(1,3,6)])
abun_snv <- left_join(abun_snv, mag_sns_change[, c(1,3,6)])
abun_snv <- left_join(abun_snv, mag_n_snv_change[, c(1,3,6)])
abun_snv <- left_join(abun_snv, mag_n_sns_change[, c(1,3,6)])
abun_snv <- abun_snv[complete.cases(abun_snv), ]

ggplot(abun_snv, aes(x = abun_change, y = snv_change, colour = EPSPS_allele))+
  geom_point()+
  facet_wrap(~Treatment)

ggplot(abun_snv, aes(x = abun_change, y = N_snv_change, colour = EPSPS_allele))+
  geom_point()+
  facet_wrap(~Treatment)

ggplot(abun_snv, aes(x = abun_change, y = sns_change, colour = EPSPS_allele))+
  geom_point()+
  facet_wrap(~Treatment)

ggplot(abun_snv, aes(x = abun_change, y = N_sns_change, colour = EPSPS_allele))+
  geom_point()+
  facet_wrap(~Treatment)

test <- glm(N_snv_change ~ total_hits, data = abun_snv)
summary(test)

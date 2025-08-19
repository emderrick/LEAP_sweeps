library(tidyverse)
library(ggplot2)

rel_abun <- read_csv("refined data files/MAG_avg_abun.csv")

ggplot(rel_abun, aes(x = `antibiotic efflux`, y = log10(Day_28_avg)))+
  geom_boxplot(aes(group = `antibiotic efflux`), alpha = 0)+
  geom_jitter(aes(colour = EPSPS_allele), width = 0.2)+
  theme_bw()+
  facet_wrap(~Treatment)

ggplot(rel_abun, aes(x = `antibiotic target alteration`, y = log10(Day_28_avg)))+
  geom_boxplot(aes(group = `antibiotic target alteration`), alpha = 0)+
  geom_jitter(aes(colour = EPSPS_allele), width = 0.2)+
  theme_bw()+
  facet_wrap(~Treatment)

ggplot(rel_abun, aes(x = `antibiotic inactivation`, y = log10(Day_28_avg)))+
  geom_boxplot(aes(group = `antibiotic inactivation`), alpha = 0)+
  geom_jitter(aes(colour = EPSPS_allele), width = 0.2)+
  theme_bw()+
  facet_wrap(~Treatment)

ggplot(rel_abun, aes(x = total_hits, y = log10(Day_28_avg)))+
  geom_boxplot(aes(group = total_hits), alpha = 0)+
  geom_jitter(aes(colour = EPSPS_allele), width = 0.2)+
  theme_bw()+
  facet_wrap(~Treatment)

ggplot(rel_abun, aes(x = EPSPS_allele, y = log10(Day_28_avg)))+
  geom_boxplot(aes(group = EPSPS_allele), alpha = 0)+
  geom_jitter(aes(colour = EPSPS_allele), width = 0.2)+
  theme_bw()+
  facet_wrap(~Treatment)

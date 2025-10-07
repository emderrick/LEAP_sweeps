library(tidyverse)
library(lme4)
#library(car)

#options(scipen=999)
setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

mag_SNVs <- read_csv("data files/T1_subsamp_SNV_summary_MAG.csv")
mag_SNVs$Treatment = as.factor(mag_SNVs$Treatment)

mag_SNVs$SNSs_Kbp <- mag_SNVs$SNSs_Mbp / 1000
mag_SNVs$SNVs_Kbp <- mag_SNVs$SNVs_Mbp / 1000

mag_SNVs$SNSs_Kbp <- round(mag_SNVs$SNSs_Kbp)
mag_SNVs$SNVs_Kbp <- round(mag_SNVs$SNVs_Kbp) 

mean_sns_kbp <- mean(mag_SNVs$SNSs_Kbp)
variance_sns_kbp <- var(mag_SNVs$SNSs_Kbp)

mean_snv_kbp <- mean(mag_SNVs$SNVs_Kbp)
variance_snv_kbp <- var(mag_SNVs$SNVs_Kbp)

# 
# ggplot(mag_SNVs, aes(x = SNVs_Kbp, y = SNVs_Mbp))+
#   geom_point()
# 
# ggplot(mag_SNVs, aes(x = SNSs_Kbp, y = SNSs_Mbp))+
#   geom_point()
# 
# ggplot(mag_SNVs, aes(x = SNVs_Mbp, y = SNSs_Mbp))+
#   geom_point()
# 
# ggplot(mag_SNVs, aes(x = SNVs_Kbp, y = SNSs_Kbp))+
#   geom_point()

snv_glmm <- glmer(SNVs_Kbp ~  Treatment * Time + (1|Name), family = poisson, data = mag_SNVs)
summary(snv_glmm)

sns_glmm <- glmer(SNSs_Kbp ~ Treatment * Time + (1|Name), family = poisson, data = mag_SNVs)
summary(sns_glmm)






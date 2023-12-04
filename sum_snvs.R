library(tidyverse)
library(dplyr)
library(ggplot2)
library(cowplot)

EPSPS_class_1 <- list("I4_MAG_00006", "L7_MAG_00028", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")
EPSPS_unclass <- list("L2_MAG_00052", "L4_MAG_00099", "L7_MAG_00020")

all_MAG_snvs <- read_csv("all_MAG_SNVs_med_Aug15.csv")
all_MAG_snvs <- subset(all_MAG_snvs, str_detect(new_name, "T2"))
all_SNV_sum <- all_MAG_snvs %>% group_by(mag, mag_length, new_name) %>% summarize(SNVs = sum(class %in% "SNV"), SNSs = sum(class %in% "SNS"))
all_SNV_sum$SNV_Mbp <- (all_SNV_sum$SNVs / all_SNV_sum$mag_length) *10^6
all_SNV_sum$SNS_Mbp <- (all_SNV_sum$SNSs / all_SNV_sum$mag_length) *10^6
write.csv(all_SNV_sum, "all_SNV_sum.csv", row.names = F)

all_SNV_sum <- subset(all_SNV_sum, select = -c(SNVs, SNSs))
all_SNV_wide <- pivot_wider(all_SNV_sum, names_from = "new_name", values_from = c(SNV_Mbp, SNS_Mbp))
all_SNV_wide <- all_SNV_wide %>% rowwise %>% mutate(control_SNV_mean = mean(c_across(contains("SNV_Mbp_Control")), na.rm = T))
all_SNV_wide <- all_SNV_wide %>% rowwise %>% mutate(GBH_SNV_mean = mean(c_across(contains("SNV_Mbp_GBH")), na.rm = T))
all_SNV_wide <- all_SNV_wide %>% rowwise %>% mutate(control_SNS_mean = mean(c_across(contains("SNS_Mbp_Control")), na.rm = T))
all_SNV_wide <- all_SNV_wide %>% rowwise %>% mutate(GBH_SNS_mean = mean(c_across(contains("SNS_Mbp_GBH")), na.rm = T))
write.csv(all_SNV_wide, "SNV_MAG_sums_table.csv", row.names = F)
small_SNV_wide <- select(all_SNV_wide, c('mag', 'mag_length', 'control_SNV_mean', 'GBH_SNV_mean', 'control_SNS_mean', 'GBH_SNS_mean'))
small_SNV_wide$EPSPS_class <- with(small_SNV_wide, ifelse(mag %in% EPSPS_class_1, "Class I", "Class II"))
small_SNV_wide$EPSPS_class <- with(small_SNV_wide, ifelse(mag %in% EPSPS_unclass, "Unclassified", EPSPS_class))
small_SNV_wide$EPSPS_class_fisher <- with(small_SNV_wide, ifelse(EPSPS_class == "Class II", "Class II", "Class I"))
small_SNV_wide$SNV_to_control <- with(small_SNV_wide, ifelse(control_SNV_mean > GBH_SNV_mean, "Decrease", "Increase"))
small_SNV_wide$SNS_to_control <- with(small_SNV_wide, ifelse(control_SNS_mean > GBH_SNS_mean, "Decrease", "Increase"))
small_SNV_wide$Sweep <- with(small_SNV_wide, ifelse(SNV_to_control == "Decrease" & SNS_to_control == "Increase", "Yes", "No"))
write.csv(small_SNV_wide, "SNV_summary_table.csv", row.names = F)

fishers_snv <- small_SNV_wide %>% group_by(EPSPS_class_fisher) %>% summarise(Increase = sum(SNV_to_control == "Increase"), Decrease = sum(SNV_to_control == "Decrease"))
fisher.test(fishers_snv)

fishers_sns <- small_SNV_wide %>% group_by(EPSPS_class_fisher) %>% summarise(Increase = sum(SNS_to_control == "Increase"), Decrease = sum(SNS_to_control == "Decrease"))
fishers_sns_test <-fisher.test(fishers_sns[2:3])

fishers_sweep <- small_SNV_wide %>% group_by(EPSPS_class_fisher) %>% summarise(Yes = sum(Sweep == "Yes"), No = sum(Sweep == "No"))
fishers_sweep_test <- fisher.test(fishers_sweep[2:3])


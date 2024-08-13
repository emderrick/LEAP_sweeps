library(tidyverse)
library(dplyr)

setwd("/Users/Emma/Documents/manuscript version/")

EPSPS_class_1 <- list("I4_MAG_00006", "L7_MAG_00028", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")
EPSPS_unclass <- list("L2_MAG_00052", "L4_MAG_00099")

all_MAG_snvs <- read_csv("all_MAG_SNVs_med_Aug12.csv")
all_MAG_snvs$SNS_count <- with(all_MAG_snvs, ifelse(final_ref_freq == 0, 1, 0))
all_MAG_snvs$SNS_count <- with(all_MAG_snvs, ifelse(is.na(SNS_count), 0, SNS_count))
all_MAG_snvs$SNV_count <- with(all_MAG_snvs, ifelse(final_ref_freq > 0 & final_ref_freq < 1, 1, 0))
all_MAG_snvs$SNV_count <- with(all_MAG_snvs, ifelse(is.na(SNV_count), 0, SNV_count))

all_SNV_sum <- all_MAG_snvs %>% group_by(mag, mag_length, new_name) %>% summarize(SNS = sum(SNS_count), SNV = sum(SNV_count))
all_SNV_sum$SNV_Mbp <- (all_SNV_sum$SNV / all_SNV_sum$mag_length) *10^6
all_SNV_sum$SNS_Mbp <- (all_SNV_sum$SNS / all_SNV_sum$mag_length) *10^6
write.csv(all_SNV_sum, "all_SNV_sum_subsamp.csv", row.names = F)

all_SNV_sum <- subset(all_SNV_sum, select = -c(SNV, SNS))
all_SNV_wide <- pivot_wider(all_SNV_sum, names_from = "new_name", values_from = c(SNV_Mbp, SNS_Mbp))
all_SNV_wide$control_SNV_mean <- rowMeans(all_SNV_wide[,c('SNV_Mbp_Control A at T2', 'SNV_Mbp_Control B at T2', 'SNV_Mbp_Control C at T2', 'SNV_Mbp_Control D at T2', 'SNV_Mbp_Control E at T2')], na.rm = T)
all_SNV_wide$GBH_SNV_mean <- rowMeans(all_SNV_wide[,c('SNV_Mbp_GBH A at T2', 'SNV_Mbp_GBH B at T2', 'SNV_Mbp_GBH C at T2', 'SNV_Mbp_GBH D at T2')], na.rm = T)
all_SNV_wide$control_SNS_mean <- rowMeans(all_SNV_wide[,c('SNS_Mbp_Control A at T2', 'SNS_Mbp_Control B at T2', 'SNS_Mbp_Control C at T2', 'SNS_Mbp_Control D at T2', 'SNS_Mbp_Control E at T2')], na.rm = T)
all_SNV_wide$GBH_SNS_mean <- rowMeans(all_SNV_wide[,c('SNS_Mbp_GBH A at T2', 'SNS_Mbp_GBH B at T2', 'SNS_Mbp_GBH C at T2', 'SNS_Mbp_GBH D at T2')], na.rm = T)
write.csv(all_SNV_wide, "SNV_MAG_sums_table_subsamp.csv", row.names = F)

small_SNV_wide <- select(all_SNV_wide, c('mag', 'mag_length', 'control_SNV_mean', 'GBH_SNV_mean', 'control_SNS_mean', 'GBH_SNS_mean'))
small_SNV_wide$EPSPS_class <- with(small_SNV_wide, ifelse(mag %in% EPSPS_class_1, "Class I", "Class II"))
small_SNV_wide$EPSPS_class <- with(small_SNV_wide, ifelse(mag %in% EPSPS_unclass, "Unclassified", EPSPS_class))
small_SNV_wide$EPSPS_class_fisher <- with(small_SNV_wide, ifelse(EPSPS_class == "Class II", "Class II", "Class I"))
small_SNV_wide$SNV_to_control <- with(small_SNV_wide, ifelse(control_SNV_mean > GBH_SNV_mean, "Decrease", "Increase"))
small_SNV_wide$SNS_to_control <- with(small_SNV_wide, ifelse(control_SNS_mean > GBH_SNS_mean, "Decrease", "Increase"))
write.csv(small_SNV_wide, "SNV_summary_table_subsamp.csv", row.names = F)

mag_cov <- read_csv("mag_coverage_subsamp.csv")
mag_cov_snv_sum <- left_join(mag_cov, all_SNV_sum)
mag_cov_snv_sum <- mag_cov_snv_sum %>% mutate(mag_name = case_when(mag == "I4_MAG_00006" ~ "Burkholderiaceae 1", mag == "I4_MAG_00065" ~ "Roseomonas_A", mag == "L2_MAG_00052" ~ "Erythrobacter",
                                                                     mag == "L3_MAG_00058" ~ "Prosthecobacter", mag == "L4_MAG_00099" ~ "Bosea sp001713455", mag == "L7_MAG_00020" ~ "Sphingorhabdus_B",
                                                                     mag == "L7_MAG_00028" ~ "Burkholderiaceae 2", mag == "L7_MAG_00043" ~ "Luteolibacter", mag == "L8_MAG_00011" ~ "Verrucomicrobiae",
                                                                     mag == "L8_MAG_00019" ~ "Flavobacteriales 1", mag == "L8_MAG_00042" ~ "Flavobacteriales 2"))
write.csv(mag_cov_snv_sum, "mag_coverage_snv_sum_subsamp.csv", row.names = F)

fishers_snv <- small_SNV_wide %>% group_by(EPSPS_class_fisher) %>% summarise(Increase = sum(SNV_to_control == "Increase"), Decrease = sum(SNV_to_control == "Decrease"))
fishers_snv_test <- fisher.test(fishers_snv[2:3])

fishers_sns <- small_SNV_wide %>% group_by(EPSPS_class_fisher) %>% summarise(Increase = sum(SNS_to_control == "Increase"), Decrease = sum(SNS_to_control == "Decrease"))
fishers_sns_test <-fisher.test(fishers_sns[2:3])

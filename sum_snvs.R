library(tidyverse)
library(dplyr)

#read in the SNV file with means
all_snv <- read_csv("all_MAG_SNVs_med_Aug15.csv")

all_SNV_sum <- all_snv %>% group_by(mag, mag_length, new_name) %>% summarize(SNVs = sum(class %in% "SNV"), SNSs = sum(class %in% "SNS"))
write.csv(all_SNV_sum, "all_SNV_sum.csv", row.names = F)
all_SNV_wide <- pivot_wider(all_SNV_sum, names_from = "new_name", values_from = c(SNVs,SNSs))
all_SNV_wide <- all_SNV_wide %>% select(-contains('T1'))
all_SNV_wide <- all_SNV_wide %>% rowwise %>% mutate(control_SNV_mean = mean(c_across(contains("SNVs_Control")), na.rm = T))
all_SNV_wide <- all_SNV_wide %>% rowwise %>% mutate(GBH_SNV_mean = mean(c_across(contains("SNVs_GBH")), na.rm = T))
all_SNV_wide <- all_SNV_wide %>% rowwise %>% mutate(control_SNS_mean = mean(c_across(contains("SNSs_Control")), na.rm = T))
all_SNV_wide <- all_SNV_wide %>% rowwise %>% mutate(GBH_SNS_mean = mean(c_across(contains("SNSs_GBH")), na.rm = T))
write.csv(all_SNV_wide, "SNV_MAG_sums_table.csv", row.names = F)
small_SNV_wide <- select(all_SNV_wide, c('mag','control_SNV_mean', 'GBH_SNV_mean', 'control_SNS_mean', 'GBH_SNS_mean'))
write.csv(small_SNV_wide, "SNV_summary_table.csv", row.names = F)

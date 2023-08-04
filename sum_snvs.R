library(tidyverse)
library(dplyr)

#read in the SNV file with means
all_snv <- read_csv("all_MAG_SNVs_med_July25.csv")

all_SNV_sum <- all_snv %>% group_by(mag, mag_length, name) %>% summarize(SNVs=sum(class %in% "SNV"), SNSs=sum(class %in% "SNS"))
write.csv(all_SNV_sum, "all_SNV_sum.csv", row.names = F)
all_SNV_wide<- pivot_wider(all_SNV_sum, names_from="name", values_from=c(SNVs,SNSs))
all_SNV_wide <- all_SNV_wide[c(1:11, 16:24)]
all_SNV_wide<- all_SNV_wide %>% rowwise %>% mutate(control_SNV_mean = mean(c_across(contains("SNVs_Control")), na.rm = T))
all_SNV_wide$control_SNV_mean <- rowMeans(all_SNV_wide(c_across(contains("SNVs_Control")), na.rm=T))
all_SNV_wide<- all_SNV_wide %>% rowwise %>% mutate(GBH_SNV_mean = mean(c_across(contains("SNVs_GBH")), na.rm = T))
all_SNV_wide<- all_SNV_wide %>% rowwise %>% mutate(control_SNS_mean = mean(c_across(contains("SNSs_Control")), na.rm = T))
all_SNV_wide<- all_SNV_wide %>% rowwise %>% mutate(GBH_SNS_mean = mean(c_across(contains("SNSs_GBH")), na.rm = T))
small_SNV_wide<- select(all_SNV_wide, c('mag','control_SNV_mean', 'GBH_SNV_mean', 'control_SNS_mean', 'GBH_SNS_mean'))
write.csv(small_SNV_wide, "SNV_summary_table.csv", row.names=F)

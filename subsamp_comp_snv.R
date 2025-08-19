library(tidyverse)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

subsamp_snv <- read_csv("subsamp data files/T1_subsamp_SNV_summary_MAG.csv")
subsamp_snv <- subsamp_snv[, c(1,2,13:16)]
colnames(subsamp_snv) <- c("mag", "Name_Time", "sub_SNS_Mbp", "sub_SNV_Mbp", "sub_N_SNS_Mbp", "sub_N_SNV_Mbp")
all_snv <- read_csv("refined data files/T1_refined_SNV_summary_MAG.csv")
all_snv <- all_snv[, c(1,2,13:19)]

comp_snv <- left_join(subsamp_snv, all_snv)

ggplot(comp_snv, aes(x = sub_SNV_Mbp, y = SNVs_Mbp))+
  geom_point()+
  xlim(0,70000)

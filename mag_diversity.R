library(tidyverse)
library(dplyr)
library(ggplot2)
library(rstatix)

options(scipen=999)
setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

mag_SNVs <- read_csv("subsamp data files/T1_subsamp_SNV_summary_MAG.csv")
mag_SNVs <- subset(mag_SNVs, med_cov >= 4)

mag_SNVs <- mag_SNVs %>% mutate(across(c(7,8,11:16), ~ replace_na(., 0)))

ggplot(mag_SNVs, aes(x = mag_coverage, y = SNVs_Mbp))+
  geom_point()
       
mag_SNV_summary <- mag_SNVs %>% group_by(Treatment, Time) %>%
  get_summary_stats(SNVs_Mbp, type = "common")

mag_SNS_summary <- mag_SNVs %>% group_by(Treatment, Time) %>%
  get_summary_stats(SNSs_Mbp, type = "common")

mag_N_SNV_summary <- mag_SNVs %>% group_by(Treatment, Time) %>%
  get_summary_stats(N_SNVs_Mbp, type = "common")

mag_N_SNS_summary <- mag_SNVs %>% group_by(Treatment, Time) %>%
  get_summary_stats(N_SNSs_Mbp, type = "common")

mag_NS_SNV_summary <- mag_SNVs  %>% group_by(Treatment, Time) %>%
  get_summary_stats(NS_ratio_SNV, type = "common")

mag_NS_SNS_summary <- mag_SNVs  %>% group_by(Treatment, Time) %>%
  get_summary_stats(NS_ratio_SNS, type = "common")

snv_5x <- dunn_test(mag_SNVs, SNVs_Mbp ~ Treatment_Time, p.adjust.method = "BH", detailed = FALSE)
sns_5x <- dunn_test(mag_SNVs, SNSs_Mbp ~ Treatment_Time, p.adjust.method = "BH", detailed = FALSE)
n_snv_5x <- dunn_test(mag_SNVs, N_SNVs_Mbp ~ Treatment_Time, p.adjust.method = "BH", detailed = FALSE)
n_sns_5x <- dunn_test(mag_SNVs, N_SNSs_Mbp ~ Treatment_Time, p.adjust.method = "BH", detailed = FALSE)
NS_SNV_5x <- dunn_test(mag_SNVs, NS_ratio_SNV ~ Treatment_Time, p.adjust.method = "BH", detailed = FALSE)
NS_SNS_5x <- dunn_test(mag_SNVs, NS_ratio_SNS ~ Treatment_Time, p.adjust.method = "BH", detailed = FALSE)

treatment_vs_SNV_boxplot <- ggplot(mag_SNVs, aes(x = Time, y = SNVs_Mbp, fill = Treatment, colour = Treatment))+
  geom_boxplot(outliers = F)+
  geom_point(aes(colour = Treatment), size = 0.25, position = position_jitterdodge(jitter.width = 0.25))+
  scale_fill_manual(values = c("white", "white"))+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  labs(x = "", y = "Polymorphic sites per Mbp")+  
  theme_bw()+
  theme(text = element_text(size = 12), axis.title.y = element_text(m= margin(10, 10, 10, 10)), axis.text = element_text(colour = "black"))
ggsave("subsamp figures/treatment_vs_SNVs_boxplot.pdf", treatment_vs_SNV_boxplot, limitsize = F, width = 5, height = 4)

treatment_vs_SNS_boxplot <- ggplot(mag_SNVs, aes(x = Time, y = SNSs_Mbp, fill = Treatment, colour = Treatment))+
  geom_boxplot(outliers = F)+
  geom_point(aes(colour = Treatment), size = 0.25, position = position_jitterdodge(jitter.width = 0.25))+
  scale_fill_manual(values = c("white", "white"))+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  labs(x = "", y = "Fixed substitutions per Mbp")+  
  theme_bw()+
  theme(text = element_text(size = 12), axis.title.y = element_text(m= margin(10, 10, 10, 10)), axis.text = element_text(colour = "black"))
ggsave("subsamp figures/treatment_vs_SNSs_boxplot.pdf", treatment_vs_SNS_boxplot, limitsize = F, width = 5, height = 4)

treatment_vs_N_SNV_boxplot <- ggplot(mag_SNVs, aes(x = Time, y = N_SNVs_Mbp, fill = Treatment, colour = Treatment))+
  geom_boxplot(outliers = F)+
  geom_point(aes(colour = Treatment), size = 0.25, position = position_jitterdodge(jitter.width = 0.25))+
  scale_fill_manual(values = c("white", "white"))+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  labs(x = "", y = "nonsyn polymorphic sites per Mbp")+  
  theme_bw()+
  theme(text = element_text(size = 12), axis.title.y = element_text(m= margin(10, 10, 10, 10)), axis.text = element_text(colour = "black"))
ggsave("subsamp figures/treatment_vs_N_SNVs_boxplot.pdf", treatment_vs_N_SNV_boxplot, limitsize = F, width = 5, height = 4)

treatment_vs_N_SNS_boxplot <- ggplot(mag_SNVs, aes(x = Time, y = N_SNSs_Mbp, fill = Treatment, colour = Treatment))+
  geom_boxplot(outliers = F)+
  geom_point(aes(colour = Treatment), size = 0.25, position = position_jitterdodge(jitter.width = 0.25))+
  scale_fill_manual(values = c("white", "white"))+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  labs(x = "", y = "nonsyn fixed substitutions per Mbp")+  
  theme_bw()+
  theme(text = element_text(size = 12), axis.title.y = element_text(m= margin(10, 10, 10, 10)), axis.text = element_text(colour = "black"))
ggsave("subsamp figures/treatment_vs_N_SNSs_boxplot.pdf", treatment_vs_N_SNS_boxplot, limitsize = F, width = 5, height = 4)

treatment_vs_NS_boxplot <- ggplot(mag_SNVs, aes(x = Time, y = NS_ratio_SNV, fill = Treatment, colour = Treatment))+
  geom_boxplot(outliers = F)+
  geom_point(aes(colour = Treatment), size = 0.25, position = position_jitterdodge(jitter.width = 0.25))+
  scale_fill_manual(values = c("white", "white"))+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  labs(x = "", y = "N:S ratio")+  
  theme_bw()+
  theme(text = element_text(size = 12), axis.title.y = element_text(m= margin(10, 10, 10, 10)), axis.text = element_text(colour = "black"))
ggsave("subsamp figures/treatment_vs_NS_SNV_boxplot.pdf", treatment_vs_NS_boxplot, limitsize = F, width = 4.75, height = 4)

treatment_vs_NS_SNS_boxplot <- ggplot(mag_SNVs, aes(x = Time, y = NS_ratio_SNS, fill = Treatment, colour = Treatment))+
  geom_boxplot(outliers = F)+
  geom_point(aes(colour = Treatment), size = 0.25, position = position_jitterdodge(jitter.width = 0.25))+
  scale_fill_manual(values = c("white", "white"))+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  labs(x = "", y = "N:S ratio")+  
  theme_bw()+
  theme(text = element_text(size = 12), axis.title.y = element_text(m= margin(10, 10, 10, 10)), axis.text = element_text(colour = "black"))
ggsave("subsamp figures/treatment_vs_NS_SNS_boxplot.pdf", treatment_vs_NS_SNS_boxplot, limitsize = F, width = 4.75, height = 4)

pond_vs_SNV_boxplot <- ggplot(mag_SNVs, aes(x = Name, y = SNVs_Mbp, fill = Time, colour = Time))+
  geom_boxplot(outliers = F)+
  geom_point(aes(colour = Time), size = 0.25, position = position_jitterdodge(jitter.width = 0.25))+
  scale_fill_manual(values = c("white", "white"))+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  labs(x = "", y = "Polymorphic sites per Mbp")+  
  theme_bw()+
  theme(text = element_text(size = 12), axis.title.y = element_text(m= margin(10, 10, 10, 10)), axis.text = element_text(colour = "black"))
ggsave("subsamp figures/pond_vs_SNV_boxplot.pdf", pond_vs_SNV_boxplot, limitsize = F, width = 10, height = 4)

pond_vs_SNS_boxplot <- ggplot(mag_SNVs, aes(x = Name, y = SNSs_Mbp, fill = Time, colour = Time))+
  geom_boxplot(outliers = F)+
  geom_point(aes(colour = Time), size = 0.25, position = position_jitterdodge(jitter.width = 0.25))+
  scale_fill_manual(values = c("white", "white"))+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  labs(x = "", y = "Fixed substiturions per Mbp")+  
  theme_bw()+
  theme(text = element_text(size = 12), axis.title.y = element_text(m= margin(10, 10, 10, 10)), axis.text = element_text(colour = "black"))
ggsave("subsamp figures/pond_vs_SNS_boxplot.pdf", pond_vs_SNS_boxplot, limitsize = F, width = 10, height = 4)

pond_vs_N_SNV_boxplot <- ggplot(mag_SNVs, aes(x = Name, y = N_SNVs_Mbp, fill = Time, colour = Time))+
  geom_boxplot(outliers = F)+
  geom_point(aes(colour = Time), size = 0.25, position = position_jitterdodge(jitter.width = 0.25))+
  scale_fill_manual(values = c("white", "white"))+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  labs(x = "", y = "nonsyn polymorphic sites per Mbp")+  
  theme_bw()+
  theme(text = element_text(size = 12), axis.title.y = element_text(m= margin(10, 10, 10, 10)), axis.text = element_text(colour = "black"))
ggsave("subsamp figures/pond_vs_N_SNV_boxplot.pdf", pond_vs_N_SNV_boxplot, limitsize = F, width = 10, height = 4)

pond_vs_N_SNS_boxplot <- ggplot(mag_SNVs, aes(x = Name, y = N_SNSs_Mbp, fill = Time, colour = Time))+
  geom_boxplot(outliers = F)+
  geom_point(aes(colour = Time), size = 0.25, position = position_jitterdodge(jitter.width = 0.25))+
  scale_fill_manual(values = c("white", "white"))+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  labs(x = "", y = "nonsyn fixed substiturions per Mbp")+  
  theme_bw()+
  theme(text = element_text(size = 12), axis.title.y = element_text(m= margin(10, 10, 10, 10)), axis.text = element_text(colour = "black"))
ggsave("subsamp figures/pond_vs_N_SNS_boxplot.pdf", pond_vs_N_SNS_boxplot, limitsize = F, width = 10, height = 4)

# treatment_vs_coverage <- ggplot(mag_SNVs, aes(x = Treatment_Time, y = med_cov, colour = Treatment))+
#   geom_point(size = 0.5, position = position_jitterdodge(jitter.width = 0.25))+
#   theme_bw()+
#   scale_colour_manual(values = c("darkgreen", "darkmagenta"))
# ggsave("refined figures/treatment_vs_coverage_T1_5x.pdf", treatment_vs_coverage, limitsize = F)
# 
# coverage_vs_diversity_point <- ggplot(mag_SNVs, aes(x = med_cov, y = SNVs_Mbp, colour = Treatment))+
#   geom_point(size = 1)+
#   theme_bw()+
#   labs(x = "MAG Coverage", y = "Polymorphic sites per Mbp")+  
#   scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
#   theme(text = element_text(size = 12), axis.title.y = element_text(m= margin(10, 10, 10, 10)), axis.text = element_text(colour = "black"))
# ggsave("refined figures/coverage_vs_diversity_T1_5x.pdf", coverage_vs_diversity_point, limitsize = F, width = 5, height = 3)
# 
# coverage_vs_SNS_point <- ggplot(mag_SNVs, aes(x = mag_coverage, y = SNSs_Mbp, colour = Treatment))+
#   geom_point(size = 1)+
#   theme_bw()+
#   labs(x = "MAG Coverage", y = "SNSs per Mbp")+  
#   scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
#   theme(text = element_text(size = 12), axis.title.y = element_text(m= margin(10, 10, 10, 10)), axis.text = element_text(colour = "black"))
# ggsave("refined figures/coverage_vs_SNS_T1_5x.pdf", coverage_vs_SNS_point, limitsize = F, width = 4, height = 4)


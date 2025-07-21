library(tidyverse)
library(ggplot2)
library(rstatix)

options(scipen=999)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

mag_SNVs <- read_csv("data files/T1_SNV_summary_MAG.csv")

mag_SNVs_10x <- subset(mag_SNVs, mag_coverage >= 10 & mag_breadth >= 0.5)
mag_diversity_summary_10x <- mag_SNVs_10x %>% group_by(Treatment, Time) %>%
  get_summary_stats(SNVs_Mbp, type = "common")

mag_coverage_summary_10x <- mag_SNVs_10x %>% group_by(Treatment, Time) %>%
  get_summary_stats(mag_coverage, type = "common")

kruskal.test(mag_SNVs_10x$SNVs_Mbp ~ mag_SNVs_10x$Treatment_Time)
snv_10x <- dunn_test(mag_SNVs_10x, SNVs_Mbp ~ Treatment_Time, p.adjust.method = "BH", detailed = FALSE)

treatment_vs_diversity_boxplot <- ggplot(mag_SNVs_10x, aes(x = Time, y = SNVs_Mbp, fill = Treatment, colour = Treatment))+
  geom_boxplot(outliers = F)+
  geom_point(aes(colour = Treatment), size = 0.25, position = position_jitterdodge(jitter.width = 0.25))+
  scale_fill_manual(values = c("white", "white"))+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  labs(x = "", y = "Polymorphic sites per Mbp")+  
  theme_bw()+
  theme(text = element_text(size = 12), axis.title.y = element_text(m= margin(10, 10, 10, 10)), axis.text = element_text(colour = "black"))
ggsave("figures/treatment_vs_diversity_boxplot_T1_10x.pdf", treatment_vs_diversity_boxplot, limitsize = F, width = 5, height = 4)

treatment_vs_coverage <- ggplot(subset(mag_SNVs_10x, mag_coverage < 1000), aes(x = Treatment_Time, y = mag_coverage, colour = Treatment))+
  geom_point(size = 0.5, position = position_jitterdodge(jitter.width = 0.25))+
  theme_bw()+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))
ggsave("figures/treatment_vs_coverage_T1_10x.pdf", treatment_vs_coverage, limitsize = F)

coverage_vs_diversity_point <- ggplot(subset(mag_SNVs_10x, mag_coverage < 400), aes(x = mag_coverage, y = SNVs_Mbp, colour = Treatment))+
  geom_point(size = 1)+
  theme_bw()+
  labs(x = "MAG Coverage", y = "Polymorphic sites per Mbp")+  
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  theme(text = element_text(size = 12), axis.title.y = element_text(m= margin(10, 10, 10, 10)), axis.text = element_text(colour = "black"))
ggsave("figures/coverage_vs_diversity_T1_10x.pdf", coverage_vs_diversity_point, limitsize = F, width = 5, height = 3)

## SNSs

T1_MAGs <- read_csv("data files/T1_mag_info.csv")
T1_MAGs$Name_Time <- paste(T1_MAGs$Name, T1_MAGs$time, sep = " ")
T1_MAGs <- subset(T1_MAGs, mag_coverage >= 10 & mag_breadth >= 0.5)
T1_MAGs <- T1_MAGs[, c("mag", "Name_Time", "SNV_count", "SNS_count")]

mag_SNSs <- read_csv("data files/T1_SNS_summary_MAG.csv")
mag_SNSs_10x <- subset(mag_SNSs, mag_coverage >= 10 & mag_breadth >= 0.5)

mag_SNSs_10x  <- full_join(mag_SNSs_10x, T1_MAGs)
mag_SNSs_10x$SNSs_Mbp <- with(mag_SNSs_10x, ifelse(is.na(n), 0, SNSs_Mbp))
mag_SNSs_10x <- mag_SNSs_10x %>% group_by(Name_Time) %>% fill(Name, Time, Treatment, Treatment_Time, .direction = "updown")
mag_SNSs_10x <- mag_SNSs_10x %>% ungroup()

mag_SNS_summary_10x <- mag_SNSs_10x %>% group_by(Treatment, Time) %>%
  get_summary_stats(SNSs_Mbp, type = "common")

kruskal.test(mag_SNSs_10x$SNSs_Mbp ~ mag_SNSs_10x$Treatment_Time)
sns_10x <- dunn_test(mag_SNSs_10x, SNSs_Mbp ~ Treatment_Time, p.adjust.method = "BH", detailed = FALSE)

treatment_vs_SNS_boxplot <- ggplot(mag_SNSs_10x, aes(x = Time, y = SNSs_Mbp, fill = Treatment, colour = Treatment))+
  geom_boxplot(outliers = F)+
  geom_point(aes(colour = Treatment), size = 0.25, position = position_jitterdodge(jitter.width = 0.25))+
  scale_fill_manual(values = c("white", "white"))+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  labs(x = "", y = "Fixed substitutions per Mbp")+  
  theme_bw()+
  theme(text = element_text(size = 12), axis.title.y = element_text(m= margin(10, 10, 10, 10)), axis.text = element_text(colour = "black"))
ggsave("figures/treatment_vs_SNSs_boxplot_T1_10x.pdf", treatment_vs_SNS_boxplot, limitsize = F, width = 5, height = 4)

coverage_vs_SNS_point <- ggplot(subset(mag_SNSs_10x, mag_coverage < 1000), aes(x = mag_coverage, y = SNSs_Mbp, colour = Treatment))+
  geom_point(size = 1)+
  theme_bw()+
  labs(x = "MAG Coverage", y = "SNSs per Mbp")+  
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  theme(text = element_text(size = 12), axis.title.y = element_text(m= margin(10, 10, 10, 10)), axis.text = element_text(colour = "black"))
ggsave("figures/coverage_vs_SNS_T1_10x.pdf", coverage_vs_SNS_point, limitsize = F, width = 4, height = 4)

#NS ratio

NS_MAG <- read_csv("data files/T1_NS_summary_MAG.csv")
NS_MAG_10x <- subset(NS_MAG, mag_coverage >= 10 & mag_breadth >= 0.5 & class == "SNV")

mag_NS_summary_10x <- NS_MAG_10x %>% group_by(Treatment, Time) %>%
  get_summary_stats(NS_ratio, type = "common")

kruskal.test(NS_MAG_10x$NS_ratio ~ NS_MAG_10x$Treatment_Time)
NS_10x <- dunn_test(NS_MAG_10x, NS_ratio ~ Treatment_Time, p.adjust.method = "BH")

treatment_vs_NS_boxplot <- ggplot(NS_MAG_10x, aes(x = Time, y = NS_ratio, fill = Treatment, colour = Treatment))+
  geom_boxplot(outliers = F)+
  geom_point(aes(colour = Treatment), size = 0.25, position = position_jitterdodge(jitter.width = 0.25))+
  scale_fill_manual(values = c("white", "white"))+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  labs(x = "", y = "N:S ratio")+  
  theme_bw()+
  theme(text = element_text(size = 12), axis.title.y = element_text(m= margin(10, 10, 10, 10)), axis.text = element_text(colour = "black"))
ggsave("figures/treatment_vs_NS_boxplot_T1_10x.pdf", treatment_vs_NS_boxplot, limitsize = F, width = 5, height = 4)


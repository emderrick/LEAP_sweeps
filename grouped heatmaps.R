library(tidyverse)
library(ggplot2)
library(viridis)
library(patchwork)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

mag_list <- c("bin.305", "bin.609", "bin.676") 

all_snv <- read_csv("data files/mag_SNV_depth_info.csv")
all_snv$time <- ifelse(grepl("1", all_snv$time), "Day 0", "Day 28")
mag_snv_avg <- read_csv("data files/all_avg_freq.csv")
all_snv <- left_join(all_snv, mag_snv_avg)
all_snv$Treatment <- ifelse(grepl("CTRL", all_snv$Name), "Control", "GBH")

all_sum <- read_csv("data files/T1_SNV_summary_MAG.csv")
all_sum <- subset(all_sum, mag %in% mag_list & mag_coverage >= 5 & mag_breadth >= 0.5)

all_SNS <- read_csv("data files/T1_SNS_summary_MAG.csv")
all_SNS <- subset(all_SNS, mag %in% mag_list & mag_coverage >= 5 & mag_breadth >= 0.5)

bin.305_snv <- subset(all_snv, mag == "bin.305")
bin.305_sum <- subset(all_sum, mag == "bin.305")
bin.305_sns <- subset(all_SNS, mag == "bin.305")

bin.305_snv_sum <- ggplot(bin.305_sum, aes(x = Name, y = SNVs_Mbp, fill = Treatment))+
  geom_col()+
  theme_bw()+
  scale_fill_manual(values = c("darkgreen", "darkmagenta"))+
  labs(x = "", y = "Polymorphic sites per Mbp", fill = "")+
  theme(text = element_text(size = 16), axis.text = element_text(colour = "black"), strip.text.x = element_text(size = 18), legend.position = "none")+
  facet_wrap(~Time, scales = "fixed")

bin.305_sns_sum <- ggplot(bin.305_sns, aes(x = Name, y = SNSs_Mbp, fill = Treatment))+
  geom_col()+
  theme_bw()+
  scale_fill_manual(values = c("darkgreen", "darkmagenta"))+
  labs(x = "", y = "Fixed subsitutions per Mbp", fill = "Treatment")+
  theme(text = element_text(size = 16), axis.text = element_text(colour = "black"), strip.text.x = element_blank())+
  facet_wrap(~Time, scales = "fixed")

bin.305_snvs <- ggplot(subset(bin.305_snv, final_ref_freq < 1 & final_ref_freq > 0), aes(x = Name, y = final_ref_freq, colour = Treatment)) +
  geom_boxplot(outliers = T,outlier.size = 0.5 )+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  ylim(0,1)+
  theme_bw()+
  labs(x = "Pond", y = "Reference frequency of SNVs")+
  theme(text = element_text(size = 16), axis.text = element_text(colour = "black"), strip.text.x = element_blank(), legend.position = "none")+
  facet_wrap(~time, scales = "fixed")

bin_305_plot <- bin.305_snv_sum /  bin.305_sns_sum / bin.305_snvs
ggsave("figures/bin.305_snvs.pdf", bin_305_plot, limitsize = F, width = 16, height = 12)

bin.609_snv <- subset(all_snv, mag == "bin.609")
bin.609_sum <- subset(all_sum, mag == "bin.609")
bin.609_sns <- subset(all_SNS, mag == "bin.609")

bin.609_snv_sum <- ggplot(bin.609_sum, aes(x = Name, y = SNVs_Mbp, fill = Treatment))+
  geom_col()+
  theme_bw()+
  scale_fill_manual(values = c("darkgreen", "darkmagenta"))+
  labs(x = "", y = "Polymorphic sites per Mbp", fill = "")+
  theme(text = element_text(size = 16), axis.text = element_text(colour = "black"), strip.text.x = element_text(size = 18), legend.position = "none")+
  facet_wrap(~Time, scales = "fixed")

bin.609_sns_sum <- ggplot(bin.609_sns, aes(x = Name, y = SNSs_Mbp, fill = Treatment))+
  geom_col()+
  theme_bw()+
  scale_fill_manual(values = c("darkgreen", "darkmagenta"))+
  labs(x = "", y = "Fixed subsitutions per Mbp", fill = "Treatment")+
  theme(text = element_text(size = 16), axis.text = element_text(colour = "black"), strip.text.x = element_blank())+
  facet_wrap(~Time, scales = "fixed")

bin.609_snvs <- ggplot(subset(bin.609_snv, final_ref_freq < 1 & final_ref_freq > 0), aes(x = Name, y = final_ref_freq, colour = Treatment)) +
  geom_boxplot(outliers = T,outlier.size = 0.5 )+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  ylim(0,1)+
  theme_bw()+
  labs(x = "Pond", y = "Reference frequency of SNVs")+
  theme(text = element_text(size = 16), axis.text = element_text(colour = "black"), strip.text.x = element_blank(), legend.position = "none")+
  facet_wrap(~time, scales = "fixed")

bin_609_plot <- bin.609_snv_sum /  bin.609_sns_sum / bin.609_snvs
ggsave("figures/bin.609_snvs.pdf", bin_609_plot, limitsize = F, width = 16, height = 13)

bin.676_snv <- subset(all_snv, mag == "bin.676")
bin.676_sum <- subset(all_sum, mag == "bin.676")
bin.676_sns <- subset(all_SNS, mag == "bin.676")

bin.676_snv_sum <- ggplot(bin.676_sum, aes(x = Name, y = SNVs_Mbp, fill = Treatment))+
  geom_col()+
  theme_bw()+
  scale_fill_manual(values = c("darkgreen", "darkmagenta"))+
  labs(x = "", y = "Polymorphic sites per Mbp", fill = "")+
  theme(text = element_text(size = 16), axis.text = element_text(colour = "black"), strip.text.x = element_text(size = 18),legend.position = "none")+
  facet_wrap(~Time, scales = "fixed")

bin.676_sns_sum <- ggplot(bin.676_sns, aes(x = Name, y = SNSs_Mbp, fill = Treatment))+
  geom_col()+
  theme_bw()+
  scale_fill_manual(values = c("darkgreen", "darkmagenta"))+
  labs(x = "", y = "Fixed subsitutions per Mbp", fill = "Treatment")+
  theme(text = element_text(size = 16), axis.text = element_text(colour = "black"), strip.text.x = element_blank())+
  facet_wrap(~Time, scales = "fixed")

bin.676_snvs <- ggplot(subset(bin.676_snv, final_ref_freq < 1 & final_ref_freq > 0), aes(x = Name, y = final_ref_freq, colour = Treatment)) +
  geom_boxplot(outliers = T,outlier.size = 0.5 )+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  ylim(0,1)+
  theme_bw()+
  labs(x = "Pond", y = "Reference frequency of SNVs")+
  theme(text = element_text(size = 16), axis.text = element_text(colour = "black"), strip.text.x = element_blank(), legend.position = "none")+
  facet_wrap(~time, scales = "fixed")

bin_676_plot <- bin.676_snv_sum /  bin.676_sns_sum / bin.676_snvs
ggsave("figures/bin.676_snvs.pdf", bin_676_plot, limitsize = F, width = 16, height = 13)

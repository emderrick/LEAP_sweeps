library(tidyverse)
library(ggpubr)

options(scipen=999)
setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

mag_SNVs <- read_csv("data files/T1_subsamp_SNV_summary_MAG.csv")
mag_SNVs$SNSs_Kbp <- mag_SNVs$SNSs_Mbp / 1000
mag_SNVs$SNVs_Kbp <- mag_SNVs$SNVs_Mbp / 1000

treatment_vs_SNV_boxplot <- ggplot(mag_SNVs, aes(x = Time, y = SNVs_Kbp, fill = Treatment, colour = Treatment))+
  geom_boxplot(outliers = F)+
  geom_point(aes(colour = Treatment), size = 0.25, position = position_jitterdodge(jitter.width = 0.25))+
  scale_fill_manual(values = c("white", "white"))+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  labs(x = "", y = "Polymorphic sites per Kbp")+  
  theme_bw()+
  theme(text = element_text(size = 12), axis.title.y = element_text(m= margin(10, 10, 10, 10)), axis.text = element_text(colour = "black"))
ggsave("figures/treatment_vs_SNVs_boxplot.pdf", treatment_vs_SNV_boxplot, units = "cm", width = 17, height = 10)

treatment_vs_SNS_boxplot <- ggplot(mag_SNVs, aes(x = Time, y = SNSs_Kbp, fill = Treatment, colour = Treatment))+
  geom_boxplot(outliers = F)+
  geom_point(aes(colour = Treatment), size = 0.25, position = position_jitterdodge(jitter.width = 0.25))+
  scale_fill_manual(values = c("white", "white"))+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  labs(x = "", y = "Fixed substitutions per Kbp")+  
  theme_bw()+
  theme(text = element_text(size = 12), axis.title.y = element_text(m= margin(10, 10, 10, 10)), axis.text = element_text(colour = "black"))
ggsave("figures/treatment_vs_SNSs_boxplot.pdf", treatment_vs_SNS_boxplot, units = "cm", width = 17, height = 10)

pond_vs_SNV_boxplot <- ggplot(mag_SNVs, aes(x = Name, y = SNVs_Kbp))+
  geom_boxplot(outliers = F, aes(colour = Treatment_Time))+
  geom_point(aes(colour = Treatment_Time), size = 0.25, position = position_jitterdodge(jitter.width = 0.25))+
  scale_fill_manual(values = c("white", "white"))+
  scale_colour_manual(values = c("#98BF64", "darkgreen", "#BE93D4", "darkmagenta"))+
  labs(x = "", y = "Polymorphic sites per Kbp", colour = "Treatment")+
  theme_classic()+
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 12, colour = "black"),
        plot.margin = unit(c(1, 1, 0, 0.25), "cm"))
ggsave("figures/pond_vs_SNV_boxplot.pdf", pond_vs_SNV_boxplot, units = "cm", width = 17, height = 10)

pond_vs_SNS_boxplot <- ggplot(mag_SNVs, aes(x = Name, y = SNSs_Kbp))+
  geom_boxplot(outliers = F, aes(colour = Treatment_Time))+
  geom_point(aes(colour = Treatment_Time), size = 0.25, position = position_jitterdodge(jitter.width = 0.25))+
  scale_fill_manual(values = c("white", "white"))+
  scale_colour_manual(values = c("#98BF64", "darkgreen", "#BE93D4", "darkmagenta"))+
  labs(x = "", y = "Fixed substiturions per Kbp", colour = "Treatment")+  
  theme_classic()+
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 12, colour = "black"),
        plot.margin = unit(c(1, 1, 0, 0.25), "cm"))
          
ggsave("figures/pond_vs_SNS_boxplot.pdf", pond_vs_SNS_boxplot, units = "cm", width = 17, height = 10)

pond_sns_snvs <- ggarrange(pond_vs_SNV_boxplot, pond_vs_SNS_boxplot, nrow = 2, common.legend = TRUE,
                           legend = "bottom", labels = c("A", "B"))

ggsave("figures/pond_vs_SNV_SNS.pdf", pond_sns_snvs, units = "cm", width = 17, height = 16)

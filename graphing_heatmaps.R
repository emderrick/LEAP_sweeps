library(tidyverse)
library(dplyr)
library(ggplot2)
library(viridis)
library(stringr)
library(ggpubr)

#read in all the SNV files
I4_MAG_00006 <- read_csv("I4_MAG_00006_SNVs.csv")
I4_MAG_00065 <- read_csv("I4_MAG_00065_SNVs.csv")
L2_MAG_00048 <- read_csv("L2_MAG_00048_SNVs.csv")
L2_MAG_00052 <- read_csv("L2_MAG_00052_SNVs.csv")
L3_MAG_00058 <- read_csv("L3_MAG_00058_SNVs.csv")
L4_MAG_00099 <- read_csv("L4_MAG_00099_SNVs.csv")
L7_MAG_00020 <- read_csv("L7_MAG_00020_SNVs.csv")
L7_MAG_00028 <- read_csv("L7_MAG_00028_SNVs.csv")
L7_MAG_00043 <- read_csv("L7_MAG_00043_SNVs.csv")
L8_MAG_00011 <- read_csv("L8_MAG_00011_SNVs.csv")
L8_MAG_00019 <- read_csv("L8_MAG_00019_SNVs.csv")
L8_MAG_00042 <- read_csv("L8_MAG_00042_SNVs.csv")

#for I4_MAG_00006
ggplot(I4_MAG_00006, aes(x = name, y = groups, fill= final_ref_freq)) +
  geom_tile()+
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(I4_MAG_00006$groups)[order(I4_MAG_00006$all_mean)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("SJAQ100 sp016735685 "),"assembled from Control B"))))+
  theme(axis.title.y = element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank(), 
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+
  guides(fill = guide_legend(reverse=TRUE))
ggsave(filename="I4_MAG_00006_SNV_95_heatmap.png", dpi = 500)

I4_MAG_00006_SNS_sum <- I4_MAG_00006 %>% group_by(name, mag_length) %>% summarize(SNVs=sum(class %in% "SNV"), SNSs=sum(class %in% "SNS"))

I4_MAG_00006_snv_frac <- ggplot(I4_MAG_00006_SNS_sum, aes(x = name, y=(SNSs/(SNSs+SNVs))))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  ylim(0, 1)+
  labs(x="Pond",  y= "fraction of SNVs dominated by a single allele")
ggsave("I4_MAG_00006_#SNV.png", limitsize = F, dpi = 500)

I4_MAG_00006_snv_tot <- ggplot(I4_MAG_00006_SNS_sum, aes(x = name, y=((SNSs+SNVs)/mag_length)*10^6))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  scale_fill_manual(values=c("#E0E0E0","#424242"))+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  ylim(0, 1)+
  labs(x="Pond",  y= "fraction of SNVs dominated by a single allele")
ggsave("I4_MAG_00006_#SNV.png", limitsize = F, dpi = 500)

ggarrange(I4_MAG_00006_heat, I4_MAG_00006_snv, ncol=1, nrow=2, align="hv", label.x="Pond",common.legend = T, legend="right")
ggsave(filename = "I4_MAG_00006.png", limitsize=F, width=16, height=18)

#for I4_MAG_00065
I4_MAG_00065_heat <- ggplot(I4_MAG_00065, aes(x = name, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(I4_MAG_00065$groups)[order(I4_MAG_00065$all_mean)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("Roseomonas sp. "),"assembled from Control B"))))+
  theme(axis.title.y = element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+
  guides(fill = guide_legend(reverse=TRUE)) 
ggsave(filename="I4_MAG_00065_SNV_95_heatmap.png", dpi = 500)

I4_MAG_00065_SNSsum <- I4_MAG_00065 %>% group_by(name) %>% summarize_at(c('number_SNSs', 'number_SNVs', 'number_divergent'), sum, na.rm=T)

I4_MAG_00065_snv <- ggplot(I4_MAG_00065_SNSsum, aes(x = name, y=(number_SNSs/number_divergent)))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  ylim(0, 1)+
  labs(x="Pond",  y= "fraction of SNVs dominated by a single allele")
ggsave("I4_MAG_00065_#SNV.png", limitsize = F, dpi = 500)

ggarrange(I4_MAG_00065_heat, I4_MAG_00065_snv, ncol=1, nrow=2, align="hv", label.x="Pond", common.legend = T, legend="right")
ggsave(filename = "I4_MAG_00065.png", limitsize=F, width=16, height=18)

#for L2_MAG_00048
L2_MAG_00048_heat <- ggplot(L2_MAG_00048, aes(x = name, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(L2_MAG_00048$groups)[order(L2_MAG_00048$all_mean)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("Flavobacterium sp. "),"assembled from GBH A"))))+
  theme(axis.title.y = element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+
  guides(fill = guide_legend(reverse=TRUE)) 
ggsave(filename="L2_MAG_00048_SNV_95_heatmap.png", dpi = 500)

L2_MAG_00048_SNSsum <- L2_MAG_00048 %>% group_by(name) %>% summarize_at(c('number_SNSs', 'number_SNVs', 'number_divergent'), sum, na.rm=T)

L2_MAG_00048_snv <- ggplot(L2_MAG_00048_SNSsum, aes(x = name, y=(number_SNSs/number_divergent)))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  ylim(0, 1)+
  labs(x="Pond",  y= "fraction of SNVs dominated by a single allele")
ggsave("L2_MAG_00048_#SNV.png", limitsize = F, dpi = 500)

ggarrange(L2_MAG_00048_heat, L2_MAG_00048_snv, ncol=1, nrow=2, align="hv", label.x="Pond", common.legend = T, legend="right")
ggsave(filename = "L2_MAG_00048.png", limitsize=F, width=16, height=18)

#for L2_MAG_00052
L2_MAG_00052_heat <- ggplot(L2_MAG_00052, aes(x = name, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(L2_MAG_00052$groups)[order(L2_MAG_00052$all_mean)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("Erythrobacter sp. "),"assembled from GBH A"))))+
  theme(axis.title.y = element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+
  guides(fill = guide_legend(reverse=TRUE)) 
ggsave(filename="L2_MAG_00052_SNV_95_heatmap.png", dpi = 500)

L2_MAG_00052_SNSsum <- L2_MAG_00052 %>% group_by(name) %>% summarize_at(c('number_SNSs', 'number_SNVs', 'number_divergent'), sum, na.rm=T)

L2_MAG_00052_snv <- ggplot(L2_MAG_00052_SNSsum, aes(x = name, y=(number_SNSs/number_divergent)))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  ylim(0, 1)+
  labs(x="Pond",  y= "fraction of SNVs dominated by a single allele")
ggsave("L2_MAG_00052_#SNV.png", limitsize = F, dpi = 500)

ggarrange(L2_MAG_00052_heat, L2_MAG_00052_snv, ncol=1, nrow=2, align="hv", label.x="Pond", common.legend = T, legend="right")
ggsave(filename = "L2_MAG_00052.png", limitsize=F, width=16, height=18)

#for L3_MAG_00058
L3_MAG_00058_heat <- ggplot(L3_MAG_00058, aes(x = name, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(L3_MAG_00058$groups)[order(L3_MAG_00058$all_mean)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("Prosthecobacter sp. "),"assembled from Control C"))))+
  theme(axis.title.y = element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+
  guides(fill = guide_legend(reverse=TRUE)) 
ggsave(filename="L3_MAG_00058_SNV_95_heatmap.png", dpi = 500)

L3_MAG_00058_SNSsum <- L3_MAG_00058 %>% group_by(name) %>% summarize_at(c('number_SNSs', 'number_SNVs', 'number_divergent'), sum, na.rm=T)

L3_MAG_00058_snv <- ggplot(L3_MAG_00058_SNSsum, aes(x = name, y=(number_SNSs/number_divergent)))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  ylim(0, 1)+
  labs(x="Pond",  y= "fraction of SNVs dominated by a single allele")
ggsave("L3_MAG_00058_#SNV.png", limitsize = F, dpi = 500)

ggarrange(L3_MAG_00058_heat, L3_MAG_00058_snv, ncol=1, nrow=2, align="hv", label.x="Pond", common.legend = T, legend="right")
ggsave(filename = "L3_MAG_00058.png", limitsize=F, width=16, height=18)

#for L4_MAG_00099
L4_MAG_00099_heat <- ggplot(L4_MAG_00099, aes(x = name, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(L4_MAG_00099$groups)[order(L4_MAG_00099$all_mean)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("Bosea sp001713455 "),"assembled from Control D"))))+
  theme(axis.title.y = element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+
  guides(fill = guide_legend(reverse=TRUE)) 
ggsave(filename="L4_MAG_00099_SNV_95_heatmap.png", dpi = 500)

L4_MAG_00099_SNSsum <- L4_MAG_00099 %>% group_by(name) %>% summarize_at(c('number_SNSs', 'number_SNVs', 'number_divergent'), sum, na.rm=T)

L4_MAG_00099_snv <- ggplot(L4_MAG_00099_SNSsum, aes(x = name, y=(number_SNSs/number_divergent)))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  ylim(0, 1)+
  labs(x="Pond",  y= "fraction of SNVs dominated by a single allele")
ggsave("L4_MAG_00099_#SNV.png", limitsize = F, dpi = 500)

ggarrange(L4_MAG_00099_heat, L4_MAG_00099_snv, ncol=1, nrow=2, align="hv", label.x="Pond", common.legend = T, legend="right")
ggsave(filename = "L4_MAG_00099.png", limitsize=F, width=16, height=18)

#for L7_MAG_00020
L7_MAG_00020_heat <- ggplot(L7_MAG_00020, aes(x = name, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(L7_MAG_00020$groups)[order(L7_MAG_00020$all_mean)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("Sphingorhabdus_B sp. "),"assembled from GBH C"))))+
  theme(axis.title.y = element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+
  guides(fill = guide_legend(reverse=TRUE)) 
ggsave(filename="L7_MAG_00020_SNV_95_heatmap.png", dpi = 500)

L7_MAG_00020_SNSsum <- L7_MAG_00020 %>% group_by(name) %>% summarize_at(c('number_SNSs', 'number_SNVs', 'number_divergent'), sum, na.rm=T)

L7_MAG_00020_snv <- ggplot(L7_MAG_00020_SNSsum, aes(x=name, y=(number_SNSs/number_divergent)))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  ylim(0, 1)+
  labs(x="Pond",  y= "fraction of SNVs dominated by a single allele")
ggsave("L7_MAG_00020_#SNV.png", limitsize = F, dpi = 500)

ggarrange(L7_MAG_00020_heat, L7_MAG_00020_snv, ncol=1, nrow=2, align="hv", label.x="Pond", common.legend = T, legend="right")
ggsave(filename = "L7_MAG_00020.png", limitsize=F, width=16, height=18)

#for L7_MAG_00028
L7_MAG_00028_heat <- ggplot(L7_MAG_00028, aes(x = name, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(L7_MAG_00028$groups)[order(L7_MAG_00028$all_mean)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("SYFN01 sp. "),"assembled from GBH C"))))+
  theme(axis.title.y = element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+
  guides(fill = guide_legend(reverse=TRUE)) 
ggsave(filename="L7_MAG_00028_SNV_95_heatmap.png", dpi = 500)

L7_MAG_00028_SNSsum <- L7_MAG_00028 %>% group_by(name) %>% summarize_at(c('number_SNSs', 'number_SNVs', 'number_divergent'), sum, na.rm=T)

L7_MAG_00028_snv <- ggplot(L7_MAG_00028_SNSsum, aes(x = name, y=(number_SNSs/number_divergent)))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  ylim(0, 1)+
  labs(x="Pond",  y= "fraction of SNVs dominated by a single allele")
ggsave("L7_MAG_00028_#SNV.png", limitsize = F, dpi = 500)

ggarrange(L7_MAG_00028_heat, L7_MAG_00028_snv, ncol=1, nrow=2, align="hv", label.x="Pond", common.legend = T, legend="right")
ggsave(filename = "L7_MAG_00028.png", limitsize=F, width=16, height=18)

#for L7_MAG_00043
L7_MAG_00043_heat <- ggplot(L7_MAG_00043, aes(x = name, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(L7_MAG_00043$groups)[order(L7_MAG_00043$all_mean)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("Luteolibacter sp. "),"assembled from GBH C"))))+
  theme(axis.title.y = element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+
  guides(fill = guide_legend(reverse=TRUE)) 
ggsave(filename="L7_MAG_00043_SNV_95_heatmap.png", dpi = 500)

L7_MAG_00043_SNSsum <- L7_MAG_00043 %>% group_by(name) %>% summarize_at(c('number_SNSs', 'number_SNVs', 'number_divergent'), sum, na.rm=T)

L7_MAG_00043_snv <- ggplot(L7_MAG_00043_SNSsum, aes(x = name, y=(number_SNSs/number_divergent)))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  ylim(0, 1)+
  labs(x="Pond",  y= "fraction of SNVs dominated by a single allele")
ggsave("L7_MAG_00043_#SNV.png", limitsize = F, dpi = 500)

ggarrange(L7_MAG_00043_heat, L7_MAG_00043_snv, ncol=1, nrow=2, align="hv", label.x="Pond", common.legend = T, legend="right")
ggsave(filename = "L7_MAG_00043.png", limitsize=F, width=16, height=18)

#for L8_MAG_00011
L8_MAG_00011_heat <- ggplot(L8_MAG_00011, aes(x = name, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(L8_MAG_00011$groups)[order(L8_MAG_00011$all_mean)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("UBA953 sp. "),"assembled from GBH D"))))+
  theme(axis.title.y = element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+
  guides(fill = guide_legend(reverse=TRUE)) 
ggsave(filename="L8_MAG_00011_SNV_95_heatmap.png", dpi = 500)

L8_MAG_00011_SNSsum <- L8_MAG_00011 %>% group_by(name) %>% summarize_at(c('number_SNSs', 'number_SNVs', 'number_divergent'), sum, na.rm=T)

L8_MAG_00011_snv <- ggplot(L8_MAG_00011_SNSsum, aes(x = name, y=(number_SNSs/number_divergent)))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  ylim(0, 1)+
  labs(x="Pond",  y= "fraction of SNVs dominated by a single allele")
ggsave("L8_MAG_00011_#SNV.png", limitsize = F, dpi = 500)

ggarrange(L8_MAG_00011_heat, L8_MAG_00011_snv, ncol=1, nrow=2, align="hv", label.x="Pond", common.legend = T, legend="right")
ggsave(filename = "L8_MAG_00011.png", limitsize=F, width=16, height=18)

#for L8_MAG_00019
L8_MAG_00019_heat <- ggplot(L8_MAG_00019, aes(x = name, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(L8_MAG_00019$groups)[order(L8_MAG_00019$all_mean)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("UA16 family "),"assembled from GBH D"))))+
  theme(axis.title.y = element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+
  guides(fill = guide_legend(reverse=TRUE)) 
ggsave(filename="L8_MAG_00019_SNV_95_heatmap.png", dpi = 500)

L8_MAG_00019_SNSsum <- L8_MAG_00019 %>% group_by(name) %>% summarize_at(c('number_SNSs', 'number_SNVs', 'number_divergent'), sum, na.rm=T)

L8_MAG_00019_snv <- ggplot(L8_MAG_00019_SNSsum, aes(x = name, y=(number_SNSs/number_divergent)))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  ylim(0, 1)+
  labs(x="Pond",  y= "fraction of SNVs dominated by a single allele")
ggsave("L8_MAG_00019_#SNV.png", limitsize = F, dpi = 500)

ggarrange(L8_MAG_00019_heat, L8_MAG_00019_snv, ncol=1, nrow=2, align="hv", label.x="Pond", common.legend = T, legend="right")
ggsave(filename = "L8_MAG_00019.png", limitsize=F, width=16, height=18)

#for L8_MAG_00042
L8_MAG_00042_heat <- ggplot(L8_MAG_00042, aes(x = name, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(L8_MAG_00042$groups)[order(L8_MAG_00042$all_mean)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("UBA4660 sp. "),"assembled from GBH D"))))+
  theme(axis.title.y = element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+
  guides(fill = guide_legend(reverse=TRUE)) 
ggsave(filename="L8_MAG_00042_SNV_95_heatmap.png", dpi = 500)

L8_MAG_00042_SNSsum <- L8_MAG_00042 %>% group_by(name) %>% summarize_at(c('class'), sum, na.rm=T)

L8_MAG_00042_snv <- ggplot(L8_MAG_00042_SNSsum, aes(x = name, y=(number_SNSs/number_divergent)))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  ylim(0, 1)+
  labs(x="Pond",  y= "fraction of SNVs dominated by a single allele")
ggsave("L8_MAG_00042_#SNV.png", limitsize = F, dpi = 500)

ggarrange(L8_MAG_00042_heat, L8_MAG_00042_snv, ncol=1, nrow=2, align="hv", label.x="Pond", common.legend = T, legend="right")
ggsave(filename = "L8_MAG_00042.png", limitsize=F, width=16, height=18)

library(tidyverse)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggpubr)
library(forcats)

all_snv <- read_csv("all_snv_Jun_5.csv")
all_sum <- read_csv("all_snv_sum.csv")
all_sum_long <- pivot_longer(all_sum, cols=contains("S"), names_to="class", values_to="divergent_sites", values_drop_na = F)

#for I4_MAG_00006
I4_MAG_00006 <- all_snv %>% subset(mag=="I4_MAG_00006")
I4_MAG_00006_heat <- ggplot(I4_MAG_00006, aes(x = name, y = groups, fill= final_ref_freq)) +
    geom_tile()+
    scale_fill_viridis(direction=-1, na.value = "white") +
    theme_classic() +
    scale_y_discrete(limits=(I4_MAG_00006$groups)[order(I4_MAG_00006$all_mean)]) +
    labs(fill="Reference Freq.", title=(expression(paste(italic("SJAQ100 sp016735685 "),"assembled from Control B"))))+
    theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank(), 
          plot.title = element_text(hjust = 0.5), text=element_text(size=15))+
    guides(fill = guide_legend(reverse=TRUE))
ggsave(filename="I4_MAG_00006_SNV_95_heatmap.png", dpi = 500)

I4_MAG_00006_snv_frac <- ggplot(subset(all_sum, mag=="I4_MAG_00006"), aes(x = name, y=SNSs/(SNSs+SNVs)))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  ylim(0, 1)+
  labs(y= "fraction of SNVs dominated by a single allele")
ggsave("I4_MAG_00006_SNV_frac.png", limitsize = F, dpi = 500)

I4_MAG_00006_snv_sum <- ggplot(subset(all_sum_long, mag=="I4_MAG_00006"), aes(x = name, y=(((divergent_sites)/mag_length)*10^6), fill=class))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  scale_fill_manual(values=c("#E0E0E0","#424242"))+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  labs(y= "SNVs / Mbp")
ggsave("I4_MAG_00006_#SNV.png", limitsize = F, dpi = 500)

#for I4_MAG_00065
I4_MAG_00065 <- all_snv %>% subset(mag=="I4_MAG_00065")
I4_MAG_00065_heat <- ggplot(I4_MAG_00065, aes(x = name, y = groups, fill= final_ref_freq)) +
  geom_tile()+
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(I4_MAG_00065$groups)[order(I4_MAG_00065$all_mean)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("Roseomonas sp. "),"assembled from Control B"))))+
  theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank(), 
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+
  guides(fill = guide_legend(reverse=TRUE))
ggsave(filename="I4_MAG_00065_SNV_95_heatmap.png", dpi = 500)

I4_MAG_00065_snv_frac <- ggplot(subset(all_sum, mag=="I4_MAG_00065"), aes(x = name, y=SNSs/(SNSs+SNVs)))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  ylim(0, 1)+
  labs(y= "fraction of SNVs dominated by a single allele")
ggsave("I4_MAG_00065_SNV_frac.png", limitsize = F, dpi = 500)

I4_MAG_00065_snv_sum <- ggplot(subset(all_sum_long, mag=="I4_MAG_00065"), aes(x = name, y=(((divergent_sites)/mag_length)*10^6), fill=class))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  scale_fill_manual(values=c("#E0E0E0","#424242"))+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  labs(y= "SNVs / Mbp")
ggsave("I4_MAG_00065_#SNV.png", limitsize = F, dpi = 500)

#for L2_MAG_00052
L2_MAG_00052 <- all_snv %>% subset(mag=="L2_MAG_00052")
L2_MAG_00052_heat <- ggplot(L2_MAG_00052, aes(x = name, y = groups, fill= final_ref_freq)) +
  geom_tile()+
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(L2_MAG_00052$groups)[order(L2_MAG_00052$all_mean)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("Erythrobacter sp. "),"assembled from GBH A"))))+
  theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank(), 
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+
  guides(fill = guide_legend(reverse=TRUE))
ggsave(filename="L2_MAG_00052_SNV_95_heatmap.png", dpi = 500)

L2_MAG_00052_snv_frac <- ggplot(subset(all_sum, mag=="L2_MAG_00052"), aes(x = name, y=SNSs/(SNSs+SNVs)))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  ylim(0, 1)+
  labs(y= "fraction of SNVs dominated by a single allele")
ggsave("L2_MAG_00052_SNV_frac.png", limitsize = F, dpi = 500)

L2_MAG_00052_snv_sum <- ggplot(subset(all_sum_long, mag=="L2_MAG_00052"), aes(x = name, y=(((divergent_sites)/mag_length)*10^6), fill=class))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  scale_fill_manual(values=c("#E0E0E0","#424242"))+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  labs(y= "SNVs / Mbp")
ggsave("L2_MAG_00052_#SNV.png", limitsize = F, dpi = 500)

#for L3_MAG_00058
L3_MAG_00058 <- all_snv %>% subset(mag=="L3_MAG_00058")
L3_MAG_00058_heat <- ggplot(L3_MAG_00058, aes(x = name, y = groups, fill= final_ref_freq)) +
  geom_tile()+
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(L3_MAG_00058$groups)[order(L3_MAG_00058$all_mean)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("Prosthecobacter sp. "),"assembled from Control C"))))+
  theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank(), 
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+
  guides(fill = guide_legend(reverse=TRUE))
ggsave(filename="L3_MAG_00058_SNV_95_heatmap.png", dpi = 500)

L3_MAG_00058_snv_frac <- ggplot(subset(all_sum, mag=="L3_MAG_00058"), aes(x = name, y=SNSs/(SNSs+SNVs)))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  ylim(0, 1)+
  labs(y= "fraction of SNVs dominated by a single allele")
ggsave("L3_MAG_00058_SNV_frac.png", limitsize = F, dpi = 500)

L3_MAG_00058_snv_sum <- ggplot(subset(all_sum_long, mag=="L3_MAG_00058"), aes(x = name, y=(((divergent_sites)/mag_length)*10^6), fill=class))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  scale_fill_manual(values=c("#E0E0E0","#424242"))+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  labs(y= "SNVs / Mbp")
ggsave("L3_MAG_00058_#SNV.png", limitsize = F, dpi = 500)

#for L4_MAG_00099
L4_MAG_00099 <- all_snv %>% subset(mag=="L4_MAG_00099")
L4_MAG_00099_heat <- ggplot(L4_MAG_00099, aes(x = name, y = groups, fill= final_ref_freq)) +
  geom_tile()+
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(L4_MAG_00099$groups)[order(L4_MAG_00099$all_mean)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("Bosea sp001713455 "),"assembled from Control D"))))+
  theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank(), 
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+
  guides(fill = guide_legend(reverse=TRUE))
ggsave(filename="L4_MAG_00099_SNV_95_heatmap.png", dpi = 500)

L4_MAG_00099_snv_frac <- ggplot(subset(all_sum, mag=="L4_MAG_00099"), aes(x = name, y=SNSs/(SNSs+SNVs)))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  ylim(0, 1)+
  labs(y= "fraction of SNVs dominated by a single allele")
ggsave("L4_MAG_00099_SNV_frac.png", limitsize = F, dpi = 500)

L4_MAG_00099_snv_sum <- ggplot(subset(all_sum_long, mag=="L4_MAG_00099"), aes(x = name, y=(((divergent_sites)/mag_length)*10^6), fill=class))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  scale_fill_manual(values=c("#E0E0E0","#424242"))+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  labs(y= "SNVs / Mbp")
ggsave("L4_MAG_00099_#SNV.png", limitsize = F, dpi = 500)

#for L7_MAG_00020
L7_MAG_00020 <- all_snv %>% subset(mag=="L7_MAG_00020")
L7_MAG_00020_heat <- ggplot(L7_MAG_00020, aes(x = name, y = groups, fill= final_ref_freq)) +
  geom_tile()+
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(L7_MAG_00020$groups)[order(L7_MAG_00020$all_mean)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("Sphingorhabdus_B sp. "),"assembled from GBH C"))))+
  theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank(), 
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+
  guides(fill = guide_legend(reverse=TRUE))
ggsave(filename="L7_MAG_00020_SNV_95_heatmap.png", dpi = 500)

L7_MAG_00020_snv_frac <- ggplot(subset(all_sum, mag=="L7_MAG_00020"), aes(x = name, y=SNSs/(SNSs+SNVs)))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  ylim(0, 1)+
  labs(y= "fraction of SNVs dominated by a single allele")
ggsave("L7_MAG_00020_SNV_frac.png", limitsize = F, dpi = 500)

L7_MAG_00020_snv_sum <- ggplot(subset(all_sum_long, mag=="L7_MAG_00020"), aes(x = name, y=(((divergent_sites)/mag_length)*10^6), fill=class))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  scale_fill_manual(values=c("#E0E0E0","#424242"))+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  labs(y= "SNVs / Mbp")
ggsave("L7_MAG_00020_#SNV.png", limitsize = F, dpi = 500)

#for L7_MAG_00028
L7_MAG_00028 <- all_snv %>% subset(mag=="L7_MAG_00028")
L7_MAG_00028_heat <- ggplot(L7_MAG_00028, aes(x = name, y = groups, fill= final_ref_freq)) +
  geom_tile()+
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(L7_MAG_00028$groups)[order(L7_MAG_00028$all_mean)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("SYFN01 sp. "),"assembled from GBH C"))))+
  theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank(), 
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+
  guides(fill = guide_legend(reverse=TRUE))
ggsave(filename="L7_MAG_00028_SNV_95_heatmap.png", dpi = 500)

L7_MAG_00028_snv_frac <- ggplot(subset(all_sum, mag=="L7_MAG_00028"), aes(x = name, y=SNSs/(SNSs+SNVs)))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  ylim(0, 1)+
  labs(y= "fraction of SNVs dominated by a single allele")
ggsave("L7_MAG_00028_SNV_frac.png", limitsize = F, dpi = 500)

L7_MAG_00028_snv_sum <- ggplot(subset(all_sum_long, mag=="L7_MAG_00028"), aes(x = name, y=(((divergent_sites)/mag_length)*10^6), fill=class))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  scale_fill_manual(values=c("#E0E0E0","#424242"))+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  labs(y= "SNVs / Mbp")
ggsave("L7_MAG_00028_#SNV.png", limitsize = F, dpi = 500)

#for L7_MAG_00043
L7_MAG_00043 <- all_snv %>% subset(mag=="L7_MAG_00043")
L7_MAG_00043_heat <- ggplot(L7_MAG_00043, aes(x = name, y = groups, fill= final_ref_freq)) +
  geom_tile()+
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(L7_MAG_00043$groups)[order(L7_MAG_00043$all_mean)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("Luteolibacter sp. "),"assembled from GBH C"))))+
  theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank(), 
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+
  guides(fill = guide_legend(reverse=TRUE))
ggsave(filename="L7_MAG_00043_SNV_95_heatmap.png", dpi = 500)

L7_MAG_00043_snv_frac <- ggplot(subset(all_sum, mag=="L7_MAG_00043"), aes(x = name, y=SNSs/(SNSs+SNVs)))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  ylim(0, 1)+
  labs(y= "fraction of SNVs dominated by a single allele")
ggsave("L7_MAG_00043_SNV_frac.png", limitsize = F, dpi = 500)

L7_MAG_00043_snv_sum <- ggplot(subset(all_sum_long, mag=="L7_MAG_00043"), aes(x = name, y=(((divergent_sites)/mag_length)*10^6), fill=class))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  scale_fill_manual(values=c("#E0E0E0","#424242"))+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  labs(y= "SNVs / Mbp")
ggsave("L7_MAG_00043_#SNV.png", limitsize = F, dpi = 500)

#for L8_MAG_00011
L8_MAG_00011 <- all_snv %>% subset(mag=="L8_MAG_00011")
L8_MAG_00011_heat <- ggplot(L8_MAG_00011, aes(x = name, y = groups, fill= final_ref_freq)) +
  geom_tile()+
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(L8_MAG_00011$groups)[order(L8_MAG_00011$all_mean)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("UBA953 sp. "),"assembled from GBH D"))))+
  theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank(), 
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+
  guides(fill = guide_legend(reverse=TRUE))
ggsave(filename="L8_MAG_00011_SNV_95_heatmap.png", dpi = 500)

L8_MAG_00011_snv_frac <- ggplot(subset(all_sum, mag=="L8_MAG_00011"), aes(x = name, y=SNSs/(SNSs+SNVs)))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  ylim(0, 1)+
  labs(y= "fraction of SNVs dominated by a single allele")
ggsave("L8_MAG_00011_SNV_frac.png", limitsize = F, dpi = 500)

L8_MAG_00011_snv_sum <- ggplot(subset(all_sum_long, mag=="L8_MAG_00011"), aes(x = name, y=(((divergent_sites)/mag_length)*10^6), fill=class))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  scale_fill_manual(values=c("#E0E0E0","#424242"))+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  labs(y= "SNVs / Mbp")
ggsave("L8_MAG_00011_#SNV.png", limitsize = F, dpi = 500)

#for L8_MAG_00019
L8_MAG_00019 <- all_snv %>% subset(mag=="L8_MAG_00019")
L8_MAG_00019_heat <- ggplot(L8_MAG_00019, aes(x = name, y = groups, fill= final_ref_freq)) +
  geom_tile()+
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(L8_MAG_00019$groups)[order(L8_MAG_00019$all_mean)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("UA16 family "),"assembled from GBH D"))))+
  theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank(), 
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+
  guides(fill = guide_legend(reverse=TRUE))
ggsave(filename="L8_MAG_00019_SNV_95_heatmap.png", dpi = 500)

L8_MAG_00019_snv_frac <- ggplot(subset(all_sum, mag=="L8_MAG_00019"), aes(x = name, y=SNSs/(SNSs+SNVs)))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  ylim(0, 1)+
  labs(y= "fraction of SNVs dominated by a single allele")
ggsave("L8_MAG_00019_SNV_frac.png", limitsize = F, dpi = 500)

L8_MAG_00019_snv_sum <- ggplot(subset(all_sum_long, mag=="L8_MAG_00019"), aes(x = name, y=(((divergent_sites)/mag_length)*10^6), fill=class))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  scale_fill_manual(values=c("#E0E0E0","#424242"))+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  labs(y= "SNVs / Mbp")
ggsave("L8_MAG_00019_#SNV.png", limitsize = F, dpi = 500)

#for L8_MAG_00042
L8_MAG_00042 <- all_snv %>% subset(mag=="L8_MAG_00042") %>% arrange(all_mean)
ggplot(L8_MAG_00042, aes(x = name, y = groups, fill= final_ref_freq)) +
  geom_tile()+
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  labs(fill="Reference Freq.", title=(expression(paste(italic("UBA4660 sp. "),"assembled from GBH D"))))+
  theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank(), 
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+
  guides(fill = guide_legend(reverse=TRUE))
ggsave(filename="L8_MAG_00042_SNV_95_heatmap.png", dpi = 500)

L8_MAG_00042_snv_frac <- ggplot(subset(all_sum, mag=="L8_MAG_00042"), aes(x = name, y=SNSs/(SNSs+SNVs)))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  ylim(0, 1)+
  labs(y= "fraction of SNVs dominated by a single allele")
ggsave("L8_MAG_00042_SNV_frac.png", limitsize = F, dpi = 500)

L8_MAG_00042_snv_sum <- ggplot(subset(all_sum_long, mag=="L8_MAG_00042"), aes(x = name, y=(((divergent_sites)/mag_length)*10^6), fill=class))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  scale_fill_manual(values=c("#E0E0E0","#424242"))+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  labs(y= "SNVs / Mbp")
ggsave("L8_MAG_00042_#SNV.png", limitsize = F, dpi = 500)

ggarrange(I4_MAG_00006_heat, I4_MAG_00065_heat, L3_MAG_00058_heat, L7_MAG_00020_heat, L8_MAG_00011_heat, L8_MAG_00019_heat,
          nrow=1, ncol=6, align="h", common.legend = T, legend="right")
ggsave("heatmaps_selection_y.png", limitsize = F, width=50, height=8)

library(tidyverse)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggpubr)

all_snv <- read_csv("all_snv_Jun_5.csv")
select_snv<- subset(all_snv, mag=="I4_MAG_00006" | mag=="I4_MAG_00065" | mag=="L3_MAG_00058" | mag=="L7_MAG_00020" | mag=="L8_MAG_00011" | mag=="L8_MAG_00019")
other_snv<- subset(all_snv, mag=="L2_MAG_00052" | mag=="L4_MAG_00099" | mag=="L7_MAG_00028" | mag=="L7_MAG_00043" | mag=="L8_MAG_00042")
all_sum <- read_csv("all_snv_sum.csv")
select_sum<- subset(all_sum, mag=="I4_MAG_00006" | mag=="I4_MAG_00065" | mag=="L3_MAG_00058" | mag=="L7_MAG_00020" | mag=="L8_MAG_00011" | mag=="L8_MAG_00019")
other_sum<- subset(all_sum, mag=="L2_MAG_00052" | mag=="L4_MAG_00099" | mag=="L7_MAG_00028" | mag=="L7_MAG_00043" | mag=="L8_MAG_00042")
all_sum_long <- pivot_longer(all_sum, cols=contains("S"), names_to="class", values_to="divergent_sites", values_drop_na = F)
select_sum_long<- subset(all_sum_long, mag=="I4_MAG_00006" | mag=="I4_MAG_00065" | mag=="L3_MAG_00058" | mag=="L7_MAG_00020" | mag=="L8_MAG_00011" | mag=="L8_MAG_00019")
other_sum_long<- subset(all_sum_long, mag=="L2_MAG_00052" | mag=="L4_MAG_00099" | mag=="L7_MAG_00028" | mag=="L7_MAG_00043" | mag=="L8_MAG_00042")

select_mag.labs<-(c(I4_MAG_00006 = "SJAQ100 sp016735685 assembled from Control B", I4_MAG_00065 = "Roseomonas sp. assembled from Control B", 
             L3_MAG_00058 ="Prosthecobacter sp. assembled from Control C", L7_MAG_00020 ="Sphingorhabdus_B sp. assembled from GBH C", 
             L8_MAG_00011 ="UBA953 sp. assembled from GBH D", L8_MAG_00019 ="UA16 family assembled from GBH D"))

other_mag.labs<-(c(L2_MAG_00052= "Erythrobacter sp. assembled from GBH A", L4_MAG_00099= "Bosea sp001713455 assembled from Control D", L7_MAG_00028= "SYFN01 sp. assembled from GBH C",
                 L7_MAG_00043= "Luteolibacter sp. assembled from GBH C", L8_MAG_00042= "UBA4660 sp. assembled from GBH D"))

select_snv_heat <- ggplot(select_snv, aes(x = name, y = reorder(groups, all_mean), fill= final_ref_freq)) +
  geom_tile()+
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank(), 
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+
  guides(fill = guide_legend(reverse=TRUE))+
  facet_wrap(~mag, nrow=1, ncol=6, scales="free", labeller=labeller(mag=select_mag.labs))
ggsave("select_snv_heat.png", limitsize=F, dpi=400, width=48, height=8)

select_snv_frac <- ggplot(select_sum, aes(x = name, y=SNSs/(SNSs+SNVs)))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15), panel.spacing=unit(1, "cm"))+
  ylim(0, 1)+
  labs(y= "fraction of SNVs dominated by a single allele")+
  facet_wrap(~mag, nrow=1, scales="free", labeller=labeller(mag=select_mag.labs))
ggsave("select_snv_frac.png", limitsize = F, width=48, height=8)

select_snv_sum <- ggplot(select_sum_long, aes(x = name, y=(((divergent_sites)/mag_length)*10^6), fill=class))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  scale_fill_manual(values=c("#E0E0E0","#424242"))+
  theme(legend.title=element_blank(),text=element_text(size=15), panel.spacing=unit(1, "cm"))+
  labs(y= "SNVs / Mbp")+
  facet_wrap(~mag, nrow=1, scales="free", labeller=labeller(mag=select_mag.labs))
ggsave("select_snv_sum.png", limitsize = F, width=48, height=8)

ggarrange(select_snv_heat, select_snv_frac, align="hv", ncol=1, nrow=2)
ggsave("heat_frac_selected.png", limitsize=F, width=48, height=16) 

ggarrange(select_snv_heat, select_snv_sum, align="hv", ncol=1, nrow=2)
ggsave("heat_sum_selected.png", limitsize=F, width=48, height=16) 

ggarrange(select_snv_heat, select_snv_sum, select_snv_frac, align="hv", ncol=1, nrow=3)
ggsave("heat_sum_frac_selected.png", limitsize=F, width=48, height=24) 


other_snv_heat <- ggplot(other_snv, aes(x = name, y = reorder(groups, all_mean), fill= final_ref_freq)) +
  geom_tile()+
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank(), 
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+
  guides(fill = guide_legend(reverse=TRUE))+
  facet_wrap(~mag, nrow=1, ncol=6, scales="free", labeller=labeller(mag=other_mag.labs))
ggsave("other_snv_heat_.png", limitsize=F, dpi=400, width=48, height=8)

other_snv_frac <- ggplot(other_sum, aes(x = name, y=SNSs/(SNSs+SNVs)))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15), panel.spacing=unit(1, "cm"))+
  ylim(0, 1)+
  labs(y= "fraction of SNVs dominated by a single allele")+
  facet_wrap(~mag, nrow=1, scales="free", labeller=labeller(mag=other_mag.labs))
ggsave("other_snv_frac.png", limitsize = F, width=48, height=8)

other_snv_sum <- ggplot(other_sum_long, aes(x = name, y=(((divergent_sites)/mag_length)*10^6), fill=class))+
  geom_bar(stat="identity")+ 
  theme_classic()+
  scale_fill_manual(values=c("#E0E0E0","#424242"))+
  theme(legend.title=element_blank(),text=element_text(size=15), panel.spacing=unit(1, "cm"))+
  labs(y= "SNVs / Mbp")+
  facet_wrap(~mag, nrow=1, scales="free", labeller=labeller(mag=other_mag.labs))
ggsave("other_snv_sum.png", limitsize = F, width=48, height=8)

ggarrange(other_snv_heat, other_snv_frac, align="hv", ncol=1, nrow=2)
ggsave("heat_frac_other.png", limitsize=F, width=48, height=16) 

ggarrange(other_snv_heat, other_snv_sum, align="hv", ncol=1, nrow=2)
ggsave("heat_sum_other.png", limitsize=F, width=48, height=16) 

ggarrange(other_snv_heat, other_snv_sum, other_snv_frac, align="hv", ncol=1, nrow=3)
ggsave("heat_sum_frac_other.png", limitsize=F, width=48, height=24) 

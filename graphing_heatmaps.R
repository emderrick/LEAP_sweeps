library(tidyverse)
library(dplyr)
library(ggplot2)
library(viridis)
library(stringr)
library(ggpubr)

#read in all the SNV files
I4_MAG_00006 <- read_csv("all_I4_MAG_00006_SNVs.csv")
I4_MAG_00065 <- read_csv("all_I4_MAG_00065_SNVs.csv")
L2_MAG_00048 <- read_csv("all_L2_MAG_00048_SNVs.csv")
L2_MAG_00052 <- read_csv("all_L2_MAG_00052_SNVs.csv")
L3_MAG_00058 <- read_csv("all_L3_MAG_00058_SNVs.csv")
L4_MAG_00099 <- read_csv("all_L4_MAG_00099_SNVs.csv")
L7_MAG_00020 <- read_csv("all_L7_MAG_00020_SNVs.csv")
L7_MAG_00028 <- read_csv("all_L7_MAG_00028_SNVs.csv")
L7_MAG_00043 <- read_csv("all_L7_MAG_00043_SNVs.csv")
L8_MAG_00011 <- read_csv("all_L8_MAG_00011_SNVs.csv")
L8_MAG_00019 <- read_csv("all_L8_MAG_00019_SNVs.csv")
L8_MAG_00042 <- read_csv("all_L8_MAG_00042_SNVs.csv")

#for I4_MAG_00006
I4_MAG_00006 <- I4_MAG_00006 %>% mutate(new_name = case_when(timepoint%>%substr(1,2) == "K1" ~ "Control A", timepoint%>%substr(1,2) == "I4" ~ "Control B",
                                                             timepoint%>%substr(1,2) == "L3" ~ "Control C", timepoint%>%substr(1,2) == "L4" ~ "Control D",
                                                             timepoint%>%substr(1,2) == "I8" ~ "Control E", timepoint%>%substr(1,2) == "L2" ~ "GBH A",
                                                             timepoint%>%substr(1,2)== "L6" ~ "GBH B", timepoint%>%substr(1,2) == "L7" ~ "GBH C",
                                                             timepoint%>%substr(1,2) == "L8" ~ "GBH D"))

I4_MAG_00006_index<-tapply(I4_MAG_00006$final_ref_freq,I4_MAG_00006$groups,FUN=median,na.rm=T)
I4_MAG_00006_index_df<-data.frame(groups=names(I4_MAG_00006_index),median=I4_MAG_00006_index)
I4_MAG_00006_median<- right_join(I4_MAG_00006_index_df, I4_MAG_00006, by=c("groups"))

I4_MAG_00006_heat <- ggplot(I4_MAG_00006_median, aes(x = new_name, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(I4_MAG_00006_median$groups)[order(I4_MAG_00006_median$median)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("SJAQ100 sp016735685 "),"assembled from Control B"))))+
  theme(axis.title.y = element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank(), 
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+
  guides(fill = guide_legend(reverse=TRUE))
ggsave(filename="I4_MAG_00006_SNV_95_heatmap.png", dpi = 500)

I4_MAG_00006$simple_class<-I4_MAG_00006$class%>%str_sub(-3,-1)

I4_MAG_00006_snv <- ggplot(I4_MAG_00006, aes(x=new_name, y=(((number_divergent/mag_length)*10^6)), fill=simple_class))+
  geom_bar(stat="identity")+ 
  scale_fill_manual(values=c("#E0E0E0","#424242"))+
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  labs(x="Pond",  y= "SNVs / Mbp")
ggsave("I4_MAG_00006_#SNV.png", limitsize = F, dpi = 500)

ggarrange(I4_MAG_00006_heat, I4_MAG_00006_snv,
          ncol=1, nrow=2,
          align="hv")
ggsave(filename = "I4_MAG_00006.png", limitsize=F, width=16, height=18)

#for I4_MAG_00065
I4_MAG_00065 <- I4_MAG_00065 %>% mutate(new_name = case_when(timepoint%>%substr(1,2) == "K1" ~ "Control A", timepoint%>%substr(1,2) == "I4" ~ "Control B",
                                                             timepoint%>%substr(1,2) == "L3" ~ "Control C", timepoint%>%substr(1,2) == "L4" ~ "Control D",
                                                             timepoint%>%substr(1,2) == "I8" ~ "Control E", timepoint%>%substr(1,2) == "L2" ~ "GBH A",
                                                             timepoint%>%substr(1,2)== "L6" ~ "GBH B", timepoint%>%substr(1,2) == "L7" ~ "GBH C",
                                                             timepoint%>%substr(1,2) == "L8" ~ "GBH D"))

I4_MAG_00065_index<-tapply(I4_MAG_00065$final_ref_freq,I4_MAG_00065$groups,FUN=median,na.rm=T)
I4_MAG_00065_index_df<-data.frame(groups=names(I4_MAG_00065_index),median=I4_MAG_00065_index)
I4_MAG_00065_median<- right_join(I4_MAG_00065_index_df, I4_MAG_00065, by=c("groups"))

I4_MAG_00065_heat <- ggplot(I4_MAG_00065_median, aes(x = new_name, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(I4_MAG_00065_median$groups)[order(I4_MAG_00065_median$median)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("Roseomonas sp. "),"assembled from Control B"))))+
  theme(axis.title.y = element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+ 
ggsave(filename="I4_MAG_00065_SNV_95_heatmap.png", dpi = 500)

I4_MAG_00065$simple_class<-I4_MAG_00065$class%>%str_sub(-3,-1)

I4_MAG_00065_snv <- ggplot(I4_MAG_00065, aes(x=new_name, y=(((number_divergent/mag_length)*10^6)), fill=simple_class))+
  geom_bar(stat="identity")+ 
  scale_fill_manual(values=c("#E0E0E0","#424242"))+
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  labs(x="Pond",  y= "SNVs / Mbp")
ggsave("I4_MAG_00065_#SNV.png", limitsize = F, dpi = 500)

ggarrange(I4_MAG_00065_heat, I4_MAG_00065_snv,
          ncol=1, nrow=2,
          align="hv")
ggsave(filename = "I4_MAG_00065.png", limitsize=F, width=16, height=18)

#for L2_MAG_00048
L2_MAG_00048 <- L2_MAG_00048 %>% mutate(new_name = case_when(timepoint%>%substr(1,2) == "K1" ~ "Control A", timepoint%>%substr(1,2) == "I4" ~ "Control B",
                                                             timepoint%>%substr(1,2) == "L3" ~ "Control C", timepoint%>%substr(1,2) == "L4" ~ "Control D",
                                                             timepoint%>%substr(1,2) == "I8" ~ "Control E", timepoint%>%substr(1,2) == "L2" ~ "GBH A",
                                                             timepoint%>%substr(1,2)== "L6" ~ "GBH B", timepoint%>%substr(1,2) == "L7" ~ "GBH C",
                                                             timepoint%>%substr(1,2) == "L8" ~ "GBH D"))

L2_MAG_00048_index<-tapply(L2_MAG_00048$final_ref_freq,L2_MAG_00048$groups,FUN=median,na.rm=T)
L2_MAG_00048_index_df<-data.frame(groups=names(L2_MAG_00048_index),median=L2_MAG_00048_index)
L2_MAG_00048_median<- right_join(L2_MAG_00048_index_df, L2_MAG_00048, by=c("groups"))

L2_MAG_00048_heat <- ggplot(L2_MAG_00048_median, aes(x = new_name, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(L2_MAG_00048_median$groups)[order(L2_MAG_00048_median$median)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("Flavobacterium sp. "),"assembled from GBH A"))))+
  theme(axis.title.y = element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+ 
ggsave(filename="L2_MAG_00048_SNV_95_heatmap.png", dpi = 500)

L2_MAG_00048$simple_class<-L2_MAG_00048$class%>%str_sub(-3,-1)

L2_MAG_00048_snv <- ggplot(L2_MAG_00048, aes(x=new_name, y=(((number_divergent/mag_length)*10^6)), fill=simple_class))+
  geom_bar(stat="identity")+ 
  scale_fill_manual(values=c("#E0E0E0","#424242"))+
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  labs(x="Pond",  y= "SNVs / Mbp")
ggsave("L2_MAG_00048_#SNV.png", limitsize = F, dpi = 500)

ggarrange(L2_MAG_00048_heat, L2_MAG_00048_snv,
          ncol=1, nrow=2,
          align="hv")
ggsave(filename = "L2_MAG_00048.png", limitsize=F, width=16, height=18)

#for L2_MAG_00052
L2_MAG_00052 <- L2_MAG_00052 %>% mutate(new_name = case_when(timepoint%>%substr(1,2) == "K1" ~ "Control A", timepoint%>%substr(1,2) == "I4" ~ "Control B",
                                                             timepoint%>%substr(1,2) == "L3" ~ "Control C", timepoint%>%substr(1,2) == "L4" ~ "Control D",
                                                             timepoint%>%substr(1,2) == "I8" ~ "Control E", timepoint%>%substr(1,2) == "L2" ~ "GBH A",
                                                             timepoint%>%substr(1,2)== "L6" ~ "GBH B", timepoint%>%substr(1,2) == "L7" ~ "GBH C",
                                                             timepoint%>%substr(1,2) == "L8" ~ "GBH D"))

L2_MAG_00052_index<-tapply(L2_MAG_00052$final_ref_freq,L2_MAG_00052$groups,FUN=median,na.rm=T)
L2_MAG_00052_index_df<-data.frame(groups=names(L2_MAG_00052_index),median=L2_MAG_00052_index)
L2_MAG_00052_median<- right_join(L2_MAG_00052_index_df, L2_MAG_00052, by=c("groups"))

L2_MAG_00052_heat <- ggplot(L2_MAG_00052_median, aes(x = new_name, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(L2_MAG_00052_median$groups)[order(L2_MAG_00052_median$median)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("Erythrobacter sp. "),"assembled from GBH A"))))+
  theme(axis.title.y = element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+ 
ggsave(filename="L2_MAG_00052_SNV_95_heatmap.png", dpi = 500)

L2_MAG_00052$simple_class<-L2_MAG_00052$class%>%str_sub(-3,-1)

L2_MAG_00052_snv <- ggplot(L2_MAG_00052, aes(x=new_name, y=(((number_divergent/mag_length)*10^6)), fill=simple_class))+
  geom_bar(stat="identity")+ 
  scale_fill_manual(values=c("#E0E0E0","#424242"))+
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  labs(x="Pond",  y= "SNVs / Mbp")
ggsave("L2_MAG_00052_#SNV.png", limitsize = F, dpi = 500)

ggarrange(L2_MAG_00052_heat, L2_MAG_00052_snv,
          ncol=1, nrow=2,
          align="hv")
ggsave(filename = "L2_MAG_00052.png", limitsize=F, width=16, height=18)

#for L3_MAG_00058
L3_MAG_00058 <- L3_MAG_00058 %>% mutate(new_name = case_when(timepoint%>%substr(1,2) == "K1" ~ "Control A", timepoint%>%substr(1,2) == "I4" ~ "Control B",
                                                             timepoint%>%substr(1,2) == "L3" ~ "Control C", timepoint%>%substr(1,2) == "L4" ~ "Control D",
                                                             timepoint%>%substr(1,2) == "I8" ~ "Control E", timepoint%>%substr(1,2) == "L2" ~ "GBH A",
                                                             timepoint%>%substr(1,2)== "L6" ~ "GBH B", timepoint%>%substr(1,2) == "L7" ~ "GBH C",
                                                             timepoint%>%substr(1,2) == "L8" ~ "GBH D"))

L3_MAG_00058_index<-tapply(L3_MAG_00058$final_ref_freq,L3_MAG_00058$groups,FUN=median,na.rm=T)
L3_MAG_00058_index_df<-data.frame(groups=names(L3_MAG_00058_index),median=L3_MAG_00058_index)
L3_MAG_00058_median<- right_join(L3_MAG_00058_index_df, L3_MAG_00058, by=c("groups"))

L3_MAG_00058_heat <- ggplot(L3_MAG_00058_median, aes(x = new_name, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(L3_MAG_00058_median$groups)[order(L3_MAG_00058_median$median)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("Prosthecobacter sp. "),"assembled from Control C"))))+
  theme(axis.title.y = element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+ 
ggsave(filename="L3_MAG_00058_SNV_95_heatmap.png", dpi = 500)

L3_MAG_00058$simple_class<-L3_MAG_00058$class%>%str_sub(-3,-1)

L3_MAG_00058_snv <- ggplot(L3_MAG_00058, aes(x=new_name, y=(((number_divergent/mag_length)*10^6)), fill=simple_class))+
  geom_bar(stat="identity")+ 
  scale_fill_manual(values=c("#E0E0E0","#424242"))+
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  labs(x="Pond",  y= "SNVs / Mbp")
ggsave("L3_MAG_00058_#SNV.png", limitsize = F, dpi = 500)

ggarrange(L3_MAG_00058_heat, L3_MAG_00058_snv,
          ncol=1, nrow=2,
          align="hv")
ggsave(filename = "L3_MAG_00058.png", limitsize=F, width=16, height=18)

#for L4_MAG_00099
L4_MAG_00099 <- L4_MAG_00099 %>% mutate(new_name = case_when(timepoint%>%substr(1,2) == "K1" ~ "Control A", timepoint%>%substr(1,2) == "I4" ~ "Control B",
                                                             timepoint%>%substr(1,2) == "L3" ~ "Control C", timepoint%>%substr(1,2) == "L4" ~ "Control D",
                                                             timepoint%>%substr(1,2) == "I8" ~ "Control E", timepoint%>%substr(1,2) == "L2" ~ "GBH A",
                                                             timepoint%>%substr(1,2)== "L6" ~ "GBH B", timepoint%>%substr(1,2) == "L7" ~ "GBH C",
                                                             timepoint%>%substr(1,2) == "L8" ~ "GBH D"))

L4_MAG_00099_index<-tapply(L4_MAG_00099$final_ref_freq,L4_MAG_00099$groups,FUN=median,na.rm=T)
L4_MAG_00099_index_df<-data.frame(groups=names(L4_MAG_00099_index),median=L4_MAG_00099_index)
L4_MAG_00099_median<- right_join(L4_MAG_00099_index_df, L4_MAG_00099, by=c("groups"))

L4_MAG_00099_heat <- ggplot(L4_MAG_00099_median, aes(x = new_name, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(L4_MAG_00099_median$groups)[order(L4_MAG_00099_median$median)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("Bosea sp001713455 "),"assembled from Control D"))))+
  theme(axis.title.y = element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+ 
ggsave(filename="L4_MAG_00099_SNV_95_heatmap.png", dpi = 500)

L4_MAG_00099$simple_class<-L4_MAG_00099$class%>%str_sub(-3,-1)

L4_MAG_00099_snv <- ggplot(L4_MAG_00099, aes(x=new_name, y=(((number_divergent/mag_length)*10^6)), fill=simple_class))+
  geom_bar(stat="identity")+ 
  scale_fill_manual(values=c("#E0E0E0","#424242"))+
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  labs(x="Pond",  y= "SNVs / Mbp")
ggsave("L4_MAG_00099_#SNV.png", limitsize = F, dpi = 500)

ggarrange(L4_MAG_00099_heat, L4_MAG_00099_snv,
          ncol=1, nrow=2,
          align="hv")
ggsave(filename = "L4_MAG_00099.png", limitsize=F, width=16, height=18)

#for L7_MAG_00020
L7_MAG_00020 <- L7_MAG_00020 %>% mutate(new_name = case_when(timepoint%>%substr(1,2) == "K1" ~ "Control A", timepoint%>%substr(1,2) == "I4" ~ "Control B",
                                                             timepoint%>%substr(1,2) == "L3" ~ "Control C", timepoint%>%substr(1,2) == "L4" ~ "Control D",
                                                             timepoint%>%substr(1,2) == "I8" ~ "Control E", timepoint%>%substr(1,2) == "L2" ~ "GBH A",
                                                             timepoint%>%substr(1,2)== "L6" ~ "GBH B", timepoint%>%substr(1,2) == "L7" ~ "GBH C",
                                                             timepoint%>%substr(1,2) == "L8" ~ "GBH D"))
L7_MAG_00020$time2<-as.numeric(L7_MAG_00020$timepoint%>%substr(9,9))+1
L7_MAG_00020$name<-paste(L7_MAG_00020$new_name, " Time ", L7_MAG_00020$time2)
L7_MAG_00020_index<-tapply(L7_MAG_00020$final_ref_freq,L7_MAG_00020$groups,FUN=median,na.rm=T)
L7_MAG_00020_index_df<-data.frame(groups=names(L7_MAG_00020_index),median=L7_MAG_00020_index)
L7_MAG_00020_median<- right_join(L7_MAG_00020_index_df, L7_MAG_00020, by=c("groups"))

L7_MAG_00020_heat <- ggplot(L7_MAG_00020_median, aes(x = name, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(L7_MAG_00020_median$groups)[order(L7_MAG_00020_median$median)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("Sphingorhabdus_B sp. "),"assembled from GBH C"))))+
  theme(axis.title.y = element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+ 
ggsave(filename="L7_MAG_00020_SNV_95_heatmap.png", dpi = 500)

L7_MAG_00020$simple_class<-L7_MAG_00020$class%>%str_sub(-3,-1)

L7_MAG_00020_snv <- ggplot(L7_MAG_00020, aes(x=name, y=(((number_divergent/mag_length)*10^6)), fill=simple_class))+
  geom_bar(stat="identity")+ 
  scale_fill_manual(values=c("#E0E0E0","#424242"))+
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  labs(x="Pond",  y= "SNVs / Mbp")
ggsave("L7_MAG_00020_#SNV.png", limitsize = F, dpi = 500)

ggarrange(L7_MAG_00020_heat, L7_MAG_00020_snv,
          ncol=1, nrow=2,
          align="hv")
ggsave(filename = "L7_MAG_00020.png", limitsize=F, width=16, height=18)

#for L7_MAG_00028
L7_MAG_00028 <- L7_MAG_00028 %>% mutate(new_name = case_when(timepoint%>%substr(1,2) == "K1" ~ "Control A", timepoint%>%substr(1,2) == "I4" ~ "Control B",
                                                             timepoint%>%substr(1,2) == "L3" ~ "Control C", timepoint%>%substr(1,2) == "L4" ~ "Control D",
                                                             timepoint%>%substr(1,2) == "I8" ~ "Control E", timepoint%>%substr(1,2) == "L2" ~ "GBH A",
                                                             timepoint%>%substr(1,2)== "L6" ~ "GBH B", timepoint%>%substr(1,2) == "L7" ~ "GBH C",
                                                             timepoint%>%substr(1,2) == "L8" ~ "GBH D"))

L7_MAG_00028_index<-tapply(L7_MAG_00028$final_ref_freq,L7_MAG_00028$groups,FUN=median,na.rm=T)
L7_MAG_00028_index_df<-data.frame(groups=names(L7_MAG_00028_index),median=L7_MAG_00028_index)
L7_MAG_00028_median<- right_join(L7_MAG_00028_index_df, L7_MAG_00028, by=c("groups"))

L7_MAG_00028_heat <- ggplot(L7_MAG_00028_median, aes(x = new_name, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(L7_MAG_00028_median$groups)[order(L7_MAG_00028_median$median)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("SYFN01 sp. "),"assembled from GBH C"))))+
  theme(axis.title.y = element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+ 
ggsave(filename="L7_MAG_00028_SNV_95_heatmap.png", dpi = 500)

L7_MAG_00028$simple_class<-L7_MAG_00028$class%>%str_sub(-3,-1)

L7_MAG_00028_snv <- ggplot(L7_MAG_00028, aes(x=new_name, y=(((number_divergent/mag_length)*10^6)), fill=simple_class))+
  geom_bar(stat="identity")+ 
  scale_fill_manual(values=c("#E0E0E0","#424242"))+
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  labs(x="Pond",  y= "SNVs / Mbp")
ggsave("L7_MAG_00028_#SNV.png", limitsize = F, dpi = 500)

ggarrange(L7_MAG_00028_heat, L7_MAG_00028_snv,
          ncol=1, nrow=2,
          align="hv")
ggsave(filename = "L7_MAG_00028.png", limitsize=F, width=16, height=18)

#for L7_MAG_00043
L7_MAG_00043 <- L7_MAG_00043 %>% mutate(new_name = case_when(timepoint%>%substr(1,2) == "K1" ~ "Control A", timepoint%>%substr(1,2) == "I4" ~ "Control B",
                                                             timepoint%>%substr(1,2) == "L3" ~ "Control C", timepoint%>%substr(1,2) == "L4" ~ "Control D",
                                                             timepoint%>%substr(1,2) == "I8" ~ "Control E", timepoint%>%substr(1,2) == "L2" ~ "GBH A",
                                                             timepoint%>%substr(1,2)== "L6" ~ "GBH B", timepoint%>%substr(1,2) == "L7" ~ "GBH C",
                                                             timepoint%>%substr(1,2) == "L8" ~ "GBH D"))

L7_MAG_00043_index<-tapply(L7_MAG_00043$final_ref_freq,L7_MAG_00043$groups,FUN=median,na.rm=T)
L7_MAG_00043_index_df<-data.frame(groups=names(L7_MAG_00043_index),median=L7_MAG_00043_index)
L7_MAG_00043_median<- right_join(L7_MAG_00043_index_df, L7_MAG_00043, by=c("groups"))

L7_MAG_00043_heat <- ggplot(L7_MAG_00043_median, aes(x = new_name, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(L7_MAG_00043_median$groups)[order(L7_MAG_00043_median$median)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("Luteolibacter sp. "),"assembled from GBH C"))))+
  theme(axis.title.y = element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+ 
ggsave(filename="L7_MAG_00043_SNV_95_heatmap.png", dpi = 500)

L7_MAG_00043$simple_class<-L7_MAG_00043$class%>%str_sub(-3,-1)

L7_MAG_00043_snv <- ggplot(L7_MAG_00043, aes(x=new_name, y=(((number_divergent/mag_length)*10^6)), fill=simple_class))+
  geom_bar(stat="identity")+ 
  scale_fill_manual(values=c("#E0E0E0","#424242"))+
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  labs(x="Pond",  y= "SNVs / Mbp")
ggsave("L7_MAG_00043_#SNV.png", limitsize = F, dpi = 500)

ggarrange(L7_MAG_00043_heat, L7_MAG_00043_snv,
          ncol=1, nrow=2,
          align="hv")
ggsave(filename = "L7_MAG_00043.png", limitsize=F, width=16, height=18)

#for L8_MAG_00011
L8_MAG_00011 <- L8_MAG_00011 %>% mutate(new_name = case_when(timepoint%>%substr(1,2) == "K1" ~ "Control A", timepoint%>%substr(1,2) == "I4" ~ "Control B",
                                                             timepoint%>%substr(1,2) == "L3" ~ "Control C", timepoint%>%substr(1,2) == "L4" ~ "Control D",
                                                             timepoint%>%substr(1,2) == "I8" ~ "Control E", timepoint%>%substr(1,2) == "L2" ~ "GBH A",
                                                             timepoint%>%substr(1,2)== "L6" ~ "GBH B", timepoint%>%substr(1,2) == "L7" ~ "GBH C",
                                                             timepoint%>%substr(1,2) == "L8" ~ "GBH D"))

L8_MAG_00011_index<-tapply(L8_MAG_00011$final_ref_freq,L8_MAG_00011$groups,FUN=median,na.rm=T)
L8_MAG_00011_index_df<-data.frame(groups=names(L8_MAG_00011_index),median=L8_MAG_00011_index)
L8_MAG_00011_median<- right_join(L8_MAG_00011_index_df, L8_MAG_00011, by=c("groups"))

L8_MAG_00011_heat <- ggplot(L8_MAG_00011_median, aes(x = new_name, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(L8_MAG_00011_median$groups)[order(L8_MAG_00011_median$median)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("UBA953 sp. "),"assembled from GBH D"))))+
  theme(axis.title.y = element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+ 
ggsave(filename="L8_MAG_00011_SNV_95_heatmap.png", dpi = 500)

L8_MAG_00011$simple_class<-L8_MAG_00011$class%>%str_sub(-3,-1)

L8_MAG_00011_snv <- ggplot(L8_MAG_00011, aes(x=new_name, y=(((number_divergent/mag_length)*10^6)), fill=simple_class))+
  geom_bar(stat="identity")+ 
  scale_fill_manual(values=c("#E0E0E0","#424242"))+
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  labs(x="Pond",  y= "SNVs / Mbp")
ggsave("L8_MAG_00011_#SNV.png", limitsize = F, dpi = 500)

ggarrange(L8_MAG_00011_heat, L8_MAG_00011_snv,
          ncol=1, nrow=2,
          align="hv")
ggsave(filename = "L8_MAG_00011.png", limitsize=F, width=16, height=18)

#for L8_MAG_00019
L8_MAG_00019 <- L8_MAG_00019 %>% mutate(new_name = case_when(timepoint%>%substr(1,2) == "K1" ~ "Control A", timepoint%>%substr(1,2) == "I4" ~ "Control B",
                                                             timepoint%>%substr(1,2) == "L3" ~ "Control C", timepoint%>%substr(1,2) == "L4" ~ "Control D",
                                                             timepoint%>%substr(1,2) == "I8" ~ "Control E", timepoint%>%substr(1,2) == "L2" ~ "GBH A",
                                                             timepoint%>%substr(1,2)== "L6" ~ "GBH B", timepoint%>%substr(1,2) == "L7" ~ "GBH C",
                                                             timepoint%>%substr(1,2) == "L8" ~ "GBH D"))

L8_MAG_00019_index<-tapply(L8_MAG_00019$final_ref_freq,L8_MAG_00019$groups,FUN=median,na.rm=T)
L8_MAG_00019_index_df<-data.frame(groups=names(L8_MAG_00019_index),median=L8_MAG_00019_index)
L8_MAG_00019_median<- right_join(L8_MAG_00019_index_df, L8_MAG_00019, by=c("groups"))

L8_MAG_00019_heat <- ggplot(L8_MAG_00019_median, aes(x = new_name, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(L8_MAG_00019_median$groups)[order(L8_MAG_00019_median$median)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("UA16 family "),"assembled from GBH D"))))+
  theme(axis.title.y = element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+ 
ggsave(filename="L8_MAG_00019_SNV_95_heatmap.png", dpi = 500)

L8_MAG_00019$simple_class<-L8_MAG_00019$class%>%str_sub(-3,-1)

L8_MAG_00019_snv <- ggplot(L8_MAG_00019, aes(x=new_name, y=(((number_divergent/mag_length)*10^6)), fill=simple_class))+
  geom_bar(stat="identity")+ 
  scale_fill_manual(values=c("#E0E0E0","#424242"))+
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  labs(x="Pond",  y= "SNVs / Mbp")
ggsave("L8_MAG_00019_#SNV.png", limitsize = F, dpi = 500)

ggarrange(L8_MAG_00019_heat, L8_MAG_00019_snv,
          ncol=1, nrow=2,
          align="hv")
ggsave(filename = "L8_MAG_00019.png", limitsize=F, width=16, height=18)

#for L8_MAG_00042
L8_MAG_00042 <- L8_MAG_00042 %>% mutate(new_name = case_when(timepoint%>%substr(1,2) == "K1" ~ "Control A", timepoint%>%substr(1,2) == "I4" ~ "Control B",
                                                             timepoint%>%substr(1,2) == "L3" ~ "Control C", timepoint%>%substr(1,2) == "L4" ~ "Control D",
                                                             timepoint%>%substr(1,2) == "I8" ~ "Control E", timepoint%>%substr(1,2) == "L2" ~ "GBH A",
                                                             timepoint%>%substr(1,2)== "L6" ~ "GBH B", timepoint%>%substr(1,2) == "L7" ~ "GBH C",
                                                             timepoint%>%substr(1,2) == "L8" ~ "GBH D"))

L8_MAG_00042_index<-tapply(L8_MAG_00042$final_ref_freq,L8_MAG_00042$groups,FUN=median,na.rm=T)
L8_MAG_00042_index_df<-data.frame(groups=names(L8_MAG_00042_index),median=L8_MAG_00042_index)
L8_MAG_00042_median<- right_join(L8_MAG_00042_index_df, L8_MAG_00042, by=c("groups"))

L8_MAG_00042_heat <- ggplot(L8_MAG_00042_median, aes(x = new_name, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  scale_y_discrete(limits=(L8_MAG_00042_median$groups)[order(L8_MAG_00042_median$median)]) +
  labs(fill="Reference Freq.", title=(expression(paste(italic("UBA4660 sp. "),"assembled from GBH D"))))+
  theme(axis.title.y = element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), text=element_text(size=15))+ 
ggsave(filename="L8_MAG_00042_SNV_95_heatmap.png", dpi = 500)

L8_MAG_00042$simple_class<-L8_MAG_00042$class%>%str_sub(-3,-1)

L8_MAG_00042_snv <- ggplot(L8_MAG_00042, aes(x=new_name, y=(((number_divergent/mag_length)*10^6)), fill=simple_class))+
  geom_bar(stat="identity")+ 
  scale_fill_manual(values=c("#E0E0E0","#424242"))+
  theme_classic()+
  theme(legend.title=element_blank(),text=element_text(size=15))+
  labs(x="Pond",  y= "SNVs / Mbp")
ggsave("L8_MAG_00042_#SNV.png", limitsize = F, dpi = 500)

ggarrange(L8_MAG_00042_heat, L8_MAG_00042_snv,
          ncol=1, nrow=2,
          align="hv")
ggsave(filename = "L8_MAG_00042.png", limitsize=F, width=16, height=18)

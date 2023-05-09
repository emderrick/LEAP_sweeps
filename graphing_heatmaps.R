library(tidyverse)
library(dplyr)
library(ggplot2)
library(viridis)
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

I4_MAG_00006_index<-tapply(I4_MAG_00006$final_ref_freq,I4_MAG_00006$groups,FUN=median,na.rm=T)
I4_MAG_00006_index_df<-data.frame(groups=names(I4_MAG_00006_index),median=I4_MAG_00006_index)
I4_MAG_00006_median<- right_join(I4_MAG_00006_index_df, I4_MAG_00006, by=c("groups"))

I4_MAG_00006_heat <- ggplot(I4_MAG_00006_median, aes(x = timepoint, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  ggtitle("I4_MAG_00006") +
  scale_y_discrete(limits=(I4_MAG_00006_median$groups)[order(I4_MAG_00006_median$median)]) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), axis.title.x = element_blank(), legend.title = element_blank(), plot.title = element_text(hjust = 0.5)) 
ggsave(filename="I4_MAG_00006_SNV_95_heatmap.png", dpi = 500)

#for I4_MAG_00065

I4_MAG_00065_index<-tapply(I4_MAG_00065$final_ref_freq,I4_MAG_00065$groups,FUN=median,na.rm=T)
I4_MAG_00065_index_df<-data.frame(groups=names(I4_MAG_00065_index),median=I4_MAG_00065_index)
I4_MAG_00065_median<- right_join(I4_MAG_00065_index_df, I4_MAG_00065, by=c("groups"))

I4_MAG_00065_heat <-ggplot(I4_MAG_00065_median, aes(x = timepoint, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  ggtitle("I4_MAG_00065") +
  scale_y_discrete(limits=(I4_MAG_00065_median$groups)[order(I4_MAG_00065_median$median)]) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), axis.title.x = element_blank(), legend.title = element_blank(), plot.title = element_text(hjust = 0.5))
ggsave(filename="I4_MAG_00065_SNV_95_heatmap.png", dpi = 500)

#for L2_MAG_00048

L2_MAG_00048_index<-tapply(L2_MAG_00048$final_ref_freq,L2_MAG_00048$groups,FUN=median,na.rm=T)
L2_MAG_00048_index_df<-data.frame(groups=names(L2_MAG_00048_index),median=L2_MAG_00048_index)
L2_MAG_00048_median<- right_join(L2_MAG_00048_index_df, L2_MAG_00048, by=c("groups"))

L2_MAG_00048_heat <-ggplot(L2_MAG_00048_median, aes(x = timepoint, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  ggtitle("L2_MAG_00048") +
  scale_y_discrete(limits=(L2_MAG_00048_median$groups)[order(L2_MAG_00048_median$median)]) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), axis.title.x = element_blank(), legend.title = element_blank(), plot.title = element_text(hjust = 0.5))
ggsave(filename="L2_MAG_00048_SNV_95_heatmap.png", dpi = 500)

#for L2_MAG_00052

L2_MAG_00052_index<-tapply(L2_MAG_00052$final_ref_freq,L2_MAG_00052$groups,FUN=median,na.rm=T)
L2_MAG_00052_index_df<-data.frame(groups=names(L2_MAG_00052_index),median=L2_MAG_00052_index)
L2_MAG_00052_median<- right_join(L2_MAG_00052_index_df, L2_MAG_00052, by=c("groups"))

L2_MAG_00052_heat <-ggplot(L2_MAG_00052_median, aes(x = timepoint, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  ggtitle("L2_MAG_00052") +
  scale_y_discrete(limits=(L2_MAG_00052_median$groups)[order(L2_MAG_00052_median$median)]) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), axis.title.x = element_blank(), legend.title = element_blank(), plot.title = element_text(hjust = 0.5))
ggsave(filename="L2_MAG_00052_SNV_95_heatmap.png", dpi = 500)

#for L3_MAG_00058

L3_MAG_00058_index<-tapply(L3_MAG_00058$final_ref_freq,L3_MAG_00058$groups,FUN=median,na.rm=T)
L3_MAG_00058_index_df<-data.frame(groups=names(L3_MAG_00058_index),median=L3_MAG_00058_index)
L3_MAG_00058_median<- right_join(L3_MAG_00058_index_df, L3_MAG_00058, by=c("groups"))

L3_MAG_00058_heat <-ggplot(L3_MAG_00058_median, aes(x = timepoint, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  ggtitle("L3_MAG_00058") +
  scale_y_discrete(limits=(L3_MAG_00058_median$groups)[order(L3_MAG_00058_median$median)]) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), axis.title.x = element_blank(), legend.title = element_blank(), plot.title = element_text(hjust = 0.5))
ggsave(filename="L3_MAG_00058_SNV_95_heatmap.png", dpi = 500)

#for L4_MAG_00099

L4_MAG_00099_index<-tapply(L4_MAG_00099$final_ref_freq,L4_MAG_00099$groups,FUN=median,na.rm=T)
L4_MAG_00099_index_df<-data.frame(groups=names(L4_MAG_00099_index),median=L4_MAG_00099_index)
L4_MAG_00099_median<- right_join(L4_MAG_00099_index_df, L4_MAG_00099, by=c("groups"))

L4_MAG_00099_heat <- ggplot(L4_MAG_00099_median, aes(x = timepoint, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  ggtitle("L4_MAG_00099") +
  scale_y_discrete(limits=(L4_MAG_00099_median$groups)[order(L4_MAG_00099_median$median)]) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), axis.title.x = element_blank(), legend.title = element_blank(), plot.title = element_text(hjust = 0.5))
ggsave(filename="L4_MAG_00099_SNV_95_heatmap.png", dpi = 500)

#for L7_MAG_00020

L7_MAG_00020_index<-tapply(L7_MAG_00020$final_ref_freq,L7_MAG_00020$groups,FUN=median,na.rm=T)
L7_MAG_00020_index_df<-data.frame(groups=names(L7_MAG_00020_index),median=L7_MAG_00020_index)
L7_MAG_00020_median<- right_join(L7_MAG_00020_index_df, L7_MAG_00020, by=c("groups"))

L7_MAG_00020_heat <- ggplot(L7_MAG_00020_median, aes(x = timepoint, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  ggtitle("L7_MAG_00020") +
  scale_y_discrete(limits=(L7_MAG_00020_median$groups)[order(L7_MAG_00020_median$median)]) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), axis.title.x = element_blank(), legend.title = element_blank(), plot.title = element_text(hjust = 0.5))
ggsave(filename="L7_MAG_00020_SNV_95_heatmap.png", dpi = 500)

#for L7_MAG_00028

L7_MAG_00028_index<-tapply(L7_MAG_00028$final_ref_freq,L7_MAG_00028$groups,FUN=median,na.rm=T)
L7_MAG_00028_index_df<-data.frame(groups=names(L7_MAG_00028_index),median=L7_MAG_00028_index)
L7_MAG_00028_median<- right_join(L7_MAG_00028_index_df, L7_MAG_00028, by=c("groups"))

L7_MAG_00028_heat <- ggplot(L7_MAG_00028_median, aes(x = timepoint, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  ggtitle("L7_MAG_00028") +
  scale_y_discrete(limits=(L7_MAG_00028_median$groups)[order(L7_MAG_00028_median$median)]) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), axis.title.x = element_blank(), legend.title = element_blank(), plot.title = element_text(hjust = 0.5))
ggsave(filename="L7_MAG_00028_SNV_95_heatmap.png", dpi = 500)

#for L7_MAG_00043

L7_MAG_00043_index<-tapply(L7_MAG_00043 $final_ref_freq,L7_MAG_00043$groups,FUN=median,na.rm=T)
L7_MAG_00043_index_df<-data.frame(groups=names(L7_MAG_00043_index),median=L7_MAG_00043_index)
L7_MAG_00043_median<- right_join(L7_MAG_00043_index_df, L7_MAG_00043 , by=c("groups"))

L7_MAG_00043_heat <- ggplot(L7_MAG_00043_median, aes(x = timepoint, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  ggtitle("L7_MAG_00043") +
  scale_y_discrete(limits=(L7_MAG_00043_median$groups)[order(L7_MAG_00043_median$median)]) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), axis.title.x = element_blank(), legend.title = element_blank(), plot.title = element_text(hjust = 0.5))
ggsave(filename="L7_MAG_00043_SNV_95_heatmap.png", dpi = 500)

#for L8_MAG_00011

L8_MAG_00011_index<-tapply(L8_MAG_00011$final_ref_freq,L8_MAG_00011$groups,FUN=median,na.rm=T)
L8_MAG_00011_index_df<-data.frame(groups=names(L8_MAG_00011_index),median=L8_MAG_00011_index)
L8_MAG_00011_median<- right_join(L8_MAG_00011_index_df, L8_MAG_00011, by=c("groups"))

L8_MAG_00011_heat <- ggplot(L8_MAG_00011_median, aes(x = timepoint, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  ggtitle("L8_MAG_00011") +
  scale_y_discrete(limits=(L8_MAG_00011_median$groups)[order(L8_MAG_00011_median$median)]) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), axis.title.x = element_blank(), legend.title = element_blank(), plot.title = element_text(hjust = 0.5))
ggsave(filename="L8_MAG_00011_SNV_95_heatmap.png", dpi = 500)

#for L8_MAG_00019

L8_MAG_00019_index<-tapply(L8_MAG_00019$final_ref_freq,L8_MAG_00019$groups,FUN=median,na.rm=T)
L8_MAG_00019_index_df<-data.frame(groups=names(L8_MAG_00019_index),median=L8_MAG_00019_index)
L8_MAG_00019_median<- right_join(L8_MAG_00019_index_df, L8_MAG_00019, by=c("groups"))

L8_MAG_00019_heat <-ggplot(L8_MAG_00019_median, aes(x = timepoint, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  ggtitle("L8_MAG_00019") +
  scale_y_discrete(limits=(L8_MAG_00019_median$groups)[order(L8_MAG_00019_median$median)]) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), axis.title.x = element_blank(), legend.title = element_blank(), plot.title = element_text(hjust = 0.5))
ggsave(filename="L8_MAG_00019_SNV_95_heatmap.png", dpi = 500)

#for L8_MAG_00042

L8_MAG_00042_index<-tapply(L8_MAG_00042$final_ref_freq,L8_MAG_00042$groups,FUN=median,na.rm=T)
L8_MAG_00042_index_df<-data.frame(groups=names(L8_MAG_00042_index),median=L8_MAG_00042_index)
L8_MAG_00042_median<- right_join(L8_MAG_00042_index_df, L8_MAG_00042, by=c("groups"))

L8_MAG_00042_heat <- ggplot(L8_MAG_00042_median, aes(x = timepoint, y = groups, fill= final_ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_classic() +
  ggtitle("L8_MAG_00042") +
  scale_y_discrete(limits=(L8_MAG_00042_median$groups)[order(L8_MAG_00042_median$median)]) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), axis.title.x = element_blank(), legend.title = element_blank(), plot.title = element_text(hjust = 0.5))
ggsave(filename="L8_MAG_00042_SNV_95_heatmap.png", dpi = 500)

#make big plot of all the heatmaps
ggarrange(I4_MAG_00006_heat, L2_MAG_00048_heat, L7_MAG_00028_heat,
          L8_MAG_00011_heat, L8_MAG_00019_heat, L8_MAG_00042_heat,
          L3_MAG_00058_heat, I4_MAG_00065_heat, L4_MAG_00099_heat,
          L2_MAG_00052_heat, L7_MAG_00043_heat, L7_MAG_00020_heat, 
          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"),
          ncol = 4, nrow = 3)
ggsave(filename="95_SNV_heat.png", limitsize = FALSE, width=30, height=20)

ggarrange(I4_MAG_00006_heat, L2_MAG_00048_heat, L7_MAG_00028_heat,
          L8_MAG_00011_heat, L8_MAG_00019_heat, L8_MAG_00042_heat,
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 3, nrow = 2)
ggsave(filename="first6_95_SNV_heat.png", limitsize = FALSE, width=12, height=10)

ggarrange(L3_MAG_00058_heat, I4_MAG_00065_heat, L4_MAG_00099_heat,
          L2_MAG_00052_heat, L7_MAG_00043_heat, L7_MAG_00020_heat, 
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 3, nrow = 2)
ggsave(filename="last_96_95_SNV_heat.png", limitsize = FALSE, width=16, height=10)
#graphing SNV bar graphs
ggplot(L7_MAG_00028, aes(x=pond.y, y=((number_SNVs/genome_length)*10^6), fill=treatment))+
  geom_bar(stat="identity") + 
  theme_classic()+
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="Pond",  y= "SNVs / Mbp", title = expression(paste(italic("Burkholderiaceae SYFN01 "),"assembled from glyphosate pond L7")))
ggsave("L7_MAG_00028_95.png", limitsize = F, dpi = 500)

L4_MAG_00099_SNV <- ggplot(L4_MAG_00099, aes(x=pond.y, y=((number_SNVs/genome_length)*10^6), fill=treatment))+
  geom_bar(stat="identity") + 
  theme_classic()+
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="Pond",  y= "SNVs / Mbp", title = expression(paste(italic("Bosea "),"assembled from control pond L4")))
ggsave("L4_MAG_00099_95.png", limitsize = F, dpi = 500)

L7_MAG_00043_SNV <- ggplot(L7_MAG_00043, aes(x=pond.y, y=((number_SNVs/genome_length)*10^6), fill=treatment))+
  geom_bar(stat="identity") + 
  theme_classic()+
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="Pond",  y= "SNVs / Mbp", title = expression(paste(italic("Luteolibacter "),"assembled from glyphosate pond L7")))
ggsave("L7_MAG_00043_95.png", limitsize = F, dpi = 500)

I8_MAG_00005_SNV <-ggplot(I8_MAG_00005, aes(x=pond.y, y=((number_SNVs/genome_length)*10^6), fill=treatment))+
  geom_bar(stat="identity") + 
  theme_classic()+
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="Pond",  y= "SNVs / Mbp", title = expression(paste(italic("Pseudomonas putida "),"assembled from control pond I8")))
ggsave("I8_MAG_00005_95.png", limitsize = F, dpi = 500)

L3_MAG_00058_SNV <-ggplot(L3_MAG_00058, aes(x=pond.y, y=((number_SNVs/genome_length)*10^6), fill=treatment))+
  geom_bar(stat="identity") + 
  theme_classic()+
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="Pond",  y= "SNVs / Mbp", title = expression(paste(italic("Prosthecobacter "),"assembled from control pond L3")))
ggsave("L3_MAG_00099_95.png", limitsize = F, dpi = 500)

L8_MAG_00019_SNV <-ggplot(L8_MAG_00019, aes(x=pond.y, y=((number_SNVs/genome_length)*10^6), fill=treatment))+
  geom_bar(stat="identity") + 
  theme_classic()+
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="Pond",  y= "SNVs / Mbp", title = expression(paste(italic("UA16 "),"assembled from glyphosate pond L8")))
ggsave("L8_MAG_00019_95.png", limitsize = F, dpi = 500)

L8_MAG_00011_SNV <- ggplot(L8_MAG_00011, aes(x=pond.y, y=((number_SNVs/genome_length)*10^6), fill=treatment))+
  geom_bar(stat="identity") + 
  theme_classic()+
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="Pond",  y= "SNVs / Mbp", title = expression(paste(italic("UBA953 "),"assembled from glyphosate pond L8")))
ggsave("L8_MAG_00011_95.png", limitsize = F, dpi = 500)

#make big plot of all the SNV bar graphs
ggarrange(L3_MAG_00058_SNV, I8_MAG_00005_SNV, L4_MAG_00099_SNV, L8_MAG_00019_SNV, L8_MAG_00011_SNV,
          L7_MAG_00043_SNV, L7_MAG_00028_SNV,
          labels = c("A", "B", "C", "D", "E", "F", "G"),
          ncol = 3, nrow = 3)
ggsave(filename="95_SNVs.png", limitsize = FALSE, width=16, height=10)
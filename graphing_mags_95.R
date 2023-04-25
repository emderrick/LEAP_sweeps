library(tidyverse)
library(ggplot2)
library(dplyr)
library(viridis)
library(stringr)
library(ggpubr)

setwd("/Users/Emma/OneDrive - McGill University/Documents/phd docs/Chapter 1/Aim 1A")

#PLOTTING REFERENCE FREQ AT EACH POSITION
#load mag snv info
mags <- read_csv("ANI_95_mag_SNVs.csv")

#filter SNVs within 100bp of beginning and end of scaffold
mags<-filter(mags, coverage >=4)
mags$pos_from_end<-mags$length-mags$position
mags<- filter(mags, position > 100)
mags<- filter(mags, pos_from_end > 100)

#assign value of 1 to everything that isn't an SNS and assign 0 to SNS
mags <- mutate(mags, number_SNVs = ifelse(class == "SNS", 0, 1))
#assign value of 1 to SNS and 0 to everythign else
mags <- mutate(mags, number_SNSs = ifelse(class == "SNS", 1, 0))
mags$number_divergent<- 1
mags$group<- paste(mags$mag, "in pond", mags$pond, "at time", mags$new_time)
mags <- filter(mags, mag != "I8_MAG_00005")
mags <- mags %>% filter(!(mag=="L3_MAG_00058" & new_time=="1"))
mags <- mags %>% filter(!(mag=="L3_MAG_00058" & new_time=="3"))
mags <- mags %>% filter(!(mag=="L4_MAG_00099" & new_time=="3"))
mags <- mags %>% filter(!(mag=="L8_MAG_00011" & new_time=="3"))
mags <- mags %>% filter(!(mag=="L7_MAG_00043" & new_time=="3"))
mags <- mags %>% filter(!(mag=="L7_MAG_00028" & new_time=="3"))
mags <- mags %>% filter(!(mag=="I4_MAG_00065" & new_time=="1"))
mags <- mags %>% filter(!(mag=="I4_MAG_00065" & new_time=="3"))
mags <- mags %>% filter(!(mag=="L7_MAG_00020" & new_time=="1" & pond=="L7"))
mags <- mags %>% filter(!(mag=="L8_MAG_00042" & new_time=="3"))
mags <- mags %>% filter(!(mag=="L2_MAG_00048" & new_time=="2" & pond=="L7"))
mags <- mags %>% filter(!(mag=="L2_MAG_00052" & new_time=="1"))

filtered_mags <- mags
write.csv(filtered_mags, "filtered_SNVs.csv", row.names=F)

#for L7_MAG_00028
L7_MAG_00028 <- filter(mags, mag=="L7_MAG_00028")
L7_MAG_00028$groups<- paste(L7_MAG_00028$scaffold, str_pad(L7_MAG_00028$position, 7, pad = "0"))
L7_MAG_00028 <- complete(L7_MAG_00028, timepoint, groups)
L7_MAG_00028_index<-tapply(L7_MAG_00028$ref_freq,L7_MAG_00028$groups,FUN=median,na.rm=T)
L7_MAG_00028_index_df<-data.frame(groups=names(L7_MAG_00028_index),median=L7_MAG_00028_index)
L7_MAG_00028_median<- right_join(L7_MAG_00028_index_df, L7_MAG_00028, by=c("groups"))

L7_MAG_00028_heat <- ggplot(L7_MAG_00028_median, aes(x = timepoint, y = groups, fill= ref_freq)) +
  geom_tile() +
  scale_fill_viridis(option="magma", na.value = "white") +
  theme_minimal() +
  ggtitle("L7_MAG_00028") +
  scale_y_discrete(limits=(L7_MAG_00028_median$groups)[order(L7_MAG_00028_median$median)]) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()) 
ggsave(filename="L7_MAG_00028_SNV_95_heatmap.png")

#for L8_MAG_00011
L8_MAG_00011 <- filter(mags, mag=="L8_MAG_00011")
L8_MAG_00011$groups<- paste(L8_MAG_00011$scaffold, str_pad(L8_MAG_00011$position, 7, pad = "0"))
L8_MAG_00011 <- complete(L8_MAG_00011, timepoint, groups)
L8_MAG_00011_index<-tapply(L8_MAG_00011$ref_freq,L8_MAG_00011$groups,FUN=median,na.rm=T)
L8_MAG_00011_index_df<-data.frame(groups=names(L8_MAG_00011_index),median=L8_MAG_00011_index)
L8_MAG_00011_median<- right_join(L8_MAG_00011_index_df, L8_MAG_00011, by=c("groups"))

L8_MAG_00011_heat <- ggplot(L8_MAG_00011_median, aes(x = timepoint, y = groups, fill= ref_freq)) +
  geom_tile() +
  scale_fill_viridis(option="magma", na.value = "white") +
  theme_minimal() +
  ggtitle("L8_MAG_00011") +
  scale_y_discrete(limits=(L8_MAG_00011_median$groups)[order(L8_MAG_00011_median$median)]) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()) 
ggsave(filename="L8_MAG_00011_SNV_95_heatmap.png")

#for L8_MAG_00019
L8_MAG_00019 <- filter(mags, mag=="L8_MAG_00019")
L8_MAG_00019$groups<- paste(L8_MAG_00019$scaffold, str_pad(L8_MAG_00019$position, 7, pad = "0"))
L8_MAG_00019 <- complete(L8_MAG_00019, timepoint, groups)
L8_MAG_00019_index<-tapply(L8_MAG_00019$ref_freq,L8_MAG_00019$groups,FUN=median,na.rm=T)
L8_MAG_00019_index_df<-data.frame(groups=names(L8_MAG_00019_index),median=L8_MAG_00019_index)
L8_MAG_00019_median<- right_join(L8_MAG_00019_index_df, L8_MAG_00019, by=c("groups"))

L8_MAG_00019_heat <-ggplot(L8_MAG_00019_median, aes(x = timepoint, y = groups, fill= ref_freq)) +
  geom_tile() +
  scale_fill_viridis(option="magma", na.value = "white") +
  theme_minimal() +
  ggtitle("L8_MAG_00019") +
  scale_y_discrete(limits=(L8_MAG_00019_median$groups)[order(L8_MAG_00019_median$median)]) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()) 
ggsave(filename="L8_MAG_00019_SNV_95_heatmap.png")

#for L3_MAG_00058
L3_MAG_00058 <- filter(mags, mag=="L3_MAG_00058")
L3_MAG_00058$groups<- paste(L3_MAG_00058$scaffold, str_pad(L3_MAG_00058$position, 7, pad = "0"))
L3_MAG_00058 <- complete(L3_MAG_00058, timepoint, groups)
L3_MAG_00058_index<-tapply(L3_MAG_00058$ref_freq,L3_MAG_00058$groups,FUN=median,na.rm=T)
L3_MAG_00058_index_df<-data.frame(groups=names(L3_MAG_00058_index),median=L3_MAG_00058_index)
L3_MAG_00058_median<- right_join(L3_MAG_00058_index_df, L3_MAG_00058, by=c("groups"))

L3_MAG_00058_heat <-ggplot(L3_MAG_00058_median, aes(x = timepoint, y = groups, fill= ref_freq)) +
  geom_tile() +
  scale_fill_viridis(option="magma", na.value = "white") +
  theme_minimal() +
  ggtitle("L3_MAG_00058") +
  scale_y_discrete(limits=(L3_MAG_00058_median$groups)[order(L3_MAG_00058_median$median)]) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()) 
ggsave(filename="L3_MAG_00058_SNV_95_heatmap.png")

#for L7_MAG_00043
L7_MAG_00043  <- filter(mags, mag=="L7_MAG_00043")
L7_MAG_00043$groups<- paste(L7_MAG_00043$scaffold, str_pad(L7_MAG_00043$position, 7, pad = "0"))
L7_MAG_00043  <- complete(L7_MAG_00043 , timepoint, groups)
L7_MAG_00043_index<-tapply(L7_MAG_00043 $ref_freq,L7_MAG_00043$groups,FUN=median,na.rm=T)
L7_MAG_00043_index_df<-data.frame(groups=names(L7_MAG_00043_index),median=L7_MAG_00043_index)
L7_MAG_00043_median<- right_join(L7_MAG_00043_index_df, L7_MAG_00043 , by=c("groups"))

L7_MAG_00043_heat <- ggplot(L7_MAG_00043_median, aes(x = timepoint, y = groups, fill= ref_freq)) +
  geom_tile() +
  scale_fill_viridis(option="magma", na.value = "white") +
  theme_minimal() +
  ggtitle("L7_MAG_00043") +
  scale_y_discrete(limits=(L7_MAG_00043_median$groups)[order(L7_MAG_00043_median$median)]) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()) 
ggsave(filename="L7_MAG_00043_SNV_95_heatmap.png")

#for L4_MAG_00099
L4_MAG_00099 <- filter(mags, mag=="L4_MAG_00099")
L4_MAG_00099$groups<- paste(L4_MAG_00099$scaffold, str_pad(L4_MAG_00099$position, 7, pad = "0"))
L4_MAG_00099 <- complete(L4_MAG_00099, timepoint, groups)
L4_MAG_00099_index<-tapply(L4_MAG_00099$ref_freq,L4_MAG_00099$groups,FUN=median,na.rm=T)
L4_MAG_00099_index_df<-data.frame(groups=names(L4_MAG_00099_index),median=L4_MAG_00099_index)
L4_MAG_00099_median<- right_join(L4_MAG_00099_index_df, L4_MAG_00099, by=c("groups"))

L4_MAG_00099_heat <- ggplot(L4_MAG_00099_median, aes(x = timepoint, y = groups, fill= ref_freq)) +
  geom_tile() +
  scale_fill_viridis(option="magma", na.value = "white") +
  theme_minimal() +
  ggtitle("L4_MAG_00099") +
  scale_y_discrete(limits=(L4_MAG_00099_median$groups)[order(L4_MAG_00099_median$median)]) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()) 
ggsave(filename="L4_MAG_00099_SNV_95_heatmap.png")

#make big plot of all the heatmaps
ggarrange(L3_MAG_00058_heat, L4_MAG_00099_heat, L8_MAG_00019_heat, L8_MAG_00011_heat,
          L7_MAG_00043_heat, L7_MAG_00028_heat,
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 3, nrow = 3)
ggsave(filename="95_SNV_heat.png", limitsize = FALSE, width=16, height=10)


#graphing SNV bar graphs
ggplot(L7_MAG_00028, aes(x=pond.y, y=((number_SNVs/genome_length)*10^6), fill=treatment))+
  geom_bar(stat="identity") + 
  theme_classic()+
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="Pond",  y= "SNVs / Mbp", title = expression(paste(italic("Burkholderiaceae SYFN01 "),"assembled from glyphosate pond L7")))
ggsave("L7_MAG_00028_95.png", limitsize = F)

L4_MAG_00099_SNV <- ggplot(L4_MAG_00099, aes(x=pond.y, y=((number_SNVs/genome_length)*10^6), fill=treatment))+
  geom_bar(stat="identity") + 
  theme_classic()+
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="Pond",  y= "SNVs / Mbp", title = expression(paste(italic("Bosea "),"assembled from control pond L4")))
ggsave("L4_MAG_00099_95.png", limitsize = F)

L7_MAG_00043_SNV <- ggplot(L7_MAG_00043, aes(x=pond.y, y=((number_SNVs/genome_length)*10^6), fill=treatment))+
  geom_bar(stat="identity") + 
  theme_classic()+
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="Pond",  y= "SNVs / Mbp", title = expression(paste(italic("Luteolibacter "),"assembled from glyphosate pond L7")))
ggsave("L7_MAG_00043_95.png", limitsize = F)

I8_MAG_00005_SNV <-ggplot(I8_MAG_00005, aes(x=pond.y, y=((number_SNVs/genome_length)*10^6), fill=treatment))+
  geom_bar(stat="identity") + 
  theme_classic()+
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="Pond",  y= "SNVs / Mbp", title = expression(paste(italic("Pseudomonas putida "),"assembled from control pond I8")))
ggsave("I8_MAG_00005_95.png", limitsize = F)

L3_MAG_00058_SNV <-ggplot(L3_MAG_00058, aes(x=pond.y, y=((number_SNVs/genome_length)*10^6), fill=treatment))+
  geom_bar(stat="identity") + 
  theme_classic()+
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="Pond",  y= "SNVs / Mbp", title = expression(paste(italic("Prosthecobacter "),"assembled from control pond L3")))
ggsave("L3_MAG_00099_95.png", limitsize = F)

L8_MAG_00019_SNV <-ggplot(L8_MAG_00019, aes(x=pond.y, y=((number_SNVs/genome_length)*10^6), fill=treatment))+
  geom_bar(stat="identity") + 
  theme_classic()+
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="Pond",  y= "SNVs / Mbp", title = expression(paste(italic("UA16 "),"assembled from glyphosate pond L8")))
ggsave("L8_MAG_00019_95.png", limitsize = F)

L8_MAG_00011_SNV <- ggplot(L8_MAG_00011, aes(x=pond.y, y=((number_SNVs/genome_length)*10^6), fill=treatment))+
  geom_bar(stat="identity") + 
  theme_classic()+
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="Pond",  y= "SNVs / Mbp", title = expression(paste(italic("UBA953 "),"assembled from glyphosate pond L8")))
ggsave("L8_MAG_00011_95.png", limitsize = F)

#make big plot of all the SNV bar graphs
ggarrange(L3_MAG_00058_SNV, I8_MAG_00005_SNV, L4_MAG_00099_SNV, L8_MAG_00019_SNV, L8_MAG_00011_SNV,
          L7_MAG_00043_SNV, L7_MAG_00028_SNV,
          labels = c("A", "B", "C", "D", "E", "F", "G"),
          ncol = 3, nrow = 3)
ggsave(filename="95_SNVs.png", limitsize = FALSE, width=16, height=10)
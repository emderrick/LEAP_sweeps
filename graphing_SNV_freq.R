library(tidyverse)
library(ggplot2)
library(dplyr)
library(viridis)
library(stringr)


setwd("/Users/Emma/OneDrive - McGill University/Documents/phd docs/Chapter 1/Aim 1A")

#create list of MAGs I'm using
mag_list<-list("L3_MAG_00058", "I8_MAG_00005", "L4_MAG_00099", "L8_MAG_00019", "L8_MAG_00011",
               "L7_MAG_00043", "L7_MAG_00028")

#PLOTTING REFERENCE FREQ AT EACH POSITION
#load mag snv info from output of combine_instrain_files.R
mag_scaf_SNV <- read_csv("mag_scaf_SNV.csv")

#fix timepoints so instead of time 0, 1, 2 it is time 1, 2, 3
mag_scaf_SNV$new_time<-mag_scaf_SNV$time.y+1

#filter to only include mags I'm using
mags <- filter(mag_scaf_SNV, mag.x %in% mag_list)

#filter SNVs within 100bp of beginning and end of scaffold
mags$pos_from_end<-mags$length-mags$position
mags<- filter(mags, position > 100)
mags<- filter(mags, pos_from_end > 100)
mags<-filter(mags, coverage >=4)

#assign value of 1 to everything that isn't an SNS and assign 0 to SNS for SNV bar plot
mags <- mutate(mags, number_SNVs = ifelse(class == "SNS", 0, 1))
mags <- mutate(mags, number_SNSs = ifelse(class == "SNS", 1, 0))
mags$number_divergent<- 1

#add new column for treatment so I can colour bar plots by treatment
mags <- mutate(mags, treatment = ifelse(pond.y == "I4" | pond.y == "K1" | pond.y == "L3" | pond.y == "L4" | pond.y == "I8", "control", "glyphosate"))
mags <- mutate(mags, treatment =ifelse(pond.y == "I8", "phosphorus", treatment))

#CLASS 1 EPSPS EXAMPLES 

#for L7_MAG_00028
L7_MAG_00028 <- filter(mags, mag.x=="L7_MAG_00028", new_time=="2")
L7_MAG_00028$groups<- paste(L7_MAG_00028$scaffold, str_pad(L7_MAG_00028$position, 7, pad = "0"))
#use complete function to fill in NA values for places where SNV is not called
L7_MAG_00028 <- complete(L7_MAG_00028, pond.y, groups)
#calculate median ref freq
L7_MAG_00028_index<-tapply(L7_MAG_00028$ref_freq,L7_MAG_00028$groups,FUN=median,na.rm=T)
#sort by median 
L7_MAG_00028_index_df<-data.frame(groups=names(L7_MAG_00028_index),median=L7_MAG_00028_index)
L7_MAG_00028_median<- right_join(L7_MAG_00028_index_df, L7_MAG_00028, by=c("groups"))

#SNV heatmap
ggplot(L7_MAG_00028_median, aes(x = pond.y, y = groups, fill= ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = "white") +
  theme_minimal() +
  ggtitle("L7_MAG_00028") +
  scale_y_discrete(limits=(L7_MAG_00028_median$groups)[order(L7_MAG_00028_median$median)]) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()) 
ggsave(filename="L7_MAG_00028_SNV_heatmap.png")

#SNV barplot
ggplot(L7_MAG_00028, aes(x=pond.y, y=((number_SNVs/genome_length)*10^6), fill=treatment))+
  geom_bar(stat="identity") + 
  theme_classic()+
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="Pond",  y= "SNVs / Mbp", title = expression(paste(italic("Burkholderiaceae SYFN01 "),"assembled from glyphosate pond L7")))
ggsave("L7_MAG_00028.png", limitsize = F)

#for L8_MAG_00011
L8_MAG_00011 <- filter(mags, mag.x=="L8_MAG_00011", new_time=="2")
L8_MAG_00011$groups<- paste(L8_MAG_00011$scaffold, str_pad(L8_MAG_00011$position, 7, pad = "0"))
L8_MAG_00011 <- complete(L8_MAG_00011, pond.y, groups)
L8_MAG_00011_index<-tapply(L8_MAG_00011$ref_freq,L8_MAG_00011$groups,FUN=median,na.rm=T)
L8_MAG_00011_index_df<-data.frame(groups=names(L8_MAG_00011_index),median=L8_MAG_00011_index)
L8_MAG_00011_median<- right_join(L8_MAG_00011_index_df, L8_MAG_00011, by=c("groups"))

ggplot(L8_MAG_00011_median, aes(x = pond.y, y = groups, fill= ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value ="white") +
  theme_minimal() +
  ggtitle("L8_MAG_00011") +
  scale_y_discrete(limits=(L8_MAG_00011_median$groups)[order(L8_MAG_00011_median$median)]) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()) 
ggsave(filename="L8_MAG_00011_SNV_heatmap.png")

ggplot(L8_MAG_00011, aes(x=pond.y, y=((number_SNVs/genome_length)*10^6), fill=treatment))+
  geom_bar(stat="identity") + 
  theme_classic()+
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="Pond",  y= "SNVs / Mbp", title = expression(paste(italic("UBA953 "),"assembled from glyphosate pond L8")))
ggsave("L8_MAG_00011.png", limitsize = F)

#for L8_MAG_00019
L8_MAG_00019 <- filter(mags, mag.x=="L8_MAG_00019", new_time=="2")
L8_MAG_00019$groups<- paste(L8_MAG_00019$scaffold, str_pad(L8_MAG_00019$position, 7, pad = "0"))
L8_MAG_00019 <- complete(L8_MAG_00019, pond.y, groups)
L8_MAG_00019_index<-tapply(L8_MAG_00019$ref_freq,L8_MAG_00019$groups,FUN=median,na.rm=T)
L8_MAG_00019_index_df<-data.frame(groups=names(L8_MAG_00019_index),median=L8_MAG_00019_index)
L8_MAG_00019_median<- right_join(L8_MAG_00019_index_df, L8_MAG_00019, by=c("groups"))

ggplot(L8_MAG_00019_median, aes(x = pond.y, y = groups, fill= ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = 'white') +
  theme_minimal() +
  ggtitle("L8_MAG_00019") +
  scale_y_discrete(limits=(L8_MAG_00019_median$groups)[order(L8_MAG_00019_median$median)]) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()) 
ggsave(filename="L8_MAG_00019_SNV_heatmap.png")

ggplot(L8_MAG_00019, aes(x=pond.y, y=((number_SNVs/genome_length)*10^6), fill=treatment))+
  geom_bar(stat="identity") + 
  theme_classic()+
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="Pond",  y= "SNVs / Mbp", title = expression(paste(italic("UA16 "),"assembled from glyphosate pond L8")))
ggsave("L8_MAG_00019.png", limitsize = F)

#CLASS 2 EXAMPLES

#for L3_MAG_00058
L3_MAG_00058 <- filter(mags, mag.x=="L3_MAG_00058", new_time=="2")
L3_MAG_00058$groups<- paste(L3_MAG_00058$scaffold, str_pad(L3_MAG_00058$position, 7, pad = "0"))
L3_MAG_00058 <- complete(L3_MAG_00058, pond.y, groups)
L3_MAG_00058_index<-tapply(L3_MAG_00058$ref_freq,L3_MAG_00058$groups,FUN=median,na.rm=T)
L3_MAG_00058_index_df<-data.frame(groups=names(L3_MAG_00058_index),median=L3_MAG_00058_index)
L3_MAG_00058_median<- right_join(L3_MAG_00058_index_df, L3_MAG_00058, by=c("groups"))

ggplot(L3_MAG_00058_median, aes(x = pond.y, y = groups, fill= ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = 'white') +
  theme_minimal() +
  ggtitle("L3_MAG_00058") +
  scale_y_discrete(limits=(L3_MAG_00058_median$groups)[order(L3_MAG_00058_median$median)]) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()) 
ggsave(filename="L3_MAG_00058_SNV_heatmap.png")

ggplot(L3_MAG_00058, aes(x=pond.y, y=((number_SNVs/genome_length)*10^6), fill=treatment))+
  geom_bar(stat="identity") + 
  theme_classic()+
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="Pond",  y= "SNVs / Mbp", title = expression(paste(italic("Prosthecobacter "),"assembled from control pond L3")))
ggsave("L3_MAG_00099.png", limitsize = F)

#for I8_MAG_00005 probably not using this example
I8_MAG_00005 <- filter(mags, mag.x=="I8_MAG_00005", new_time=="3")
I8_MAG_00005$groups<- paste(I8_MAG_00005$scaffold, str_pad(I8_MAG_00005$position, 7, pad = "0"))
I8_MAG_00005 <- complete(I8_MAG_00005, pond.y, groups)
I8_MAG_00005_index<-tapply(I8_MAG_00005$ref_freq,I8_MAG_00005$groups,FUN=median,na.rm=T)
I8_MAG_00005_index_df<-data.frame(groups=names(I8_MAG_00005_index),median=I8_MAG_00005_index)
I8_MAG_00005_median<- right_join(I8_MAG_00005_index_df, I8_MAG_00005, by=c("groups"))

ggplot(I8_MAG_00005_median, aes(x = pond.y, y = groups, fill= ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = 'white') +
  theme_minimal() +
  ggtitle("I8_MAG_00005") +
  scale_y_discrete(limits=(I8_MAG_00005_median$groups)[order(I8_MAG_00005_median$median)]) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()) 
ggsave(filename="I8_MAG_00005_SNV_heatmap.png")

ggplot(I8_MAG_00005, aes(x=pond.y, y=((number_SNVs/genome_length)*10^6), fill=treatment))+
  geom_bar(stat="identity") + 
  theme_classic()+
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="Pond",  y= "SNVs / Mbp", title = expression(paste(italic("Pseudomonas putida "),"assembled from control pond I8")))
ggsave("I8_MAG_00005.png", limitsize = F)

#for L7_MAG_00043
L7_MAG_00043  <- filter(mags, mag.y=="L7_MAG_00043",new_time=="2")
L7_MAG_00043 $groups<- paste(L7_MAG_00043$scaffold, str_pad(L7_MAG_00043$position, 7, pad = "0"))
L7_MAG_00043  <- complete(L7_MAG_00043 , pond.y, groups)
L7_MAG_00043_index<-tapply(L7_MAG_00043 $ref_freq,L7_MAG_00043$groups,FUN=median,na.rm=T)
L7_MAG_00043_index_df<-data.frame(groups=names(L7_MAG_00043_index),median=L7_MAG_00043_index)
L7_MAG_00043_median<- right_join(L7_MAG_00043_index_df, L7_MAG_00043 , by=c("groups"))

ggplot(L7_MAG_00043_median, aes(x = pond.y, y = groups, fill= ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = 'white') +
  theme_minimal() +
  ggtitle("L7_MAG_00043") +
  scale_y_discrete(limits=(L7_MAG_00043_median$groups)[order(L7_MAG_00043_median$median)]) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()) 
ggsave(filename="L7_MAG_00043_SNV_heatmap.png")

ggplot(L7_MAG_00043, aes(x=pond.y, y=((number_SNVs/genome_length)*10^6), fill=treatment))+
  geom_bar(stat="identity") + 
  theme_classic()+
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="Pond",  y= "SNVs / Mbp", title = expression(paste(italic("Luteolibacter "),"assembled from glyphosate pond L7")))
ggsave("L7_MAG_00043.png", limitsize = F)

#UNCLASSIFIED EPSPS and probably not using this example
#for L4_MAG_00099
L4_MAG_00099 <- filter(mags, mag.x=="L4_MAG_00099", new_time=="2")
L4_MAG_00099$groups<- paste(L4_MAG_00099$scaffold, str_pad(L4_MAG_00099$position, 7, pad = "0"))
L4_MAG_00099 <- complete(L4_MAG_00099, pond.y, groups)
L4_MAG_00099_index<-tapply(L4_MAG_00099$ref_freq,L4_MAG_00099$groups,FUN=median,na.rm=T)
L4_MAG_00099_index_df<-data.frame(groups=names(L4_MAG_00099_index),median=L4_MAG_00099_index)
L4_MAG_00099_median<- right_join(L4_MAG_00099_index_df, L4_MAG_00099, by=c("groups"))

ggplot(L4_MAG_00099_median, aes(x = pond.y, y = groups, fill= ref_freq)) +
  geom_tile() +
  scale_fill_viridis(direction=-1, na.value = 'white') +
  theme_minimal() +
  ggtitle("L4_MAG_00099") +
  scale_y_discrete(limits=(L4_MAG_00099_median$groups)[order(L4_MAG_00099_median$median)]) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()) 
ggsave(filename="old_L4_MAG_00099_SNV_heatmap.png")

ggplot(L4_MAG_00099, aes(x=pond.y, y=((number_SNVs/genome_length)*10^6), fill=treatment))+
  geom_bar(stat="identity") + 
  theme_classic()+
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="Pond",  y= "SNVs / Mbp", title = expression(paste(italic("Bosea "),"assembled from control pond L4")))
ggsave("L4_MAG_00099.png", limitsize = F)

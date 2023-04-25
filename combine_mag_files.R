library(tidyverse)
library(ggplot2)
library(dplyr)
library(viridis)
library(stringr)

setwd("/Users/Emma/OneDrive - McGill University/Documents/phd docs/Chapter 1/Aim 1A")

mag_list<-list("L3_MAG_00058", "I8_MAG_00005", "L4_MAG_00099", "L8_MAG_00019", "L8_MAG_00011",
               "L7_MAG_00043", "L7_MAG_00028", "I4_MAG_00006", "I4_MAG_00065", "L2_MAG_00052",
               "L7_MAG_00020", "L8_MAG_00042", "L2_MAG_00048")

#list genome files from instrain output in directory
genome_files<-list.files("95_profiles/",recursive = T, pattern=".*genome_info.tsv",full.names = T)

#create an empty dataframe
all_mags<-data.frame()

#add genome files into dataframe and add column that is the name of the file
for(i in 1:length(genome_files)){
  pond_time_mags<-read.table(genome_files[i],sep="\t",header=T)
  timepoint<-gsub(".*profile_output/", "", genome_files[i]) %>% substr(1,9)
  pond_time_mags<-cbind(pond_time_mags,timepoint=rep(timepoint,nrow(pond_time_mags)))
  all_mags<-rbind(all_mags,pond_time_mags)
}

#adding new columns with info I need to big file
all_mags$pond<-all_mags$timepoint%>%substr(1,2)
all_mags$time<-all_mags$timepoint%>%substr(9,9)
#fix timepoints so instead of 0, 1, 2 it is 1, 2, 3
all_mags$time<- as.numeric(all_mags$time)
all_mags$new_time<-all_mags$time +1
#add treatment
all_mags <- mutate(all_mags, treatment = ifelse(pond == "I4" | pond == "K1" | pond == "L3" | pond == "L4" | pond == "I8", "control", "glyphosate"))
all_mags <- mutate(all_mags, treatment =ifelse(pond == "I8", "phosphorus", treatment))
all_mags <- all_mags %>% rename("mag_length"="length")
all_mags <- all_mags %>% rename("mag"="genome")
all_mags <- all_mags %>% rename("mag_coverage"="coverage")
write.csv(all_mags, "ANI_95_all_mags.csv", row.names = F)

#will use this file for graphing summary figures
summary_mags <- filter(all_mags, mag %in% mag_list)
write.csv(summary_mags, "ANI_95_mags.csv", row.names = F)

scaffold_files<-list.files("95_profiles/",recursive = T, pattern=".*scaffold_info.tsv",full.names = T)
all_scaffolds<-data.frame()

for(i in 1:length(scaffold_files)){
  pond_time_scaffolds<-read.table(scaffold_files[i],sep="\t",header=T)
  timepoint<-gsub(".*profile_output/", "", scaffold_files[i]) %>% substr(1,9)
  pond_time_scaffolds<-cbind(pond_time_scaffolds,timepoint=rep(timepoint,nrow(pond_time_scaffolds)))
  all_scaffolds<-rbind(all_scaffolds,pond_time_scaffolds)
} 

all_scaffolds$mag<- all_scaffolds$scaffold%>%substr(1,12)

#Merge the scaffold dataframe with the dataframe containing the MAG genome information
#all_mags = y all_scaffolds = x
mag_scaffolds <-all_scaffolds %>% left_join(all_mags, by=c("mag", "timepoint"))
write.csv(mag_scaffolds, "ANI_95_all_scaffolds.csv", row.names = F)

SNV_files<-list.files("95_profiles/",recursive = T, pattern=".*SNVs.tsv",full.names = T)
all_SNVs<-data.frame()

for(i in 1:length(SNV_files)){
  pond_time_SNV<-read.table(SNV_files[i],sep="\t",header=T)
  timepoint<-gsub(".*profile_output/", "", SNV_files[i]) %>% substr(1,9)
  pond_time_SNV<-cbind(pond_time_SNV,timepoint=rep(timepoint,nrow(pond_time_SNV)))
  all_SNVs<-rbind(all_SNVs,pond_time_SNV)
}  

#merge SNV file with scaffold file 
mag_scaf_SNV<- right_join(all_SNVs, mag_scaffolds, by=c("timepoint","scaffold"))
#add a column that says the MAG in what pond
mag_scaf_SNV$group<- paste(mag_scaf_SNV$mag, "in pond", mag_scaf_SNV$pond)

#save so I don't have to keep re-running this
write.csv(mag_scaf_SNV, "ANI_95_all_SNVs.csv", row.names = F)

#filter to only include candidate mags
mags <- filter(mag_scaf_SNV, mag %in% mag_list)

#save version with just MAGs of interest
write.csv(mags, "ANI_95_mag_SNVs.csv", row.names = F)


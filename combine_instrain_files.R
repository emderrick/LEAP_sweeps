library(tidyverse)
library(ggplot2)
library(dplyr)
library(viridis)
library(stringr)

setwd("/Users/Emma/OneDrive - McGill University/Documents/phd docs/Chapter 1/Aim 1A")

#list genome files from instrain output in directory
genome_files<-list.files("instrain_ANI99/",recursive = T, pattern=".*genome_info.tsv",full.names = T)

#create an empty dataframe
all_mags<-data.frame()

#add genome files into dataframe and add column that is the name of the file
for(i in 1:length(genome_files)){
  pond_time_mags<-read.table(genome_files[i],sep="\t",header=T)
  timepoint<-gsub(".*profile_output/", "", genome_files[i]) %>% substr(1,9)
  pond_time_mags<-cbind(pond_time_mags,timepoint=rep(timepoint,nrow(pond_time_mags)))
  all_mags<-rbind(all_mags,pond_time_mags)
}

#adding new columns with info I need
all_mags$genome_length<- all_mags$length
all_mags$mag<- all_mags$genome

#make a new dataframe with only 3 columns I need to add to the scaffold file
mag_length <- all_mags[, c("mag", "genome_length", "timepoint", "coverage")]
mag_length <- mag_length %>% rename("mag_coverage"="coverage")

scaffold_files<-list.files("instrain_ANI99/",recursive = T, pattern=".*scaffold_info.tsv",full.names = T)
all_scaffolds<-data.frame()

for(i in 1:length(scaffold_files)){
  pond_time_scaffolds<-read.table(scaffold_files[i],sep="\t",header=T)
  timepoint<-gsub(".*profile_output/", "", scaffold_files[i]) %>% substr(1,9)
  pond_time_scaffolds<-cbind(pond_time_scaffolds,timepoint=rep(timepoint,nrow(pond_time_scaffolds)))
  all_scaffolds<-rbind(all_scaffolds,pond_time_scaffolds)
} 

all_scaffolds$mag<-all_scaffolds$scaffold%>%substr(1,12)
all_scaffolds$pond<-all_scaffolds$timepoint%>%substr(1,2)
all_scaffolds$time<-all_scaffolds$timepoint%>%substr(9,9)

#Merge the scaffold dataframe with the dataframe containing the MAG genome information
mag_scaffolds <-all_scaffolds %>% left_join(mag_length, by=c("mag", "timepoint"))
write.csv(mag_scaffolds, "mag_scaffolds.csv", row.names = F)

SNV_files<-list.files("instrain_ANI99/",recursive = T, pattern=".*SNVs.tsv",full.names = T)
all_SNVs<-data.frame()

for(i in 1:length(SNV_files)){
  pond_time_SNV<-read.table(SNV_files[i],sep="\t",header=T)
  timepoint<-gsub(".*profile_output/", "", SNV_files[i]) %>% substr(1,9)
  pond_time_SNV<-cbind(pond_time_SNV,timepoint=rep(timepoint,nrow(pond_time_SNV)))
  all_SNVs<-rbind(all_SNVs,pond_time_SNV)
}  

all_SNVs$mag<-all_SNVs$scaffold%>%substr(1,12)
all_SNVs$pond<-all_SNVs$timepoint%>%substr(1,2)
all_SNVs$time<-all_SNVs$timepoint%>%substr(9,9)

#merge SNV file with scaffold file 
mag_scaf_SNV<- all_SNVs %>% right_join(mag_scaffolds, by=c("timepoint","scaffold"))

#save so I don't have to keep re-running this
write.csv(mag_scaf_SNV, "mag_scaf_SNV.csv", row.names = F)

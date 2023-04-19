library(tidyverse)
library(ggplot2)
library(dplyr)
library(viridis)
library(stringr)

setwd("/Users/Emma/OneDrive - McGill University/Documents/phd docs/Chapter 1/Aim 1A")

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

#adding new columns with info I need
all_mags$genome_length<- all_mags$length
all_mags$mag<- all_mags$genome

#make a new dataframe with only 3 columns I need to add to the scaffold file
mag_length <- all_mags[, c("mag", "genome_length", "timepoint", "coverage")]
mag_length <- mag_length %>% rename("mag_coverage"="coverage")

scaffold_files<-list.files("95_profiles/",recursive = T, pattern=".*scaffold_info.tsv",full.names = T)
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
write.csv(mag_scaffolds, "ANI_95_mag_scaffolds.csv", row.names = F)

SNV_files<-list.files("95_profiles/",recursive = T, pattern=".*SNVs.tsv",full.names = T)
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

#add a column that says the MAG in what pond
mag_scaf_SNV$group<- paste(mag_scaf_SNV$mag.x, "in pond", mag_scaf_SNV$pond.x)

#fix timepoints so instead of 0, 1, 2 it is 1, 2, 3
mag_scaf_SNV$time.y <- as.numeric(mag_scaf_SNV$time.y)
mag_scaf_SNV$new_time<-mag_scaf_SNV$time.y +1

#list of mags with potential
mag_list<-list("L3_MAG_00058", "I8_MAG_00005", "L4_MAG_00099", "L8_MAG_00019", "L8_MAG_00011",
               "L7_MAG_00043", "L7_MAG_00028")

#filter to only include candidate mags
mags <- filter(mag_scaf_SNV, mag.x %in% mag_list)

#filter SNVs within 100bp of beginning and end of scaffold
mags$pos_from_end<-mags$length-mags$position
mags<- filter(mags, position > 100)
mags<- filter(mags, pos_from_end > 100)
mags<-filter(mags, coverage >=4)

#assign value of 1 to everything that isn't an SNS and assign 0 to SNS
mags <- mutate(mags, number_SNVs = ifelse(class == "SNS", 0, 1))
#assign value of 1 to SNS and 0 to everythign else
mags <- mutate(mags, number_SNSs = ifelse(class == "SNS", 1, 0))
mags$number_divergent<- 1

#add new column for treatment
mags <- mutate(mags, treatment = ifelse(pond.y == "I4" | pond.y == "K1" | pond.y == "L3" | pond.y == "L4" | pond.y == "I8", "control", "glyphosate"))
mags <- mutate(mags, treatment =ifelse(pond.y == "I8", "phosphorus", treatment))

#save so I don't have to keep re-running this
write.csv(mags, "ANI_95_mags.csv", row.names = F)

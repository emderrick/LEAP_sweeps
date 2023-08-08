library(tidyverse)
library(ggplot2)
library(dplyr)
library(viridis)
library(stringr)

mag_list<-list("L3_MAG_00058", "L4_MAG_00099", "L8_MAG_00019", "L8_MAG_00011", "L7_MAG_00043", "L7_MAG_00028",
               "I4_MAG_00006", "I4_MAG_00065", "L2_MAG_00052", "L7_MAG_00020", "L8_MAG_00042")

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

all_mags$time<-as.numeric(all_mags$timepoint%>%substr(9,9))
all_mags$new_time<-all_mags$time +1
all_mags$pond<-all_mags$timepoint%>%substr(1,2)

all_mags <- mutate(all_mags, treatment = ifelse(pond == "I4" | pond == "K1" | pond == "L3" | pond == "L4" | pond == "I8", "control", "glyphosate"))
all_mags <- mutate(all_mags, treatment =ifelse(pond == "I8", "phosphorus", treatment))
all_mags <- all_mags %>% rename("mag_length"="length")
all_mags <- all_mags %>% rename("mag"="genome")
all_mags <- all_mags %>% rename("mag_coverage"="coverage")

all_mags <- all_mags %>% mutate(name = case_when(timepoint%>%substr(1,2) == "K1" ~ "Control A", timepoint%>%substr(1,2) == "I4" ~ "Control B",
                                timepoint%>%substr(1,2) == "L3" ~ "Control C", timepoint%>%substr(1,2) == "L4" ~ "Control D",
                                timepoint%>%substr(1,2) == "I8" ~ "Control E", timepoint%>%substr(1,2) == "L2" ~ "GBH A",
                                timepoint%>%substr(1,2)== "L6" ~ "GBH B", timepoint%>%substr(1,2) == "L7" ~ "GBH C",
                                timepoint%>%substr(1,2) == "L8" ~ "GBH D"))
all_mags$new_name <- paste(paste(all_mags$name, "at T"), all_mags$new_time, sep="")

write.csv(all_mags, "ANI_95_all_mags.csv", row.names = F)

#filter out mags and times I don't need
mag_cov <- subset(all_mags, mag %in% mag_list)
mag_cov <- mag_cov %>% subset(!(mag=="L3_MAG_00058" & new_time=="1"))
mag_cov <- mag_cov %>% subset(!(mag=="L3_MAG_00058" & new_time=="3"))
mag_cov <- mag_cov %>% subset(!(mag=="L4_MAG_00099" & new_time=="3"))
mag_cov <- mag_cov %>% subset(!(mag=="L8_MAG_00011" & new_time=="3"))
mag_cov <- mag_cov %>% subset(!(mag=="L7_MAG_00043" & new_time=="3"))
mag_cov <- mag_cov %>% subset(!(mag=="L7_MAG_00028" & new_time=="3"))
mag_cov <- mag_cov %>% subset(!(mag=="I4_MAG_00065" & new_time=="1"))
mag_cov <- mag_cov %>% subset(!(mag=="I4_MAG_00065" & new_time=="3"))
mag_cov <- mag_cov %>% subset(!(mag=="L8_MAG_00042" & new_time=="3"))
mag_cov <- mag_cov %>% subset(!(mag=="L7_MAG_00020" & new_time=="1" & pond=="L7"))
mag_cov <- mag_cov %>% subset(!(mag=="L2_MAG_00052" & new_time=="1"))
mag_cov <- mag_cov[, c(1,2,3,37)]
write.csv(mag_cov, "candidate_mag_coverage.csv", row.names = F)


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
mag_scaf_SNV<- left_join(all_SNVs, mag_scaffolds, by=c("timepoint","scaffold"))
write.csv(mag_scaf_SNV, "ANI_95_all_SNVs.csv", row.names = F)

#filter to only include candidate mags and times I'm interested in
mags <- subset(mag_scaf_SNV, mag %in% mag_list)
write.csv(mags, "ANI_95_mag_SNVs.csv", row.names = F)
mags <- mags %>% subset(!(mag=="L3_MAG_00058" & new_time=="1"))
mags <- mags %>% subset(!(mag=="L3_MAG_00058" & new_time=="3"))
mags <- mags %>% subset(!(mag=="L4_MAG_00099" & new_time=="3"))
mags <- mags %>% subset(!(mag=="L8_MAG_00011" & new_time=="3"))
mags <- mags %>% subset(!(mag=="L7_MAG_00043" & new_time=="3"))
mags <- mags %>% subset(!(mag=="L7_MAG_00028" & new_time=="3"))
mags <- mags %>% subset(!(mag=="I4_MAG_00065" & new_time=="1"))
mags <- mags %>% subset(!(mag=="I4_MAG_00065" & new_time=="3"))
mags <- mags %>% subset(!(mag=="L8_MAG_00042" & new_time=="3"))
mags <- mags %>% subset(!(mag=="L2_MAG_00052" & new_time=="1"))
mags <- mags %>% subset(!(mag=="L7_MAG_00020" & new_time=="1" & pond=="L7"))

#filter SNVs within 100bp of beginning and end of scaffold
mags$pos_from_end<-mags$length-mags$position
mags<- subset(mags, position > 100)
mags<- subset(mags, pos_from_end > 100)

#assign value of 1 to everything that isn't an SNS and assign 0 to SNS
mags <- mutate(mags, number_SNVs = ifelse(class == "SNS", 0, 1))
#assign value of 1 to SNS and 0 to everythign else
mags <- mutate(mags, number_SNSs = ifelse(class == "SNS", 1, 0))
mags$number_divergent<- 1
mags$full_group <- paste(mags$mag, mags$new_name)
mags <- mags %>% subset(!(position_coverage > mag_coverage*3)) 
mags <- mags %>% subset(!(position_coverage < mag_coverage/3))
mags <- subset(mags, allele_count <=2)
write.csv(mags, "filtered_ANI_95_mag_SNVs.csv", row.names=F)

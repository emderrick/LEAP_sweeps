library(tidyverse)

mag_list<-list("L3_MAG_00058", "L4_MAG_00099", "L8_MAG_00019", "L8_MAG_00011", "L7_MAG_00043", "L7_MAG_00028",
               "I4_MAG_00006", "I4_MAG_00065", "L2_MAG_00052", "L7_MAG_00020", "L8_MAG_00042")

#list genome files from instrain output in directory
genome_files<-list.files("subsampled_instrain/",recursive = T, pattern=".*genome_info.tsv",full.names = T)

#create an empty dataframe
all_mags<-data.frame()

#add genome files into dataframe and add column that is the name of the file
for(i in 1:length(genome_files)){
  pond_time_mags<-read.table(genome_files[i],sep="\t",header=T)
  timepoint<-gsub(".*instrain_output/", "", genome_files[i]) %>% substr(9,17)
  pond_time_mags<-cbind(pond_time_mags,timepoint=rep(timepoint,nrow(pond_time_mags)))
  all_mags<-rbind(all_mags,pond_time_mags)
}

all_mags$time <- all_mags$timepoint %>% substr(9,9)
all_mags$new_time <- as.numeric(all_mags$time) + 1
all_mags$pond <- all_mags$timepoint %>% substr(1,2)

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
write.csv(all_mags, "all_mags_subsamp.csv", row.names = F)

mag_cov <- all_mags[, c(1,2,3,37)]
write.csv(mag_cov, "mag_coverage_subsamp.csv", row.names = F)

scaffold_files<-list.files("subsampled_instrain/",recursive = T, pattern=".*scaffold_info.tsv",full.names = T)
all_scaffolds<-data.frame()

for(i in 1:length(scaffold_files)){
  pond_time_scaffolds<-read.table(scaffold_files[i],sep="\t",header=T)
  timepoint<-gsub(".*instrain_output/", "", scaffold_files[i]) %>% substr(9,17)
  pond_time_scaffolds<-cbind(pond_time_scaffolds,timepoint=rep(timepoint,nrow(pond_time_scaffolds)))
  all_scaffolds<-rbind(all_scaffolds,pond_time_scaffolds)
} 

all_scaffolds$mag<- all_scaffolds$scaffold %>% substr(1,12)

mag_scaffolds <- all_scaffolds %>% left_join(all_mags, by = c("mag", "timepoint"))
write.csv(mag_scaffolds, "all_scaffolds_subsamp.csv", row.names = F)

SNV_files<-list.files("subsampled_instrain/", recursive = T, pattern=".*SNVs.tsv",full.names = T)
all_SNVs<-data.frame()

for(i in 1:length(SNV_files)){
  pond_time_SNV<-read.table(SNV_files[i],sep="\t",header=T)
  timepoint<-gsub(".*instrain_output/", "", SNV_files[i]) %>% substr(9,17)
  pond_time_SNV<-cbind(pond_time_SNV,timepoint=rep(timepoint,nrow(pond_time_SNV)))
  all_SNVs<-rbind(all_SNVs,pond_time_SNV)
}  

mag_scaf_SNV <- left_join(all_SNVs, mag_scaffolds, by=c("timepoint","scaffold"))
write.csv(mag_scaf_SNV, "all_SNVs_subsamp.csv", row.names = F)

mag_scaf_SNV$pos_from_end <- mag_scaf_SNV$length - mag_scaf_SNV$position
mag_scaf_SNV <- subset(mag_scaf_SNV, position > 100)
mag_scaf_SNV <- subset(mag_scaf_SNV, pos_from_end > 100)

mag_scaf_SNV$number_SNVs <- with(mag_scaf_SNV, ifelse(class == "SNV", 1, 0))
mag_scaf_SNV$number_SNSs <- with(mag_scaf_SNV, ifelse(class == "SNS", 1, 0))

mag_scaf_SNV$number_divergent <- 1
mag_scaf_SNV$full_group <- paste(mag_scaf_SNV$mag, mag_scaf_SNV$new_name)
mag_scaf_SNV <- mag_scaf_SNV %>% subset(!(position_coverage > mag_coverage*3)) 
mag_scaf_SNV <- mag_scaf_SNV %>% subset(!(position_coverage < mag_coverage/3))
mag_scaf_SNV <- subset(mag_scaf_SNV, allele_count <= 2)
write.csv(mag_scaf_SNV, "filtered_mag_SNVs_subsamp.csv", row.names=F)

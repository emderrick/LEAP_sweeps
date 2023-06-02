library(tidyverse)
library(dplyr)
library(tidyr)

L3_MAG_00058 <- read_csv("all_L3_MAG_00058_SNVs.csv")
L4_MAG_00099 <- read_csv("all_L4_MAG_00099_SNVs.csv")
L8_MAG_00019 <- read_csv("all_L8_MAG_00019_SNVs.csv")
L8_MAG_00011 <- read_csv("all_L8_MAG_00011_SNVs.csv")
L7_MAG_00043 <- read_csv("all_L7_MAG_00043_SNVs.csv")
L7_MAG_00028 <- read_csv("all_L7_MAG_00028_SNVs.csv")
I4_MAG_00006 <- read_csv("all_I4_MAG_00006_SNVs.csv")
I4_MAG_00065 <- read_csv("all_I4_MAG_00065_SNVs.csv")
L7_MAG_00020 <- read_csv("all_L7_MAG_00020_SNVs.csv")
L8_MAG_00042 <- read_csv("all_L8_MAG_00042_SNVs.csv")
L2_MAG_00052 <- read_csv("all_L2_MAG_00052_SNVs.csv")

all_MAG_SNVs<-rbind(L3_MAG_00058, L4_MAG_00099, L8_MAG_00019, L8_MAG_00011, L7_MAG_00043, L7_MAG_00028, I4_MAG_00006, I4_MAG_00065, L7_MAG_00020, L8_MAG_00042, L2_MAG_00052) 
all_MAG_SNVs <- all_MAG_SNVs %>% mutate(new_name = case_when(timepoint%>%substr(1,2) == "K1" ~ "Control A", timepoint%>%substr(1,2) == "I4" ~ "Control B",
                                                             timepoint%>%substr(1,2) == "L3" ~ "Control C", timepoint%>%substr(1,2) == "L4" ~ "Control D",
                                                             timepoint%>%substr(1,2) == "I8" ~ "Control E", timepoint%>%substr(1,2) == "L2" ~ "GBH A",
                                                             timepoint%>%substr(1,2)== "L6" ~ "GBH B", timepoint%>%substr(1,2) == "L7" ~ "GBH C",
                                                             timepoint%>%substr(1,2) == "L8" ~ "GBH D"))
all_MAG_SNVs$time2<-as.numeric(all_MAG_SNVs$timepoint%>%substr(9,9))+1
all_MAG_SNVs$name<-paste(all_MAG_SNVs$new_name, sep=" at T", all_MAG_SNVs$time2)
all_MAG_SNVs$mag<-all_MAG_SNVs$groups%>%substr(1,12)

small_MAG_SNVs<- all_MAG_SNVs %>% select(c('groups', 'mag', 'name', 'final_ref_freq'))
gene_MAG_SNVs<- all_MAG_SNVs %>% select(c('gene', 'groups'))

gene_dict<-distinct(gene_MAG_SNVs)
gene_dict<-gene_dict[complete.cases(gene_dict),]
small_MAG_SNVs<- left_join(small_MAG_SNVs,gene_dict, by="groups")
small_MAG_SNVs_hor <- spread(small_MAG_SNVs, key=name, value=final_ref_freq)
write.csv(small_MAG_SNVs_hor, "small_MAG_SNVs_hor.csv", row.names = F)

small_MAG_SNVs_hor <- read_csv("small_MAG_SNVs_hor.csv")
small_MAG_SNVs_hor<- small_MAG_SNVs_hor %>% rowwise %>% mutate(all_mean = mean(c_across(where(is.numeric)), na.rm=TRUE))
small_MAG_SNVs_hor<- small_MAG_SNVs_hor %>% rowwise %>% mutate(control_mean = mean(c_across((starts_with("Control"))), na.rm=TRUE))
small_MAG_SNVs_hor<- small_MAG_SNVs_hor %>% rowwise %>% mutate(GBH_mean = mean(c_across((starts_with("GBH"))), na.rm=TRUE))
write.csv(small_MAG_SNVs_hor, "small_MAG_SNVs_hor_means.csv", row.names = F)
small_MAG_SNVs_hor$abs_val<- abs(small_MAG_SNVs_hor$control_mean - small_MAG_SNVs_hor$GBH_mean)
threshold_snvs<- subset(small_MAG_SNVs_hor, abs_val >= 0.5)

#get all gene coordinates
#list genome files from instrain output in directory
gene_files<-list.files("95_profiles/",recursive = T, pattern=".*gene_info.tsv",full.names = T)

#create an empty dataframe
all_genes<-data.frame()

#add genome files into dataframe and add column that is the name of the file
for(i in 1:length(gene_files)){
  pond_time_genes<-read.table(gene_files[i],sep="\t",header=T)
  timepoint<-gsub(".*profile_output/", "", gene_files[i]) %>% substr(1,9)
  pond_time_genes<-cbind(pond_time_genes,timepoint=rep(timepoint,nrow(pond_time_genes)))
  all_genes<-rbind(all_genes,pond_time_genes)
}

all_gene_coord <- subset(all_genes, select=c(gene, start, end))
all_gene_coord <- all_gene_coord %>% distinct(gene, .keep_all = T)
all_gene_coord$new_start<-all_gene_coord$start+1
all_gene_coord$new_end<-all_gene_coord$end+1

only_genes <- subset(threshold_snvs, gene!="NA", select=c(gene))
unique_genes <- distinct(only_genes)
gene_locations<- left_join(unique_genes, all_gene_coord, by=c("gene"))

write_tsv(unique_genes, file="all_genes.tsv")
write_tsv(gene_locations, file="gene_locations.tsv")

#graphing

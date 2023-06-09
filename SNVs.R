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
write.csv(all_MAG_SNVs, "all_MAG_SNVs.csv", row.names=F)
all_MAG_SNVs <- read_csv("all_MAG_SNVs.csv")

small_MAG_SNVs<- all_MAG_SNVs %>% select(c('groups', 'gene', 'mag', 'mag_length', 'name', 'final_ref_freq'))
small_MAG_SNVs <- small_MAG_SNVs %>% group_by(groups) %>% fill(gene, mag_length, .direction="updown")
small_MAG_SNVs_hor <- spread(small_MAG_SNVs, key=name, value=final_ref_freq)
small_MAG_SNVs_hor$all_mean <-rowMeans(small_MAG_SNVs_hor[c(6,7,9,11,12,14,15,16,17)], na.rm=T)
small_MAG_SNVs_hor$control_mean <-rowMeans(small_MAG_SNVs_hor[c(6,7,9,11,12)], na.rm=T)
small_MAG_SNVs_hor$GBH_mean <-rowMeans(small_MAG_SNVs_hor[c(14,15,16,17)], na.rm=T)
small_MAG_SNVs_hor$abs_val<- abs(small_MAG_SNVs_hor$control_mean - small_MAG_SNVs_hor$GBH_mean)
write.csv(small_MAG_SNVs_hor, "small_MAG_SNVs_hor_means.csv", row.names = F)

small_MAG_SNVs_hor <- read_csv("small_MAG_SNVs_hor_means.csv")
threshold_snvs<- subset(small_MAG_SNVs_hor, abs_val >= 0.5)
strict <- threshold_snvs
strict$control <- with(strict, ifelse(control_mean <= 0.5, 'low', 'high'))
strict$control_ref<- with(strict, ifelse(control=="high", apply(strict[c(6,7,9,11,12)], 1, min, na.rm=T), apply(strict[c(6,7,9,11,12)], 1, max, na.rm=T)))
strict$GBH_ref <- with(strict, ifelse(control=="high", apply(strict[c(14,15,16,17)], 1, max, na.rm=T), apply(strict[c(14,15,16,17)], 1, min, na.rm=T)))
only_strict <- strict %>% subset((control=="high" & control_ref > GBH_ref) | (control=="low" & control_ref < GBH_ref))
not_strict <- strict %>% subset((control=="high" & control_ref <= GBH_ref) | (control=="low" & control_ref >= GBH_ref))

#EXTRACT GENE POSITIONS I'M INTERESTED IN
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

#get all gene coordinates
all_gene_coord <- subset(all_genes, select=c(gene, start, end))
all_gene_coord <- all_gene_coord %>% distinct(gene, .keep_all = T)
all_gene_coord$new_start<-all_gene_coord$start+1
all_gene_coord$new_end<-all_gene_coord$end+1
#get gene positions that pass threshold
only_genes <- subset(only_strict, gene!="NA", select=c(gene))
only_genes_sum <- only_genes %>% count(gene, name="snvs_in_gene")
gene_locations<- left_join(only_genes_sum, all_gene_coord, by=c("gene"))
#save file of gene list and position
write_tsv(unique_genes, file="all_genes.tsv")
write_tsv(gene_locations, file="gene_locations.tsv")

#SAVE NEW FILES WITH MEANS
class_MAG_SNVs<-all_MAG_SNVs %>% select(c('groups','name','class'))
class_MAG_SNVs$class<-class_MAG_SNVs$class%>%str_sub(-3,-1)

L3_MAG_00058_SNVs <- subset(small_MAG_SNVs_hor, mag=="L3_MAG_00058") %>% 
  select(c("groups", "mag", "mag_length", "gene", "Control C at T2", "Control D at T2","GBH C at T2","GBH D at T2", "all_mean", "control_mean", "GBH_mean", "abs_val")) %>%
  pivot_longer(cols=contains("at"),  names_to="name", values_to="final_ref_freq", values_drop_na=F) %>%
  left_join(class_MAG_SNVs, by=c("groups", "name"))
write.csv(L3_MAG_00058_SNVs, "L3_MAG_00058_SNVs.csv", row.names = F)

L4_MAG_00099_SNVs <- subset(small_MAG_SNVs_hor, mag=="L4_MAG_00099") %>% 
  select(c("groups", "mag", "mag_length","gene", "Control D at T2", "GBH A at T2", "GBH B at T2","GBH C at T2","GBH D at T2", "all_mean", "control_mean", "GBH_mean", "abs_val")) %>%
  pivot_longer(cols=contains("at"),  names_to="name", values_to="final_ref_freq", values_drop_na=F) %>% 
  left_join(class_MAG_SNVs, by=c("groups", "name"))
write.csv(L4_MAG_00099_SNVs, "L4_MAG_00099_SNVs.csv", row.names = F)

L8_MAG_00019_SNVs <- subset(small_MAG_SNVs_hor, mag=="L8_MAG_00019") %>% 
  select(c("groups", "mag", "mag_length","gene", "Control E at T2", "GBH A at T2","GBH B at T2","GBH D at T2", "all_mean", "control_mean", "GBH_mean", "abs_val")) %>%
  pivot_longer(cols=contains("at"),  names_to="name", values_to="final_ref_freq", values_drop_na=F) %>% 
  left_join(class_MAG_SNVs, by=c("groups", "name"))
write.csv(L8_MAG_00019_SNVs, "L8_MAG_00019_SNVs.csv", row.names = F)

L8_MAG_00011_SNVs <- subset(small_MAG_SNVs_hor, mag=="L8_MAG_00011") %>% 
  select(c("groups", "mag", "mag_length","gene", "Control E at T2", "GBH A at T2","GBH D at T2", "all_mean", "control_mean", "GBH_mean", "abs_val")) %>%
  pivot_longer(cols=contains("at"),  names_to="name", values_to="final_ref_freq", values_drop_na=F) %>% 
  left_join(class_MAG_SNVs, by=c("groups", "name")) 
write.csv(L8_MAG_00011_SNVs, "L8_MAG_00011_SNVs.csv", row.names = F)

L7_MAG_00043_SNVs <- subset(small_MAG_SNVs_hor, mag=="L7_MAG_00043") %>% 
  select(c("groups", "mag", "mag_length","gene", "Control D at T2", "GBH B at T2","GBH C at T2", "all_mean", "control_mean", "GBH_mean", "abs_val")) %>%
  pivot_longer(cols=contains("at"),  names_to="name", values_to="final_ref_freq", values_drop_na=F) %>% 
  left_join(class_MAG_SNVs, by=c("groups", "name"))
write.csv(L7_MAG_00043_SNVs, "L7_MAG_00043_SNVs.csv", row.names = F)

L7_MAG_00028_SNVs <- subset(small_MAG_SNVs_hor, mag=="L7_MAG_00028") %>% 
  select(c("groups", "mag", "mag_length","gene", "Control E at T2", "GBH A at T2","GBH B at T2","GBH C at T2", "all_mean", "control_mean", "GBH_mean", "abs_val")) %>%
  pivot_longer(cols=contains("at"),  names_to="name", values_to="final_ref_freq", values_drop_na=F) %>% 
  left_join(class_MAG_SNVs, by=c("groups", "name"))
write.csv(L7_MAG_00028_SNVs, "L7_MAG_00028_SNVs.csv", row.names = F)

I4_MAG_00006_SNVs <- subset(small_MAG_SNVs_hor, mag=="I4_MAG_00006") %>% 
  select(c("groups", "mag", "mag_length","gene", "Control B at T2", "Control E at T2","GBH A at T2","GBH D at T2", "all_mean", "control_mean", "GBH_mean", "abs_val")) %>%
  pivot_longer(cols=contains("at"),  names_to="name", values_to="final_ref_freq", values_drop_na=F) %>% 
  left_join(class_MAG_SNVs, by=c("groups", "name"))
write.csv(I4_MAG_00006_SNVs, "I4_MAG_00006_SNVs.csv", row.names = F)

I4_MAG_00065_SNVs <- subset(small_MAG_SNVs_hor, mag=="I4_MAG_00065") %>% 
  select(c("groups", "mag", "mag_length","gene", "Control A at T2","Control B at T2","Control C at T2", "Control D at T2","Control E at T2","GBH B at T2", "all_mean", "control_mean", "GBH_mean", "abs_val")) %>%
  pivot_longer(cols=contains("at"),  names_to="name", values_to="final_ref_freq", values_drop_na=F) %>% 
  left_join(class_MAG_SNVs, by=c("groups", "name"))
write.csv(I4_MAG_00065_SNVs, "I4_MAG_00065_SNVs.csv", row.names = F)

L7_MAG_00020_SNVs <- subset(small_MAG_SNVs_hor, mag=="L7_MAG_00020") %>% 
  select(c("groups", "mag", "mag_length","gene", "Control A at T1","Control A at T2","Control C at T1", "Control C at T2","Control D at T1",
           "Control D at T2","GBH A at T1", "GBH A at T2","all_mean", "control_mean", "GBH_mean", "abs_val")) %>%
  pivot_longer(cols=contains("at"),  names_to="name", values_to="final_ref_freq", values_drop_na=F) %>% 
  left_join(class_MAG_SNVs, by=c("groups", "name"))
write.csv(L7_MAG_00020_SNVs, "L7_MAG_00020_SNVs.csv", row.names = F)

L8_MAG_00042_SNVs <- subset(small_MAG_SNVs_hor, mag=="L8_MAG_00042") %>% 
  select(c("groups", "mag", "mag_length","gene", "Control A at T2", "Control C at T2","Control D at T2","Control E at T2","GBH D at T2", "all_mean", "control_mean", "GBH_mean", "abs_val")) %>%
  pivot_longer(cols=contains("at"),  names_to="name", values_to="final_ref_freq", values_drop_na=F) %>% 
  left_join(class_MAG_SNVs, by=c("groups", "name"))
write.csv(L8_MAG_00042_SNVs, "L8_MAG_00042_SNVs.csv", row.names = F)

L2_MAG_00052_SNVs <- subset(small_MAG_SNVs_hor, mag=="L2_MAG_00052") %>% 
  select(c("groups", "mag", "mag_length","gene", "Control A at T2", "Control B at T2","Control D at T2","Control E at T2","GBH A at T2", "all_mean", "control_mean", "GBH_mean", "abs_val")) %>%
  pivot_longer(cols=contains("at"),  names_to="name", values_to="final_ref_freq", values_drop_na=F) %>% 
  left_join(class_MAG_SNVs, by=c("groups", "name"))
write.csv(L2_MAG_00052_SNVs, "L2_MAG_00052_SNVs.csv", row.names = F)


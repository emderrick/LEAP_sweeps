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
write.csv(all_MAG_SNVs, "all_MAG_SNVs.csv", row.names=F)

small_MAG_SNVs<- all_MAG_SNVs %>% select(c('scaffold.x', 'position.x', 'length', 'groups', 'gene', 'mag', 'mag_length', 'name', 'final_ref_freq'))
small_MAG_SNVs <- small_MAG_SNVs %>% group_by(groups) %>% fill(gene, mag_length, length, .direction="updown")
small_MAG_SNVs_hor <- spread(small_MAG_SNVs, key=name, value=final_ref_freq)
small_MAG_SNVs_hor$all_mean <-rowMeans(small_MAG_SNVs_hor[c(9,10,12,14,15,17,18,19,20)], na.rm=T)
small_MAG_SNVs_hor$control_mean <-rowMeans(small_MAG_SNVs_hor[c(9,10,12,14,15)], na.rm=T)
small_MAG_SNVs_hor$GBH_mean <-rowMeans(small_MAG_SNVs_hor[c(17,18,19,20)], na.rm=T)
small_MAG_SNVs_hor$abs_val<- abs(small_MAG_SNVs_hor$control_mean - small_MAG_SNVs_hor$GBH_mean)
write.csv(small_MAG_SNVs_hor, "small_MAG_SNVs_hor_means.csv", row.names = F)

threshold_snvs<- subset(small_MAG_SNVs_hor, abs_val >= 0.5)
threshold_snvs$control <- with(threshold_snvs, ifelse(control_mean < 0.5, 'low', 'high'))
threshold_snvs$control <- with(threshold_snvs, ifelse(control_mean == 0.5, 0.5, control))
threshold_snvs$control_ref<- with(threshold_snvs, ifelse(control=="high", apply(threshold_snvs[c(9,10,12,14,15)], 1, min, na.rm=T), apply(threshold_snvs[c(9,10,12,14,15)], 1, max, na.rm=T)))
threshold_snvs$control_ref<- with(threshold_snvs, ifelse(control==0.5, ifelse(GBH_mean >= 0.5, apply(threshold_snvs[c(9,10,12,14,15)], 1, max, na.rm=T), apply(threshold_snvs[c(9,10,12,14,15)], 1, min, na.rm=T)), control_ref))
threshold_snvs$GBH_ref <- with(threshold_snvs, ifelse(control=="high", apply(threshold_snvs[c(17,18,19,20)], 1, max, na.rm=T), apply(threshold_snvs[c(17,18,19,20)], 1, min, na.rm=T)))
threshold_snvs$GBH_ref<- with(threshold_snvs, ifelse(control==0.5, ifelse(GBH_mean >= 0.5, apply(threshold_snvs[c(17,18,19,20)], 1, min, na.rm=T), apply(threshold_snvs[c(17,18,19,20)], 1, max, na.rm=T)), GBH_ref))
threshold_snvs$pass <- with(threshold_snvs, ifelse(((control=="high" & GBH_ref < control_ref) | (control=="low" & GBH_ref > control_ref) | (control== 0.5 & GBH_mean >= 0.5 & GBH_ref > control_ref) | (control== 0.5 & GBH_mean < 0.5 & GBH_ref < control_ref)), "yes", "no"))

strict <- threshold_snvs %>% subset(pass =="yes")
not_strict <- threshold_snvs %>% subset(pass=="no")
write.csv(threshold_snvs, "threshold_snvs.csv", row.names = F)

small_mag_hor_pass <- small_MAG_SNVs_hor %>% left_join(threshold_snvs[c(4,28)], by=c('groups'))
write.csv(small_mag_hor_pass, "small_mag_hor_pass.csv", row.names = F)

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
genes_sum <- strict %>% group_by(scaffold.x) %>% count(gene, name="snvs_in_gene")
gene_locations<- left_join(genes_sum, all_gene_coord, by=c("gene"))
#get list of snvs in non-gene positons
no_gene <- subset(strict, is.na(gene))

write.csv(all_gene_coord, file="all_genes.csv", row.names = F)
write.csv(gene_locations, file="gene_locations.csv", row.names = F)
write.csv(no_gene, file="snvs_no_gene.csv", row.names= F)

#SAVE NEW FILES WITH MEANS
class_MAG_SNVs<-all_MAG_SNVs %>% select(c('groups','name', 10:13, 15:21, 23:24, 26 ))
class_MAG_SNVs$class<-class_MAG_SNVs$class%>%str_sub(-3,-1)

L3_MAG_00058_SNVs <- subset(small_MAG_SNVs_hor, mag=="L3_MAG_00058") %>% 
  select(c(1:7, "Control C at T2", "Control D at T2","GBH C at T2","GBH D at T2", 21:24)) %>%
  pivot_longer(cols=contains("at"),  names_to="name", values_to="final_ref_freq", values_drop_na=F) %>%
  left_join(class_MAG_SNVs, by=c("groups", "name", "final_ref_freq")) %>%
  subset(scaffold.x !="NA")
write.csv(L3_MAG_00058_SNVs, "L3_MAG_00058_SNVs.csv", row.names = F)

L4_MAG_00099_SNVs <- subset(small_MAG_SNVs_hor, mag=="L4_MAG_00099") %>% 
  select(c(1:7, "Control D at T2", "GBH A at T2", "GBH B at T2","GBH C at T2","GBH D at T2", 21:24)) %>%
  pivot_longer(cols=contains("at"),  names_to="name", values_to="final_ref_freq", values_drop_na=F) %>% 
  left_join(class_MAG_SNVs, by=c("groups", "name", "final_ref_freq")) %>%
  subset(scaffold.x !="NA")
write.csv(L4_MAG_00099_SNVs, "L4_MAG_00099_SNVs.csv", row.names = F)

L8_MAG_00019_SNVs <- subset(small_MAG_SNVs_hor, mag=="L8_MAG_00019") %>% 
  select(c(1:7, "Control E at T2", "GBH A at T2","GBH B at T2","GBH D at T2", 21:24)) %>%
  pivot_longer(cols=contains("at"),  names_to="name", values_to="final_ref_freq", values_drop_na=F) %>% 
  left_join(class_MAG_SNVs, by=c("groups", "name", "final_ref_freq")) %>% 
  subset(scaffold.x !="NA")
write.csv(L8_MAG_00019_SNVs, "L8_MAG_00019_SNVs.csv", row.names = F)

L8_MAG_00011_SNVs <- subset(small_MAG_SNVs_hor, mag=="L8_MAG_00011") %>% 
  select(c(1:7, "Control E at T2", "GBH A at T2","GBH D at T2", 21:24)) %>%
  pivot_longer(cols=contains("at"),  names_to="name", values_to="final_ref_freq", values_drop_na=F) %>% 
  left_join(class_MAG_SNVs, by=c("groups", "name", "final_ref_freq")) %>% 
  subset(scaffold.x !="NA")
write.csv(L8_MAG_00011_SNVs, "L8_MAG_00011_SNVs.csv", row.names = F)

L7_MAG_00043_SNVs <- subset(small_MAG_SNVs_hor, mag=="L7_MAG_00043") %>% 
  select(c(1:7, "Control D at T2", "GBH B at T2","GBH C at T2", 21:24)) %>%
  pivot_longer(cols=contains("at"),  names_to="name", values_to="final_ref_freq", values_drop_na=F) %>% 
  left_join(class_MAG_SNVs, by=c("groups", "name", "final_ref_freq")) %>%
  subset(scaffold.x !="NA")
write.csv(L7_MAG_00043_SNVs, "L7_MAG_00043_SNVs.csv", row.names = F)

L7_MAG_00028_SNVs <- subset(small_MAG_SNVs_hor, mag=="L7_MAG_00028") %>% 
  select(c(1:7, "Control E at T2", "GBH A at T2","GBH B at T2","GBH C at T2", 21:24)) %>%
  pivot_longer(cols=contains("at"),  names_to="name", values_to="final_ref_freq", values_drop_na=F) %>% 
  left_join(class_MAG_SNVs, by=c("groups", "name", "final_ref_freq")) %>%
  subset(scaffold.x !="NA")
write.csv(L7_MAG_00028_SNVs, "L7_MAG_00028_SNVs.csv", row.names = F)

I4_MAG_00006_SNVs <- subset(small_MAG_SNVs_hor, mag=="I4_MAG_00006") %>% 
  select(c(1:7, "Control B at T2", "Control E at T2","GBH A at T2","GBH D at T2", 21:24)) %>%
  pivot_longer(cols=contains("at"),  names_to="name", values_to="final_ref_freq", values_drop_na=F) %>% 
  left_join(class_MAG_SNVs, by=c("groups", "name", "final_ref_freq")) %>%
  subset(scaffold.x !="NA")
write.csv(I4_MAG_00006_SNVs, "I4_MAG_00006_SNVs.csv", row.names = F)

I4_MAG_00065_SNVs <- subset(small_MAG_SNVs_hor, mag=="I4_MAG_00065") %>% 
  select(c(1:7, "Control A at T2","Control B at T2","Control C at T2", "Control D at T2","Control E at T2","GBH B at T2", 21:24)) %>%
  pivot_longer(cols=contains("at"),  names_to="name", values_to="final_ref_freq", values_drop_na=F) %>% 
  left_join(class_MAG_SNVs, by=c("groups", "name", "final_ref_freq")) %>%
  subset(scaffold.x !="NA")
write.csv(I4_MAG_00065_SNVs, "I4_MAG_00065_SNVs.csv", row.names = F)

L7_MAG_00020_SNVs <- subset(small_MAG_SNVs_hor, mag=="L7_MAG_00020") %>% 
  select(c(1:7, "Control A at T1","Control A at T2","Control C at T1", "Control C at T2","Control D at T1", "Control D at T2","GBH A at T1", "GBH A at T2", 21:24)) %>%
  pivot_longer(cols=contains("at"),  names_to="name", values_to="final_ref_freq", values_drop_na=F) %>% 
  left_join(class_MAG_SNVs, by=c("groups", "name", "final_ref_freq")) %>%
  subset(scaffold.x !="NA")
write.csv(L7_MAG_00020_SNVs, "L7_MAG_00020_SNVs.csv", row.names = F)

L8_MAG_00042_SNVs <- subset(small_MAG_SNVs_hor, mag=="L8_MAG_00042") %>% 
  select(c(1:7, "Control A at T2", "Control C at T2","Control D at T2","Control E at T2","GBH D at T2", 21:24)) %>%
  pivot_longer(cols=contains("at"),  names_to="name", values_to="final_ref_freq", values_drop_na=F) %>% 
  left_join(class_MAG_SNVs, by=c("groups", "name", "final_ref_freq")) %>%
  subset(scaffold.x !="NA")
write.csv(L8_MAG_00042_SNVs, "L8_MAG_00042_SNVs.csv", row.names = F)

L2_MAG_00052_SNVs <- subset(small_MAG_SNVs_hor, mag=="L2_MAG_00052") %>% 
  select(c(1:7, "Control A at T2", "Control B at T2","Control D at T2","Control E at T2","GBH A at T2", 21:24)) %>%
  pivot_longer(cols=contains("at"),  names_to="name", values_to="final_ref_freq", values_drop_na=F) %>% 
  left_join(class_MAG_SNVs, by=c("groups", "name", "final_ref_freq")) %>%
  subset(scaffold.x !="NA")
write.csv(L2_MAG_00052_SNVs, "L2_MAG_00052_SNVs.csv", row.names = F)

all_finished_SNVs<-rbind(L3_MAG_00058_SNVs, L4_MAG_00099_SNVs, L8_MAG_00019_SNVs, L8_MAG_00011_SNVs, L7_MAG_00043_SNVs,
                         L7_MAG_00028_SNVs, I4_MAG_00006_SNVs, I4_MAG_00065_SNVs, L7_MAG_00020_SNVs, L8_MAG_00042_SNVs, L2_MAG_00052_SNVs) 

all_finished_SNVs <- all_finished_SNVs %>% left_join(threshold_snvs[c(1:5, 28)], by=c("scaffold.x", "position.x", "length", "groups", "gene"))

write.csv(all_finished_SNVs, "all_MAG_SNVs_med_July25.csv", row.names=F)

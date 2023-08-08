library(tidyverse)
library(dplyr)
library(tidyr)

mag_list <- list("I4_MAG_00006", "I4_MAG_00065", "L2_MAG_00052", "L3_MAG_00058", "L4_MAG_00099",
                 "L7_MAG_00020", "L7_MAG_00028", "L7_MAG_00043", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")

ponds_list <- list(I4_MAG_00006 = list("Control B at T2", "Control E at T2", "GBH A at T2", "GBH D at T2"), 
                   I4_MAG_00065 = list("Control A at T2", "Control B at T2", "Control C at T2", "Control D at T2", "Control E at T2", "GBH B at T2"), 
                   L2_MAG_00052 = list("Control A at T2", "Control B at T2", "Control D at T2", "Control E at T2", "GBH A at T2"), 
                   L3_MAG_00058 = list("Control C at T2", "Control D at T2", "GBH C at T2", "GBH D at T2"), 
                   L4_MAG_00099 = list("Control D at T2", "GBH A at T2", "GBH B at T2", "GBH C at T2", "GBH D at T2"), 
                   L7_MAG_00020 = list("Control A at T1", "Control A at T2", "Control C at T1", "Control C at T2", "Control D at T1", "Control D at T2", "GBH A at T1", "GBH A at T2"), 
                   L7_MAG_00028 = list("Control E at T2", "GBH A at T2", "GBH B at T2", "GBH C at T2"), 
                   L7_MAG_00043 = list("Control D at T2", "GBH B at T2", "GBH C at T2"), 
                   L8_MAG_00011 = list("Control E at T2", "GBH A at T2", "GBH D at T2"), 
                   L8_MAG_00019 = list("Control E at T2", "GBH A at T2", "GBH B at T2", "GBH D at T2"), 
                   L8_MAG_00042 = list("Control A at T2", "Control C at T2", "Control D at T2", "Control E at T2", "GBH D at T2"))

all_MAG_SNVs <- read_csv("all_MAG_SNVs.csv")

small_MAG_SNVs<- all_MAG_SNVs[c('scaffold', 'position', 'length', 'groups', 'gene', 'mag', 'mag_length', 'new_name', 'final_ref_freq')]
small_MAG_SNVs_hor <- spread(small_MAG_SNVs, key=new_name, value=final_ref_freq)
small_MAG_SNVs_hor$all_mean <-rowMeans(small_MAG_SNVs_hor[c(8:20)], na.rm=T)
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

small_mag_hor_pass <- left_join(small_MAG_SNVs_hor, threshold_snvs[,c(4,28)], by=c('groups'))
write.csv(small_mag_hor_pass, "small_mag_hor_pass.csv", row.names = F)

#new snv files with means
class_MAG_SNVs<-all_MAG_SNVs %>% select(c('groups','new_name', 9:14, 20:33))
class_MAG_SNVs$class<-class_MAG_SNVs$class%>%str_sub(-3,-1)

all_finished_SNVs <- data.frame()

for(MAG in mag_list){
  MAG_index<-which(names(ponds_list)==MAG)
  MAG_SNVs <- subset(small_mag_hor_pass, mag==MAG)
  MAG_SNVs <- MAG_SNVs[, c(1:7, which((names(MAG_SNVs) %in% ponds_list[[MAG]])==T), 21:25)] %>%
  pivot_longer(cols=contains("at"),  names_to="new_name", values_to="final_ref_freq", values_drop_na=F) %>%
  left_join(class_MAG_SNVs, by = join_by(length, groups, mag, mag_length, new_name, final_ref_freq))
  write.csv(MAG_SNVs, paste(MAG, "_SNVs.csv", sep=""), row.names=F)
  all_finished_SNVs <- rbind(all_finished_SNVs, MAG_SNVs)
  print(paste("done", MAG))
}

write.csv(all_finished_SNVs, "all_MAG_SNVs_med_July25.csv", row.names=F)


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
all_gene_coord <- all_genes[c('gene', 'start', 'end')]
all_gene_coord <- all_gene_coord %>% distinct(gene, .keep_all = T)
all_gene_coord$new_start<-all_gene_coord$start+1
all_gene_coord$new_end<-all_gene_coord$end+1
#get gene positions that pass threshold
genes_sum <- strict %>% count(scaffold, gene, name="snvs_in_gene")
gene_locations<- left_join(genes_sum, all_gene_coord, by=c("gene"))
#get list of snvs in non-gene positons
no_gene <- subset(strict, is.na(gene))

write.csv(all_gene_coord, file="all_genes.csv", row.names = F)
write.csv(gene_locations, file="gene_locations.csv", row.names = F)
write.csv(no_gene, file="snvs_no_gene.csv", row.names= F)


# for sweeps
sweep_hor <- spread(sweep, key=new_name, value=final_ref_freq)
sweep_hor$all_mean <-rowMeans(sweep_hor[,c(8:20)], na.rm=T)
sweep_hor$control_mean <-rowMeans(sweep_hor[,c(9,10,12,14,15)], na.rm=T)
sweep_hor$GBH_mean <-rowMeans(sweep_hor[,c(17,18,19,20)], na.rm=T)
sweep_hor$abs_val<- abs(sweep_hor$control_mean - sweep_hor$GBH_mean)

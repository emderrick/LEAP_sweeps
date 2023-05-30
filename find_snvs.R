library(tidyverse)
library(dplyr)

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

L3_MAG_00058 <- read_csv("all_L3_MAG_00058_SNVs.csv")
#L3 and L4 are controls, L7 and L8 are gly
L3_MAG_00058_L3_1 <- filter(L3_MAG_00058, (pond=="L3" & final_ref_freq >= 0.5))
L3_MAG_00058_L4_1 <- filter(L3_MAG_00058, (pond=="L4" & final_ref_freq >= 0.5))
L3_MAG_00058_L7_1 <- filter(L3_MAG_00058, (pond=="L7" & final_ref_freq < 0.5))
L3_MAG_00058_L8_1 <- filter(L3_MAG_00058, (pond=="L8" & final_ref_freq < 0.5))
L3_MAG_00058_50_1 <- (list(L3_MAG_00058_L3_1,L3_MAG_00058_L4_1,L3_MAG_00058_L7_1,L3_MAG_00058_L8_1) %>% reduce(inner_join, by='groups'))

L3_MAG_00058_L3_2 <- filter(L3_MAG_00058, (pond=="L3" & final_ref_freq < 0.5))
L3_MAG_00058_L4_2 <- filter(L3_MAG_00058, (pond=="L4" & final_ref_freq < 0.5))
L3_MAG_00058_L7_2 <- filter(L3_MAG_00058, (pond=="L7" & final_ref_freq >= 0.5))
L3_MAG_00058_L8_2 <- filter(L3_MAG_00058, (pond=="L8" & final_ref_freq >= 0.5))
L3_MAG_00058_50_2 <- (list(L3_MAG_00058_L3_2,L3_MAG_00058_L4_2,L3_MAG_00058_L7_2,L3_MAG_00058_L8_2) %>% reduce(inner_join, by='groups'))
L3_MAG_00058_50_all <- bind_rows(L3_MAG_00058_50_1, L3_MAG_00058_50_2)
write.csv(L3_MAG_00058_50_all, "L3_MAG_00058_50.csv", row.names=F)

L4_MAG_00099 <- read_csv("all_L4_MAG_00099_SNVs.csv")
#L4 is a control, L2, L6, L7, and L8 are gly
L4_MAG_00099_L4_1 <- filter(L4_MAG_00099, (pond=="L4" & final_ref_freq >= 0.5))
L4_MAG_00099_L2_1 <- filter(L4_MAG_00099, (pond=="L2" & final_ref_freq < 0.5))
L4_MAG_00099_L6_1 <- filter(L4_MAG_00099, (pond=="L6" & final_ref_freq < 0.5))
L4_MAG_00099_L7_1 <- filter(L4_MAG_00099, (pond=="L7" & final_ref_freq < 0.5))
L4_MAG_00099_L8_1 <- filter(L4_MAG_00099, (pond=="L8" & final_ref_freq < 0.5))
L4_MAG_00099_50_1 <- (list(L4_MAG_00099_L4_1,L4_MAG_00099_L2_1,L4_MAG_00099_L6_1,L4_MAG_00099_L7_1,L4_MAG_00099_L8_1) %>% reduce(inner_join, by='groups'))

L4_MAG_00099_L4_2 <- filter(L4_MAG_00099, (pond=="L4" & final_ref_freq < 0.5))
L4_MAG_00099_L2_2 <- filter(L4_MAG_00099, (pond=="L2" & final_ref_freq >= 0.5))
L4_MAG_00099_L6_2 <- filter(L4_MAG_00099, (pond=="L6" & final_ref_freq >= 0.5))
L4_MAG_00099_L7_2 <- filter(L4_MAG_00099, (pond=="L7" & final_ref_freq >= 0.5))
L4_MAG_00099_L8_2 <- filter(L4_MAG_00099, (pond=="L8" & final_ref_freq >= 0.5))
L4_MAG_00099_50_2 <- (list(L4_MAG_00099_L4_2,L4_MAG_00099_L2_2,L4_MAG_00099_L6_2,L4_MAG_00099_L7_2,L4_MAG_00099_L8_2) %>% reduce(inner_join, by='groups'))
L4_MAG_00099_50_all <- bind_rows(L4_MAG_00099_50_1, L4_MAG_00099_50_2)
write.csv(L4_MAG_00099_50_all, "L4_MAG_00099_50.csv", row.names=F)

L8_MAG_00019 <- read_csv("all_L8_MAG_00019_SNVs.csv")
#I8 is a control, L2, L6 and L8 are gly
L8_MAG_00019_I8_1 <- filter(L8_MAG_00019, (pond=="I8" & final_ref_freq < 0.5))
L8_MAG_00019_L2_1 <- filter(L8_MAG_00019, (pond=="L2" & final_ref_freq >= 0.5))
L8_MAG_00019_L6_1 <- filter(L8_MAG_00019, (pond=="L6" & final_ref_freq >= 0.5))
L8_MAG_00019_L8_1 <- filter(L8_MAG_00019, (pond=="L8" & final_ref_freq >= 0.5))
L8_MAG_00019_50_1 <- (list(L8_MAG_00019_I8_1,L8_MAG_00019_L2_1,L8_MAG_00019_L6_1,L8_MAG_00019_L8_1) %>% reduce(inner_join, by='groups'))

L8_MAG_00019_I8_2 <- filter(L8_MAG_00019, (pond=="I8" & final_ref_freq >= 0.5))
L8_MAG_00019_L2_2 <- filter(L8_MAG_00019, (pond=="L2" & final_ref_freq < 0.5))
L8_MAG_00019_L6_2 <- filter(L8_MAG_00019, (pond=="L6" & final_ref_freq < 0.5))
L8_MAG_00019_L8_2 <- filter(L8_MAG_00019, (pond=="L8" & final_ref_freq < 0.5))
L8_MAG_00019_50_2 <- (list(L8_MAG_00019_I8_2,L8_MAG_00019_L2_2,L8_MAG_00019_L6_2,L8_MAG_00019_L8_2) %>% reduce(inner_join, by='groups'))
L8_MAG_00019_50_all <- bind_rows(L8_MAG_00019_50_1, L8_MAG_00019_50_2)
write.csv(L8_MAG_00019_50_all, "L8_MAG_00019_50.csv", row.names=F)

L8_MAG_00011 <- read_csv("all_L8_MAG_00011_SNVs.csv")
#I8 is a control, L2, and L8 are gly
L8_MAG_00011_I8_1 <- filter(L8_MAG_00011, (pond=="I8" & final_ref_freq < 0.5))
L8_MAG_00011_L2_1 <- filter(L8_MAG_00011, (pond=="L2" & final_ref_freq >= 0.5))
L8_MAG_00011_L8_1 <- filter(L8_MAG_00011, (pond=="L8" & final_ref_freq >= 0.5))
L8_MAG_00011_50_1 <- (list(L8_MAG_00011_I8_1,L8_MAG_00011_L2_1,L8_MAG_00011_L8_1) %>% reduce(inner_join, by='groups'))

L8_MAG_00011_I8_2 <- filter(L8_MAG_00011, (pond=="I8" & final_ref_freq >= 0.5))
L8_MAG_00011_L2_2 <- filter(L8_MAG_00011, (pond=="L2" & final_ref_freq < 0.5))
L8_MAG_00011_L8_2 <- filter(L8_MAG_00011, (pond=="L8" & final_ref_freq < 0.5))
L8_MAG_00011_50_2 <- (list(L8_MAG_00011_I8_2,L8_MAG_00011_L2_2,L8_MAG_00011_L8_2) %>% reduce(inner_join, by='groups'))
L8_MAG_00011_50_all <- bind_rows(L8_MAG_00011_50_1, L8_MAG_00011_50_2)
write.csv(L8_MAG_00011_50_all, "L8_MAG_00011_50.csv", row.names=F)

L7_MAG_00043 <- read_csv("all_L7_MAG_00043_SNVs.csv")
#L4 is a control, L6 and L7 are gly
L7_MAG_00043_L4_1 <- filter(L7_MAG_00043, (pond=="L4" & final_ref_freq < 0.5))
L7_MAG_00043_L6_1 <- filter(L7_MAG_00043, (pond=="L6" & final_ref_freq >= 0.5))
L7_MAG_00043_L7_1 <- filter(L7_MAG_00043, (pond=="L7" & final_ref_freq >= 0.5))
L7_MAG_00043_50_1 <- (list(L7_MAG_00043_L4_1,L7_MAG_00043_L6_1,L7_MAG_00043_L7_1) %>% reduce(inner_join, by='groups'))

L7_MAG_00043_L4_2 <- filter(L7_MAG_00043, (pond=="L4" & final_ref_freq >= 0.5))
L7_MAG_00043_L6_2 <- filter(L7_MAG_00043, (pond=="L6" & final_ref_freq < 0.5))
L7_MAG_00043_L7_2 <- filter(L7_MAG_00043, (pond=="L7" & final_ref_freq < 0.5))
L7_MAG_00043_50_2 <- (list(L7_MAG_00043_L4_2,L7_MAG_00043_L6_2,L7_MAG_00043_L7_2) %>% reduce(inner_join, by='groups'))
L7_MAG_00043_50_all <- bind_rows(L7_MAG_00043_50_1, L7_MAG_00043_50_2)
write.csv(L7_MAG_00043_50_all, "L7_MAG_00043_50.csv", row.names=F)

L7_MAG_00028 <- read_csv("all_L7_MAG_00028_SNVs.csv")
#I8 is a control, L2, L6, and L7 are gly
L7_MAG_00028_I8_1 <- filter(L7_MAG_00028, (pond=="I8" & final_ref_freq < 0.5))
L7_MAG_00028_L2_1 <- filter(L7_MAG_00028, (pond=="L2" & final_ref_freq >= 0.5))
L7_MAG_00028_L6_1 <- filter(L7_MAG_00028, (pond=="L6" & final_ref_freq >= 0.5))
L7_MAG_00028_L7_1 <- filter(L7_MAG_00028, (pond=="L7" & final_ref_freq >= 0.5))
L7_MAG_00028_50_1 <- (list(L7_MAG_00028_I8_1,L7_MAG_00028_L2_1,L7_MAG_00028_L6_1,L7_MAG_00028_L7_1) %>% reduce(inner_join, by='groups'))

L7_MAG_00028_I8_2 <- filter(L7_MAG_00028, (pond=="I8" & final_ref_freq >= 0.5))
L7_MAG_00028_L2_2 <- filter(L7_MAG_00028, (pond=="L2" & final_ref_freq < 0.5))
L7_MAG_00028_L6_2 <- filter(L7_MAG_00028, (pond=="L6" & final_ref_freq < 0.5))
L7_MAG_00028_L7_2 <- filter(L7_MAG_00028, (pond=="L7" & final_ref_freq < 0.5))
L7_MAG_00028_50_2 <- (list(L7_MAG_00028_I8_2,L7_MAG_00028_L2_2,L7_MAG_00028_L6_2,L7_MAG_00028_L7_2) %>% reduce(inner_join, by='groups'))
L7_MAG_00028_50_all <- bind_rows(L7_MAG_00028_50_1, L7_MAG_00028_50_2)
write.csv(L7_MAG_00028_50_all, "L7_MAG_00028_50.csv", row.names=F)

I4_MAG_00006 <- read_csv("all_I4_MAG_00006_SNVs.csv")
#I4 and I8 are controls, L2 and L8 are gly
I4_MAG_00006_I4_1 <- filter(I4_MAG_00006, pond=="I4" & final_ref_freq >= 0.5)
I4_MAG_00006_I8_1 <- filter(I4_MAG_00006, pond=="I8" & final_ref_freq >= 0.5)
I4_MAG_00006_L2_1 <- filter(I4_MAG_00006, pond=="L2" & final_ref_freq < 0.5) 
I4_MAG_00006_L8_1 <- filter(I4_MAG_00006, pond=="L8" & final_ref_freq < 0.5)
I4_MAG_00006_50_1 <- (list(I4_MAG_00006_I4_1, I4_MAG_00006_I8_1, I4_MAG_00006_L2_1, I4_MAG_00006_L8_1) %>% reduce(inner_join, by='groups'))

I4_MAG_00006_I4_2 <- filter(I4_MAG_00006, pond=="I4" & final_ref_freq < 0.5)
I4_MAG_00006_I8_2 <- filter(I4_MAG_00006, pond=="I8" & final_ref_freq < 0.5)
I4_MAG_00006_L2_2 <- filter(I4_MAG_00006, pond=="L2" & final_ref_freq >= 0.5) 
I4_MAG_00006_L8_2 <- filter(I4_MAG_00006, pond=="L8" & final_ref_freq >= 0.5)
I4_MAG_00006_50_2 <- (list(I4_MAG_00006_I4_2, I4_MAG_00006_I8_2, I4_MAG_00006_L2_2, I4_MAG_00006_L8_2) %>% reduce(inner_join, by='groups'))
I4_MAG_00006_50_all <- bind_rows(I4_MAG_00006_50_1, I4_MAG_00006_50_2)
write.csv(I4_MAG_00006_50_all, "I4_MAG_00006_50.csv", row.names=F)

I4_MAG_00065 <- read_csv("all_I4_MAG_00065_SNVs.csv")
#I4, I8, K1, L3, L4 are controls, L6 is gly
I4_MAG_00065_I4_1 <- filter(I4_MAG_00065, pond=="I4" & final_ref_freq >= 0.5)
I4_MAG_00065_I8_1 <- filter(I4_MAG_00065, pond=="I8" & final_ref_freq >= 0.5)
I4_MAG_00065_K1_1 <- filter(I4_MAG_00065, pond=="K1" & final_ref_freq >= 0.5)
I4_MAG_00065_L3_1 <- filter(I4_MAG_00065, pond=="L3" & final_ref_freq >= 0.5)
I4_MAG_00065_L4_1 <- filter(I4_MAG_00065, pond=="L4" & final_ref_freq >= 0.5) 
I4_MAG_00065_L6_1 <- filter(I4_MAG_00065, pond=="L6" & final_ref_freq <= 0.5)
I4_MAG_00065_50_1 <- (list(I4_MAG_00065_I4_1, I4_MAG_00065_I8_1, I4_MAG_00065_K1_1, I4_MAG_00065_L3_1, I4_MAG_00065_L4_1, I4_MAG_00065_L6_1) %>% reduce(inner_join, by='groups'))

I4_MAG_00065_I4_2 <- filter(I4_MAG_00065, pond=="I4" & final_ref_freq <= 0.5)
I4_MAG_00065_I8_2 <- filter(I4_MAG_00065, pond=="I8" & final_ref_freq <= 0.5)
I4_MAG_00065_K1_2 <- filter(I4_MAG_00065, pond=="K1" & final_ref_freq <= 0.5)
I4_MAG_00065_L3_2 <- filter(I4_MAG_00065, pond=="L3" & final_ref_freq <= 0.5)
I4_MAG_00065_L4_2 <- filter(I4_MAG_00065, pond=="L4" & final_ref_freq <= 0.5) 
I4_MAG_00065_L6_2 <- filter(I4_MAG_00065, pond=="L6" & final_ref_freq >= 0.5)
I4_MAG_00065_50_2 <- (list(I4_MAG_00065_I4_2, I4_MAG_00065_I8_2, I4_MAG_00065_K1_2, I4_MAG_00065_L3_2, I4_MAG_00065_L4_2, I4_MAG_00065_L6_2) %>% reduce(inner_join, by='groups'))
I4_MAG_00065_50_all <- bind_rows(I4_MAG_00065_50_1, I4_MAG_00065_50_2)
write.csv(I4_MAG_00065_50_all, "I4_MAG_00065_50.csv", row.names=F)

L7_MAG_00020 <- read_csv("all_L7_MAG_00020_SNVs.csv")
#K1, L3, and L4 are a control, L2 is gly
L7_MAG_00020_K1_1 <- filter(L7_MAG_00020, (pond=="K1" & final_ref_freq <= 0.5))
L7_MAG_00020_L3_1 <- filter(L7_MAG_00020, (pond=="L3" & final_ref_freq <= 0.5))
L7_MAG_00020_L4_1 <- filter(L7_MAG_00020, (pond=="L4" & final_ref_freq <= 0.5))
L7_MAG_00020_L2_1 <- filter(L7_MAG_00020, (pond=="L2" & final_ref_freq >= 0.5))
L7_MAG_00020_50_1 <- (list(L7_MAG_00020_K1_1,L7_MAG_00020_L3_1,L7_MAG_00020_L4_1,L7_MAG_00020_L2_1) %>% reduce(inner_join, by='groups'))

L7_MAG_00020_K1_2 <- filter(L7_MAG_00020, (pond=="K1" & final_ref_freq >= 0.5))
L7_MAG_00020_L3_2 <- filter(L7_MAG_00020, (pond=="L3" & final_ref_freq >= 0.5))
L7_MAG_00020_L4_2 <- filter(L7_MAG_00020, (pond=="L4" & final_ref_freq >= 0.5))
L7_MAG_00020_L2_2 <- filter(L7_MAG_00020, (pond=="L2" & final_ref_freq <= 0.5))
L7_MAG_00020_50_2 <- (list(L7_MAG_00020_K1_2,L7_MAG_00020_L3_2,L7_MAG_00020_L4_2,L7_MAG_00020_L2_2) %>% reduce(inner_join, by='groups'))
L7_MAG_00020_50_all <- bind_rows(L7_MAG_00020_50_1, L7_MAG_00020_50_2)
write.csv(L7_MAG_00020_50_all, "L7_MAG_00020_50.csv", row.names=F)

L8_MAG_00042 <- read_csv("all_L8_MAG_00042_SNVs.csv")
#I8, K1, L3, and L4 are a control, L8 is gly
L8_MAG_00042_I8_1 <- filter(L8_MAG_00042, (pond=="I8" & final_ref_freq <= 0.5))
L8_MAG_00042_K1_1 <- filter(L8_MAG_00042, (pond=="K1" & final_ref_freq <= 0.5))
L8_MAG_00042_L3_1 <- filter(L8_MAG_00042, (pond=="L3" & final_ref_freq <= 0.5))
L8_MAG_00042_L4_1 <- filter(L8_MAG_00042, (pond=="L4" & final_ref_freq <= 0.5))
L8_MAG_00042_L8_1 <- filter(L8_MAG_00042, (pond=="L8" & final_ref_freq >= 0.5))
L8_MAG_00042_50_1 <- (list(L8_MAG_00042_I8_1,L8_MAG_00042_K1_1,L8_MAG_00042_L3_1,L8_MAG_00042_L4_1,L8_MAG_00042_L8_1) %>% reduce(inner_join, by='groups'))

L8_MAG_00042_I8_2 <- filter(L8_MAG_00042, (pond=="I8" & final_ref_freq >= 0.5))
L8_MAG_00042_K1_2 <- filter(L8_MAG_00042, (pond=="K1" & final_ref_freq >= 0.5))
L8_MAG_00042_L3_2 <- filter(L8_MAG_00042, (pond=="L3" & final_ref_freq >= 0.5))
L8_MAG_00042_L4_2 <- filter(L8_MAG_00042, (pond=="L4" & final_ref_freq >= 0.5))
L8_MAG_00042_L8_2 <- filter(L8_MAG_00042, (pond=="L8" & final_ref_freq <= 0.5))
L8_MAG_00042_50_2 <- (list(L8_MAG_00042_I8_2,L8_MAG_00042_K1_2,L8_MAG_00042_L3_2,L8_MAG_00042_L4_2,L8_MAG_00042_L8_2) %>% reduce(inner_join, by='groups'))
L8_MAG_00042_50_all <- bind_rows(L8_MAG_00042_50_1, L8_MAG_00042_50_2)
write.csv(L8_MAG_00042_50_all, "L8_MAG_00042_50.csv", row.names=F)

L2_MAG_00052 <- read_csv("all_L2_MAG_00052_SNVs.csv")
#I4, I8, K1, and L4 are a control, L2, is gly
L2_MAG_00052_I4_1 <- filter(L2_MAG_00052, (pond=="I4" & final_ref_freq <= 0.5))
L2_MAG_00052_I8_1 <- filter(L2_MAG_00052, (pond=="I8" & final_ref_freq <= 0.5))
L2_MAG_00052_K1_1 <- filter(L2_MAG_00052, (pond=="K1" & final_ref_freq <= 0.5))
L2_MAG_00052_L4_1 <- filter(L2_MAG_00052, (pond=="L4" & final_ref_freq <= 0.5))
L2_MAG_00052_L2_1 <- filter(L2_MAG_00052, (pond=="L2" & final_ref_freq >= 0.5))
L2_MAG_00052_50_1 <- (list(L2_MAG_00052_I4_1,L2_MAG_00052_I8_1,L2_MAG_00052_K1_1,L2_MAG_00052_L4_1,L2_MAG_00052_L2_1) %>% reduce(inner_join, by='groups'))

L2_MAG_00052_I4_2 <- filter(L2_MAG_00052, (pond=="I4" & final_ref_freq >= 0.5))
L2_MAG_00052_I8_2 <- filter(L2_MAG_00052, (pond=="I8" & final_ref_freq >= 0.5))
L2_MAG_00052_K1_2 <- filter(L2_MAG_00052, (pond=="K1" & final_ref_freq >= 0.5))
L2_MAG_00052_L4_2 <- filter(L2_MAG_00052, (pond=="L4" & final_ref_freq >= 0.5))
L2_MAG_00052_L2_2 <- filter(L2_MAG_00052, (pond=="L2" & final_ref_freq <= 0.5))
L2_MAG_00052_50_2 <- (list(L2_MAG_00052_I4_2,L2_MAG_00052_I8_2,L2_MAG_00052_K1_2,L2_MAG_00052_L4_2,L2_MAG_00052_L2_2) %>% reduce(inner_join, by='groups'))
L2_MAG_00052_50_all <- bind_rows(L2_MAG_00052_50_1, L2_MAG_00052_50_2)
write.csv(L2_MAG_00052_50_all, "L2_MAG_00052_50.csv", row.names=F)

I4_MAG_00006_genes<- subset(I4_MAG_00006_50_all, select=c(mag.x, groups, gene.x))
L7_MAG_00028_genes<- subset(L7_MAG_00028_50_all, select=c(mag.x, groups, gene.x))
L8_MAG_00011_genes<- subset(L8_MAG_00011_50_all, select=c(mag.x, groups, gene.x))
L8_MAG_00019_genes<- subset(L8_MAG_00019_50_all, select=c(mag.x, groups, gene.x))
L8_MAG_00042_genes<- subset(L8_MAG_00042_50_all, select=c(mag.x, groups, gene.x))
L3_MAG_00058_genes<- subset(L3_MAG_00058_50_all, select=c(mag.x, groups, gene.x))
I4_MAG_00065_genes<- subset(I4_MAG_00065_50_all, select=c(mag.x, groups, gene.x))
L4_MAG_00099_genes<- subset(L4_MAG_00099_50_all, select=c(mag.x, groups, gene.x))
L2_MAG_00052_genes<- subset(L2_MAG_00052_50_all, select=c(mag.x, groups, gene.x))
L7_MAG_00043_genes<- subset(L7_MAG_00043_50_all, select=c(mag.x, groups, gene.x))
L7_MAG_00020_genes<- subset(L7_MAG_00020_50_all, select=c(mag.x, groups, gene.x))

genes <- rbind(I4_MAG_00006_genes, L7_MAG_00028_genes, L8_MAG_00011_genes, L8_MAG_00019_genes,
                   L8_MAG_00042_genes, L3_MAG_00058_genes, I4_MAG_00065_genes, L4_MAG_00099_genes,
                   L2_MAG_00052_genes, L7_MAG_00043_genes, L7_MAG_00020_genes)
genes <- rename(genes, gene = gene.x, mag=mag.x)

all_gene_coord <- subset(all_genes, select=c(gene, start, end))
all_gene_coord <- all_gene_coord %>% distinct(gene, .keep_all = T)
all_gene_coord$new_start<-all_gene_coord$start+1
all_gene_coord$new_end<-all_gene_coord$end+1
only_genes <- subset(genes, gene!="NA", select=c(gene))
unique_genes <- distinct(only_genes)
gene_locations<- left_join(unique_genes, all_gene_coord, by=c("gene"))

write_tsv(genes, file="all_genes.tsv")
write_tsv(gene_locations, file="gene_locations.tsv")

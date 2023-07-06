library(tidyverse)
library(dplyr)

gene_files<-list.files("full_bakta_output",recursive = T, pattern=".tsv",full.names = T)

all_genes<-data.frame()
for(i in 1:length(gene_files)){
  genes<-read_tsv(gene_files[i], skip=5)
  mag<-gsub(".*bakta_output/", "", gene_files[i]) %>% substr(1,12)
  genes<-cbind(genes,mag=rep(mag,nrow(genes)))
  all_genes<-rbind(all_genes,genes)
}
#add in the direction
colnames(all_genes)[1]="Sequence_ID"
all_genes$contig <- all_genes$Sequence_ID %>% substr(8,10)
gene_locations<- read_csv("gene_locations.csv")
gene_locations$contig <- gene_locations$gene %>% substr(23,25) %>% str_remove("^0+")
gene_locations$mag <- gene_locations$gene %>% substr(1,12)
gene_matches<- left_join(gene_locations,all_genes, by=c("mag", "contig", "new_start"="Start", "new_end"="Stop"))

gene_no_NA <- gene_matches %>% drop_na(Gene)

parallel_genes <- gene_no_NA %>% count(Gene, mag)
parallel_gene_mag <- parallel_genes %>% count(Gene)

I4_MAG_00006_top_genes <- subset(gene_no_NA, mag=="I4_MAG_00006")
I4_MAG_00006_top_genes <- I4_MAG_00006_top_genes[c(2, 8, 13, 14)]

I4_MAG_00065_top_genes <- subset(gene_no_NA, mag=="I4_MAG_00065")
I4_MAG_00065_top_genes <- I4_MAG_00065_top_genes[c(2, 8, 13, 14)]

L2_MAG_00052_top_genes <- subset(gene_no_NA, mag=="L2_MAG_00052")
L2_MAG_00052_top_genes <- L2_MAG_00052_top_genes[c(2, 8, 13, 14)]

L3_MAG_00058_top_genes <- subset(gene_no_NA, mag=="L3_MAG_00058")
L3_MAG_00058_top_genes <- L3_MAG_00058_top_genes[c(2, 8, 13, 14)]

L4_MAG_00099_top_genes <- subset(gene_no_NA, mag=="L4_MAG_00099")
L4_MAG_00099_top_genes <- L4_MAG_00099_top_genes[c(2, 8, 13, 14)]

L7_MAG_00020_top_genes <- subset(gene_no_NA, mag=="L7_MAG_00020")
L7_MAG_00020_top_genes <- L7_MAG_00020_top_genes[c(2, 8, 13, 14)]

L7_MAG_00028_top_genes <- subset(gene_no_NA, mag=="L7_MAG_00028")
L7_MAG_00028_top_genes <- L7_MAG_00028_top_genes[c(2, 8, 13, 14)]

L7_MAG_00043_top_genes <- subset(gene_no_NA, mag=="L7_MAG_00043")
L7_MAG_00043_top_genes <- L7_MAG_00043_top_genes[c(2, 8, 13, 14)]

L8_MAG_00011_top_genes <- subset(gene_no_NA, mag=="L8_MAG_00011")
L8_MAG_00011_top_genes <- L8_MAG_00011_top_genes[c(2, 8, 13, 14)]

L8_MAG_00019_top_genes <- subset(gene_no_NA, mag=="L8_MAG_00019")
L8_MAG_00019_top_genes <- L8_MAG_00019_top_genes[c(2, 8, 13, 14)]

L8_MAG_00042_top_genes <- subset(gene_no_NA, mag=="L8_MAG_00042")
L8_MAG_00042_top_genes <- L8_MAG_00042_top_genes[c(2, 8, 13, 14)]



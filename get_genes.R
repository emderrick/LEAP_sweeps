library(tidyverse)
library(dplyr)

gene_files<-list.files("bakta_output",recursive = T, pattern=".*tsv",full.names = T)
all_genes<-data.frame()
for(i in 1:length(gene_files)){
  genes<-read_tsv(gene_files[i], skip=5)
  mag<-gsub(".*bakta_output/", "", gene_files[i]) %>% substr(1,12)
  genes<-cbind(genes,mag=rep(mag,nrow(genes)))
  all_genes<-rbind(all_genes,genes)
}

colnames(all_genes)[1]="Sequence_ID"
all_genes$contig <- all_genes$Sequence_ID %>% substr(8,10)
gene_locations<- read_tsv("gene_locations.tsv")
gene_locations$contig <- gene_locations$gene %>% substr(23,25) %>% str_remove("^0+")
gene_locations$mag <- gene_locations$gene %>% substr(1,12)
gene_matches<- left_join(gene_locations,all_genes, by=c("mag", "contig", "new_start"="Start", "new_end"="Stop"))
gene_cog <- filter(gene_matches, str_detect(DbXrefs, "COG"))
gene_cog$cog_class<- gene_cog$DbXrefs %>% substr(1,19)

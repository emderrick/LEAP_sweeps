library(tidyverse)
library(dplyr)

threshold_snvs <- read_csv("threshold_snvs.csv")
threshold_snvs <- subset(threshold_snvs, pass=="yes")

gene_files <- list.files("95_profiles/",recursive = T, pattern=".*gene_info.tsv",full.names = T)
all_genes <- data.frame()

for(i in 1:length(gene_files)){
  pond_time_genes <- read.table(gene_files[i],sep="\t",header=T)
  timepoint <- gsub(".*profile_output/", "", gene_files[i]) %>% substr(1,9)
  pond_time_genes <- cbind(pond_time_genes,timepoint=rep(timepoint,nrow(pond_time_genes)))
  all_genes <- rbind(all_genes,pond_time_genes)
}

all_genes$mag <- all_genes$gene %>% substr(1,12)
test <- subset(all_genes, mag == "L8_MAG_00011")
#get all gene coordinates
all_gene_coord <- all_genes[c('gene', 'start', 'end')]
all_gene_coord <- all_gene_coord %>% distinct(gene, .keep_all = T)
all_gene_coord$new_start <- all_gene_coord$start+1
all_gene_coord$new_end <- all_gene_coord$end+1

bakta_files <- list.files("full_bakta_output",recursive = T, pattern=".tsv",full.names = T)
all_bakta <- data.frame()

for(i in 1:length(bakta_files)){
  genes <- read_tsv(bakta_files[i], skip=5)
  mag <- gsub(".*bakta_output/", "", bakta_files[i]) %>% substr(1,12)
  genes <- cbind(genes,mag=rep(mag,nrow(genes)))
  all_bakta <- rbind(all_bakta,genes)
}

parevol_genes <- read_csv("not_strict_MAG_significant_genes.csv")
parevol_genes <- left_join(parevol_genes, all_gene_coord)
parevol_gene_match <- inner_join(parevol_genes, all_bakta, by = c("mag", "contig", "new_start"="Start", "new_end"="Stop"))

colnames(all_bakta)[1]="Sequence_ID"
all_bakta$contig <- all_bakta$Sequence_ID %>% substr(8,10) 
all_gene_coord$contig <- all_gene_coord$gene %>% substr(23,25) %>% str_remove("^0+")
all_gene_coord$mag <- all_gene_coord$gene %>% substr(1,12)
gene_matches <- inner_join(all_gene_coord, all_bakta, by = c("mag", "contig", "new_start"="Start", "new_end"="Stop"))





#get gene positions that pass threshold
genes_sum <- threshold_snvs %>% count(scaffold, gene, name="snvs_in_gene")
gene_locations <- left_join(genes_sum, all_gene_coord, by=c("gene"))
#get list of snvs in non-gene positons
no_gene <- subset(threshold_snvs, is.na(gene))

write.csv(all_gene_coord, "all_genes.csv", row.names = F)
write.csv(gene_locations, "gene_locations.csv", row.names = F)
write.csv(no_gene, "snvs_no_gene.csv", row.names= F)

parallel_genes <- gene_no_NA %>% count(Gene, mag) %>% count(Gene)

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

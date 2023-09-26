library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggVennDiagram)
library(viridis)
library(cowplot)

mag_list <- list("I4_MAG_00006", "I4_MAG_00065", "L2_MAG_00052", "L3_MAG_00058", "L4_MAG_00099",
                 "L7_MAG_00020", "L7_MAG_00028", "L7_MAG_00043", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")

gene_files <- list.files("95_profiles/",recursive = T, pattern=".*gene_info.tsv",full.names = T)
all_genes <- data.frame()

for(i in 1:length(gene_files)){
  pond_time_genes <- read.table(gene_files[i],sep="\t",header=T)
  timepoint <- gsub(".*profile_output/", "", gene_files[i]) %>% substr(1,9)
  pond_time_genes <- cbind(pond_time_genes,timepoint=rep(timepoint,nrow(pond_time_genes)))
  all_genes <- rbind(all_genes,pond_time_genes)
}

bakta_files <- list.files("MAG bakta output/",recursive = T, pattern=".tsv",full.names = T)
all_bakta <- data.frame()

for(i in 1:length(bakta_files)){
  genes <- read_tsv(bakta_files[i], skip=5)
  mag <- gsub(".*bakta_output/", "", bakta_files[i]) %>% substr(1,12)
  genes <- cbind(genes,mag=rep(mag,nrow(genes)))
  all_bakta <- rbind(all_bakta,genes)
}

prokka_files <- list.files("prokka_output/",recursive = T, pattern=".tsv",full.names = T)
all_prokka <- data.frame()

for(i in 1:length(prokka_files)){
  genes <- read_tsv(prokka_files[i])
  mag <- prokka_files[i] %>% substr(16,27)
  genes <- cbind(genes, mag = rep(mag, nrow(genes)))
  all_prokka <- rbind(all_prokka, genes)
}

#get all gene coordinates
all_gene_coord <- all_genes[c('gene', 'start', 'end')]
all_gene_coord <- subset(all_gene_coord, mag %in% mag_list)
all_gene_coord <- all_gene_coord %>% distinct(gene, .keep_all = T)
all_gene_coord$start <- all_gene_coord$start+1
all_gene_coord$end <- all_gene_coord$end+1
all_gene_coord$scaffold <- all_gene_coord$gene %>% substr(1,25)
all_gene_coord$mag <- all_gene_coord$gene %>% substr(1,12)

colnames(all_bakta)[1]="scaffold"
bakta_gene_coord <- left_join(all_bakta, all_gene_coord, by = c("mag", "scaffold", "Start" = "start", "Stop" = "end"))
bakta_gene_coord$COG_bakta<- str_extract(bakta_gene_coord$DbXrefs, "COG\\d{4}")

prokka_gene_COG <- all_prokka[, c("gene", "COG")]
prokka_gene_COG <- distinct(prokka_gene_COG) %>% na.omit()

bakta_with_COG <- left_join(bakta_gene_coord, prokka_gene_COG, by = c("Gene" = "gene"))

all_background_genes <- bakta_with_COG[, c("mag", "scaffold", "gene", "Gene", "Product", "COG_bakta", "COG")] %>% subset(is.na(gene) == F)
write.csv(all_background_genes, "all_background_genes.csv", row.names = F)

not_strict_parevol_genes <- read_csv("not_strict_MAG_significant_genes.csv")
not_strict_significant_genes <- left_join(not_strict_parevol_genes, all_background_genes)
write.csv(not_strict_significant_genes, "significant_genes_not_strict.csv", row.names = F)

strict_parevol_genes <- read_csv("MAG_significant_genes.csv")
strict_significant_genes <- left_join(strict_parevol_genes, all_background_genes)
write.csv(strict_significant_genes, "significant_genes_strict.csv", row.names = F)

threshold_snvs <- read_csv("threshold_snvs.csv")
threshold_snvs <- subset(threshold_snvs, pass=="yes")
threshold_snvs$snv_count <- 1
threshold_snvs_sum <- threshold_snvs %>% group_by(mag, scaffold, gene) %>% summarize(snvs_in_gene = sum(snv_count))
threshold_significant_genes <- left_join(threshold_snvs_sum, all_background_genes)
write.csv(threshold_significant_genes, "threshold_significant_genes.csv", row.names = F)

threshold_snvs_sum  <- na.omit(threshold_snvs_sum)
threshold_genes <- as.vector(threshold_snvs_sum$gene)
parevol_strict_genes <- as.vector(strict_parevol_genes$gene)
parevol_not_strict_genes <- as.vector(not_strict_parevol_genes$gene)

all_sig_genes <- list(threshold_genes, parevol_strict_genes, parevol_not_strict_genes)

gene_VD <- ggVennDiagram(all_sig_genes, 
              category.names = c("Allele Freq", "Strict", "Loose"), 
              label = c("count"), edge_size = 0) +
              scale_x_continuous(expand = expansion(mult = .2))+
              scale_fill_gradient('Genes', high = "purple4", low = "thistle")

save_plot("gene_venndiagram.jpeg", gene_VD)

not_strict_parevol_genes$count <- 1
not_strict_parevol_sum <- not_strict_parevol_genes %>% group_by(mag) %>% summarize(genes = sum(count))

strict_parevol_genes$count <- 1
strict_parevol_sum <- strict_parevol_genes %>% group_by(mag) %>% summarize(genes = sum(count))

threshold_snvs_sum$count <- 1
threshold_snvs_gene_sum  <- threshold_snvs_sum  %>% group_by(mag) %>% summarize(genes = sum(count))

mag_list2 <- list("I4_MAG_00006", "I4_MAG_00065", "L2_MAG_00052", "L3_MAG_00058", "L4_MAG_00099",
                  "L7_MAG_00020", "L7_MAG_00028", "L7_MAG_00043", "L8_MAG_00019")

MAG_overlap <- data.frame()
for(MAG in mag_list2){
  MAG_parevol_genes <- subset(not_strict_parevol_genes, mag == MAG)
  MAG_threshold_genes <- subset(threshold_snvs_sum, mag == MAG)
  MAG_parevol_gene_vector <- as.vector(MAG_parevol_genes$gene)
  MAG_threshold_gene_vector <- as.vector(MAG_threshold_genes$gene)
  overlaping_genes <- as.data.frame(intersect(MAG_parevol_genes$gene, MAG_threshold_genes$gene))
  overlaping_genes$mag <- MAG
  MAG_overlap <- rbind(MAG_overlap, overlaping_genes)
}

MAG_overlap$count <- 1
overlap_summary <- MAG_overlap %>% group_by(mag) %>% summarise(total_overlap = sum(count))

# threshold_multiple_matches <- as.data.frame(threshold_significant_genes$gene[duplicated(threshold_significant_genes$gene)]) %>% na.omit()
# #get gene positions that pass threshold
# genes_sum <- threshold_snvs %>% count(scaffold, gene, name="snvs_in_gene")
# gene_locations <- left_join(genes_sum, all_gene_coord, by=c("gene"))
# #get list of snvs in non-gene positons
# no_gene <- subset(threshold_snvs, is.na(gene))
# 
# write.csv(all_gene_coord, "all_genes.csv", row.names = F)
# write.csv(gene_locations, "gene_locations.csv", row.names = F)
# write.csv(no_gene, "snvs_no_gene.csv", row.names= F)
# 
# parallel_genes <- gene_no_NA %>% count(Gene, mag) %>% count(Gene)


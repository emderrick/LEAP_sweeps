library(tidyverse)

COG_gene_to_category <- read.table("cog-20.to_category.tsv", header = FALSE, sep = "\t")

gene_cov_sig <- read_csv("gene_cov_significant.csv")
gene_cov_sig_cat <- left_join(gene_cov_sig, COG_gene_to_category, by = c("COG_ID" = "V1"))
loose_parevol_genes <- read_csv("significant_genes_not_strict.csv")
loose_parevol_genes_cat <- left_join(loose_parevol_genes, COG_gene_to_category, by = c("COG_ID" = "V1"))
strict_parevol_genes <- read_csv("significant_genes_strict.csv")
strict_parevol_genes_cat <- left_join(strict_parevol_genes, COG_gene_to_category, by = c("COG_ID" = "V1"))
allele_freq_sig <- read_csv("threshold_significant_genes.csv")
allele_freq_sig_cat <- left_join(allele_freq_sig, COG_gene_to_category, by = c("COG_ID" = "V1"))

gene_cov_sig_E <- subset(gene_cov_sig_cat, V2 == "E")
loose_parevol_genes_E <- subset(loose_parevol_genes_cat, V2 == "E")
strict_parevol_genes_E <- subset(strict_parevol_genes_cat, V2 == "E")
allele_freq_sig_E <- subset(allele_freq_sig_cat, V2 == "E")


allele_freq_genes <- as.list(allele_freq_sig$gene)
parevol_strict_genes <- as.list(strict_parevol_genes$gene)
parevol_loose_genes <- as.list(loose_parevol_genes$gene)
gene_cov_genes <- as.list(gene_cov_sig$gene)

gene_overlap <- as.data.frame(intersect(allele_freq_sig$gene, strict_parevol_genes$gene))
colnames(gene_overlap)[1]="gene"
gene_overlap$mag <- gene_overlap$gene %>% substr(1,12)
gene_overlap$count = 1
gene_overlap_sum <- gene_overlap %>% group_by(mag) %>% summarize(overlap = sum(count))

ge
#gene_overlap_E <- Reduce(intersect, list(loose_parevol_genes_E$gene, strict_parevol_genes_E$gene, allele_freq_sig_E$gene))

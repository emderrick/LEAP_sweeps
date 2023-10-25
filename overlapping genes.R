library(tidyverse)

gene_cov_sig <- read_csv("gene_cov_significant.csv")
loose_parevol_genes <- read_csv("significant_genes_not_strict.csv")
strict_parevol_genes <- read_csv("significant_genes_strict.csv")
allele_freq_sig <- read_csv("threshold_significant_genes.csv")

allele_freq_genes <- as.list(allele_freq_sig$gene)
parevol_strict_genes <- as.list(strict_parevol_genes$gene)
parevol_loose_genes <- as.list(loose_parevol_genes$gene)
gene_cov_genes <- as.list(gene_cov_sig$gene)

gene_overlap <- Reduce(intersect, list(allele_freq_genes, parevol_loose_genes, parevol_strict_genes, gene_cov_genes))

gene_cov_sig_K <- subset(gene_cov_sig, COG_category == "K")
loose_parevol_genes_K <- subset(loose_parevol_genes, COG_category == "K")
strict_parevol_genes_K <- subset(strict_parevol_genes, COG_category == "K")
allele_freq_sig_K <- subset(allele_freq_sig, COG_category == "K")

gene_overlap_K <- Reduce(intersect, list(loose_parevol_genes_K$gene, strict_parevol_genes_K$gene, allele_freq_sig_K$gene))

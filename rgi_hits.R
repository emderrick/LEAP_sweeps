library(tidyverse)
setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

arg_files <- list.files("refined data files/refined_MAG_RGI/strict_RGI/", recursive = T, pattern = ".*RGI.txt", full.names = T)

arg_genes <- data.frame()
for(i in 1:length(arg_files)){
  mag_args <- read_tsv(arg_files[i])
  mag_args$mag <- arg_files[i] %>% substr(48,58)
  arg_genes <- rbind(arg_genes, mag_args)
}

arg_genes$mag  <- arg_genes$mag  %>% str_remove("_[A-Z]")
arg_genes <- arg_genes[, c(29,16)]
arg_gene_sum <- arg_genes %>% group_by(mag, `Resistance Mechanism`) %>% count()
arg_gene_sum <- pivot_wider(arg_gene_sum, names_from = `Resistance Mechanism`, values_from = n)
arg_gene_sum[is.na(arg_gene_sum)] <- 0
arg_gene_sum$total_hits <- rowSums(arg_gene_sum[, c(2:4)])

write.csv(arg_gene_sum, "refined data files/arg_hits.csv", row.names = F)

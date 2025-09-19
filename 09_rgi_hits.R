library(tidyverse)
setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

mag_list <- c("MAG_00097_1", "MAG_00103_1", "MAG_00110_1", "MAG_00179_1",
              "MAG_00194_1", "MAG_00197_1", "MAG_00201_1", "MAG_00674_1")

arg_files <- list.files("data files/refined_MAG_RGI/", recursive = T, pattern = ".*RGI.txt", full.names = T)

arg_genes <- data.frame()
for(i in 1:length(arg_files)){
  mag_args <- read_tsv(arg_files[i])
  mag_args$mag <- arg_files[i] %>% substr(29,39)
  arg_genes <- rbind(arg_genes, mag_args)
}

arg_genes$mag  <- arg_genes$mag  %>% str_remove("_[A-Z]")
mag_args <- subset(arg_genes, mag %in% mag_list)
mag_args <- mag_args[, c(2,9,15:17,29)]
colnames(mag_args)[1] <- "gene"
write.csv(mag_args, "data files/mag_arg_hits.csv", row.names = F)


arg_genes <- arg_genes[, c(29,16)]
arg_gene_sum <- arg_genes %>% group_by(mag, `Resistance Mechanism`) %>% count()
arg_gene_sum <- pivot_wider(arg_gene_sum, names_from = `Resistance Mechanism`, values_from = n)
arg_gene_sum[is.na(arg_gene_sum)] <- 0
arg_gene_sum$total_hits <- rowSums(arg_gene_sum[, c(2:4)])
write.csv(arg_gene_sum, "data files/arg_hits.csv", row.names = F)


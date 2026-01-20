library(tidyverse)
setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")


nonsyn_sig_genes <- read_csv("data files/nonsyn_allele_shifts_07_enriched_genes.csv")
syn_sig_genes <- read_csv("data files/syn_allele_shifts_07_enriched_genes.csv")

aro_nonsyn <- subset(nonsyn_sig_genes, str_detect(Preferred_name, pattern = "aro"))
colnames(aro_nonsyn)[3] <- "Nonsynonymous allele shifts"

aro_syn <- subset(syn_sig_genes, str_detect(Preferred_name, pattern = "aro"))
colnames(aro_syn)[3] <- "Synonymous allele shifts"

aro_genes <- full_join(aro_syn, aro_nonsyn)
aro_genes <- aro_genes[, c(1,6,4,8,3,9)]
aro_genes[is.na(aro_genes)] <- 0
colnames(aro_genes)[1] <- "MAG"
colnames(aro_genes)[2] <- "Gene"

mag_names <- read_csv("data files/mag_comp_contam.csv")
aro_genes <- left_join(aro_genes, mag_names)
aro_genes <- aro_genes[, c(9,1:6)]
write.csv(aro_genes, "data files/aro_gene_table.csv", row.names = F)

nonsyn_sig_genes$operon <- nonsyn_sig_genes$Preferred_name %>% substr(1,3)
nonsyn_operon_sum <- nonsyn_sig_genes %>% group_by(operon) %>% count()

syn_sig_genes$operon <- syn_sig_genes$Preferred_name %>% substr(1,3)
syn_operon_sum <- syn_sig_genes %>% group_by(operon) %>% count()



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

all_genes <- read_csv("data files/cog_background_snv_freq_genes.csv")
all_genes$operon <- all_genes$Preferred_name %>% substr(1,3)
all_genes <- subset(all_genes, !(operon == "-"))
all_genes$operon <- tolower(all_genes$operon)
all_operon_sum <- all_genes %>% group_by(operon) %>% count()
colnames(all_operon_sum)[2] <- "total"

operons <- left_join(all_operon_sum, nonsyn_operon_sum)
operons$ratio <- operons$n / operons$total
colnames(operons)[3] <- "Nonsynonymous allele shifts"
operons <- subset(operons, !(nchar(operon) != 3))
operons <- subset(operons, str_detect(operon, "\\d") == F)
operons[is.na(operons)] <- 0
write.csv(operons, "data files/operon_gene_sum.csv", row.names = F)

ggplot(operons, aes(x = total, y = `Nonsynonymous allele shifts`))+
  geom_point()+
  theme_bw()+
  xlim(0,120)+
  ylim(0,120)+
  geom_abline(slope = 1)

operons <- subset(operons, ratio > 0)
ggplot(operons, aes(x = reorder(operon, ratio), y = ratio))+
  geom_point()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 4))



#operons <- subset(operons, total > 5)



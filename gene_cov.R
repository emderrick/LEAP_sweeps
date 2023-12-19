library(tidyverse)
library(dplyr)

mag_list <- list("I4_MAG_00006", "I4_MAG_00065", "L2_MAG_00052", "L3_MAG_00058", "L4_MAG_00099",
                 "L7_MAG_00020", "L7_MAG_00028", "L7_MAG_00043", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")

mag_labs <- (c(I4_MAG_00006 = "SJAQ100 sp016735685", I4_MAG_00065 = "Roseomonas sp.", L2_MAG_00052 = "Erythrobacter sp.",
               L3_MAG_00058 = "Prosthecobacter sp.", L4_MAG_00099 = "Bosea sp001713455", L7_MAG_00020 = "Sphingorhabdus_B sp.",
               L7_MAG_00028 = "SYFN01 sp.", L7_MAG_00043 = "Luteolibacter sp.",L8_MAG_00011 = "UBA953 sp.",
               L8_MAG_00019 = "UA16", L8_MAG_00042 = "UBA4660 sp."))

gene_files <- list.files("subsampled_instrain/", recursive = T, pattern = ".*gene_info.tsv", full.names = T)
all_genes <- data.frame()
for(i in 1:length(gene_files)){
  pond_time_genes <- read.table(gene_files[i],sep = "\t",header = T)
  timepoint <- gsub(".*instrain_output/", "", gene_files[i]) %>% substr(9, 17)
  pond_time_genes <- cbind(pond_time_genes, timepoint = rep(timepoint, nrow(pond_time_genes)))
  all_genes <- rbind(all_genes, pond_time_genes)
}

all_genes$mag <- all_genes$gene %>% substr(1,12)
all_genes <- all_genes[, c(1:4, 21:22)]

all_mags <- read_csv("all_mags_subsamp.csv")
all_mags <- all_mags[, c("mag", "mag_coverage", "timepoint", "name", "new_time")]
all_genes_mag <- left_join(all_genes, all_mags, by = c("mag", "timepoint"))
all_genes_mag <- subset(all_genes_mag, new_time == 2)
all_genes_mag <- all_genes_mag %>% group_by(mag) %>% complete(gene, timepoint, fill = list(coverage = 0))
all_genes_mag <- all_genes_mag %>% group_by(mag, timepoint) %>% fill(mag_coverage, new_time, name, .direction = c("updown"))
all_genes_mag <- all_genes_mag %>% group_by(gene) %>% fill(scaffold, gene_length, .direction = c("updown"))
write.csv(all_genes_mag, "MAG_gene_info_subsamp.csv", row.names = F)

all_genes_mag$rel_cov <- all_genes_mag$coverage / all_genes_mag$mag_coverage
all_genes_mag <- all_genes_mag[, c("gene", "name", "rel_cov")]
all_genes_wide <- pivot_wider(all_genes_mag, names_from = "name", values_from = "rel_cov")
all_genes_wide$mean <- rowMeans(all_genes_wide[, c('Control A', 'Control B', 'Control C', 'Control D', 'Control E', 'GBH A', 'GBH B', 'GBH C', 'GBH D')], na.rm = T)

all_genes_long <- pivot_longer(all_genes_wide, cols = c('Control A', 'Control B', 'Control C', 'Control D', 'Control E', 'GBH A', 'GBH B', 'GBH C', 'GBH D'), names_to = "name", values_to = "rel_cov")
all_genes_long <- na.omit(all_genes_long)
all_genes_long$mag <- all_genes_long$gene %>% substr(1,12)
write.csv(all_genes_long, "gene_rel_cov_subsamp.csv",row.names = F)

all_genes_wide$control_mean <- rowMeans(all_genes_wide[, c('Control A', 'Control B', 'Control C', 'Control D', 'Control E')], na.rm = T)
all_genes_wide$GBH_mean <- rowMeans(all_genes_wide[,c('GBH A', 'GBH B', 'GBH C', 'GBH D')], na.rm = T)
all_genes_wide$cov_dif <- all_genes_wide$GBH_mean - all_genes_wide$control_mean
all_genes_wide$abs_val <- abs(all_genes_wide$cov_dif)
all_genes_wide <- all_genes_wide[-which(all_genes_wide[2:10] > 3, arr.ind = TRUE)[, 1], ]
write.csv(all_genes_wide, "gene_rel_cov_wide_subsamp.csv",row.names = F)

gene_changes <- subset(all_genes_wide, abs_val > 0.5)
gene_changes$control <- with(gene_changes, ifelse(control_mean < GBH_mean, 'low', 'high'))
gene_changes$control_ref <- with(gene_changes, ifelse(control == "high",
                                                      apply(gene_changes[,c('Control A', 'Control B', 'Control C', 'Control D', 'Control E')], 1, min, na.rm = T),
                                                      apply(gene_changes[,c('Control A', 'Control B', 'Control C', 'Control D', 'Control E')], 1, max, na.rm = T)))

gene_changes$GBH_ref <- with(gene_changes, ifelse(control == "high",
                                                  apply(gene_changes[,c('GBH A', 'GBH B', 'GBH C', 'GBH D')], 1, max, na.rm = T),
                                                  apply(gene_changes[,c('GBH A', 'GBH B', 'GBH C', 'GBH D')], 1, min, na.rm = T)))

gene_changes$pass <- with(gene_changes, ifelse(((control == "high" & GBH_ref < control_ref) | (control == "low" & GBH_ref > control_ref)), "yes", "no"))
gene_changes_pass <- subset(gene_changes, pass == "yes")

gene_changes_pass$mag <- gene_changes_pass$gene %>% substr(1,12)
write.csv(gene_changes_pass, "gene_coverage_sig_genes_subsamp.csv", row.names = F)

gene_decrease <- subset(gene_changes_pass, cov_dif < -0.5)
gene_increase <- subset(gene_changes_pass, cov_dif > 0.5)

gene_decrease_sum <- gene_decrease %>% group_by(mag) %>% summarize(cov_decrease_total = sum(pass == "yes"))
gene_increase_sum <- gene_increase %>% group_by(mag) %>% summarize(cov_increase_total = sum(pass == "yes"))

gene_sum <- left_join(gene_increase_sum, gene_decrease_sum)
write.csv(gene_sum, "gene_cov_change_sum_subsamp.csv", row.names = F)

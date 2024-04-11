library(tidyverse)
library(dplyr)

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
all_mags <- all_mags[, c("mag", "mag_coverage", "timepoint", "new_name", "new_time")]
all_genes_mag <- left_join(all_genes, all_mags, by = c("mag", "timepoint"))
all_genes_mag <- subset(all_genes_mag, new_name != "GBH C at T1")
all_genes_mag <- all_genes_mag %>% group_by(mag) %>% complete(gene, timepoint, fill = list(coverage = 0))
all_genes_mag <- all_genes_mag %>% group_by(mag, timepoint) %>% fill(mag_coverage, new_time, new_name, .direction = c("updown"))
all_genes_mag <- all_genes_mag %>% group_by(gene) %>% fill(scaffold, gene_length, .direction = c("updown"))
write.csv(all_genes_mag, "MAG_gene_info_subsamp.csv", row.names = F)

all_genes_mag$rel_cov <- all_genes_mag$coverage / all_genes_mag$mag_coverage
all_genes_mag <- all_genes_mag[, c("gene", "new_name", "rel_cov")]
all_genes_wide <- pivot_wider(all_genes_mag, names_from = "new_name", values_from = "rel_cov")

all_genes_wide$control_T2_mean <- rowMeans(all_genes_wide[, c('Control A at T2', 'Control B at T2', 'Control C at T2', 'Control D at T2', 'Control E at T2')], na.rm = T)
all_genes_wide$GBH_T2_mean <- rowMeans(all_genes_wide[,c('GBH A at T2', 'GBH B at T2', 'GBH C at T2', 'GBH D at T2')], na.rm = T)
all_genes_wide$T2_cov_dif <- all_genes_wide$GBH_T2_mean - all_genes_wide$control_T2_mean
all_genes_wide$T2_abs_val <- abs(all_genes_wide$T2_cov_dif)
all_genes_wide <- all_genes_wide[-which(all_genes_wide[2:10] > 3, arr.ind = TRUE)[, 1], ]
all_genes_wide$mag <- all_genes_wide$gene %>% substr(1,12)
write.csv(all_genes_wide, "gene_rel_cov_wide_subsamp.csv",row.names = F)

gene_changes <- subset(all_genes_wide, T2_abs_val > 0.5)
gene_changes$control <- with(gene_changes, ifelse(control_T2_mean < GBH_T2_mean, 'low', 'high'))
gene_changes$control_ref <- with(gene_changes, ifelse(control == "high",
                                                      apply(gene_changes[,c('Control A at T2', 'Control B at T2', 'Control C at T2', 'Control D at T2', 'Control E at T2')], 1, min, na.rm = T),
                                                      apply(gene_changes[,c('Control A at T2', 'Control B at T2', 'Control C at T2', 'Control D at T2', 'Control E at T2')], 1, max, na.rm = T)))

gene_changes$GBH_ref <- with(gene_changes, ifelse(control == "high",
                                                  apply(gene_changes[,c('GBH A at T2', 'GBH B at T2', 'GBH C at T2', 'GBH D at T2')], 1, max, na.rm = T),
                                                  apply(gene_changes[,c('GBH A at T2', 'GBH B at T2', 'GBH C at T2', 'GBH D at T2')], 1, min, na.rm = T)))

gene_changes$pass <- with(gene_changes, ifelse(((control == "high" & GBH_ref < control_ref) | (control == "low" & GBH_ref > control_ref)), "yes", "no"))
gene_changes_pass <- subset(gene_changes, pass == "yes")
write.csv(gene_changes_pass, "gene_coverage_sig_genes_subsamp.csv", row.names = F)

gene_decrease <- subset(gene_changes_pass, T2_cov_dif < -0.5)
gene_increase <- subset(gene_changes_pass, T2_cov_dif > 0.5)

gene_decrease_sum <- gene_decrease %>% group_by(mag) %>% summarize(cov_decrease_total = sum(pass == "yes"))
gene_increase_sum <- gene_increase %>% group_by(mag) %>% summarize(cov_increase_total = sum(pass == "yes"))
gene_sum <- left_join(gene_increase_sum, gene_decrease_sum)
write.csv(gene_sum, "gene_cov_change_sum_subsamp.csv", row.names = F)

all_genes_wide_L7 <- subset(all_genes_wide, mag == "L7_MAG_00020")
all_genes_wide_L7 <- all_genes_wide_L7[, c(19, 1, 4, 6:8, 11:14)]
all_genes_long_L7 <- pivot_longer(all_genes_wide_L7, cols = c(3:10), names_to = "new_name", values_to = "rel_cov")
all_genes_long_L7$time <- all_genes_long_L7$new_name %>% str_sub(-1)
all_genes_long_L7$pond <- all_genes_long_L7$new_name %>% str_sub(end = -7)
all_genes_long_L7 <- subset(all_genes_long_L7, select = -c(new_name))
all_genes_wide_L7 <- all_genes_long_L7 %>% group_by(pond) %>% pivot_wider(names_from = time, values_from = rel_cov)
all_genes_wide_L7$cov_dif <- all_genes_wide_L7$'2' - all_genes_wide_L7$'1'
all_genes_wide_L7$abs_val <- abs(all_genes_wide_L7$cov_dif)
all_genes_wide_L7 <- all_genes_wide_L7[-which(all_genes_wide_L7[4:5] > 3, arr.ind = TRUE)[, 1], ]
all_genes_long_L7 <- pivot_longer(all_genes_wide_L7, cols = c('1', '2'), names_to = "new_time", values_to = "rel_cov")
all_genes_long_L7$gene_pond <- paste(all_genes_long_L7$gene, all_genes_long_L7$pond, sep = "")
write.csv(all_genes_long_L7, "gene_cov_change_L7.csv", row.names = F)

gene_changes_L7 <- subset(all_genes_wide_L7, abs_val > 0.5)
gene_changes_L7_sum <- gene_changes_L7 %>% group_by(pond) %>% summarize(cov_decrease_total = sum(cov_dif < -0.5), cov_increase_total = sum(cov_dif > 0.5))
write.csv(gene_changes_L7_sum, "gene_cov_change_sum_L7.csv", row.names = F)

gene_decrease_L7 <- subset(gene_decrease, mag == "L7_MAG_00020")
gene_increase_L7 <- subset(gene_increase, mag == "L7_MAG_00020")
gene_changes_increase_L7 <- subset(gene_changes_L7, cov_dif > 0.5 & pond == "GBH A")
gene_changes_decrease_L7 <- subset(gene_changes_L7, cov_dif < -0.5 & pond == "Control D")
L7_increase_overlap <- as.data.frame(intersect(gene_changes_increase_L7$gene, gene_increase_L7$gene)) 
L7_decrease_overlap <- as.data.frame(intersect(gene_changes_decrease_L7$gene, gene_decrease_L7$gene)) 

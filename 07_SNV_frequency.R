library(tidyverse)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

mag_list <- c("MAG_00097_1", "MAG_00110_1", "MAG_00179_1",
              "MAG_00194_1", "MAG_00197_1", "MAG_00201_1","MAG_00674_1")

all_mag_SNVs <- read_csv("data files/MAG_SNV_depth_info.csv")
all_mag_SNVs <- subset(all_mag_SNVs, mag %in% mag_list)
all_mag_SNVs <- subset(all_mag_SNVs, mag_coverage >= 5 & mag_breadth >= 0.5)

all_mag_SNVs$time <- with(all_mag_SNVs, ifelse(time == "1", "Day 0", "Day 28"))
all_mag_SNVs$pond_time <- paste(all_mag_SNVs$time, all_mag_SNVs$Name, sep = " ")
mag_snvs <- all_mag_SNVs[, c("gene", "pond_time", "mag", "group", "mutation_type", "final_ref_freq")]
postions_keep <- mag_snvs %>% group_by(group) %>% count(mutation_type)
postions_keep <- spread(postions_keep, mutation_type, n)
postions_keep$na_count <- rowSums(is.na(postions_keep[, c("I", "M", "N", "S")]))
postions_keep <- subset(postions_keep, na_count >= 3)
mag_snvs <- subset(mag_snvs, group %in% postions_keep$group)
mag_snvs <- mag_snvs %>% group_by(group) %>% fill(mutation_type, .direction = "updown") %>% ungroup()
mag_snvs_wide <- spread(mag_snvs, pond_time, final_ref_freq)

mag_snvs_wide$CTRL_0 <-rowMeans(mag_snvs_wide[, grep("Day 0 CTRL", colnames(mag_snvs_wide))], na.rm = T)
mag_snvs_wide$GBH_0 <- rowMeans(mag_snvs_wide[, grep("Day 0 GBH", colnames(mag_snvs_wide))], na.rm = T)
mag_snvs_wide$CTRL_28 <- rowMeans(mag_snvs_wide[, grep("Day 28 CTRL", colnames(mag_snvs_wide))], na.rm = T)
mag_snvs_wide$GBH_28 <- rowMeans(mag_snvs_wide[, grep("Day 28 GBH", colnames(mag_snvs_wide))], na.rm = T)
mag_snvs_wide <- mag_snvs_wide %>% mutate_all(function(x) ifelse(is.nan(x), NA, x))
mag_snvs_wide <- mag_snvs_wide[, c("mag", "gene", "group", "mutation_type", "CTRL_0", "GBH_0", "CTRL_28", "GBH_28")]
write.csv(mag_snvs_wide, "data files/all_snv_frequency.csv", row.names = F)

snv_frequency <- mag_snvs_wide[complete.cases(mag_snvs_wide), ]
snv_frequency$GBH_change <- snv_frequency$GBH_28 - snv_frequency$GBH_0
snv_frequency$CTRL_change <- snv_frequency$CTRL_28 - snv_frequency$CTRL_0
snv_frequency$GBH_CTRL_change <- snv_frequency$GBH_change - snv_frequency$CTRL_change
snv_frequency$GBH_CTRL_T1 <- snv_frequency$GBH_0 - snv_frequency$CTRL_0
snv_frequency$GBH_CTRL_T2 <- snv_frequency$GBH_28 - snv_frequency$CTRL_28

#snv_frequency$GBH_CTRL_change_abs <- abs(snv_frequency$GBH_CTRL_change)
#snv_frequency$GBH_change_abs <- abs(snv_frequency$GBH_change)
#snv_frequency$CTRL_change_abs <- abs(snv_frequency$CTRL_change)
snv_frequency$GBH_CTRL_T1_abs <- abs(snv_frequency$GBH_CTRL_T1)
#snv_frequency$GBH_CTRL_T2_abs <- abs(snv_frequency$GBH_CTRL_T2)

snv_frequency$case <- with(snv_frequency, ifelse((GBH_CTRL_T1_abs <= 0.3 & GBH_change <= -0.7 & GBH_CTRL_T2 <= -0.7), "case_a", "no"))
snv_frequency$case <- with(snv_frequency, ifelse((GBH_CTRL_T1_abs <= 0.3 & CTRL_change <= -0.7 & GBH_CTRL_T2 >= 0.7), "case_b", case))
write.csv(snv_frequency, "data files/snv_frequency_changes_07.csv", row.names = F)


mag_genes <- read_csv("data files/MAG_gene_info.csv")
mag_genes$rel_cov <- mag_genes$gene_coverage / mag_genes$mag_coverage
mag_genes <- subset(mag_genes, rel_cov <= 1.2 & rel_cov >= 0.6)

eggnog_genes <- read_tsv("data files/eggnog_genes.emapper.annotations", skip = 4)
colnames(eggnog_genes)[1]="gene"
eggnog_genes$COG_ID<- str_extract(eggnog_genes$eggNOG_OGs, "COG\\d{4}")
eggnog_genes <- eggnog_genes[, c("gene", "COG_ID", "Description", "Preferred_name")]
eggnog_genes$scaffold <- eggnog_genes$gene %>% str_extract("[^_]*_[^_]*")
background_cog_snv <- subset(eggnog_genes, gene %in% mag_genes$gene)
background_cog_snv <- subset(background_cog_snv, is.na(COG_ID) == F)
write.csv(background_cog_snv, "data files/cog_background_snv_freq_genes.csv", row.names = F)

sig_a_changes <- subset(snv_frequency, case == "case_a")
sig_a_nonsyn_changes <- subset(sig_a_changes, mutation_type == "N")
sig_a_syn_changes <- subset(sig_a_changes, mutation_type == "S")

sig_a_snvs_sum <- sig_a_changes %>% group_by(mag, gene) %>% count()
sig_a_snvs_sum <- subset(sig_a_snvs_sum, !str_detect(gene, ","))

sig_a_nonsyn_sum <- sig_a_nonsyn_changes %>% group_by(mag, gene) %>% count()
sig_a_syn_sum <- sig_a_syn_changes %>% group_by(mag, gene) %>% count()

significant_a_genes <- left_join(sig_a_snvs_sum, background_cog_snv) 
write.csv(significant_a_genes, "data files/allele_shifts_a_significant_genes_07.csv", row.names = F)
sum_significant_a_genes <- significant_a_genes %>% group_by(mag) %>% count()

significant_a_nonsyn_genes <- left_join(sig_a_nonsyn_sum, background_cog_snv) 
write.csv(significant_a_nonsyn_genes, "data files/nonsyn_allele_shifts_a_significant_genes_07.csv", row.names = F)

significant_a_syn_genes <- left_join(sig_a_syn_sum, background_cog_snv) 
write.csv(significant_a_syn_genes, "data files/syn_allele_shifts_a_significant_genes_07.csv", row.names = F)

sig_b_changes <- subset(snv_frequency, case == "case_b")
sig_b_nonsyn_changes <- subset(sig_b_changes, mutation_type == "N")
sig_b_syn_changes <- subset(sig_b_changes, mutation_type == "S")

sig_b_snvs_sum <- sig_b_changes %>% group_by(mag, gene) %>% count()
sig_b_snvs_sum <- subset(sig_b_snvs_sum, !str_detect(gene, ","))

sig_b_nonsyn_sum <- sig_b_nonsyn_changes %>% group_by(mag, gene) %>% count()
sig_b_syn_sum <- sig_b_syn_changes %>% group_by(mag, gene) %>% count()

significant_b_genes <- left_join(sig_b_snvs_sum, background_cog_snv) 
write.csv(significant_b_genes, "data files/allele_shifts_b_significant_genes_07.csv", row.names = F)

significant_b_nonsyn_genes <- left_join(sig_b_nonsyn_sum, background_cog_snv) 
write.csv(significant_b_nonsyn_genes, "data files/nonsyn_allele_shifts_b_significant_genes_07.csv", row.names = F)

significant_b_syn_genes <- left_join(sig_b_syn_sum, background_cog_snv)
write.csv(significant_b_syn_genes, "data files/syn_allele_shifts_b_significant_genes_07.csv", row.names = F)


# snv_frequency_plot <- pivot_longer(snv_frequency, cols = c(5:8), values_to = "avg_freq", names_to = "Treatment_Time")
# snv_frequency_plot$Treatment <- snv_frequency_plot$Treatment_Time %>% substr(1,4)
# snv_frequency_plot$Treatment <- snv_frequency_plot$Treatment %>% str_remove("_")
# snv_frequency_plot$Treatment_group <- paste(snv_frequency_plot$Treatment, snv_frequency_plot$group)
# snv_frequency_plot$Time <- snv_frequency_plot$Treatment_Time %>% str_sub(-2) %>% str_remove("_")
# snv_frequency_plot$Time <- paste("Day ", snv_frequency_plot$Time, sep = "")

# snv_frequency_plot_a <- subset(snv_frequency_plot, case == "case_a")
# plot_a <- ggplot(snv_frequency_plot_a, aes(x = Time, y = avg_freq, colour = Treatment))+
#   geom_smooth(aes(group = Treatment), method = "lm")+
#   #geom_point()+
#   #geom_line(aes(group = Treatment_group))+
#   scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
#   scale_x_discrete(expand = c(0, 0))+
#   ylim(0,1)+
#   labs(y = "Reference Frequency")+
#   theme_bw()
# ggsave("figures/snv_freq_case_a_plot.pdf", plot_a, units = "cm", width = 10, height = 6)
# 
# snv_frequency_plot_b <- subset(snv_frequency_plot, case == "case_b")
# plot_b <- ggplot(snv_frequency_plot_b, aes(x = Time, y = avg_freq, colour = Treatment))+
#   geom_smooth(aes(group = Treatment), method = "lm")+
#   #geom_point()+
#   #geom_line(aes(group = Treatment_group))+
#   scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
#   scale_x_discrete(expand = c(0, 0))+
#   ylim(0,1)+
#   labs(y = "Reference Frequency")+
#   theme_bw()
# ggsave("figures/snv_freq_case_b_plot.pdf", plot_b, units = "cm", width = 10, height = 6)


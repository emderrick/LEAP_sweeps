library(tidyverse)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)


# for Emma data
mag_list <- list("I4_MAG_00006", "I4_MAG_00065", "L2_MAG_00052", "L3_MAG_00058", "L4_MAG_00099",
                 "L7_MAG_00020", "L7_MAG_00028", "L7_MAG_00043", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")

mag_labs <- (c(I4_MAG_00006 = "SJAQ100 sp016735685", I4_MAG_00065 = "Roseomonas sp.", L2_MAG_00052 = "Erythrobacter sp.",
               L3_MAG_00058 = "Prosthecobacter sp.", L4_MAG_00099 = "Bosea sp001713455", L7_MAG_00020 = "Sphingorhabdus_B sp.",
               L7_MAG_00028 = "SYFN01 sp.", L7_MAG_00043 = "Luteolibacter sp.",L8_MAG_00011 = "UBA953 sp.",
               L8_MAG_00019 = "UA16", L8_MAG_00042 = "UBA4660 sp."))


gene_files <- list.files("95_profiles/",recursive = T, pattern=".*gene_info.tsv",full.names = T)
all_genes <- data.frame()
for(i in 1:length(gene_files)){
  pond_time_genes <- read.table(gene_files[i],sep="\t",header=T)
  timepoint <- gsub(".*profile_output/", "", gene_files[i]) %>% substr(1,9)
  pond_time_genes <- cbind(pond_time_genes,timepoint=rep(timepoint,nrow(pond_time_genes)))
  all_genes <- rbind(all_genes,pond_time_genes)
}

all_genes$mag <- all_genes$gene %>% substr(1,12)
all_genes <- subset(all_genes, mag %in% mag_list)
all_genes <- all_genes[, c(1:4, 21:22)]

all_mags <- read_csv("ANI_95_all_mags.csv")
all_mags <- subset(all_mags, mag %in% mag_list)
all_mags <- all_mags[, c("mag", "mag_coverage", "timepoint", "name", "new_time")]
all_genes_mag <- left_join(all_genes, all_mags, by = c("mag", "timepoint"))
all_genes_mag <- subset(all_genes_mag, new_time == 2)
all_genes_mag <- all_genes_mag %>% group_by(mag) %>% complete(gene, timepoint, fill = list(coverage = 0))
all_genes_mag <- all_genes_mag %>% group_by(mag, timepoint) %>% fill(mag_coverage, new_time, name, .direction = c("updown"))
all_genes_mag$rel_cov <- all_genes_mag$coverage / all_genes_mag$mag_coverage
all_genes_graph <- all_genes_mag[, c("gene", "name", "rel_cov")]
all_genes_wide <- pivot_wider(all_genes_graph, names_from = "name", values_from = "rel_cov")
all_genes_wide$mean <- rowMeans(all_genes_wide[,c('Control A', 'Control B', 'Control C', 'Control D', 'Control E', 'GBH A', 'GBH B', 'GBH C', 'GBH D')], na.rm = T)

all_genes_long <- pivot_longer(all_genes_wide, cols = c('Control A', 'Control B', 'Control C', 'Control D', 'Control E', 'GBH A', 'GBH B', 'GBH C', 'GBH D'), names_to = "name", values_to = "rel_cov")
all_genes_long <- na.omit(all_genes_long)
all_genes_long$mag <- all_genes_long$gene %>% substr(1,12)
write.csv(all_genes_long, "gene_rel_cov.csv",row.names = F)

all_genes_wide$control_mean <- rowMeans(all_genes_wide[, c('Control A', 'Control B', 'Control C', 'Control D', 'Control E')], na.rm = T)
all_genes_wide$GBH_mean <- rowMeans(all_genes_wide[,c('GBH A', 'GBH B', 'GBH C', 'GBH D')], na.rm = T)
all_genes_wide$cov_dif <- all_genes_wide$GBH_mean - all_genes_wide$control_mean
all_genes_wide$abs_val <- abs(all_genes_wide$cov_dif)
all_genes_wide <- all_genes_wide[-which(all_genes_wide[2:10] > 3, arr.ind = TRUE)[, 1], ]
write.csv(all_genes_wide, "gene_rel_cov_wide.csv",row.names = F)

all_genes_line <- pivot_longer(all_genes_wide, cols = c('control_mean', 'GBH_mean'), names_to = "treatment", values_to = "average_copies")
all_genes_line$mag <- all_genes_line$gene %>% substr(1,12)
all_genes_line$mag_order = factor(all_genes_line$mag, levels=c('I4_MAG_00006', 'I4_MAG_00065', 'L3_MAG_00058', 'L7_MAG_00020', 'L8_MAG_00011', 'L8_MAG_00019',
                                                             'L2_MAG_00052', 'L4_MAG_00099',  'L7_MAG_00028', 'L7_MAG_00043', 'L8_MAG_00042'))

gene_copy <- ggplot(all_genes_line %>% arrange(abs_val), aes(x = treatment, y = average_copies, colour = cov_dif, group = reorder(gene, abs_val)))+ 
  geom_line()+
  geom_point()+
  #scale_colour_gradientn(colours = c("#CC0000", "#FF3333", "#FF7070", "#FFADAD", "#e7dada","#F5F5F5","#EBEEFF", "#C2CDFF", "#708AFF", "#3358FF","#0022E0"), limits = c(-3,3))+
  scale_colour_gradientn(colours = brewer.pal(11, "RdBu"), limits = c(-3,3))+
  labs(y ="Gene Copies", colour = NULL) +
  theme_classic()+
  theme(text = element_text(size = 30, colour = 'black'),
        axis.text = element_text(colour = "black"),
        axis.title = element_text(size = 35), 
        axis.ticks.length = unit(.4, "cm"),
        axis.ticks = element_line(linewidth = 1.5, colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.line = element_line(linewidth = 1),
        strip.text.x.top = element_text(size = 30, margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), face = "bold"),
        strip.background = element_rect(linewidth = 2),
        legend.position = ("bottom"),
        legend.margin =  margin(t = 60, r = 0, b = 0, l = 0),
        legend.key.size = unit(2, 'cm'),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 35),
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.spacing.y = unit(5, "lines"))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag_order, nrow = 2, ncol = 6, scales = "free", labeller = labeller(mag_order = mag_labs))
 
save_plot("gene_copy_dif.jpeg", gene_copy, base_height = 7.5, base_width = 6, ncol = 6, nrow = 2, dpi = 300, limitsize = F)


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
write.csv(gene_changes_pass, "gene_coverage_sig_genes.csv", row.names = F)
gene_decrease <- subset(gene_changes_pass, cov_dif < -0.5)
gene_increase <- subset(gene_changes_pass, cov_dif > 0.5)

background_cog <- read_csv("cog_background_genes.csv")
gene_cov_significant <- left_join(gene_changes_pass[, c("gene", "abs_val")], background_cog)
gene_cov_sig_increase <- left_join(gene_increase[, c("gene", "abs_val")], background_cog)
gene_cov_sig_decrease <- left_join(gene_decrease[, c("gene", "abs_val")], background_cog)
write.csv(gene_cov_significant, "gene_cov_significant.csv", row.names = F)
write.csv(gene_cov_sig_increase, "gene_cov_sig_increase.csv", row.names = F)
write.csv(gene_cov_sig_decrease, "gene_cov_sig_decrease.csv", row.names = F)


library(tidyverse)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

all_genes_wide <- read_csv("gene_rel_cov_wide_subsamp.csv")
all_genes_line <- pivot_longer(all_genes_wide, cols = c('control_mean', 'GBH_mean'), names_to = "treatment", values_to = "average_copies")
all_genes_line$mag = factor(all_genes_line$gene %>% substr(1,12), levels=c('L2_MAG_00052', 'L4_MAG_00099', 'L7_MAG_00020', 'L7_MAG_00028', 'L7_MAG_00043', 'L8_MAG_00042', 'I4_MAG_00006', 'I4_MAG_00065', 'L3_MAG_00058', 'L8_MAG_00011', 'L8_MAG_00019'))
                                                                      

mag_labs <- (c(I4_MAG_00006 = "SJAQ100 sp016735685", I4_MAG_00065 = "Roseomonas sp.", L2_MAG_00052 = "Erythrobacter sp.",
               L3_MAG_00058 = "Prosthecobacter sp.", L4_MAG_00099 = "Bosea sp001713455", L7_MAG_00020 = "Sphingorhabdus_B sp.",
               L7_MAG_00028 = "SYFN01 sp.", L7_MAG_00043 = "Luteolibacter sp.",L8_MAG_00011 = "UBA953 sp.",
               L8_MAG_00019 = "UA16", L8_MAG_00042 = "UBA4660 sp."))

gene_copy <- ggplot(all_genes_line, aes(x = treatment, y = average_copies, colour = cov_dif, group = reorder(gene, abs_val)))+ 
  geom_line()+
  geom_point()+
  scale_colour_gradientn(colours = brewer.pal(11, "RdBu"), limits = c(-3,3))+
  labs(y ="Gene Copies", colour = "Gene Copy Diff.") +
  theme_classic()+
  theme(text = element_text(size = 16, colour = 'black'),
        axis.text = element_text(size = 16, colour = "black"),
        axis.title = element_text(size = 16), 
        axis.ticks.length = unit(.4, "cm"),
        axis.ticks = element_line(linewidth = 1, colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.line = element_line(linewidth = 1),
        strip.text.x.top = element_text(size = 18, margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"), face = "bold"),
        strip.background = element_rect(linewidth = 2),
        legend.position = ("bottom"),
        legend.margin =  margin(t = 30, r = 0, b = 0, l = 0),
        legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.spacing.y = unit(3, "lines"))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 3))+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag, nrow = 2, ncol = 6, scales = "free", as.table = F, labeller = labeller(mag = mag_labs))
save_plot("gene_copy_dif.jpeg", gene_copy, base_height = 5, base_width = 4.5, ncol = 6, nrow = 2, dpi = 300, limitsize = F)


#looking at L7_MAG_00020 through time
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
all_genes_L7 <- subset(all_genes_mag, mag == "L7_MAG_00020")

all_genes_L7 <- all_genes_L7 %>% complete(gene, timepoint, fill = list(coverage = 0))
all_genes_L7 <- all_genes_L7 %>% group_by(timepoint) %>% fill(mag_coverage, new_time, name, .direction = c("updown"))
all_genes_L7$rel_cov <- all_genes_L7$coverage / all_genes_L7$mag_coverage
all_genes_L7 <- all_genes_L7[, c("gene", "name", "new_time", "rel_cov")]
all_genes_L7_wide <- pivot_wider(all_genes_L7, names_from = "new_time", values_from = "rel_cov")
all_genes_L7_wide$cov_dif <- all_genes_L7_wide$'2' - all_genes_L7_wide$'1'
all_genes_L7_wide$abs_val <- abs(all_genes_L7_wide$cov_dif)
all_genes_L7_wide <- all_genes_L7_wide[-which(all_genes_L7_wide[3:4] > 3, arr.ind = TRUE)[, 1], ]
all_genes_L7_long <- pivot_longer(all_genes_L7_wide, cols = c('1', '2'), names_to = "new_time", values_to = "rel_cov")
all_genes_L7_long$gene_pond <- paste(all_genes_L7_long$gene, all_genes_L7_long$name, sep = "")

gene_copy_L7 <- ggplot(all_genes_L7_long, aes(x = new_time, y = rel_cov, colour = cov_dif, group = reorder(gene_pond, abs_val)))+ 
  geom_line()+
  scale_colour_gradientn(colours = brewer.pal(11, "RdBu"), limits = c(-3,3))+
  labs(y ="Gene Copies", colour = "Gene Copy Diff.") +
  theme_classic()+
  theme(text = element_text(size = 14, colour = 'black'),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 14), 
        axis.ticks.length = unit(.4, "cm"),
        axis.ticks = element_line(linewidth = 1, colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 10, r = 20, b = 0, l = 0)),
        axis.line = element_line(linewidth = 1),
        strip.text.x.top = element_text(size = 16, margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"), face = "bold"),
        strip.background = element_rect(linewidth = 2),
        legend.margin =  margin(t = 0, r = 0, b = 0, l = 20),
        legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.margin = unit(c(1, 1, 1, 1), "cm"))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 3))+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~name, nrow = 1, ncol = 4, scales = "free")
save_plot("gene_copy_dif_L7.jpeg", gene_copy_L7, base_height = 4, base_width = 4.5, ncol = 4, nrow = 1, dpi = 300, limitsize = F)

library(tidyverse)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# for Emma data
mag_list <- list("I4_MAG_00006", "I4_MAG_00065", "L2_MAG_00052", "L3_MAG_00058", "L4_MAG_00099",
                 "L7_MAG_00020", "L7_MAG_00028", "L7_MAG_00043", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")

mag_labs <- (c(I4_MAG_00006 = "SJAQ100 sp016735685 assembled from Control B", I4_MAG_00065 = "Roseomonas sp. assembled from Control B", L2_MAG_00052 = "Erythrobacter sp. assembled from GBH A",
               L3_MAG_00058 = "Prosthecobacter sp. assembled from Control C", L4_MAG_00099 = "Bosea sp001713455 assembled from Control D", L7_MAG_00020 = "Sphingorhabdus_B sp. assembled from GBH C",
               L7_MAG_00028 = "SYFN01 sp. assembled from GBH C", L7_MAG_00043 = "Luteolibacter sp. assembled from GBH C",L8_MAG_00011 = "UBA953 sp. assembled from GBH D",
               L8_MAG_00019 = "UA16 family assembled from GBH D", L8_MAG_00042 = "UBA4660 sp. assembled from GBH D"))


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
all_genes_mag <- subset(all_genes_mag, coverage <= 3 * mag_coverage)
all_genes_mag$rel_cov <- all_genes_mag$coverage / all_genes_mag$mag_coverage
all_genes_graph <- all_genes_mag[, c("gene", "name", "rel_cov")]
all_genes_wide <- pivot_wider(all_genes_graph, names_from = "name", values_from = "rel_cov")
all_genes_wide$mean <- rowMeans(all_genes_wide[,c('Control A', 'Control B', 'Control C', 'Control D', 'Control E', 'GBH A', 'GBH B', 'GBH C', 'GBH D')], na.rm = T)
all_genes_wide <- subset(all_genes_wide, mean != 0) #one gene has a mean of 0 because the only pond with coverage of that gene was excluded because the coverage was > 3x the mean, so that gene isn't in any ponds 
all_genes_long <- pivot_longer(all_genes_wide, cols = c('Control A', 'Control B', 'Control C', 'Control D', 'Control E', 'GBH A', 'GBH B', 'GBH C', 'GBH D'), names_to = "name", values_to = "rel_cov")
all_genes_long <- na.omit(all_genes_long)
all_genes_long$mag <- all_genes_long$gene %>% substr(1,12)
write.csv(all_genes_long, "gene_rel_cov.csv",row.names = F)

gene_heatmap <- ggplot(all_genes_long, aes(x = name, y = reorder(gene, mean), fill = rel_cov)) +
  geom_tile()+
  scale_fill_gradientn(colours = brewer.pal(11, "PiYG"), limits = c(0,3))+
  theme_classic() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), 
        text = element_text(size = 17, face = "bold"), strip.text.x.top = element_text(size = 20))+
  labs(legend = "Gene Copies")+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag, nrow = 4, ncol = 3, scales = "free", labeller = labeller(mag = mag_labs))
ggsave("gene_heatmap.png", limitsize = F, dpi = 200, width = 32, height = 32)

all_genes_wide$control_mean <- rowMeans(all_genes_wide[, c('Control A', 'Control B', 'Control C', 'Control D', 'Control E')], na.rm = T)
all_genes_wide$GBH_mean <- rowMeans(all_genes_wide[,c('GBH A', 'GBH B', 'GBH C', 'GBH D')], na.rm = T)
all_genes_wide$abs_val <- abs(all_genes_wide$control_mean - all_genes_wide$GBH_mean)
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

library(tidyverse)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(ggh4x)

sens_mags <- list("I4_MAG_00006", "L7_MAG_00028", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")
res_mags <- list("I4_MAG_00065", "L3_MAG_00058", "L7_MAG_00020", "L7_MAG_00043")
unclass_mags <- list("L2_MAG_00052", "L4_MAG_00099")

mag_labs <- (c(I4_MAG_00006 = "SJAQ100 sp016735685", I4_MAG_00065 = "Roseomonas sp.", L2_MAG_00052 = "Erythrobacter sp.",
               L3_MAG_00058 = "Prosthecobacter sp.", L4_MAG_00099 = "Bosea sp001713455", L7_MAG_00020 = "Sphingorhabdus_B sp.",
               L7_MAG_00028 = "SYFN01 sp.", L7_MAG_00043 = "Luteolibacter sp.",L8_MAG_00011 = "UBA953 sp.",
               L8_MAG_00019 = "UA16", L8_MAG_00042 = "UBA4660 sp."))

all_genes_wide <- read_csv("gene_rel_cov_wide_subsamp.csv")
all_genes_line <- pivot_longer(all_genes_wide, cols = c('control_T2_mean', 'GBH_T2_mean'), names_to = "treatment", values_to = "average_copies")

all_genes_sens <- subset(all_genes_line, mag %in% sens_mags)
all_genes_res <- subset(all_genes_line, mag %in% res_mags)
all_genes_unclass <- subset(all_genes_line, mag %in% unclass_mags)

gene_copy_sens <- ggplot(all_genes_sens, aes(x = treatment, y = average_copies, colour = T2_cov_dif, group = reorder(gene, T2_abs_val)))+ 
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
        legend.position = ("none"),
        panel.spacing.y = unit(3, "lines"))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 3))+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag, nrow = 1, ncol = 5, scales = "free", labeller = labeller(mag = mag_labs))+
  force_panelsizes(cols = unit(4, "in"))

save_plot("gene_copy_dif_sens.jpeg", gene_copy_sens, base_height = 4, base_width = 24, dpi = 300, limitsize = F)

gene_copy_res <- ggplot(all_genes_res, aes(x = treatment, y = average_copies, colour = T2_cov_dif, group = reorder(gene, T2_abs_val)))+ 
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
        panel.spacing.y = unit(3, "lines"))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 3))+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag, nrow = 1, ncol = 4, scales = "free", labeller = labeller(mag = mag_labs))+
  force_panelsizes(cols = unit(4, "in"))
save_plot("gene_copy_dif_res.jpeg", gene_copy_res, base_height = 4, base_width = 22, dpi = 300, limitsize = F)

gene_copy_unclass <- ggplot(all_genes_unclass, aes(x = treatment, y = average_copies, colour = T2_cov_dif, group = reorder(gene, T2_abs_val)))+ 
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
        legend.margin =  margin(t = 0, r = 20, b = 0, l = 0),
        legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        panel.spacing.y = unit(3, "lines"))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 3))+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~mag, nrow = 1, ncol = 2, scales = "free", labeller = labeller(mag = mag_labs))+
  force_panelsizes(cols = unit(4, "in"))
save_plot("gene_copy_dif_unclass.jpeg", gene_copy_unclass, base_height = 4, base_width = 18, dpi = 300, limitsize = F)

#looking at L7_MAG_00020 through time
all_genes_long_L7 <- read_csv("gene_cov_change_L7.csv")

gene_copy_L7 <- ggplot(all_genes_long_L7, aes(x = new_time, y = rel_cov, colour = cov_dif, group = reorder(gene_pond, abs_val)))+ 
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
        panel.spacing.y = unit(3, "lines"))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 3))+
  scale_x_discrete(expand = c(0, 0))+
  facet_wrap(~pond, nrow = 1, ncol = 4, scales = "free", labeller = labeller(mag = mag_labs))+
  force_panelsizes(cols = unit(4, "in"))
save_plot("gene_copy_dif_L7.jpeg", gene_copy_L7, base_height = 4, base_width = 22, dpi = 300, limitsize = F)

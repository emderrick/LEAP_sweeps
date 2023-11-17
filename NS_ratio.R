library(tidyverse)
library(dplyr)
library(ggplot2)
library(cowplot)

mag_labs <- (c(I4_MAG_00006 = "SJAQ100 sp016735685", I4_MAG_00065 = "Roseomonas sp.", L2_MAG_00052 = "Erythrobacter sp.",
               L3_MAG_00058 = "Prosthecobacter sp.", L4_MAG_00099 = "Bosea sp001713455", L7_MAG_00020 = "Sphingorhabdus_B sp.",
               L7_MAG_00028 = "SYFN01 sp.", L7_MAG_00043 = "Luteolibacter sp.",L8_MAG_00011 = "UBA953 sp.",
               L8_MAG_00019 = "UA16", L8_MAG_00042 = "UBA4660 sp."))

select_mags <- list("I4_MAG_00006", "I4_MAG_00065", "L3_MAG_00058", "L7_MAG_00020", "L8_MAG_00011", "L8_MAG_00019")
EPSPS_class_1 <- list("I4_MAG_00006", "L7_MAG_00028", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")
EPSPS_unclass <- list("L2_MAG_00052", "L4_MAG_00099", "L7_MAG_00020")

all_MAG_snvs <- read_csv("all_MAG_SNVs_med_Aug15.csv")
all_MAG_snvs <- subset(all_MAG_snvs, str_detect(new_name, "T2"))
MAG_snvs <- all_MAG_snvs[, c('mag', 'gene', 'mutation_type', 'new_name', 'mag_length')] %>% subset(is.na(mutation_type) == F)

MAG_NS <- MAG_snvs %>% group_by(mag, mag_length, new_name) %>% summarize(N = sum(mutation_type=="N"), S = sum(mutation_type=="S"), I = sum(mutation_type=="I"), M = sum(mutation_type=="M"))
MAG_NS$total <- MAG_NS$N + MAG_NS$S + MAG_NS$I + MAG_NS$M
MAG_NS$graph_name <- str_sub(MAG_NS$new_name, end = -6) 
MAG_NS$graph_name <- gsub('Control', 'CTL', MAG_NS$graph_name)

# MAG_NS_long <- pivot_longer(MAG_NS, cols = c("N", "S", "I", "M"), names_to = "mutation_type", values_to = "count")
# 
# ggplot(MAG_NS_long, aes(x = graph_name, y = count/total, fill = mutation_type))+
#   geom_bar(stat = "identity")+
#   theme_classic()+
#   labs(y = "Fraction of SNVs", x = "pond", legend = "Mutation Type")+
#   facet_wrap(~mag, nrow = 4, ncol = 4, scales="free", labeller = labeller(mag = mag_labs))
# save_plot("NS_MAG_plot.jpeg", NS_plot, ncol = 4, nrow = 4, dpi = 300)

MAG_NS$NS_ratio <-  MAG_NS$N / MAG_NS$S
MAG_NS$SNVs_MBp <- (MAG_NS$total / MAG_NS$mag_length) * 1000000
MAG_NS$treatment <- str_sub(MAG_NS$new_name, end = -8)
MAG_NS_wide <- MAG_NS[, c('mag', 'new_name', 'SNVs_MBp')] %>% pivot_wider(names_from = new_name, values_from = SNVs_MBp)
MAG_NS_wide <- MAG_NS_wide %>% rowwise() %>% mutate(max_SNVs_mag = max(c_across(contains("T")), na.rm = T))
MAG_NS_max <- pivot_longer(MAG_NS_wide, cols = contains("T2"), names_to = 'new_name', values_to = 'SNVs_MBp') %>% subset(is.na(SNVs_MBp) == F) 
MAG_NS_max <- right_join(MAG_NS, MAG_NS_max, by = c('mag', 'new_name', 'SNVs_MBp'))
MAG_NS_max$sweep <- with(MAG_NS_max, ifelse(mag %in% select_mags, "Yes", "No"))
MAG_NS_plot <- subset(MAG_NS_max, SNVs_MBp == max_SNVs_mag)
MAG_NS_plot$EPSPS_class <- with(MAG_NS_plot, ifelse(mag %in% EPSPS_class_1, "Class I - Sensitive", "Class II - Resistant"))
MAG_NS_plot$EPSPS_class <- with(MAG_NS_plot, ifelse(mag %in% EPSPS_unclass, "Unclassified", EPSPS_class))

NS_SNV_EPSPS <- ggplot(MAG_NS_plot, aes(x = NS_ratio, y = SNVs_MBp, colour = EPSPS_class))+
  geom_point(aes(shape=sweep), size=2.5)+
  scale_colour_manual(values = c("#0A9396", "#EE9B00", "#AE2012"))+
  labs(y = "Variable sites in population per Mbp", x = "N:S ratio", colour = "EPSPS Class", shape = "Potential Sweep")+
  theme_classic()+
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(colour = "black", size = 8), 
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        legend.position = c(0.85, 0.65))
save_plot("NS_SNV_max_EPSPS_plot.jpeg", NS_SNV_EPSPS, dpi = 300, base_height = 4, base_width = 6.7)


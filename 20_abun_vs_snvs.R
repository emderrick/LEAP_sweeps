library(tidyverse)
library(rstatix)
library(ggpubr)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")
# 
# all_MAG_SNVs <- read_csv("data files/snvs_for_depth.csv")
# all_MAG_SNVs$time <- as.character(all_MAG_SNVs$time)
# 
# sample_names <- read_csv("data files/chapter_1_sample_names.csv")
# limited_depth <- list.files("data files", recursive = T, pattern = "MAG_SNV_depth_info_", full.names = T)
# 
# all_depth <- data_frame()
# for(i in 1:length(limited_depth)){
#   sample_depth <- read_csv(limited_depth[i])
#   all_depth <- rbind(all_depth, sample_depth)
# }
# 
# colnames(all_depth)[3] <- "depth"
# all_depth <- left_join(all_depth, sample_names)
# all_depth$time <- as.character(all_depth$Pond %>% substr(4,4))
# all_depth <- all_depth[, c(1:3,5,6)]
#write.csv(all_depth, "data files/all_mag_depth_info.csv", row.names = F)
# 
#all_depth <- read_csv("data files/all_mag_depth_info.csv")
#
# N_SNVs <- subset(all_MAG_SNVs, mutation_type == "N")
# all_depth$time <- as.character(all_depth$time)
# N_depth <- subset(all_depth, group %in% N_SNVs$group)
# N_SNV_depth <- full_join(N_SNVs, N_depth)
# N_SNV_depth$new_ref_freq <- with(N_SNV_depth, ifelse(depth >= 5, 1, NA))
# N_SNV_depth$final_ref_freq <- with(N_SNV_depth, ifelse(is.na(ref_freq), new_ref_freq, ref_freq))
# N_SNV_depth <- N_SNV_depth %>% group_by(group) %>% fill(mag, .direction = "updown")
# N_SNV_depth <- N_SNV_depth %>% group_by(Sample) %>% fill(Name, time, .direction = "updown")
# N_SNV_depth <- N_SNV_depth %>% group_by(Sample, mag) %>% fill(mag_pond, time, .direction = "updown")
# N_SNV_depth <- N_SNV_depth %>% ungroup()
# N_SNV_depth <- N_SNV_depth[, c("mag", "Name", "time", "group", "final_ref_freq")]
# write.csv(N_SNV_depth, "data files/nonsyn_snv_depth.csv", row.names = F)
# 
# N_SNV_depth <- read_csv("data files/nonsyn_snv_depth.csv")
# snv_freq_wide <- pivot_wider(N_SNV_depth, names_from = "time", values_from = "final_ref_freq")
# snv_freq_wide <- snv_freq_wide[complete.cases(snv_freq_wide), ]
# snv_freq_wide$freq_dif <- snv_freq_wide$`1` - snv_freq_wide$`3`
# write.csv(snv_freq_wide, "data files/snv_freq_wide.csv", row.names = F)

snv_freq_wide <- read_csv("data files/snv_freq_wide.csv")
rel_abun <- read_csv("data files/MAG_rel_abun_change.csv")
snv_freq_wide <- left_join(snv_freq_wide, rel_abun[, c(1:3,6)])
snv_freq_wide <- subset(snv_freq_wide, !(`1` == "1" & `3` == "1"))
snv_freq_wide$freq_dif_abs <- abs(snv_freq_wide$freq_dif)
snv_freq_wide$mag_pond <- paste(snv_freq_wide$mag, snv_freq_wide$Name)
snv_freq_wide$plot_abun_change <- round(snv_freq_wide$abun_change, 1)

abun_freq_GBH <- subset(snv_freq_wide, Treatment == "GBH")
abun_freq_CTRL <- subset(snv_freq_wide, Treatment == "Control")

plot_labels <- data.frame(mag = c(unique(abun_freq_GBH$mag)),
                          species = c("28-YEA-48 sp027484655", "Erythrobacter sp027486035", "Sphingorhabdus_B sp027488785",
                          "Sphingorhabdus_B sp021298455", "Erythrobacter sp.", "Aquariibacter lacus",
                          "JAIXWP01 sp027491075", "Prosthecobacter sp.", "CAMBQF01 sp.", "Brevundimonas sp027487935"))

abun_freq_GBH <- left_join(abun_freq_GBH, plot_labels)      
abun_freq_GBH$species_pond <- paste(abun_freq_GBH$species, abun_freq_GBH$Name)
abun_freq_GBH$order_species <- reorder(abun_freq_GBH$species_pond, abun_freq_GBH$abun_change)


GBH_freq_abun <- ggplot(abun_freq_GBH, aes(x = plot_abun_change, y = freq_dif_abs, fill = order_species))+
  geom_boxplot(outlier.size = 0.5)+
  geom_smooth(aes(group = Treatment), method = "glm", colour = "red", size = 0.5)+
  labs(x = "MAG abundance change", y = "nonsyn-SNV frequency change", fill = "", title = "MAGs from GBH ponds")+
  scale_fill_brewer(palette = "Paired")+
  theme_classic()+
  theme(legend.position = "bottom",
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 10, colour = "black"),
        legend.text = element_text(size = 7, colour = "black"),
        title = element_text(size = 10, colour = "black"),
        legend.key.size = unit(0.4, "cm"))+
  guides(fill = guide_legend(ncol = 3), )

ggsave("figures/GBH_snv_freq_abun.pdf", GBH_freq_abun, units = "cm", width = 17, height = 15)

CTRL_freq_abun <- ggplot(abun_freq_CTRL, aes(x = plot_abun_change, y = freq_dif_abs))+
  geom_boxplot(aes(group = mag_pond),outlier.size = 0.5)+
  geom_smooth(aes(group = Treatment), method = "glm", colour = "red", size = 0.5)+
  labs(x = "MAG abundance change", y = "nonsyn-SNV frequency change", title = "MAGs from control ponds")+
  theme_classic()+
  theme(legend.position = "None",
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 10, colour = "black"),
        title = element_text(size = 10, colour = "black"))

ggsave("figures/CTRL_snv_freq_abun.pdf", CTRL_freq_abun, units = "cm", width = 17, height = 14)

ctrl_gbh <- ggarrange(GBH_freq_abun, CTRL_freq_abun, common.legend = F, legend = "bottom", ncol = 1, heights = c(1.5,1),labels = c("A", "B"))
ggsave("figures/gbh_ctrl_freq_abun.pdf", ctrl_gbh, units = "cm", width = 17, height = 19)


gbh_med_freq_abun <- abun_freq_GBH %>% group_by(mag_pond) %>% summarise(freq_median = median(freq_dif_abs),
                                                                        abun_change = median(abun_change))

ctrl_med_freq_abun <- abun_freq_CTRL %>% group_by(mag_pond) %>% summarise(freq_median = median(freq_dif_abs),
                                                                          abun_change = median(abun_change))

gbh_freq_abun <- glm(freq_median ~ abun_change, data = gbh_med_freq_abun)
summary(gbh_freq_abun)

ctrl_freq_abun <- glm(freq_median ~ abun_change, data = ctrl_med_freq_abun)
summary(ctrl_freq_abun)



# ggplot(abun_snv_all, aes(x = abun_change, y = freq_dif, colour = mag_pond))+
#   geom_point()+
#   geom_smooth(aes(group = Treatment), method = "glm", colour = "red", size = 0.5)+
#   scale_colour_manual(values = c("darkmagenta", "darkgreen", "grey80"))+
#   labs(x = "MAG abundance change", y = "SNV abundance change")+
#   theme_classic()+
#   theme(axis.text = element_text(size = 8, colour = "black"), legend.position = "none")+
#   facet_wrap(~Treatment, scales = "free")
# 
# 
# snv_abun <- glm(N_snv_change ~ abun_change, data = abun_snv_all)
# summary(snv_abun)
# 
# abun_snv_GBH <- subset(abun_snv_all, Treatment == "GBH")
# abun_snv_CTRL <- subset(abun_snv_all, Treatment == "Control")
# 
# gbh_snv_abun <- glm(N_snv_change ~ abun_change, data = abun_snv_GBH)
# summary(gbh_snv_abun)
# 
# ctrl_snv_abun <- glm(N_snv_change ~ abun_change, data = abun_snv_CTRL)
# summary(ctrl_snv_abun)
# 
# ggplot(abun_snv_all, aes(x = abun_change, y = N_snv_change, colour = EPSPS_allele))+
#   geom_point()+
#   geom_smooth(aes(group = Treatment), method = "glm", colour = "red", size = 0.5)+
#   scale_colour_manual(values = c("darkmagenta", "darkgreen", "grey80"))+
#   labs(x = "MAG abundance change", y = "SNV abundance change")+
#   theme_classic()+
#   theme(axis.text = element_text(size = 8, colour = "black"), legend.position = "none")+
#   facet_wrap(~Treatment, scales = "free")




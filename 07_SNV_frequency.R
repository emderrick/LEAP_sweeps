library(tidyverse)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

mag_list <- c("MAG_00097_1", "MAG_00103_1","MAG_00110_1", "MAG_00179_1",
              "MAG_00194_1", "MAG_00197_1", "MAG_00201_1","MAG_00674_1")

all_mag_SNVs <- read_csv("data files/MAG_SNV_depth_info.csv")
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
write.csv(snv_frequency, "data files/snv_frequency_changes.csv", row.names = F)

snv_frequency_plot <- pivot_longer(snv_frequency, cols = c(5:8), values_to = "avg_freq", names_to = "Treatment_Time")
snv_frequency_plot$Treatment <- snv_frequency_plot$Treatment_Time %>% substr(1,4)
snv_frequency_plot$Treatment <- snv_frequency_plot$Treatment %>% str_remove("_")
snv_frequency_plot$Treatment_group <- paste(snv_frequency_plot$Treatment, snv_frequency_plot$group)
snv_frequency_plot$Time <- snv_frequency_plot$Treatment_Time %>% str_sub(-2) %>% str_remove("_")
snv_frequency_plot$Time <- paste("Day ", snv_frequency_plot$Time, sep = "")

snv_frequency_plot_a <- subset(snv_frequency_plot, case == "case_a")
plot_a <- ggplot(snv_frequency_plot_a, aes(x = Time, y = avg_freq, colour = Treatment))+
  geom_smooth(aes(group = Treatment), method = "lm")+
  #geom_point()+
  #geom_line(aes(group = Treatment_group))+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  scale_x_discrete(expand = c(0, 0))+
  ylim(0,1)+
  labs(y = "Reference Frequency")+
  theme_bw()
ggsave("figures/snv_freq_case_a_plot.pdf", plot_a, units = "cm", width = 10, height = 6)

snv_frequency_plot_b <- subset(snv_frequency_plot, case == "case_b")
plot_b <- ggplot(snv_frequency_plot_b, aes(x = Time, y = avg_freq, colour = Treatment))+
  geom_smooth(aes(group = Treatment), method = "lm")+
  #geom_point()+
  #geom_line(aes(group = Treatment_group))+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  scale_x_discrete(expand = c(0, 0))+
  ylim(0,1)+
  labs(y = "Reference Frequency")+
  theme_bw()
ggsave("figures/snv_freq_case_b_plot.pdf", plot_b, units = "cm", width = 10, height = 6)


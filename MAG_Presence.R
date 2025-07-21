library(tidyverse)
library(ggplot2)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

mag_pres <- read_csv("data files/T1_SNV_summary_MAG.csv")
mag_pres$coverage <- as.integer(mag_pres$mag_coverage)
mag_pres_5x <- subset(mag_pres, mag_coverage >= 5 & mag_breadth >= 0.5)

mag_pres_plot <- ggplot(mag_pres_5x, aes(x = Time, y = fct_reorder(Name, desc(Treatment)), colour = Treatment))+
  geom_text(aes(label = coverage), size = 6)+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text = element_text(size = 12),
        legend.position = "none")+
  facet_wrap(~mag)
ggsave("figures/T1_MAG_pres.pdf", mag_pres_plot, limitsize = F, width=48, height=48)

mag_pres_5x <- mag_pres_5x[order(mag_pres_5x$Name_Time),]
mag_pres_5x <- mag_pres_5x[order(mag_pres_5x$Time),]
present_MAGs <- mag_pres_5x[, c(1,2,13)]
present_MAGs_wide <- pivot_wider(present_MAGs, names_from = Name_Time, values_from = coverage)
present_MAGs_wide$T1_CTRL_present <- 5 - rowSums(is.na(present_MAGs_wide[2:6]))
present_MAGs_wide$T1_GBH_present <- 4 - rowSums(is.na(present_MAGs_wide[7:10]))
present_MAGs_wide$T3_CTRL_present <- 5 - rowSums(is.na(present_MAGs_wide[11:15]))
present_MAGs_wide$T3_GBH_present <- 4 - rowSums(is.na(present_MAGs_wide[16:19]))

good_examples <- subset(present_MAGs_wide, (T3_CTRL_present >= 1 & T3_GBH_present >= 1 & T1_CTRL_present >= 1 & T1_GBH_present >= 1))
good_MAGs <- pivot_longer(good_examples, cols = c(2:19), names_to = "Sample", values_to = "Coverage")
good_MAGs$Name <- good_MAGs$Sample %>% str_sub(end = -2)
good_MAGs$Time <- ifelse(grepl("1", good_MAGs$Sample), "Day 0", "Day 28")
good_MAGs$Treatment <- ifelse(grepl("CTRL", good_MAGs$Name), "Control", "GBH")
good_MAGs$plot_order <- order(good_MAGs$Name, good_MAGs$Name)
good_MAGs$plot_order <- order(good_MAGs$plot_order, good_MAGs$Treatment)
  
good_mag_pres <- ggplot(good_MAGs, aes(x = Time, y = fct_reorder(Name, plot_order, .desc = T), colour = Treatment))+
  geom_text(aes(label = Coverage), size = 5)+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text = element_text(size = 12),
        legend.position = "none")+
  facet_wrap(~mag)
ggsave("figures/T1_one_one_MAG_pres.pdf", good_mag_pres, width = 8, height = 8)





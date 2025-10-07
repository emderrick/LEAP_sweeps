library(tidyverse)
library(ggplot2)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

mag_pres <- read_csv("data files/T1_mag_info.csv")
mag_pres$Treatment <- ifelse(grepl("CTRL", mag_pres$Name), "Control", "GBH")
mag_pres$Name_Time <- paste(mag_pres$Name, mag_pres$time, sep = " ")
mag_pres$coverage <- as.integer(mag_pres$mag_coverage)

mag_pres <- subset(mag_pres, mag_breadth >= 0.7)

mag_pres_plot <- ggplot(mag_pres, aes(x = as.character(time), y = fct_reorder(Name, desc(Treatment)), colour = Treatment))+
  geom_text(aes(label = coverage), size = 6)+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text = element_text(size = 12),
        legend.position = "none")+
  facet_wrap(~mag)
ggsave("figures/T1_MAG_pres.pdf", mag_pres_plot, limitsize = F, width=48, height=48)

mag_pres <- mag_pres[order(mag_pres$Name),]
mag_pres <- mag_pres[order(mag_pres$time),]
present_MAGs <- mag_pres[, c("mag", "Name_Time", "coverage")]
present_MAGs <- subset(present_MAGs, coverage >= 5)
present_MAGs_wide <- pivot_wider(present_MAGs, names_from = Name_Time, values_from = coverage)
present_MAGs_wide$T1_CTRL_present <- 4 - rowSums(is.na(present_MAGs_wide[2:5]))
present_MAGs_wide$T1_GBH_present <- 4 - rowSums(is.na(present_MAGs_wide[6:9]))
present_MAGs_wide$T3_CTRL_present <- 4 - rowSums(is.na(present_MAGs_wide[10:13]))
present_MAGs_wide$T3_GBH_present <- 4 - rowSums(is.na(present_MAGs_wide[14:17]))

good_examples <- subset(present_MAGs_wide, (T3_CTRL_present >= 1 & T3_GBH_present >= 1 & T1_CTRL_present >= 1 & T1_GBH_present >= 1))

mag_list <- good_examples$mag

mag_plot <- subset(mag_pres, mag %in% mag_list)
mag_plot <- complete(mag_plot, mag, Name, time)
mag_plot$coverage <- with(mag_plot, ifelse(is.na(coverage), "<1", coverage))
mag_plot$Time <- ifelse(mag_plot$time == "1", "Day 0", "Day 28")
mag_plot$Treatment <- ifelse(grepl("CTRL", mag_plot$Name), "Control", "GBH")

mag_plot$plot_order <- order(mag_plot$Name, mag_plot$Name)
mag_plot$plot_order <- order(mag_plot$plot_order, mag_plot$Treatment)

good_mag_pres <- ggplot(mag_plot, aes(x = Time, y = fct_reorder(Name, plot_order, .desc = T)))+
  geom_point(aes(colour = Treatment), shape = 21, size = 8)+
  geom_text(aes(label = coverage, colour = Treatment), size = 3)+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  theme_classic()+
  theme(text = element_text(size = 10), axis.text = element_text(colour = "black"),
        axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none")+
  facet_wrap(~mag)
ggsave("figures/good_T1_MAG_pres.pdf", good_mag_pres, units = "cm", width = 17, height = 20)





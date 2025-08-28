library(tidyverse)
library(ggplot2)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

mag_pres <- read_csv("data files/T1_mag_info.csv")
mag_pres$Treatment <- ifelse(grepl("CTRL", mag_pres$Name), "Control", "GBH")
mag_pres$Name_Time <- paste(mag_pres$Name, mag_pres$time, sep = " ")
mag_pres$coverage <- round(mag_pres$mag_coverage)

mag_pres <- subset(mag_pres, mag_coverage >= 5 & mag_breadth >= 0.5)


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
present_MAGs_wide <- pivot_wider(present_MAGs, names_from = Name_Time, values_from = coverage)
present_MAGs_wide$T1_CTRL_present <- 4 - rowSums(is.na(present_MAGs_wide[2:5]))
present_MAGs_wide$T1_GBH_present <- 4 - rowSums(is.na(present_MAGs_wide[6:9]))
present_MAGs_wide$T3_CTRL_present <- 4 - rowSums(is.na(present_MAGs_wide[10:13]))
present_MAGs_wide$T3_GBH_present <- 4 - rowSums(is.na(present_MAGs_wide[14:17]))

good_examples <- subset(present_MAGs_wide, (T3_CTRL_present >= 1 & T3_GBH_present >= 1 & T1_CTRL_present >= 1 & T1_GBH_present >= 1))
good_MAGs <- pivot_longer(good_examples, cols = c(2:17), names_to = "Sample", values_to = "coverage")
good_MAGs$Name <- good_MAGs$Sample %>% str_sub(end = -2)
good_MAGs$Time <- ifelse(grepl("1", good_MAGs$Sample), "Day 0", "Day 28")
good_MAGs$Treatment <- ifelse(grepl("CTRL", good_MAGs$Name), "Control", "GBH")
good_MAGs$plot_order <- order(good_MAGs$Name, good_MAGs$Name)
good_MAGs$plot_order <- order(good_MAGs$plot_order, good_MAGs$Treatment)
  
good_mag_pres <- ggplot(good_MAGs, aes(x = Time, y = fct_reorder(Name, plot_order, .desc = T), colour = Treatment))+
  geom_text(aes(label = coverage), size = 5)+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text = element_text(size = 12),
        legend.position = "none")+
  facet_wrap(~mag)
ggsave("figures/good_T1_MAG_pres.pdf", good_mag_pres, width = 8, height = 8)





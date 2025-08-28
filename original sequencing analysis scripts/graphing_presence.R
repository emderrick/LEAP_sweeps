library(tidyverse)
library(ggplot2)
library(dplyr)


#load mag info
summary_mags <- read_csv("ANI_95_all_mags.csv")
completed_mags <- complete(summary_mags, mag, timepoint)

mag_pres <- ggplot(completed_mags, aes(x = new_time, y=pond, colour=treatment))+
  geom_point(size=4)+
  scale_colour_manual(values = c("control"="#F56B5CFF", "phosphorus"="#CC3F71FF", "glyphosate"="#3F0F72FF"))+
  theme_classic(base_size=4)+
  scale_x_discrete(limits=c("1", "2", "3"))+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))+
  facet_wrap(~mag, nrow=18, ncol=18)
ggsave("all_mag_pres.png", limitsize = F, width=48, height=48)


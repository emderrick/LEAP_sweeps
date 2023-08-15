library(tidyverse)
library(ggplot2)
library(dplyr)
library(viridis)
library(stringr)
library(cowplot)

mag_list <- list("I4_MAG_00006", "I4_MAG_00065", "L2_MAG_00052", "L3_MAG_00058", "L4_MAG_00099",
                 "L7_MAG_00020", "L7_MAG_00028", "L7_MAG_00043", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")
              
#load mag snv info
all_snv <- read_csv("filtered_ANI_95_mag_SNVs.csv")
scaffold_info <- all_snv[, c('length', 'scaffold', 'coverage', 'new_name', 'mag')] %>% distinct()
mag_cov<- all_snv %>% group_by(scaffold, new_name) %>% summarize(SNV_SNS_tot = sum(number_divergent))
mag_cov <- left_join(mag_cov, scaffold_info, by = c("scaffold", "new_name"))

all_MAG_cov <-  ggplot(mag_cov, aes(x = coverage, y=(SNV_SNS_tot/length)*10^6, colour=new_name)) + 
  geom_point()+
  scale_colour_viridis(discrete = T, option= "magma")+
  labs(y="SNVs / MBp", x="Coverage (x)") +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.title=element_blank())+
  facet_wrap(~mag, nrow = 4, ncol = 4, scales="free")

save_plot("MAG_cov_SNV.jpeg", all_MAG_cov, ncol = 4, nrow = 4, dpi = 300)


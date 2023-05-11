library(tidyverse)
library(dplyr)

setwd("/Users/Emma/OneDrive - McGill University/Documents/phd docs/Chapter 1/Aim 1A")

#load mag snv info
mags <- read_csv("ANI_95_mag_SNVs.csv")

#filter SNVs within 100bp of beginning and end of scaffold
mags<-filter(mags, coverage >=4)
mags$pos_from_end<-mags$length-mags$position
mags<- filter(mags, position > 100)
mags<- filter(mags, pos_from_end > 100)

#assign value of 1 to everything that isn't an SNS and assign 0 to SNS
mags <- mutate(mags, number_SNVs = ifelse(class == "SNS", 0, 1))
#assign value of 1 to SNS and 0 to everythign else
mags <- mutate(mags, number_SNSs = ifelse(class == "SNS", 1, 0))
mags$number_divergent<- 1
mags$group<- paste(mags$mag, "in pond", mags$pond, "at time", mags$new_time)
write.csv(filtered_mags, "filtered_ANI_95_mag_SNVs.csv", row.names=F)

#filter out timepoints I won't use for heatmaps
mags <- filter(mags, mag != "I8_MAG_00005")
mags <- mags %>% filter(!(mag=="L3_MAG_00058" & new_time=="1"))
mags <- mags %>% filter(!(mag=="L3_MAG_00058" & new_time=="3"))
mags <- mags %>% filter(!(mag=="L4_MAG_00099" & new_time=="3"))
mags <- mags %>% filter(!(mag=="L8_MAG_00011" & new_time=="3"))
mags <- mags %>% filter(!(mag=="L7_MAG_00043" & new_time=="3"))
mags <- mags %>% filter(!(mag=="L7_MAG_00028" & new_time=="3"))
mags <- mags %>% filter(!(mag=="I4_MAG_00065" & new_time=="1"))
mags <- mags %>% filter(!(mag=="I4_MAG_00065" & new_time=="3"))
mags <- mags %>% filter(!(mag=="L8_MAG_00042" & new_time=="3"))
mags <- mags %>% filter(!(mag=="L2_MAG_00048" & new_time=="2" & pond=="L7"))
mags <- mags %>% filter(!(mag=="L2_MAG_00052" & new_time=="1"))
filtered_mags <- mags

write.csv(filtered_mags, "filtered_TP_SNVs.csv", row.names=F)

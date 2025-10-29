library(tidyverse)
library(ggplot2)
options(scipen=999)
setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

sample_names <- read_csv("data files/chapter_1_sample_names.csv")
EPSPS <- read_csv("data files/EPSPS_MAG_Classification.csv")
EPSPS <- EPSPS[, c("mag", "Classification")]
args <- read_csv("data files/arg_hits.csv")
rel_abun <- read_tsv("data files/T1_refined_coverM.tsv")
rel_abun$Sample <- rel_abun$Sample %>% substr(1,12)
rel_abun <- left_join(rel_abun, sample_names)
colnames(rel_abun)[2] <- "mag"
colnames(rel_abun)[3] <- "Relative_Abundance"
rel_abun <- subset(rel_abun, !(mag == "unmapped"))
rel_abun$Time <- rel_abun$Pond %>% substr(4,4)
rel_abun$Treatment <- ifelse(grepl("CTRL", rel_abun$Name), "Control", "GBH")
rel_abun <- rel_abun[, c(2,10:12,3,5)]
rel_abun <- subset(rel_abun, !(Name == "CTRL E"))
write.csv(rel_abun, "data files/MAG_rel_abun.csv", row.names = F)

rel_abun <- rel_abun[, c(1:5)]
rel_abun_change <- pivot_wider(rel_abun, names_from = Time, values_from = Relative_Abundance)
rel_abun_change$abun_change <- rel_abun_change$`3` - rel_abun_change$`1`
rel_abun_change <- left_join(rel_abun_change, EPSPS)
rel_abun_change <- left_join(rel_abun_change, args)
rel_abun_change$Classification <- with(rel_abun_change, ifelse(is.na(Classification), "Missing", Classification))
rel_abun_change[is.na(rel_abun_change)] <- 0
rel_abun_change$EPSPS_allele <- with(rel_abun_change, ifelse(Classification == "Class I", "Sensitive", "Resistant"))
rel_abun_change$EPSPS_allele <- with(rel_abun_change, ifelse(Classification == "Missing" | Classification == "Unclassified", "Unclassified", EPSPS_allele))
write.csv(rel_abun_change, "data files/MAG_rel_abun_change.csv", row.names = F)

rel_abun_AVG <- rel_abun_change %>% group_by(mag, Treatment, Classification, EPSPS_allele, `antibiotic efflux`, `antibiotic target alteration`, `antibiotic inactivation`, total_hits) %>% summarise(Day_28_avg = mean(`3`), Day_0_avg = mean(`1`))
write.csv(rel_abun_AVG, "data files/MAG_avg_abun.csv", row.names = F)



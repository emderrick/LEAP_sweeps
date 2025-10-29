library(tidyverse)
setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

EPSPS <- read_csv("data files/EPSPS classification.csv")
EPSPS$mag <- EPSPS$Name %>% substr(1,11)
EPSPS$mag  <- EPSPS$mag  %>% str_remove("_[A-Z]")
EPSPS$motifs  <- EPSPS$ClassIII  %>% substr(7,7)
EPSPS$max_classifation <- as.numeric(apply(EPSPS[, c(2:5)], 1, max, na.rm=TRUE))
EPSPS$Classification <- with(EPSPS, ifelse(max_classifation >= 0.9, "Classified", "Unclassified"))
EPSPS$Classification <- with(EPSPS, ifelse(Classification == "Classified"
                                           & (max_classifation == ClassIalpha
                                           | max_classifation == ClassIbeta),
                                           "Class I", Classification))
EPSPS$Classification <- with(EPSPS, ifelse(Classification == "Classified" 
                                           & max_classifation == ClassII,
                                           "Class II", Classification))
EPSPS$Classification <- with(EPSPS, ifelse(motifs >= 1, "Class III", Classification))
EPSPS <- EPSPS[, c(7, 2:6, 8:10)]
write.csv(EPSPS, "data files/EPSPS_MAG_Classification.csv", row.names = F)

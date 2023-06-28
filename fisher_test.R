library(tidyverse)
library(dplyr)

epsps_sns <- data.frame("increase" = c(3,2),
                        "decrease" = c(2,1),
                        row.names = c("Class I", "Class II"))

fish_test_sns<- fisher.test(epsps_sns, alternative = "less")

epsps_snv <- data.frame("increase" = c(1,1),
                        "decrease" = c(4,2),
                        row.names = c("Class I", "Class II"))

fish_test_snv<- fisher.test(epsps_snv)

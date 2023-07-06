library(tidyverse)
library(dplyr)

epsps_sns <- data.frame("increase" = c(5,2),
                        "decrease" = c(3,1),
                        row.names = c("Class I", "Class II"))

fish_test_sns<- fisher.test(epsps_sns, alternative = "greater")
                            

epsps_snv <- data.frame("increase" = c(4,1),
                        "decrease" = c(4,2),
                        row.names = c("Class I", "Class II"))

fish_test_snv<- fisher.test(epsps_snv, alternative = "less")

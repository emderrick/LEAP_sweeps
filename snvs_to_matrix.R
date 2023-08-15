library(tidyverse)
library(dplyr)

all_MAG_snvs <- read_csv("all_MAG_SNVs_med_Aug15.csv")

mag_list <- list("I4_MAG_00006", "I4_MAG_00065", "L2_MAG_00052", "L3_MAG_00058", "L4_MAG_00099",
                 "L7_MAG_00020", "L7_MAG_00028", "L7_MAG_00043", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")

for(MAG in mag_list){
  MAG_snvs <- subset(all_MAG_snvs, mag==MAG & number_divergent==1 & is.na(gene) == F)
  MAG_snvs$name <- sub(" ", "_", MAG_snvs$name)
  gene_pop_matrix <- dcast(MAG_snvs, name~gene, value.var = "number_divergent", sum)
  write.csv(gene_pop_matrix, paste(MAG, "_gene_pop_matrix.csv", sep=""), row.names = F)
}

library(tidyverse)
library(dplyr)

L3_MAG_00058 <- read_csv("all_L3_MAG_00058_SNVs.csv")
#L3 and L4 are controls, L7 and L8 are gly
L3_MAG_00058_L3 <- filter(L3_MAG_00058, (pond=="L3" & final_ref_freq >= 0.5))
L3_MAG_00058_L4 <- filter(L3_MAG_00058, (pond=="L4" & final_ref_freq >= 0.5))
L3_MAG_00058_L7 <- filter(L3_MAG_00058, (pond=="L7" & final_ref_freq <= 0.5))
L3_MAG_00058_L8 <- filter(L3_MAG_00058, (pond=="L8" & final_ref_freq <= 0.5))
L3_MAG_00058_50 <- (list(L3_MAG_00058_L3,L3_MAG_00058_L4,L3_MAG_00058_L7,L3_MAG_00058_L8) %>% reduce(inner_join, by='groups'))
write.csv(L3_MAG_00058_50, "L3_MAG_00058_50.csv", row.names=F)

L4_MAG_00099 <- read_csv("all_L4_MAG_00099_SNVs.csv")
#L4 is a control, L2, L6, L7, and L8 are gly
L4_MAG_00099_L4 <- filter(L4_MAG_00099, (pond=="L4" & final_ref_freq >= 0.5))
L4_MAG_00099_L2 <- filter(L4_MAG_00099, (pond=="L2" & final_ref_freq <= 0.5))
L4_MAG_00099_L6 <- filter(L4_MAG_00099, (pond=="L6" & final_ref_freq <= 0.5))
L4_MAG_00099_L7 <- filter(L4_MAG_00099, (pond=="L7" & final_ref_freq <= 0.5))
L4_MAG_00099_L8 <- filter(L4_MAG_00099, (pond=="L8" & final_ref_freq <= 0.5))
L4_MAG_00099_50 <- (list(L4_MAG_00099_L4,L4_MAG_00099_L2,L4_MAG_00099_L6,L4_MAG_00099_L7,L4_MAG_00099_L8) %>% reduce(inner_join, by='groups'))
write.csv(L4_MAG_00099_50, "L4_MAG_00099_50.csv", row.names=F)

L8_MAG_00019 <- read_csv("all_L8_MAG_00019_SNVs.csv")
#I8 is a control, L2, L6 and L8 are gly
L8_MAG_00019_I8 <- filter(L8_MAG_00019, (pond=="I8" & final_ref_freq <= 0.5))
L8_MAG_00019_L2 <- filter(L8_MAG_00019, (pond=="L2" & final_ref_freq >= 0.5))
L8_MAG_00019_L6 <- filter(L8_MAG_00019, (pond=="L6" & final_ref_freq >= 0.5))
L8_MAG_00019_L8 <- filter(L8_MAG_00019, (pond=="L8" & final_ref_freq >= 0.5))
L8_MAG_00019_50 <- (list(L8_MAG_00019_I8,L8_MAG_00019_L2,L8_MAG_00019_L6,L8_MAG_00019_L8) %>% reduce(inner_join, by='groups'))
write.csv(L8_MAG_00019_50, "L8_MAG_00019_50.csv", row.names=F)

L8_MAG_00011 <- read_csv("all_L8_MAG_00011_SNVs.csv")
#I8 is a control, L2, and L8 are gly
L8_MAG_00011_I8 <- filter(L8_MAG_00011, (pond=="I8" & final_ref_freq <= 0.5))
L8_MAG_00011_L2 <- filter(L8_MAG_00011, (pond=="L2" & final_ref_freq >= 0.5))
L8_MAG_00011_L8 <- filter(L8_MAG_00011, (pond=="L8" & final_ref_freq >= 0.5))
L8_MAG_00011_50 <- (list(L8_MAG_00011_I8,L8_MAG_00011_L2,L8_MAG_00011_L8) %>% reduce(inner_join, by='groups'))
write.csv(L8_MAG_00011_50, "L8_MAG_00011_50.csv", row.names=F)

L7_MAG_00043 <- read_csv("all_L7_MAG_00043_SNVs.csv")
#L4 is a control, L6 and L7 are gly
L7_MAG_00043_L4 <- filter(L7_MAG_00043, (pond=="L4" & final_ref_freq <= 0.5))
L7_MAG_00043_L6 <- filter(L7_MAG_00043, (pond=="L6" & final_ref_freq >= 0.5))
L7_MAG_00043_L7 <- filter(L7_MAG_00043, (pond=="L7" & final_ref_freq >= 0.5))
L7_MAG_00043_50 <- (list(L7_MAG_00043_L4,L7_MAG_00043_L6,L7_MAG_00043_L7) %>% reduce(inner_join, by='groups'))
write.csv(L7_MAG_00043_50, "L7_MAG_00043_50.csv", row.names=F)

L7_MAG_00028 <- read_csv("all_L7_MAG_00028_SNVs.csv")
#I8 is a control, L2, L6, and L7 are gly
L7_MAG_00028_I8 <- filter(L7_MAG_00028, (pond=="I8" & final_ref_freq <= 0.5))
L7_MAG_00028_L2 <- filter(L7_MAG_00028, (pond=="L2" & final_ref_freq >= 0.5))
L7_MAG_00028_L6 <- filter(L7_MAG_00028, (pond=="L6" & final_ref_freq >= 0.5))
L7_MAG_00028_L7 <- filter(L7_MAG_00028, (pond=="L7" & final_ref_freq >= 0.5))
L7_MAG_00028_50 <- (list(L7_MAG_00028_I8,L7_MAG_00028_L2,L7_MAG_00028_L6,L7_MAG_00028_L7) %>% reduce(inner_join, by='groups'))
write.csv(L7_MAG_00028_50, "L7_MAG_00028_50.csv", row.names=F)

I4_MAG_00006 <- read_csv("all_I4_MAG_00006_SNVs.csv")
#I4 and I8 are controls, L2 and L8 are gly
I4_MAG_00006_I4 <- filter(I4_MAG_00006, pond=="I4" & final_ref_freq >= 0.5)
I4_MAG_00006_I8 <- filter(I4_MAG_00006, pond=="I8" & final_ref_freq >= 0.5)
I4_MAG_00006_L2 <- filter(I4_MAG_00006, pond=="L2" & final_ref_freq <= 0.5) 
I4_MAG_00006_L8 <- filter(I4_MAG_00006, pond=="L8" & final_ref_freq <= 0.5)
I4_MAG_00006_50 <- (list(I4_MAG_00006_I4, I4_MAG_00006_I8, I4_MAG_00006_L2, I4_MAG_00006_L8 ) %>% reduce(inner_join, by='groups'))
write.csv(I4_MAG_00006_50, "I4_MAG_00006_50.csv", row.names=F)

I4_MAG_00065 <- read_csv("all_I4_MAG_00065_SNVs.csv")
#I4, I8, K1, L3, L4 are controls, L6 is gly
I4_MAG_00065_I4 <- filter(I4_MAG_00065, pond=="I4" & final_ref_freq >= 0.5)
I4_MAG_00065_I8 <- filter(I4_MAG_00065, pond=="I8" & final_ref_freq >= 0.5)
I4_MAG_00065_K1 <- filter(I4_MAG_00065, pond=="K1" & final_ref_freq >= 0.5)
I4_MAG_00065_L3 <- filter(I4_MAG_00065, pond=="L3" & final_ref_freq >= 0.5)
I4_MAG_00065_L4 <- filter(I4_MAG_00065, pond=="L4" & final_ref_freq >= 0.5) 
I4_MAG_00065_L6 <- filter(I4_MAG_00065, pond=="L6" & final_ref_freq <= 0.5)
I4_MAG_00065_50 <- (list(I4_MAG_00065_I4, I4_MAG_00065_I8, I4_MAG_00065_K1, I4_MAG_00065_L3, I4_MAG_00065_L4, I4_MAG_00065_L6) %>% reduce(inner_join, by='groups'))
write.csv(I4_MAG_00065_50, "I4_MAG_00065_50.csv", row.names=F)

L7_MAG_00020 <- read_csv("all_L7_MAG_00020_SNVs.csv")
#K1, L3, and L4 are a control, L2 is gly
L7_MAG_00020_K1 <- filter(L7_MAG_00020, (pond=="K1" & final_ref_freq <= 0.5))
L7_MAG_00020_L3 <- filter(L7_MAG_00020, (pond=="L3" & final_ref_freq <= 0.5))
L7_MAG_00020_L4 <- filter(L7_MAG_00020, (pond=="L4" & final_ref_freq <= 0.5))
L7_MAG_00020_L2 <- filter(L7_MAG_00020, (pond=="L2" & final_ref_freq >= 0.5))
L7_MAG_00020_50 <- (list(L7_MAG_00020_K1,L7_MAG_00020_L3,L7_MAG_00020_L4,L7_MAG_00020_L2) %>% reduce(inner_join, by='groups'))
write.csv(L7_MAG_00020_50, "L7_MAG_00020_50.csv", row.names=F)

L8_MAG_00042 <- read_csv("all_L8_MAG_00042_SNVs.csv")
#I8, K1, L3, and L4 are a control, L8 is gly
L8_MAG_00042_I8 <- filter(L8_MAG_00042, (pond=="I8" & final_ref_freq <= 0.5))
L8_MAG_00042_K1 <- filter(L8_MAG_00042, (pond=="K1" & final_ref_freq <= 0.5))
L8_MAG_00042_L3 <- filter(L8_MAG_00042, (pond=="L3" & final_ref_freq <= 0.5))
L8_MAG_00042_L4 <- filter(L8_MAG_00042, (pond=="L4" & final_ref_freq <= 0.5))
L8_MAG_00042_L8 <- filter(L8_MAG_00042, (pond=="L8" & final_ref_freq >= 0.5))
L8_MAG_00042_50 <- (list(L8_MAG_00042_I8,L8_MAG_00042_K1,L8_MAG_00042_L3,L8_MAG_00042_L4,L8_MAG_00042_L8) %>% reduce(inner_join, by='groups'))
write.csv(L8_MAG_00042_50, "L8_MAG_00042_50.csv", row.names=F)

L2_MAG_00052 <- read_csv("all_L2_MAG_00052_SNVs.csv")
#I4, I8, K1, and L4 are a control, L2, is gly
L2_MAG_00052_I4 <- filter(L2_MAG_00052, (pond=="I4" & final_ref_freq <= 0.5))
L2_MAG_00052_I8 <- filter(L2_MAG_00052, (pond=="I8" & final_ref_freq <= 0.5))
L2_MAG_00052_K1 <- filter(L2_MAG_00052, (pond=="K1" & final_ref_freq <= 0.5))
L2_MAG_00052_L4 <- filter(L2_MAG_00052, (pond=="L4" & final_ref_freq <= 0.5))
L2_MAG_00052_L2 <- filter(L2_MAG_00052, (pond=="L2" & final_ref_freq >= 0.5))
L2_MAG_00052_50 <- (list(L2_MAG_00052_I4,L2_MAG_00052_I8,L2_MAG_00052_K1,L2_MAG_00052_L4,L2_MAG_00052_L2) %>% reduce(inner_join, by='groups'))
write.csv(L2_MAG_00052_50, "L2_MAG_00052_50.csv", row.names=F)

I4_MAG_00006_genes<- subset(I4_MAG_00006_50, select=c(mag.x, groups, gene.x))
L7_MAG_00028_genes<- subset(L7_MAG_00028_50, select=c(mag.x, groups, gene.x))
L8_MAG_00011_genes<- subset(L8_MAG_00011_50, select=c(mag.x, groups, gene.x))
L8_MAG_00019_genes<- subset(L8_MAG_00019_50, select=c(mag.x, groups, gene.x))
L8_MAG_00042_genes<- subset(L8_MAG_00042_50, select=c(mag.x, groups, gene.x))
L3_MAG_00058_genes<- subset(L3_MAG_00058_50, select=c(mag.x, groups, gene.x))
I4_MAG_00065_genes<- subset(I4_MAG_00065_50, select=c(mag.x, groups, gene.x))
L4_MAG_00099_genes<- subset(L4_MAG_00099_50, select=c(mag.x, groups, gene.x))
L2_MAG_00052_genes<- subset(L2_MAG_00052_50, select=c(mag.x, groups, gene.x))
L7_MAG_00043_genes<- subset(L7_MAG_00043_50, select=c(mag.x, groups, gene.x))
L7_MAG_00020_genes<- subset(L7_MAG_00020_50, select=c(mag.x, groups, gene.x))

all_genes <- rbind(I4_MAG_00006_genes, L7_MAG_00028_genes, L8_MAG_00011_genes, L8_MAG_00019_genes,
                   L8_MAG_00042_genes, L3_MAG_00058_genes, I4_MAG_00065_genes, L4_MAG_00099_genes,
                   L2_MAG_00052_genes, L7_MAG_00043_genes, L7_MAG_00020_genes)
all_genes <- rename(all_genes, gene = gene.x, mag=mag.x)
only_genes <- subset(all_genes, gene!="NA", select=c(gene))
write.csv(all_genes, file="all_genes.csv", row.names=F)
write.csv(only_genes, file="only_genes.csv", row.names=F)

library(tidyverse)
library(dplyr)

filtered_SNVs <- read_csv("filtered_ANI_95_mag_SNVs.csv")

#have to do each MAG individually or I run out of memory

I4_MAG_00006 <- filter(filtered_SNVs, mag=="I4_MAG_00006")
I4_MAG_00006$groups<- paste(I4_MAG_00006$scaffold, str_pad(I4_MAG_00006$position, 7, pad = "0"))
I4_MAG_00006 <- complete(I4_MAG_00006, timepoint, groups)

I4_MAG_00006_depth <- read.table("~/Documents/I4_MAG_00006_depth.txt", sep="\t", header=F)
I4_MAG_00006_depth <- I4_MAG_00006_depth %>% rename("scaffold"="V1", "position"="V2", "samtools_depth"="V3", "timepoint"="V4")
I4_MAG_00006_depth$groups<- paste(I4_MAG_00006_depth$scaffold, str_pad(I4_MAG_00006_depth$position, 7, pad = "0"))
I4_MAG_00006_depth <- mutate(I4_MAG_00006_depth, new_ref_freq = ifelse(samtools_depth >= 5, 1, NA))

all_I4_MAG_00006 <- right_join(I4_MAG_00006_depth, I4_MAG_00006, by=c("timepoint", "groups"))
all_I4_MAG_00006 <- mutate(all_I4_MAG_00006, final_ref_freq= ifelse(is.na(ref_freq), new_ref_freq, ref_freq))
all_I4_MAG_00006 = all_I4_MAG_00006[,c(1:14, 88, 15:87)]
write.csv(all_I4_MAG_00006, "all_I4_MAG_00006_SNVs.csv", row.names=F)

I4_MAG_00065 <- filter(filtered_SNVs, mag=="I4_MAG_00065")
I4_MAG_00065$groups<- paste(I4_MAG_00065$scaffold, str_pad(I4_MAG_00065$position, 7, pad = "0"))
I4_MAG_00065 <- complete(I4_MAG_00065, timepoint, groups)

I4_MAG_00065_depth <- read.table("~/Documents/I4_MAG_00065_depth.txt", sep="\t", header=F)
I4_MAG_00065_depth <- I4_MAG_00065_depth %>% rename("scaffold"="V1", "position"="V2", "samtools_depth"="V3", "timepoint"="V4")
I4_MAG_00065_depth$groups<- paste(I4_MAG_00065_depth$scaffold, str_pad(I4_MAG_00065_depth$position, 7, pad = "0"))
I4_MAG_00065_depth <- mutate(I4_MAG_00065_depth, new_ref_freq = ifelse(samtools_depth >= 5, 1, NA))

all_I4_MAG_00065 <- right_join(I4_MAG_00065_depth, I4_MAG_00065, by=c("timepoint", "groups"))
all_I4_MAG_00065 <- mutate(all_I4_MAG_00065, final_ref_freq= ifelse(is.na(ref_freq), new_ref_freq, ref_freq))
all_I4_MAG_00065 = all_I4_MAG_00065[,c(1:14, 88, 15:87)]
write.csv(all_I4_MAG_00065, "all_I4_MAG_00065_SNVs.csv", row.names=F)

L2_MAG_00052 <- filter(filtered_SNVs, mag=="L2_MAG_00052")
L2_MAG_00052$groups<- paste(L2_MAG_00052$scaffold, str_pad(L2_MAG_00052$position, 7, pad = "0"))
L2_MAG_00052 <- complete(L2_MAG_00052, timepoint, groups)

L2_MAG_00052_depth <- read.table("~/Documents/L2_MAG_00052_depth.txt", sep="\t", header=F)
L2_MAG_00052_depth <- L2_MAG_00052_depth %>% rename("scaffold"="V1", "position"="V2", "samtools_depth"="V3", "timepoint"="V4")
L2_MAG_00052_depth$groups<- paste(L2_MAG_00052_depth$scaffold, str_pad(L2_MAG_00052_depth$position, 7, pad = "0"))
L2_MAG_00052_depth <- mutate(L2_MAG_00052_depth, new_ref_freq = ifelse(samtools_depth >= 5, 1, NA))

all_L2_MAG_00052 <- right_join(L2_MAG_00052_depth, L2_MAG_00052, by=c("timepoint", "groups"))
all_L2_MAG_00052 <- mutate(all_L2_MAG_00052, final_ref_freq= ifelse(is.na(ref_freq), new_ref_freq, ref_freq))
all_L2_MAG_00052 = all_L2_MAG_00052[,c(1:14, 88, 15:87)]
write.csv(all_L2_MAG_00052, "all_L2_MAG_00052_SNVs.csv", row.names=F)


L3_MAG_00058 <- filter(filtered_SNVs, mag=="L3_MAG_00058")
L3_MAG_00058$groups<- paste(L3_MAG_00058$scaffold, str_pad(L3_MAG_00058$position, 7, pad = "0"))
L3_MAG_00058 <- complete(L3_MAG_00058, timepoint, groups)

L3_MAG_00058_depth <- read.table("~/Documents/L3_MAG_00058_depth.txt", sep="\t", header=F)
L3_MAG_00058_depth <- L3_MAG_00058_depth %>% rename("scaffold"="V1", "position"="V2", "samtools_depth"="V3", "timepoint"="V4")
L3_MAG_00058_depth$groups<- paste(L3_MAG_00058_depth$scaffold, str_pad(L3_MAG_00058_depth$position, 7, pad = "0"))
L3_MAG_00058_depth <- mutate(L3_MAG_00058_depth, new_ref_freq = ifelse(samtools_depth >= 5, 1, NA))

all_L3_MAG_00058 <- right_join(L3_MAG_00058_depth, L3_MAG_00058, by=c("timepoint", "groups"))
all_L3_MAG_00058 <- mutate(all_L3_MAG_00058, final_ref_freq= ifelse(is.na(ref_freq), new_ref_freq, ref_freq))
all_L3_MAG_00058 = all_L3_MAG_00058[,c(1:14, 88, 15:87)]
write.csv(all_L3_MAG_00058, "all_L3_MAG_00058_SNVs.csv", row.names=F)


L4_MAG_00099 <- filter(filtered_SNVs, mag=="L4_MAG_00099")
L4_MAG_00099$groups<- paste(L4_MAG_00099$scaffold, str_pad(L4_MAG_00099$position, 7, pad = "0"))
L4_MAG_00099 <- complete(L4_MAG_00099, timepoint, groups)

L4_MAG_00099_depth <- read.table("~/Documents/L4_MAG_00099_depth.txt", sep="\t", header=F)
L4_MAG_00099_depth <- L4_MAG_00099_depth %>% rename("scaffold"="V1", "position"="V2", "samtools_depth"="V3", "timepoint"="V4")
L4_MAG_00099_depth$groups<- paste(L4_MAG_00099_depth$scaffold, str_pad(L4_MAG_00099_depth$position, 7, pad = "0"))
L4_MAG_00099_depth <- mutate(L4_MAG_00099_depth, new_ref_freq = ifelse(samtools_depth >= 5, 1, NA))

all_L4_MAG_00099 <- right_join(L4_MAG_00099_depth, L4_MAG_00099, by=c("timepoint", "groups"))
all_L4_MAG_00099 <- mutate(all_L4_MAG_00099, final_ref_freq= ifelse(is.na(ref_freq), new_ref_freq, ref_freq))
all_L4_MAG_00099 = all_L4_MAG_00099[,c(1:14, 88, 15:87)]
write.csv(all_L4_MAG_00099, "all_L4_MAG_00099_SNVs.csv", row.names=F)


L7_MAG_00020 <- filter(filtered_SNVs, mag=="L7_MAG_00020")
L7_MAG_00020$groups<- paste(L7_MAG_00020$scaffold, str_pad(L7_MAG_00020$position, 7, pad = "0"))
L7_MAG_00020 <- complete(L7_MAG_00020, timepoint, groups)

L7_MAG_00020_depth <- read.table("~/Documents/L7_MAG_00020_depth.txt", sep="\t", header=F)
L7_MAG_00020_depth <- L7_MAG_00020_depth %>% rename("scaffold"="V1", "position"="V2", "samtools_depth"="V3", "timepoint"="V4")
L7_MAG_00020_depth$groups<- paste(L7_MAG_00020_depth$scaffold, str_pad(L7_MAG_00020_depth$position, 7, pad = "0"))
L7_MAG_00020_depth <- mutate(L7_MAG_00020_depth, new_ref_freq = ifelse(samtools_depth >= 5, 1, NA))

all_L7_MAG_00020 <- right_join(L7_MAG_00020_depth, L7_MAG_00020, by=c("timepoint", "groups"))
all_L7_MAG_00020 <- mutate(all_L7_MAG_00020, final_ref_freq= ifelse(is.na(ref_freq), new_ref_freq, ref_freq))
all_L7_MAG_00020 = all_L7_MAG_00020[,c(1:14, 88, 15:87)]
write.csv(all_L7_MAG_00020, "all_L7_MAG_00020_SNVs.csv", row.names=F)


L7_MAG_00028 <- filter(filtered_SNVs, mag=="L7_MAG_00028")
L7_MAG_00028$groups<- paste(L7_MAG_00028$scaffold, str_pad(L7_MAG_00028$position, 7, pad = "0"))
L7_MAG_00028 <- complete(L7_MAG_00028, timepoint, groups)

L7_MAG_00028_depth <- read.table("~/Documents/L7_MAG_00028_depth.txt", sep="\t", header=F)
L7_MAG_00028_depth <- L7_MAG_00028_depth %>% rename("scaffold"="V1", "position"="V2", "samtools_depth"="V3", "timepoint"="V4")
L7_MAG_00028_depth$groups<- paste(L7_MAG_00028_depth$scaffold, str_pad(L7_MAG_00028_depth$position, 7, pad = "0"))
L7_MAG_00028_depth <- mutate(L7_MAG_00028_depth, new_ref_freq = ifelse(samtools_depth >= 5, 1, NA))

all_L7_MAG_00028 <- right_join(L7_MAG_00028_depth, L7_MAG_00028, by=c("timepoint", "groups"))
all_L7_MAG_00028 <- mutate(all_L7_MAG_00028, final_ref_freq= ifelse(is.na(ref_freq), new_ref_freq, ref_freq))
all_L7_MAG_00028 = all_L7_MAG_00028[,c(1:14, 88, 15:87)]
write.csv(all_L7_MAG_00028, "all_L7_MAG_00028_SNVs.csv", row.names=F)


L7_MAG_00043 <- filter(filtered_SNVs, mag=="L7_MAG_00043")
L7_MAG_00043$groups<- paste(L7_MAG_00043$scaffold, str_pad(L7_MAG_00043$position, 7, pad = "0"))
L7_MAG_00043 <- complete(L7_MAG_00043, timepoint, groups)

L7_MAG_00043_depth <- read.table("~/Documents/L7_MAG_00043_depth.txt", sep="\t", header=F)
L7_MAG_00043_depth <- L7_MAG_00043_depth %>% rename("scaffold"="V1", "position"="V2", "samtools_depth"="V3", "timepoint"="V4")
L7_MAG_00043_depth$groups<- paste(L7_MAG_00043_depth$scaffold, str_pad(L7_MAG_00043_depth$position, 7, pad = "0"))
L7_MAG_00043_depth <- mutate(L7_MAG_00043_depth, new_ref_freq = ifelse(samtools_depth >= 5, 1, NA))

all_L7_MAG_00043 <- right_join(L7_MAG_00043_depth, L7_MAG_00043, by=c("timepoint", "groups"))
all_L7_MAG_00043 <- mutate(all_L7_MAG_00043, final_ref_freq= ifelse(is.na(ref_freq), new_ref_freq, ref_freq))
all_L7_MAG_00043 = all_L7_MAG_00043[,c(1:14, 88, 15:87)]
write.csv(all_L7_MAG_00043, "all_L7_MAG_00043_SNVs.csv", row.names=F)


L8_MAG_00011 <- filter(filtered_SNVs, mag=="L8_MAG_00011")
L8_MAG_00011$groups<- paste(L8_MAG_00011$scaffold, str_pad(L8_MAG_00011$position, 7, pad = "0"))
L8_MAG_00011 <- complete(L8_MAG_00011, timepoint, groups)

L8_MAG_00011_depth <- read.table("~/Documents/L8_MAG_00011_depth.txt", sep="\t", header=F)
L8_MAG_00011_depth <- L8_MAG_00011_depth %>% rename("scaffold"="V1", "position"="V2", "samtools_depth"="V3", "timepoint"="V4")
L8_MAG_00011_depth$groups<- paste(L8_MAG_00011_depth$scaffold, str_pad(L8_MAG_00011_depth$position, 7, pad = "0"))
L8_MAG_00011_depth <- mutate(L8_MAG_00011_depth, new_ref_freq = ifelse(samtools_depth >= 5, 1, NA))

all_L8_MAG_00011 <- right_join(L8_MAG_00011_depth, L8_MAG_00011, by=c("timepoint", "groups"))
all_L8_MAG_00011 <- mutate(all_L8_MAG_00011, final_ref_freq= ifelse(is.na(ref_freq), new_ref_freq, ref_freq))
all_L8_MAG_00011 = all_L8_MAG_00011[,c(1:14, 88, 15:87)]
write.csv(all_L8_MAG_00011, "all_L8_MAG_00011_SNVs.csv", row.names=F)


L8_MAG_00019 <- filter(filtered_SNVs, mag=="L8_MAG_00019")
L8_MAG_00019$groups<- paste(L8_MAG_00019$scaffold, str_pad(L8_MAG_00019$position, 7, pad = "0"))
L8_MAG_00019 <- complete(L8_MAG_00019, timepoint, groups)

L8_MAG_00019_depth <- read.table("~/Documents/L8_MAG_00019_depth.txt", sep="\t", header=F)
L8_MAG_00019_depth <- L8_MAG_00019_depth %>% rename("scaffold"="V1", "position"="V2", "samtools_depth"="V3", "timepoint"="V4")
L8_MAG_00019_depth$groups<- paste(L8_MAG_00019_depth$scaffold, str_pad(L8_MAG_00019_depth$position, 7, pad = "0"))
L8_MAG_00019_depth <- mutate(L8_MAG_00019_depth, new_ref_freq = ifelse(samtools_depth >= 5, 1, NA))

all_L8_MAG_00019 <- right_join(L8_MAG_00019_depth, L8_MAG_00019, by=c("timepoint", "groups"))
all_L8_MAG_00019 <- mutate(all_L8_MAG_00019, final_ref_freq= ifelse(is.na(ref_freq), new_ref_freq, ref_freq))
all_L8_MAG_00019 = all_L8_MAG_00019[,c(1:14, 88, 15:87)]
write.csv(all_L8_MAG_00019, "all_L8_MAG_00019_SNVs.csv", row.names=F)


L8_MAG_00042 <- filter(filtered_SNVs, mag=="L8_MAG_00042")
L8_MAG_00042$groups<- paste(L8_MAG_00042$scaffold, str_pad(L8_MAG_00042$position, 7, pad = "0"))
L8_MAG_00042 <- complete(L8_MAG_00042, timepoint, groups)

L8_MAG_00042_depth <- read.table("~/Documents/L8_MAG_00042_depth.txt", sep="\t", header=F)
L8_MAG_00042_depth <- L8_MAG_00042_depth %>% rename("scaffold"="V1", "position"="V2", "samtools_depth"="V3", "timepoint"="V4")
L8_MAG_00042_depth$groups<- paste(L8_MAG_00042_depth$scaffold, str_pad(L8_MAG_00042_depth$position, 7, pad = "0"))
L8_MAG_00042_depth <- mutate(L8_MAG_00042_depth, new_ref_freq = ifelse(samtools_depth >= 5, 1, NA))

all_L8_MAG_00042 <- right_join(L8_MAG_00042_depth, L8_MAG_00042, by=c("timepoint", "groups"))
all_L8_MAG_00042 <- mutate(all_L8_MAG_00042, final_ref_freq= ifelse(is.na(ref_freq), new_ref_freq, ref_freq))
all_L8_MAG_00042 = all_L8_MAG_00042[,c(1:14, 88, 15:87)]
write.csv(all_L8_MAG_00042, "all_L8_MAG_00042_SNVs.csv", row.names=F)












library(tidyverse)
library(ggplot2)
library(cowplot)
library(reporter)
library(tidyverse)
library(svglite)


setwd("/Users/Emma/Documents/manuscript version/")

genome_files <- list.files("95_profiles//", recursive = T, pattern = ".*genome_info.tsv", full.names = T)
all_mags <- data.frame()
for(i in 1:length(genome_files)){
  pond_time_mags <- read.table(genome_files[i], sep="\t", header = T)
  timepoint <- gsub(".*instrain_output/", "", genome_files[i]) %>% substr(15,23)
  pond_time_mags <- cbind(pond_time_mags, timepoint = rep(timepoint, nrow(pond_time_mags)))
  all_mags <- rbind(all_mags, pond_time_mags)
}

all_mags$time <- all_mags$timepoint %>% substr(9,9)
all_mags$new_time <- as.numeric(all_mags$time) + 1
all_mags$pond <- all_mags$timepoint %>% substr(1,2)
all_mags <- mutate(all_mags, treatment = ifelse(pond == "I4" | pond == "K1" | pond == "L3" | pond == "L4" | pond == "I8", "control", "glyphosate"))
all_mags <- mutate(all_mags, treatment = ifelse(pond == "I8", "phosphorus", treatment))
all_mags <- all_mags %>% rename("mag_length" = "length")
all_mags <- all_mags %>% rename("mag" = "genome")
all_mags <- all_mags %>% rename("mag_coverage" = "coverage")
all_mags <- all_mags %>% mutate(name = case_when(timepoint%>%substr(1,2) == "K1" ~ "Control A", timepoint%>%substr(1,2) == "I4" ~ "Control B",
                                                 timepoint%>%substr(1,2) == "L3" ~ "Control C", timepoint%>%substr(1,2) == "L4" ~ "Control D",
                                                 timepoint%>%substr(1,2) == "I8" ~ "Control E", timepoint%>%substr(1,2) == "L2" ~ "GBH A",
                                                 timepoint%>%substr(1,2) == "L6" ~ "GBH B", timepoint%>%substr(1,2) == "L7" ~ "GBH C",
                                                 timepoint%>%substr(1,2) == "L8" ~ "GBH D"))

all_mags$new_name <- paste(paste(all_mags$name, "at T"), all_mags$new_time, sep = "")
#write.csv(all_mags, "all_mags_not_subsamp.csv", row.names = F)

scaffold_files <- list.files("95_profiles//", recursive = T, pattern = ".*scaffold_info.tsv", full.names = T)
all_scaffolds <- data.frame()

for(i in 1:length(scaffold_files)){
  pond_time_scaffolds <- read.table(scaffold_files[i], sep = "\t", header = T)
  timepoint <- gsub(".*instrain_output/", "", scaffold_files[i]) %>% substr(15,23)
  pond_time_scaffolds <- cbind(pond_time_scaffolds, timepoint = rep(timepoint, nrow(pond_time_scaffolds)))
  all_scaffolds <- rbind(all_scaffolds, pond_time_scaffolds)
} 

all_scaffolds$mag <- all_scaffolds$scaffold %>% substr(1,12)

mag_scaffolds <- all_scaffolds %>% left_join(all_mags, by = c("mag", "timepoint"))

SNV_files <- list.files("95_profiles//", recursive = T, pattern = ".*SNVs.tsv", full.names = T)
all_SNVs <- data.frame()
for(i in 1:length(SNV_files)){
  pond_time_SNV <- read.table(SNV_files[i], sep = "\t", header = T)
  timepoint <- gsub(".*instrain_output/", "", SNV_files[i]) %>% substr(15,23)
  pond_time_SNV <- cbind(pond_time_SNV, timepoint = rep(timepoint, nrow(pond_time_SNV)))
  all_SNVs <- rbind(all_SNVs, pond_time_SNV)
}  

mag_scaf_SNV <- left_join(all_SNVs, mag_scaffolds, by = c("timepoint","scaffold"))

mag_list <- c('I4_MAG_00006','L7_MAG_00028','L8_MAG_00011', 'L8_MAG_00019', 'L8_MAG_00042', 'I4_MAG_00065',
              'L3_MAG_00058', 'L7_MAG_00020',  'L7_MAG_00043', 'L2_MAG_00052', 'L4_MAG_00099')

mag_scaf_SNV <- subset(mag_scaf_SNV, mag %in% mag_list)

#write.csv(mag_scaf_SNV, "all_SNVs_not_subsamp.csv", row.names = F)
mag_scaf_SNV <- mag_scaf_SNV[, c("scaffold", "timepoint", "length", "coverage", "ref_freq", "position", "position_coverage", "mag_coverage", "mag", "new_time", "pond", "treatment", "name", "new_name")]

mag_scaf_SNV$pos_from_end <- mag_scaf_SNV$length - mag_scaf_SNV$position
mag_scaf_SNV <- subset(mag_scaf_SNV, position > 100)
mag_scaf_SNV <- subset(mag_scaf_SNV, pos_from_end > 100)
mag_scaf_SNV$full_group <- paste(mag_scaf_SNV$mag, mag_scaf_SNV$new_name)
mag_scaf_SNV <- mag_scaf_SNV %>% subset(!(position_coverage > mag_coverage*3)) 
mag_scaf_SNV <- mag_scaf_SNV %>% subset(!(position_coverage < mag_coverage/3))
#write.csv(mag_scaf_SNV, "filtered_mag_SNVs_not_subsamp.csv", row.names = F)


mag_labs <- c(I4_MAG_00006 = "Burkholderiaceae 1", I4_MAG_00065 = "Roseomonas_A", L2_MAG_00052 = "Erythrobacter", 
              L3_MAG_00058 = "Prosthecobacter", L4_MAG_00099 = "Bosea sp001713455", L7_MAG_00020 = "Sphingorhabdus_B",
              L7_MAG_00028 = "Burkholderiaceae 2", L7_MAG_00043 = "Luteolibacter", L8_MAG_00011 = "Verrucomicrobiae", 
              L8_MAG_00019 = "Flavobacteriales 1", L8_MAG_00042 = "Flavobacteriales 2")

all_snv <- mag_scaf_SNV
all_snv <- subset(all_snv, new_time == 2)
all_snv <- subset(all_snv, ref_freq < 1 & ref_freq > 0)
all_snv$one <- 1

snv_sum <- all_snv %>% group_by(mag, scaffold, new_name, name, length, coverage) %>% summarise(count = sum(one))


snv_sum$mag_order = factor(snv_sum$mag, levels=c('I4_MAG_00006','L7_MAG_00028','L8_MAG_00011', 'L8_MAG_00019', 'L8_MAG_00042','I4_MAG_00065',
                                                 'L3_MAG_00058', 'L7_MAG_00020',  'L7_MAG_00043', 'L2_MAG_00052', 'L4_MAG_00099'))

all_MAG_scaf_cov <- ggplot(snv_sum, aes(x = coverage, y = log10((count/length)*10^6), colour = name)) + 
  geom_point(size = 1)+
  scale_colour_manual(values = c("#3b0043","#5d1c66","#7f3888","#c270ce", "#D9A5E0","#FFEA94","#ffd633","#ffc71f","#fdb721"))+
  #scale_colour_manual(values = c("#002C3D", "#005F73", "#0A9396", "#94D2BD", "#EE9B00", "#CA6702", "#BB3E03", "#AE2012", "#9B2226"))+
  labs(y = paste("log", {subsc('10')}, " SNVs / Mbp"), x="Coverage (X)", colour= "Pond") +
  theme_classic()+
  theme(text = element_text(size = 12, colour = 'black'),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        legend.position = "bottom",
        strip.text.x.top = element_text(face = "bold"))+
  scale_y_continuous(limits=c(0,5))+
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 10)))+
  facet_wrap(~mag_order, ncol = 6, scales = "free", labeller = labeller(mag_order = mag_labs))

save_plot("MAG_scaf_cov_SNV_sub.svg", all_MAG_scaf_cov, base_height = 3.5, base_width = 8.5)

library(tidyverse)
library(ggplot2)
library(viridis)
library(RColorBrewer)

setwd("/Users/Emma/Documents/manuscript version/")

genome_files <- list.files("95_profiles//", recursive = T, pattern = ".*genome_info.tsv", full.names = T)
all_mags <- data.frame()
for(i in 1:length(genome_files)){
  pond_time_mags <- read.table(genome_files[i], sep="\t", header = T)
  timepoint <- gsub(".*instrain_output/", "", genome_files[i]) %>% substr(15,23)
  pond_time_mags <- cbind(pond_time_mags, timepoint = rep(timepoint, nrow(pond_time_mags)))
  all_mags <- rbind(all_mags, pond_time_mags)
}

all_mags <- all_mags[, c("genome", "length")] 
colnames(all_mags)[1] <- "mag"
all_mags <- all_mags %>% unique()

coverm_files <- list.files("coverm_subsamp//", recursive = T, pattern = ".*_coverm", full.names = T)
all_coverm <- data.frame()
for(i in 1:length(coverm_files)){
  time_coverm <- read.table(coverm_files[i], sep="\t", header = T)
  colnames(time_coverm) <- c("mag", "cover_m_rel_abun", "mean", "trimmed_mean", "covered_bases", "read_count", "reads_per_base")
  timepoint <- coverm_files[i] %>% substr(26,34)
  time_coverm <- cbind(time_coverm, timepoint = rep(timepoint, nrow(time_coverm)))
  all_coverm <- rbind(all_coverm, time_coverm)
}

all_coverm <- left_join(all_coverm, all_mags)
all_coverm <- subset(all_coverm, mag != "L7_MAG_00037")
all_coverm <- subset(all_coverm, mag != "unmapped")
all_coverm$pond <- all_coverm$timepoint %>% substr(1,2)
all_coverm <- all_coverm %>% mutate(name = case_when(pond == "K1" ~ "Control A", pond == "I4" ~ "Control B",
                                                     pond == "L3" ~ "Control C", pond == "L4" ~ "Control D",
                                                     pond == "I8" ~ "Control E", pond == "L2" ~ "GBH A",
                                                     pond == "L6" ~ "GBH B", pond == "L7" ~ "GBH C",
                                                     pond == "L8" ~ "GBH D"))

all_coverm$time <- as.numeric(all_coverm$timepoint %>% substr(9,9)) + 1
all_coverm$new_name <- paste(paste(all_coverm$name, "at T"), all_coverm$time, sep = "")
all_coverm <- mutate(all_coverm, treatment = ifelse(pond == "I4" | pond == "K1" | pond == "L3" | pond == "L4" | pond == "I8", "Control", "GBH"))
all_coverm <- mutate(all_coverm, treatment = ifelse(pond == "I8", "Phosphorus Control", treatment))
all_coverm <- mutate(all_coverm, treatment_graph = ifelse(pond == "I4" | pond == "K1" | pond == "L3" | pond == "L4" | pond == "I8", "Control", "GBH"))
all_coverm <- mutate(all_coverm, length = ifelse(mag == "L2_MAG_00050", 983613, length))
all_coverm <- mutate(all_coverm, length = ifelse(mag == "L3_MAG_00067", 1212569, length))
all_coverm$breadth <- all_coverm$covered_bases/all_coverm$length

coverm_50 <- subset(all_coverm, mean > 1 & breadth >= 0.5)
coverm_50$one <- 1
richness <- coverm_50 %>% group_by(timepoint, name, new_name, treatment, treatment_graph, time) %>% summarize(species = sum(one))

pond_richness <- ggplot(richness, aes(y = species, x = treatment_graph, colour = treatment_graph, fill = treatment_graph))+
  geom_boxplot()+
  scale_colour_manual(values = c("#6b033e", "#ffcb30"))+
  scale_fill_manual(values = c("#BE9EAF", "#FDF2D3"))+
  theme_classic()+
  theme(legend.position = "none",
        axis.text = element_text(colour = "black", size = 12),
        axis.title.x = element_text(vjust = -1, size = 12),
        axis.title.y = element_text(vjust = 2, size = 12))+
  labs(y = "Species Richness (# MAGs)", x = "Treatment", colour = "Pond")+
  scale_y_continuous(limits = c(0, 160), expand = expansion(mult = c(0, 0.05)))
  
save_plot("pond_richness.jpeg", pond_richness, base_width = 3, base_height = 4)

GBH_species <- subset(richness, treatment_graph == "GBH")$species
control_species <- subset(richness, treatment_graph == "Control")$species
richness_stats <- wilcox.test(GBH_species, control_species)

snv_info <- read_csv("all_SNV_sum_subsamp.csv")
snv_info <- subset(snv_info, select = c("mag", "new_name", "SNV_Mbp"))
snv_richness <- inner_join(snv_info, richness, by = c("new_name"))
snv_richness <- snv_richness %>% mutate(mag_name = case_when(mag == "I4_MAG_00006" ~ "Burkholderiaceae 1", mag == "L7_MAG_00028" ~ "Burkholderiaceae 2", mag == "L8_MAG_00011" ~ "Verrucomicrobiae",
                                                            mag == "L8_MAG_00019" ~ "Flavobacteriales 1", mag == "L8_MAG_00042" ~ "Flavobacteriales 2", mag == "I4_MAG_00065" ~ "Roseomonas_A",
                                                            mag == "L3_MAG_00058" ~ "Prosthecobacter", mag == "L7_MAG_00020" ~ "Sphingorhabdus_B", mag == "L7_MAG_00043" ~ "Luteolibacter",
                                                            mag == "L2_MAG_00052" ~ "Erythrobacter", mag == "L4_MAG_00099" ~ "Bosea sp001713455"))

snv_richness$mag_order = factor(snv_richness$mag_name, levels=c('Burkholderiaceae 1', 'Burkholderiaceae 2', 'Verrucomicrobiae', 'Flavobacteriales 1', 'Flavobacteriales 2',
                                                'Roseomonas_A', 'Prosthecobacter', 'Sphingorhabdus_B',  'Luteolibacter',
                                                'Erythrobacter', 'Bosea sp001713455'))

snv_richness_plot <- ggplot(snv_richness, aes(x = species, y = SNV_Mbp, colour = mag_order))+
  geom_point(size = 3, aes(shape = factor(treatment)))+
  geom_smooth(method="lm",se=FALSE,fullrange=T,color="black", size = 0.3, aes(group = mag))+
  labs(x = "Species Richness (# MAGs)", y = "SNVs per Mbp", colour = "MAG", shape = "Treatment")+
  theme_classic()+
  scale_colour_brewer(palette = "Set3")+
  theme(axis.text = element_text(colour = "black", size = 12),
        axis.title.x = element_text(vjust = -1, size = 12),
        axis.title.y = element_text(vjust = 2, size = 12),
        legend.key.size = unit(0.25, "cm"),
        legend.key.spacing.y = unit(0.1, "cm"))+
  scale_y_continuous(limits = c(-500, 25000))

save_plot("SNV_richness.jpeg", snv_richness_plot, base_width = 8, base_height = 4)

summary(lm(SNV_Mbp ~ species, data = snv_richness))

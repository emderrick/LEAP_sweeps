library(tidyverse)
library(dplyr)
library(ggplot2)
library(cowplot)

sens_mags <- list("I4_MAG_00006", "L7_MAG_00028", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")
res_mags <- list("I4_MAG_00065", "L3_MAG_00058", "L7_MAG_00020", "L7_MAG_00043")
mag_labs <- c(I4_MAG_00006 = "Burkholderiaceae 1", I4_MAG_00065 = "Roseomonas_A", L2_MAG_00052 = "Erythrobacter", 
              L3_MAG_00058 = "Prosthecobacter", L4_MAG_00099 = "Bosea sp001713455", L7_MAG_00020 = "Sphingorhabdus_B",
              L7_MAG_00028 = "Burkholderiaceae 2", L7_MAG_00043 = "Luteolibacter", L8_MAG_00011 = "Verrucomicrobiae", 
              L8_MAG_00019 = "Flavobacteriales 1", L8_MAG_00042 = "Flavobacteriales 2")

mapped_reads <- read.table("mapped_reads_part1.txt", sep=",", header = T)
mag_info <- read_csv("all_mags_subsamp.csv")
mag_info <- subset(mag_info, select = c("mag", "mag_length", "pond", "name")) %>% unique()
mapped_reads <- as.data.frame(apply(mapped_reads, 2, function(x) gsub('\\s+', '', x)))
mapped_reads$pulse <- mapped_reads$pulse %>% substr(1,9)
mapped_reads$pond <- mapped_reads$pulse %>% substr(1,2)
mapped_reads <- right_join(mapped_reads, mag_info)
mapped_reads$EPSPS_class <- with(mapped_reads, ifelse(mag %in% sens_mags, "Sensitive", "Unclassified"))
mapped_reads$EPSPS_class <- with(mapped_reads, ifelse(mag %in% res_mags, "Resistant", EPSPS_class))
mapped_reads <- mutate(mapped_reads, treatment = ifelse(pond == "I4" | pond == "K1" | pond == "L3" | pond == "L4" | pond == "I8", "Control", "GBH"))
mapped_reads <- mutate(mapped_reads, treatment = ifelse(pond == "I8", "phosphorus", treatment))

mapped_reads$paired_mapped <- as.numeric(mapped_reads$mapped)/2
mapped_reads$relative_abundance <- mapped_reads$paired_mapped / as.numeric(mapped_reads$total)
mapped_reads$RA_genome <- (mapped_reads$relative_abundance / mapped_reads$mag_length) * 10^6
mapped_reads$time <- mapped_reads$pulse %>% substr(9,9) %>% as.numeric() + 1

mapped_reads$mag_order = factor(mapped_reads$mag, levels=c('I4_MAG_00006', 'L7_MAG_00028', 'L8_MAG_00011', 'L8_MAG_00019', 'L8_MAG_00042',
                                                            'I4_MAG_00065', 'L3_MAG_00058', 'L7_MAG_00020', 'L7_MAG_00043',
                                                            'L2_MAG_00052', 'L4_MAG_00099'))
mapped_reads_no3 <- subset(mapped_reads, !(time == "3"))

RA_T2 <- subset(mapped_reads, time == "2")
RA_T2 <- subset(RA_T2, select = c("mag", "name", "relative_abundance", "EPSPS_class", "treatment"))
write.csv(RA_T2, "RA_T2.csv", row.names = F)

RA_time <- ggplot(mapped_reads_no3, aes(x = time, y = RA_genome, group = pond, colour = treatment))+
  geom_line()+
  theme_classic()+
  theme(text = element_text(size = 10),
        axis.line.x.bottom = element_line(linewidth = 0.85),
        axis.line.y.left = element_line(linewidth = 0.85),
        axis.title.y = element_text(size = 13), 
        axis.text = element_text(colour = "black"),
        strip.text.x.top = element_text(size = 12, margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"), face = "bold"),
        strip.background = element_rect(linewidth = 2, colour = "black"))+
  labs(y = "Relative Abundance", x = "time point", colour = "Treatment")+
  scale_x_continuous(breaks = c(1, 2))+
  facet_wrap(~mag_order, scales = "free", labeller = labeller(mag_order = mag_labs))

save_plot("relative_abundance.jpeg", RA_time, dpi = 300, base_width = 12, base_height = 6)



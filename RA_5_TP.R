library(tidyverse)
library(dplyr)
library(ggplot2)
library(cowplot)
library(svglite)


setwd("/Users/emma/Documents/manuscript version/")

sens_mags <- list("I4_MAG_00006", "L7_MAG_00028", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")
res_mags <- list("I4_MAG_00065", "L3_MAG_00058", "L7_MAG_00020", "L7_MAG_00043")
mag_labs <- c(I4_MAG_00006 = "Burkholderiaceae 1", I4_MAG_00065 = "Roseomonas_A", L2_MAG_00052 = "Erythrobacter", 
              L3_MAG_00058 = "Prosthecobacter", L4_MAG_00099 = "Bosea sp001713455", L7_MAG_00020 = "Sphingorhabdus_B",
              L7_MAG_00028 = "Burkholderiaceae 2", L7_MAG_00043 = "Luteolibacter", L8_MAG_00011 = "Verrucomicrobiae", 
              L8_MAG_00019 = "Flavobacteriales 1", L8_MAG_00042 = "Flavobacteriales 2")


total_reads <- read_csv("read_stats.csv")
total_reads <- as.data.frame(apply(total_reads, 2, function(x) gsub('\\s+', '', x)))
total_reads$timepoint <- total_reads$file %>% substr(1,4)
total_reads <- total_reads[, c(9, 4)] 
total_reads_sum <- total_reads %>% unique()

mapped_reads <- read.table("reads_all_5.txt", sep = ",", header = F)
mapped_reads <- as.data.frame(apply(mapped_reads, 2, function(x) gsub('\\s+', '', x)))
mapped_reads$timepoint <- mapped_reads$V2 %>% substr(1,4)
mapped_reads <- mapped_reads[, c(1, 3:4)]

read_info <- left_join(mapped_reads, total_reads_sum)
read_info <- read_info %>% rename("mag" = "V1", "mapped" = "V3")
read_info$pond <- read_info$timepoint %>% substr(1,2)
read_info$time <- read_info$timepoint %>% substr(4,4)
read_info <- read_info %>% mutate(pulse = case_when(time == 1 ~ 0, time == 2 ~ 1, time == 3 ~ 1, time == 4 ~ 2, time == 5 ~ 2))
read_info <- read_info[, c(1,2,4,5,7)]
read_info_sum <- read_info %>% group_by(mag, pond, pulse) %>% summarize(mapped_reads = sum(as.numeric(mapped)), all_reads = sum(as.numeric(num_seqs)))


mag_info <- read_csv("all_mags_subsamp.csv")
mag_info <- subset(mag_info, select = c("mag", "name", "pond", "mag_length")) %>% unique()

read_mag_info <- left_join(read_info_sum, mag_info, by = c("mag", "pond"))
read_mag_info <- read_mag_info %>% group_by(mag) %>% fill(mag_length, .direction = "updown")
read_mag_info <- read_mag_info %>% group_by(pond) %>% fill(name, .direction = "updown")

read_mag_info$EPSPS_class <- with(read_mag_info, ifelse(mag %in% sens_mags, "Sensitive", "Unclassified"))
read_mag_info$EPSPS_class <- with(read_mag_info, ifelse(mag %in% res_mags, "Resistant", EPSPS_class))
read_mag_info <- mutate(read_mag_info, treatment = ifelse(pond == "I4" | pond == "K1" | pond == "L3" | pond == "L4" | pond == "I8", "Control", "GBH"))

read_mag_info$relative_abundance <- read_mag_info$mapped_reads / read_mag_info$all_reads
read_mag_info$RA_genome <- (read_mag_info$relative_abundance / read_mag_info$mag_length) * 10^6
read_mag_info$coverage <- (read_mag_info$mapped_reads * 150) / read_mag_info$mag_length


read_mag_info$mag_order = factor(read_mag_info$mag, levels=c('I4_MAG_00006', 'L7_MAG_00028', 'L8_MAG_00011', 'L8_MAG_00019', 'L8_MAG_00042',
                                                     'I4_MAG_00065', 'L3_MAG_00058', 'L7_MAG_00020', 'L7_MAG_00043',
                                                     'L2_MAG_00052', 'L4_MAG_00099'))

read_mag_info <- subset(read_mag_info, (mag == "I4_MAG_00006" & (name == "Control B" | name == "Control E" | name == "GBH A" | name == "GBH D") |
                                  mag == "I4_MAG_00065" & (name == "Control A" | name == "Control B" | name == "Control C" | name == "Control D" | name == "Control E" | name == "GBH B") |
                                  mag == "L2_MAG_00052" & (name == "Control A" | name == "Control B" | name == "Control D" | name == "Control E" | name == "GBH A") |
                                  mag == "L3_MAG_00058" & (name == "Control C" | name == "Control D" | name == "GBH C" | name == "GBH D") |
                                  mag == "L4_MAG_00099" & (name == "Control D" | name == "GBH A" | name == "GBH B" | name == "GBH C" | name == "GBH D") |
                                  mag == "L7_MAG_00020" & (name == "Control A" | name == "Control C" | name == "Control D" |name == "GBH A" |name == "GBH C") |
                                  mag == "L7_MAG_00028" & (name == "Control E" | name == "GBH A" | name == "GBH B" | name == "GBH C") |
                                  mag == "L7_MAG_00043" & (name == "Control D" | name == "GBH B" | name == "GBH C") |
                                  mag == "L8_MAG_00011" & (name == "Control E" | name == "GBH A" | name == "GBH D") |
                                  mag == "L8_MAG_00019" & (name == "Control E" | name == "GBH A" | name == "GBH B" | name == "GBH D") |
                                  mag == "L8_MAG_00042" & (name == "Control A" | name == "Control C" | name == "Control D" | name == "Control E" | name == "GBH D")))

both_pulse_RA <- ggplot(subset(read_mag_info, mag != "L4_MAG_00099" & mag != "L2_MAG_00052"), aes(x = pulse, y = RA_genome, group = pond, colour = treatment))+
  geom_line()+
  theme_classic()+
  scale_colour_manual(values = c("#6b033e", "#ffcb30"))+
  theme(text = element_text(size = 12),
        axis.line.x.bottom = element_line(linewidth = 0.85),
        axis.line.y.left = element_line(linewidth = 0.85),
        axis.title.y = element_text(size = 13), 
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        strip.text.x.top = element_text(size = 12, margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"), face = "bold"),
        strip.background = element_rect(linewidth = 2, colour = "black"))+
  labs(y = "Relative Abundance", x = "Time Point", colour = "Treatment")+
  scale_x_continuous(breaks = c(0,1,2), labels = c('TP 1', 'TP 2', 'TP 3'))+
  facet_wrap(~mag_order, scales = "free_y", ncol = 5, labeller = labeller(mag_order = mag_labs))

save_plot("relative_abundance_allpulse.svg", both_pulse_RA, base_width = 12, base_height = 4)


pulse_1_RA <- ggplot(subset(read_mag_info, pulse != 2 & mag != "L4_MAG_00099" & mag != "L2_MAG_00052"), aes(x = pulse, y = RA_genome, group = name, colour = name))+
  geom_line(size = 1)+
  theme_classic()+
  scale_colour_manual(values = c("#3b0043","#5d1c66","#7f3888","#c270ce", "#D9A5E0","#FFEA94","#ffd633","#ffc71f","#fdb721"))+
  theme(text = element_text(size = 10),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        legend.position = "bottom",
        strip.text.x.top = element_text(face = "bold"))+
  labs(y = "Relative Abundance", x = "Time Point", colour = "Treatment")+
  scale_x_continuous(breaks = c(0,1), labels = c('TP 1', 'TP 2'))+
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 5)))+
  facet_wrap(~mag_order, scales = "free", ncol = 5, labeller = labeller(mag_order = mag_labs))

save_plot("relative_abundance_pulse1.svg", pulse_1_RA, base_height = 3.5, base_width = 8.5)


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


total_reads <- read_csv("read_stats.csv")
total_reads <- as.data.frame(apply(total_reads, 2, function(x) gsub('\\s+', '', x)))
total_reads$pond <- total_reads$file %>% substr(1,2)
total_reads$time <- total_reads$file %>% substr(4,4)
total_reads <- total_reads %>% mutate(pulse = case_when(time == 1 ~ 0, time == 2 ~ 1, time == 3 ~ 1, time == 4 ~ 2, time == 5 ~ 2))
total_reads$timepoint <- paste(total_reads$pond, "_pulse", total_reads$pulse, sep = "")
total_reads <- total_reads[, c(12, 4)]
total_reads_sum <- total_reads %>% group_by(timepoint) %>% summarize(reads = sum(as.numeric(num_seqs)))

mapped_reads <- read.table("mapped_reads.txt", sep = ",", header = T)
mapped_reads <- as.data.frame(apply(mapped_reads, 2, function(x) gsub('\\s+', '', x)))
mapped_reads$pulse <- mapped_reads$pulse %>% substr(1,9)

read_info <- left_join(mapped_reads, total_reads_sum, by = c("pulse" = "timepoint"))
mag_info <- read_csv("all_mags_subsamp.csv")
mag_info <- subset(mag_info, select = c("mag", "mag_length", "name", "timepoint")) %>% unique()

read_info$pond <- read_info$pulse %>% substr(1,2)
read_info <- left_join(read_info, mag_info, by = c("pulse" = "timepoint", "MAG" = "mag"))
read_info$EPSPS_class <- with(read_info, ifelse(MAG %in% sens_mags, "Sensitive", "Unclassified"))
read_info$EPSPS_class <- with(read_info, ifelse(MAG %in% res_mags, "Resistant", EPSPS_class))
read_info <- mutate(read_info, treatment = ifelse(pond == "I4" | pond == "K1" | pond == "L3" | pond == "L4" | pond == "I8", "Control", "GBH"))
read_info <- mutate(read_info, treatment = ifelse(pond == "I8", "Phosphorus", treatment))
read_info <- read_info %>% group_by(MAG) %>% fill(mag_length, .direction = "updown")
read_info <- read_info %>% group_by(pond) %>% fill(name, .direction = "updown")

read_info$mapped_reads <- as.numeric(read_info$mapped_reads)
read_info$paired_total <- read_info$reads/2
read_info$relative_abundance <- read_info$mapped_reads / read_info$paired_total
read_info$RA_genome <- (read_info$relative_abundance / read_info$mag_length) * 10^6
read_info$time <- read_info$pulse %>% substr(9,9) %>% as.numeric() + 1


read_info$mag_order = factor(read_info$MAG, levels=c('I4_MAG_00006', 'L7_MAG_00028', 'L8_MAG_00011', 'L8_MAG_00019', 'L8_MAG_00042',
                                                     'I4_MAG_00065', 'L3_MAG_00058', 'L7_MAG_00020', 'L7_MAG_00043',
                                                     'L2_MAG_00052', 'L4_MAG_00099'))

read_info <- subset(read_info, (MAG == "I4_MAG_00006" & (name == "Control B" | name == "Control E" | name == "GBH A" | name == "GBH D") |
                                        MAG == "I4_MAG_00065" & (name == "Control A" | name == "Control B" | name == "Control C" | name == "Control D" | name == "Control E" | name == "GBH B") |
                                        MAG == "L2_MAG_00052" & (name == "Control A" | name == "Control B" | name == "Control D" | name == "Control E" | name == "GBH A") |
                                        MAG == "L3_MAG_00058" & (name == "Control C" | name == "Control D" | name == "GBH C" | name == "GBH D") |
                                        MAG == "L4_MAG_00099" & (name == "Control D" | name == "GBH A" | name == "GBH B" | name == "GBH C" | name == "GBH D") |
                                        MAG == "L7_MAG_00020" & (name == "Control A" | name == "Control C" | name == "Control D" |name == "GBH A") |
                                        MAG == "L7_MAG_00028" & (name == "Control E" | name == "GBH A" | name == "GBH B" | name == "GBH C") |
                                        MAG == "L7_MAG_00043" & (name == "Control D" | name == "GBH B" | name == "GBH C") |
                                        MAG == "L8_MAG_00011" & (name == "Control E" | name == "GBH A" | name == "GBH D") |
                                        MAG == "L8_MAG_00019" & (name == "Control E" | name == "GBH A" | name == "GBH B" | name == "GBH D") |
                                        MAG == "L8_MAG_00042" & (name == "Control A" | name == "Control C" | name == "Control D" | name == "Control E" | name == "GBH D")))

read_info_no3 <- subset(read_info, !(time == "3"))

RA_time <- ggplot(read_info, aes(x = time, y = RA_genome, group = pond, colour = treatment))+
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

save_plot("relative_abundance_only_select_all_time.jpeg", RA_time, dpi = 300, base_width = 12, base_height = 6)


RA_T2 <- subset(read_info, time == "2")
RA_T2 <- subset(RA_T2, select = c("mag", "name", "relative_abundance", "EPSPS_class", "treatment"))
write.csv(RA_T2, "RA_T2.csv", row.names = F)
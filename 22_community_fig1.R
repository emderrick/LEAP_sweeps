library(tidyverse)
library(ape)
library(vegan)
library(ggpubr)
library(lme4)
library(RColorBrewer)

options(scipen=999)
setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

sample_info <- read_csv("data files/chapter_1_sample_names.csv")
sample_info$Treatment <- ifelse(grepl("CTRL", sample_info$Name), "Control", "GBH")
sample_info$Time <- ifelse(grepl("_1", sample_info$Pond), "Day 0", "Day 28")
sample_info$Treatment_Time <- paste(sample_info$Treatment, sample_info$Time, sep = " ")
sample_info <- subset(sample_info, !(Name == "CTRL E"))

T1_samples <- subset(sample_info, Time == "Day 0")
T1_species_beta <- read_tsv("data files/beta_div_species_day_0.txt", skip = 8)
rownames(T1_species_beta) <- T1_species_beta$x
T1_species_beta <- T1_species_beta[, c(2:9)]
T1_species_beta <- t(T1_species_beta)
T1_species_bray <- as.dist(T1_species_beta)

T1_species_pcoa <- pcoa(T1_species_bray, correction="none", rn=NULL)
T1_relative_eig_species <- T1_species_pcoa$values$Relative_eig
T1_species_pcoa_coords <- as.data.frame(T1_species_pcoa$vectors)
T1_species_pcoa_coords <- cbind(T1_species_pcoa_coords, T1_samples)

T1_species <- ggplot(T1_species_pcoa_coords, aes(x = Axis.1, y = Axis.2, colour = Treatment)) +
  geom_point(size = 2)+
  scale_colour_manual(values = c("#98BF64", "#BE93D4"))+
  labs(title = "Day 0 Bray-Curtis Dissimilarity",
       x = paste(paste("PC 1", signif((T1_relative_eig_species[1] * 100), 4), sep = "  "),"%", sep = ""),
       y = paste(paste("PC 2", signif((T1_relative_eig_species[2] * 100), 4), sep = "  "),"%", sep = ""),
       colour = "Treatment")+
  xlim(-0.35, 0.35)+
  ylim(-0.35, 0.35)+
  theme_classic()+
  theme(title = element_text(size = 10),
        text = element_text(size = 4),
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 10, colour = "black"),
        legend.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 10, colour = "black"))+
  annotate("text", x = -0.26, y = -0.28, label = "italic(p) == 0.4426", size = 3, parse = TRUE)+
  annotate("text", x = -0.265, y = -0.33, label = "R^2 == 0.127", size = 3, parse = TRUE)

T2_samples <- subset(sample_info, Time == "Day 28")
T2_species_beta <- read_tsv("data files/beta_div_species_day_28.txt", skip = 8)
rownames(T2_species_beta) <- T2_species_beta$x
T2_species_beta <- T2_species_beta[, c(2:9)]
T2_species_beta <- t(T2_species_beta)
T2_species_bray <- as.dist(T2_species_beta)

T2_species_pcoa <- pcoa(T2_species_bray, correction="none", rn=NULL)
T2_relative_eig_species <- T2_species_pcoa$values$Relative_eig
T2_species_pcoa_coords <- as.data.frame(T2_species_pcoa$vectors)
T2_species_pcoa_coords <- cbind(T2_species_pcoa_coords, T2_samples)

T2_species <- ggplot(T2_species_pcoa_coords, aes(x = Axis.1, y = Axis.2, colour = Treatment)) +
  geom_point(size = 2)+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  labs(title = "Day 28 Bray-Curtis Dissimilarity",
       x = paste(paste("PC 1", signif((T2_relative_eig_species[1] * 100), 4), sep = "  "),"%", sep = ""),
       y = paste(paste("PC 2", signif((T2_relative_eig_species[2] * 100), 4), sep = "  "),"%", sep = ""),
       colour = "Treatment")+
  xlim(-0.35, 0.35)+
  ylim(-0.35, 0.35)+
  theme_classic()+
  theme(title = element_text(size = 10),
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 10, colour = "black"),
        legend.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 10, colour = "black"))+
  annotate("text", x = -0.27, y = -0.28, label = "italic(p) == 0.0240", size = 3, parse = TRUE)+
  annotate("text", x = -0.255, y = -0.33, label = "R^2 == 0.5337", size = 3, parse = TRUE)

species_pcoa <- ggarrange(T1_species, T2_species, nrow = 1, common.legend = F, legend = "bottom")
ggsave("figures/beta_div_species.pdf", species_pcoa, units = "cm", width = 17, height = 8)

### mags

mags <- read_csv("data files/MAG_rel_abun.csv")
mags$Name_Time <- paste(mags$Name, mags$Time, sep = " ")
mags$Time <- ifelse(grepl("1", mags$Name_Time), "Day 0", "Day 28")
mag_info <- read_csv("data files/T1_mag_info.csv")
mag_info$Time <- ifelse(mag_info$time == "1", "Day 0", "Day 28")
mags <- left_join(mags, mag_info[, c(1,4,7,10)])
mags <- mags %>% group_by(mag) %>% fill(mag_length, .direction = "updown") %>% ungroup()
mags$percent_cov <- mags$`Covered Bases`/mags$mag_length

pond_community <- subset(mags, Relative_Abundance >= 0.01 & percent_cov >= 0.1)
pond_community <- pond_community %>% group_by(Name, Time, Treatment) %>% count() %>% ungroup()
pond_community$Treatment_Time <- paste(pond_community$Treatment, pond_community$Time, sep = " ")

mag_community <- ggplot(pond_community, aes(x = Treatment, group = Treatment_Time, colour = Treatment_Time, fill = Treatment_Time, y = n))+
  geom_boxplot()+
  geom_point()+
  scale_fill_manual(values = c("#98BF64", "darkgreen", "#BE93D4", "darkmagenta"))+
  scale_colour_manual(values = c("darkgreen", "#0B4F02", "darkmagenta", "#5E0069"))+
  theme_classic()+
  theme(strip.text.x = element_text(size = 12),
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 10, colour = "black"),
        plot.margin = unit(c(0.5, 0.25, 0.5, 0.25), "cm"),
        legend.margin=margin(-10, 0, 0, 0),
        legend.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 10, colour = "black"),
        legend.position = "bottom")+
  labs(y = "Number of Detected MAGs", x = "", fill = "", colour = "")+
  ylim(0,200)+
  facet_wrap(~Time)
ggsave("figures/mag_community.pdf", mag_community, units = "cm", width = 17, height = 7)

mags_glmm <- glmer(n ~ Treatment * Time + (1|Name), family = poisson, data = pond_community)
summary(mags_glmm)

#bar plot

kraken_files <- list.files("data files/kraken_class", recursive = T, pattern = ".bracken", full.names = T)

all_kraken <- data_frame()
for(i in 1:length(kraken_files)){
  kraken_output <- read_tsv(kraken_files[i])
  kraken_output$Sample <- kraken_files[i] %>% substr(25,36)
  all_kraken <- rbind(all_kraken, kraken_output)
}

write.csv(all_kraken, "data files/all_kraken_class.csv", row.names = F)
all_kraken <- left_join(all_kraken, sample_info)
all_kraken$Time <- all_kraken$Pond %>% substr(4,4)
all_kraken$Time <- ifelse(all_kraken$Time == "1", "Day 0", "Day 28")

sample_list <- unique(all_kraken$Sample)

kraken_simple <- data.frame()
for(i in sample_list){
  sample_abun <- subset(all_kraken, Sample == i)
  sample_abun$other <- sum(sample_abun$fraction_total_reads[sample_abun$fraction_total_reads < 0.01])
  sample_abun <- sample_abun %>% add_row(Sample = i, fraction_total_reads = unique(sample_abun$other),
                                         name = "Other", Name = unique(sample_abun$Name), Time = unique(sample_abun$Time))
  kraken_simple <- rbind(kraken_simple, sample_abun)
}

kraken_simple <- subset(kraken_simple, fraction_total_reads >= 0.01 | name == "Other")

kraken_plot <- ggplot(kraken_simple, aes(x = Name, y = fraction_total_reads, fill = fct_reorder(name, fraction_total_reads)))+
  geom_col(position = "stack")+
  scale_fill_manual(values = c("#3e2137","#70377f","#17434b","#1f0e1c",
                               "#f5edba","#d9bd66", "#e4943a","#9a6348","#d79b7d",
                               "#584563","#8c8fae","#34859d","#7ec4c1",
                               "#647d34","#c0c741","#9d303b","#d26471"))+
  labs(x = "Pond", y = "Relative Abundance", fill = "Class")+
  theme_classic()+
  theme(legend.position = "bottom", legend.key.size = unit(0.25, "cm"),
        legend.key.spacing = unit(0.05, "cm"), legend.title.position = "top",
        legend.margin=margin(-10,0,0,0),
        strip.text.x = element_text(size = 12),
        axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title = element_text(size = 10, colour = "black"),
        legend.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 10, colour = "black"))+
  guides(fill = guide_legend(reverse=T))+
  facet_wrap(~Time, scales = "fixed")
ggsave("figures/class_community.pdf", kraken_plot, units = "cm", width = 18, height = 13)  


#figure 1
mags_pcoa <- ggarrange(mag_community, ggarrange(T1_species, T2_species, common.legend = TRUE, legend = "none"), kraken_plot,
                       nrow = 3, heights = c(1.15, 1, 1.5), common.legend = F, labels = c("A", "B", "C"))
ggsave("figures/mags_pcoa_bar.pdf", mags_pcoa, units = "cm", width = 17, height = 22)



# class_beta <- read_tsv("data files/beta_div_class_day_0.txt", skip = 16)
# rownames(class_beta) <- class_beta$x
# class_beta <- class_beta[, c(2:17)]
# class_beta <- t(class_beta)
# class_bray <- as.dist(class_beta)
# 
# adonis2(class_bray ~ as.factor(Treatment_Time), data = sample_info, permutations = 999, by = "onedf")
# adonis2(class_bray ~ as.factor(Treatment) + as.factor(Time) + as.factor(Treatment) : as.factor(Time), data = sample_info, permutations = 999, by = "onedf")
# 
# class_pcoa <- pcoa(class_bray, correction="none", rn=NULL)
# relative_eig_class <- class_pcoa$values$Relative_eig
# class_pcoa_coords <- as.data.frame(class_pcoa$vectors)
# class_pcoa_coords <- cbind(class_pcoa_coords, sample_info)
# 
# class_plot <- ggplot(class_pcoa_coords, aes(x = Axis.1, y = Axis.2, colour = Treatment_Time)) +
#   geom_point(size = 2)+
#   scale_colour_manual(values = c("#98BF64", "darkgreen", "#BE93D4", "darkmagenta"))+
#   labs(title = "Class Bray-Curtis Dissimilarity",
#        x = paste(paste("PC 1", signif((relative_eig_class[1] * 100), 4), sep = "  "),"%", sep = ""),
#        y = paste(paste("PC 2", signif((relative_eig_class[2] * 100), 4), sep = "  "),"%", sep = ""),
#        colour = "Treatment")+
#   theme_classic()+
#   theme(title = element_text(size = 10),
#         axis.text = element_text(size = 12, colour = "black"),
#         axis.title = element_text(size = 12, colour = "black"),
#         legend.text = element_text(size = 12, colour = "black"),
#         legend.title = element_text(size = 12, colour = "black"))


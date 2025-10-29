library(tidyverse)
library(ape)
library(vegan)
library(ggpubr)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

sample_info <- read_csv("data files/chapter_1_sample_names.csv")
sample_info$Treatment <- ifelse(grepl("CTRL", sample_info$Name), "Control", "GBH")
sample_info$Time <- ifelse(grepl("_1", sample_info$Pond), "Day 0", "Day 28")
sample_info$Treatment_Time <- paste(sample_info$Treatment, sample_info$Time, sep = " ")
sample_info <- subset(sample_info, !(Name == "CTRL E"))
T1_samples <- subset(sample_info, Time == "Day 0")
T2_samples <- subset(sample_info, Time == "Day 28")

day_0_beta_div_files <- list.files("data files", pattern = "day_0.txt", full.names = T)
day_28_beta_div_files <- list.files("data files", pattern = "day_28.txt", full.names = T)

pvals <- data.frame()
for(i in 1: length(day_0_beta_div_files)){
  kraken_beta <- read_tsv(day_0_beta_div_files[i], skip = 8)
  rownames(kraken_beta) <- kraken_beta$x
  kraken_beta  <- kraken_beta[, c(2:9)]
  kraken_beta  <- t(kraken_beta)
  kraken_bray <- as.dist(kraken_beta)
  permanova <- adonis2(kraken_bray ~ Treatment, data = T1_samples, permutations = 1000)
  level <- day_0_beta_div_files[i] %>% str_remove("data files/")
  level <- level %>% str_remove("beta_div_")
  level <- level %>% str_remove(".txt")
  result <- c(level, permanova$R2[1], permanova$`Pr(>F)`[1])
  pvals <- rbind(pvals, result)
}

for(i in 1: length(day_28_beta_div_files)){
  kraken_beta <- read_tsv(day_28_beta_div_files[i], skip = 8)
  rownames(kraken_beta) <- kraken_beta$x
  kraken_beta  <- kraken_beta[, c(2:9)]
  kraken_beta  <- t(kraken_beta)
  kraken_bray <- as.dist(kraken_beta)
  permanova <- adonis2(kraken_bray ~ Treatment, data = T2_samples, permutations = 1000)
  level <- day_28_beta_div_files[i] %>% str_remove("data files/")
  level <- level %>% str_remove("beta_div_")
  level <- level %>% str_remove(".txt")
  result <- c(level, permanova$R2[1], permanova$`Pr(>F)`[1])
  pvals <- rbind(pvals, result)

}

colnames(pvals) <- c("Level", "R2", "p_value")
pvals$padj <- p.adjust(pvals$p_value, method="BH")
#write.csv(pvals, "data files/permanovas.csv", row.names = F)

for(i in 1: length(day_0_beta_div_files)){
  level <- day_0_beta_div_files[i] %>% str_remove("data files/")
  level <- level %>% str_remove("beta_div_")
  level <- level %>% str_remove("_day_0.txt")
  
  kraken_beta <- read_tsv(day_0_beta_div_files[i], skip = 8)
  rownames(kraken_beta) <- kraken_beta$x
  kraken_beta  <- kraken_beta[, c(2:9)]
  kraken_beta  <- t(kraken_beta)
  kraken_bray <- as.dist(kraken_beta)
  
  kraken_pcoa <- pcoa(kraken_bray, correction="none", rn=NULL)
  relative_eig <- kraken_pcoa$values$Relative_eig
  pcoa_coords <- as.data.frame(kraken_pcoa$vectors)
  pcoa_coords <- cbind(pcoa_coords, T1_samples)

  pcoa_plot <- ggplot(pcoa_coords, aes(x = Axis.1, y = Axis.2, colour = Treatment_Time)) +
    geom_point(size = 2)+
    scale_colour_manual(values = c("#98BF64", "#BE93D4"))+
    labs(title = paste("Day 0", level, "level Bray-Curtis Dissimilarity", sep = " "),
         x = paste("PC 1", signif((relative_eig[1] * 100), 4), "%", sep = " "),
         y = paste("PC 2", signif((relative_eig[2] * 100), 4), "%", sep = " "),
         colour = "Treatment")+
    theme_classic()+
    theme(title = element_text(size = 10),
          axis.text = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 12, colour = "black"),
          legend.text = element_text(size = 12, colour = "black"),
          legend.title = element_text(size = 12, colour = "black"),
          legend.position = "bottom")

  ggsave(paste("figures/beta_div_pcoa_", level, "_day_0.pdf", sep = ""), pcoa_plot, limitsize = F, width = 6, height = 4)
}


for(i in 1: length(day_28_beta_div_files)){
  level <- day_28_beta_div_files[i] %>% str_remove("data files/")
  level <- level %>% str_remove("beta_div_")
  level <- level %>% str_remove("_day_28.txt")
  
  kraken_beta <- read_tsv(day_28_beta_div_files[i], skip = 8)
  rownames(kraken_beta) <- kraken_beta$x
  kraken_beta  <- kraken_beta[, c(2:9)]
  kraken_beta  <- t(kraken_beta)
  kraken_bray <- as.dist(kraken_beta)
  
  kraken_pcoa <- pcoa(kraken_bray, correction="none", rn=NULL)
  relative_eig <- kraken_pcoa$values$Relative_eig
  pcoa_coords <- as.data.frame(kraken_pcoa$vectors)
  pcoa_coords <- cbind(pcoa_coords, T2_samples)
  
  pcoa_plot <- ggplot(pcoa_coords, aes(x = Axis.1, y = Axis.2, colour = Treatment_Time)) +
    geom_point(size = 2)+
    scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
    labs(title = paste("Day 28", level, "level Bray-Curtis Dissimilarity", sep = " "),
         x = paste("PC 1", signif((relative_eig[1] * 100), 4), "%", sep = " "),
         y = paste("PC 2", signif((relative_eig[2] * 100), 4), "%", sep = " "),
         colour = "Treatment")+
    theme_classic()+
    theme(title = element_text(size = 10),
          axis.text = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 12, colour = "black"),
          legend.text = element_text(size = 12, colour = "black"),
          legend.title = element_text(size = 12, colour = "black"),
          legend.position = "bottom")
  
  ggsave(paste("figures/beta_div_pcoa_", level, "_day_28.pdf", sep = ""), pcoa_plot, limitsize = F, width = 6, height = 4)
}

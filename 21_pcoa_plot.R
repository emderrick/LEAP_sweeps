library(tidyverse)
library(ape)
library(ggpubr)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

#beta_div_files <- list.files("data files", pattern = "beta_div.csv", full.names = T)
sample_info <- read_csv("data files/chapter_1_sample_names.csv")
sample_info$Treatment <- ifelse(grepl("CTRL", sample_info$Name), "Control", "GBH")
sample_info$Time <- ifelse(grepl("_1", sample_info$Pond_Time), "Day 0", "Day 28")
sample_info$Treatment_Time <- paste(sample_info$Treatment, sample_info$Time, sep = " ")

class_beta <- read_csv("data files/class_beta_div.csv")
class_beta <- sapply(class_beta, as.numeric)
class_pcoa <- pcoa(class_beta, correction="none", rn=NULL)
relative_eig_class <- class_pcoa$values$Relative_eig

class_pcoa_coords <- as.data.frame(class_pcoa$vectors)
class_pcoa_coords <- cbind(class_pcoa_coords, sample_info)

class_plot <- ggplot(class_pcoa_coords, aes(x = Axis.1, y = Axis.2, colour = Treatment_Time)) +
  geom_point(size = 2)+
  scale_colour_manual(values = c("#98BF64", "darkgreen", "#BE93D4", "darkmagenta"))+
  labs(title = "Class level Bray-Curtis Dissimilarity",
       x = paste("PC 1", signif((relative_eig_class[1] * 100), 4), "%", sep = " "),
       y = paste("PC 2", signif((relative_eig_class[2] * 100), 4), "%", sep = " "),
       colour = "Treatment")+
  theme_classic()+
  theme(axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 12, colour = "black"),
        legend.position = "bottom")

species_beta <- read_csv("data files/species_beta_div.csv")
species_beta <- sapply(species_beta, as.numeric)
species_pcoa <- pcoa(species_beta, correction="none", rn=NULL)
relative_eig_species <- species_pcoa$values$Relative_eig

species_pcoa_coords <- as.data.frame(species_pcoa$vectors)
species_pcoa_coords <- cbind(species_pcoa_coords, sample_info)

species_plot <- ggplot(species_pcoa_coords, aes(x = Axis.1, y = Axis.2, colour = Treatment_Time)) +
  geom_point(size = 2)+
  scale_colour_manual(values = c("#98BF64", "darkgreen", "#BE93D4", "darkmagenta"))+
  labs(title = "Species level Bray-Curtis Dissimilarity",
       x = paste("PC 1", signif((relative_eig_species[1] * 100), 4), "%", sep = " "),
       y = paste("PC 2", signif((relative_eig_species[2] * 100), 4), "%", sep = " "),
       colour = "Treatment")+
  theme_classic()+
  theme(axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 12, colour = "black"),
        legend.position = "bottom")

class_species_pcoa_plot <- ggarrange(class_plot, species_plot, nrow=1, common.legend = TRUE, legend="bottom")
ggsave("figures/beta_div_pcoa_class_species.pdf", class_species_pcoa_plot, limitsize = F, width = 8.5, height = 3.5)


# for(i in 1: length(beta_div_files)){
#   level <- beta_div_files[i] %>% str_remove("data files/")
#   level <- level %>% str_remove("_beta_div.csv")
#   level <- level %>% str_to_title(level)
#   
#   kraken_beta <- read_csv(beta_div_files[i])
#   kraken_beta <- sapply(kraken_beta, as.numeric)
#   kraken_pcoa <- pcoa(kraken_beta, correction="none", rn=NULL)
#   relative_eig <- kraken_pcoa$values$Relative_eig
#   
#   pcoa_coords <- as.data.frame(kraken_pcoa$vectors)
#   pcoa_coords <- cbind(pcoa_coords, sample_info)
#   
#   pcoa_plot <- ggplot(pcoa_coords, aes(x = Axis.1, y = Axis.2, colour = Treatment_Time)) +
#     geom_point(size = 2)+
#     scale_colour_manual(values = c("#98BF64", "darkgreen", "#BE93D4", "darkmagenta"))+
#     labs(title = paste(level, "level Bray-Curtis Dissimilarity", sep = " "),
#          x = paste("PC 1", signif((relative_eig[1] * 100), 4), "%", sep = " "),
#          y = paste("PC 2", signif((relative_eig[2] * 100), 4), "%", sep = " "),
#          colour = "Treatment")+
#     theme_classic()+
#     theme(axis.text = element_text(size = 12, colour = "black"),
#           axis.title = element_text(size = 12, colour = "black"),
#           legend.position = "bottom")
#     
#   ggsave(paste("figures/beta_div_pcoa_", level, ".pdf", sep = ""), pcoa_plot, limitsize = F, width = 6, height = 4)
# }







library(tidyverse)
library(ggpubr)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

mag_list <- c("MAG_00097_1", "MAG_00110_1", "MAG_00179_1", "MAG_00194_1", "MAG_00197_1", "MAG_00201_1","MAG_00674_1")

all_samples <- read_csv("data files/chapter_1_sample_names.csv")
all_samples$Time <- ifelse(grepl("_1", all_samples$Pond), "Day 0", "Day 28")
all_samples <- subset(all_samples, !(Name == "CTRL E"))

all_sum <- read_csv("data files/T1_subsamp_SNV_summary_MAG.csv")
all_sum <- subset(all_sum, mag %in% mag_list)
all_sum$Name_Time <- paste(all_sum$Name, all_sum$Time, sep = " ")

plot_labels <- data.frame(mag = c(unique(all_sum$mag)),
                          species = c("Sphingorhabdus_B sp021298455", "Erythrobacter sp027486035", "JAIXWP01 sp027491075",
                                      "Sphingorhabdus_B sp027488785", "Prosthecobacter sp.", "Brevundimonas sp027487935",
                                      "CAMBQF01 sp."))
all_sum <- left_join(all_sum, plot_labels)

snv_mag <- ggplot(all_sum, aes(x = Name, y = SNVs_Mbp, fill = Treatment))+
  geom_col()+
  theme_bw()+
  scale_fill_manual(values = c("darkgreen", "darkmagenta"))+
  theme(axis.text = element_text(colour = "black", size = 7),
        axis.title = element_blank(),
        title = element_text(size = 7),
        strip.text.x = element_text(size = 7),
        strip.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "none")+
  labs(title = "Polymorphic sites per Mbp")+
  facet_grid(species ~ Time, scales = "free", labeller = labeller(species = label_wrap_gen(width = 18)))

sns_mag <- ggplot(all_sum, aes(x = Name, y = SNSs_Mbp, fill = Treatment))+
  geom_col()+
  theme_bw()+
  scale_fill_manual(values = c("darkgreen", "darkmagenta"))+
  theme(axis.text = element_text(colour = "black", size = 7),
        axis.title = element_blank(),
        title = element_text(size = 7),
        strip.text = element_text(size = 7),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "none")+
  labs(title = "Fixed substitutions per Mbp")+
  facet_grid(species ~ Time, scales = "free", labeller = labeller(species = label_wrap_gen(width = 20)))

mag_plot <- ggarrange(snv_mag, sns_mag, common.legend = T, legend = "none", ncol = 2, labels = c("A", "B"))
ggsave("figures/mag_sns_snv.pdf", mag_plot, units = "cm", width = 17, height = 20)

all_sum %>% group_by(species, mag, Time) %>% count()


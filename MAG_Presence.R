library(tidyverse)
setwd("/Users/emma/Documents/")

T1_bin_stats <- read_csv("bin_stats.csv")
colnames(T1_bin_stats)[1] <- "Bin"
T1_bin_stats$Bin <- str_remove(T1_bin_stats$Bin, ".fa")
bins_checkm <- read_tsv("MAGs_checkM.tsv")
colnames(bins_checkm)[1] <- "Bin"
T1_bins <- full_join(T1_bin_stats, bins_checkm)

good_bins <- subset(T1_bin_stats, Completeness >= 70 & Contamination < 10)
ggplot(T1_bin_stats, aes(y = Contamination, x = Completeness))+
  geom_point()


MAG_abun <- read_tsv("T1_MAGs_coverM.tsv")
sample_names <- read_csv("chapter_1_sample_names.csv")
MAG_abun$Sample <- MAG_abun$Sample %>% substr(1,12)
MAG_abun <- left_join(MAG_abun, sample_names)

MAG_abun$Pond <- MAG_abun$Pond_Time %>% substr(1,2)
MAG_abun$Time <- MAG_abun$Pond_Time %>% substr(4,4)
MAG_abun$present <- with(MAG_abun, ifelse(Mean >= 5 & `Covered Fraction` >= 0.7, 1, 0))
MAG_abun$Treatment <- ifelse(grepl("CTRL", MAG_abun$Name), "1", "2")
MAG_abun$Coverage <- as.integer(MAG_abun$Mean)

mag_pres <- ggplot(subset(MAG_abun, present == 1), aes(x = Time, y = fct_reorder(Name, desc(Treatment)), colour = Treatment))+
  geom_text(aes(label=Coverage), size = 8)+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text = element_text(size = 12),
        legend.position = "none")+
  facet_wrap(~Genome)
ggsave("T1_MAG_pres_0.7.pdf", limitsize = F, width=48, height=48)


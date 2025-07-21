library(tidyverse)
library(ggplot2)
library(RColorBrewer)
options(scipen=999)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

sample_names <- read_csv("data files/chapter_1_sample_names.csv")

metaphlan <- read_tsv("data files/QC_merged_metaphlan_abundance_table.tsv", skip = 1)
metaphlan_long <- pivot_longer(metaphlan, cols = c(2:19), names_to = "Sample", values_to = "Abundance")
metaphlan_long$Sample <- metaphlan_long$Sample %>% substr(1,12)
metaphlan_long <- left_join(metaphlan_long, sample_names)
metaphlan_long$Time <- metaphlan_long$Pond_Time %>% substr(4,4)
metaphlan_long$Name_Time <- paste(metaphlan_long$Name, metaphlan_long$Time, sep = "_")
metaphlan_long$Time <- with(metaphlan_long, ifelse(Time == "1", "Day 0", "Day 7"))
sample_list <- unique(metaphlan_long$Name_Time)

metaphlan_rel <- data.frame()
for(i in sample_list){
  sample_abun <- subset(metaphlan_long, Name_Time == i)
  sample_abun$bacterial_abundance <- sample_abun$Abundance[sample_abun$clade_name == "k__Bacteria"]
  metaphlan_rel <- rbind(metaphlan_rel, sample_abun)
}

metaphlan_rel$relative_abundance <- metaphlan_rel$Abundance / metaphlan_rel$bacterial_abundance

metaphlan_bacteria <- subset(metaphlan_rel, !(str_detect(clade_name, "p__")))
metaphlan_phyla <- subset(metaphlan_rel, str_detect(clade_name, "p__") & !(str_detect(clade_name, "s__")) & !(str_detect(clade_name, "g__")) & !(str_detect(clade_name, "o__")) & !(str_detect(clade_name, "c__")))
metaphlan_phyla$clade_name <- metaphlan_phyla$clade_name %>% str_sub(16)

metaphlan_phyla_simple <- data.frame()
for(i in sample_list){
  sample_abun <- subset(metaphlan_phyla, Name_Time == i)
  sample_abun$other <- sum(sample_abun$relative_abundance[sample_abun$relative_abundance < 0.01])
  sample_abun <- sample_abun %>% add_row(Name_Time = i, relative_abundance = unique(sample_abun$other),
                                         clade_name = "Other", Name = unique(sample_abun$Name), Time = unique(sample_abun$Time))
  metaphlan_phyla_simple <- rbind(metaphlan_phyla_simple, sample_abun)
}

metaphlan_phyla_simple <- subset(metaphlan_phyla_simple, relative_abundance >= 0.01)

metaphlan_phyla_plot <- ggplot(metaphlan_phyla_simple, aes(x = Name, y = relative_abundance, fill = reorder(clade_name, relative_abundance)))+
  geom_col(position = "stack")+
  labs(x = "Pond", y = "Relative Abundance", fill = "Phylum")+  
  scale_x_discrete(expand = c(0, 0))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  scale_fill_brewer(palette = "Paired", direction = -1)+
  theme_bw()+
  facet_wrap(~Time)

ggsave("figures/QC_metaphlan_phyla.pdf", metaphlan_phyla_plot, limitsize = F, width = 16, height = 7)

metaphlan_class <- subset(metaphlan_rel, str_detect(clade_name, "c__") & !(str_detect(clade_name, "s__")) & !(str_detect(clade_name, "g__")) & !(str_detect(clade_name, "o__")))

sample_list <- unique(metaphlan_class$Name_Time)
metaphlan_class_simple <- data.frame()

for(i in sample_list){
  sample_abun <- subset(metaphlan_class, Name_Time == i)
  sample_abun$other <- sum(sample_abun$relative_abundance[sample_abun$relative_abundance <= 0.05])
  sample_abun <- sample_abun %>% add_row(Name_Time = i, relative_abundance = unique(sample_abun$other),
                                         clade_name = "Other", Name = unique(sample_abun$Name), Time = unique(sample_abun$Time))
  metaphlan_class_simple <- rbind(metaphlan_class_simple, sample_abun)
}

metaphlan_class_simple <- subset(metaphlan_class_simple, relative_abundance > 0.05)

metaphlan_class_plot <- ggplot(metaphlan_class_simple, aes(x = Name, y = relative_abundance, fill = clade_name))+
  geom_col(position = "stack")+
  labs(x = "Pond", y = "Relative Abundance", fill = "Class")+  
  scale_x_discrete(expand = c(0, 0))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  theme_bw()+
  theme(panel.spacing = unit(0.5,"cm"))+
  facet_wrap(~Time)
ggsave("figures/QC_metaphlan_class.pdf", metaphlan_class_plot, limitsize = F, width = 24, height = 10)



metaphlan_species <- subset(metaphlan_rel, str_detect(clade_name, "s__"))
metaphlan_genera <- subset(metaphlan_rel, str_detect(clade_name, "g__") & !(str_detect(clade_name, "s__")))
metaphlan_order <- subset(metaphlan_rel, str_detect(clade_name, "o__") & !(str_detect(clade_name, "s__")) & !(str_detect(clade_name, "g__")))

metaphlan_order_plot <- ggplot(metaphlan_order, aes(x = Name, y = relative_abundance, fill = clade_name))+
  geom_col(position = "stack")+
  theme(legend.position = "bottom")+
  facet_wrap(~Time)
ggsave("figures/QC_metaphlan_order.pdf", metaphlan_order_plot, limitsize = F, width = 36, height = 24)

metaphlan_genera_plot <- ggplot(metaphlan_genera, aes(x = Name, y = relative_abundance, fill = clade_name))+
  geom_col(position = "stack")+
  theme(legend.position = "bottom")+
  facet_wrap(~Time)
ggsave("figures/QC_metaphlan_genera.pdf", metaphlan_genera_plot, limitsize = F, width = 36, height = 24)

metaphlan_species_plot <- ggplot(metaphlan_species, aes(x = Name, y = relative_abundance, fill = clade_name))+
  geom_col(position = "stack")+
  theme(legend.position = "none")+
  facet_wrap(~Time)
ggsave("figures/QC_metaphlan_species.pdf", metaphlan_species_plot, limitsize = F, width = 16, height = 7)

### bracken

bracken_list <- list.files(path = "data files/phylum_bracken/", pattern = "phylum.bracken",full.names = T)

bracken_phyla <- data.frame()
for(bracken_file in bracken_list){
  bracken_abun <- read_tsv(bracken_file)
  bracken_abun$Sample <- bracken_file %>% substr(28,39)
  bracken_phyla <- rbind(bracken_phyla, bracken_abun)
}

bracken_phyla <- left_join(bracken_phyla, sample_names)
bracken_phyla$Time <-bracken_phyla$Pond_Time %>% substr(4,4)
bracken_phyla$Name_Time <- paste(bracken_phyla$Name, bracken_phyla$Time, sep = "_")

sample_list <- unique(bracken_phyla$Name_Time)

bracken_phyla_simple <- data.frame()
for(i in sample_list){
  sample_abun <- subset(bracken_phyla, Name_Time == i)
  sample_abun$other <- sum(sample_abun$fraction_total_reads[sample_abun$fraction_total_reads < 0.01])
  sample_abun <- sample_abun %>% add_row(Name_Time = i, fraction_total_reads = unique(sample_abun$other),
                                         name = "Other", Name = unique(sample_abun$Name), Time = unique(sample_abun$Time))
  bracken_phyla_simple <- rbind(bracken_phyla_simple, sample_abun)
}

bracken_phyla_simple <- subset(bracken_phyla_simple, fraction_total_reads >= 0.01 | name == "Other")

bracken_phyla_plot <- ggplot(bracken_phyla_simple , aes(x = Name, y = fraction_total_reads, fill = reorder(name, fraction_total_reads)))+
  geom_col(position = "stack")+
  labs(x = "Pond", y = "Relative Abundance", fill = "Phylum")+  
  scale_fill_brewer(palette = "Paired", direction = -1)+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  theme_bw()+
  facet_wrap(~Time)

ggsave("figures/QC_bracken_phyla.pdf", bracken_phyla_plot, limitsize = F, width = 16, height = 7)


bracken_class_list <- list.files(path = "data files/class_bracken/", pattern = "class.bracken", full.names = T)

bracken_class <- data.frame()
for(bracken_file in bracken_class_list){
  bracken_abun <- read_tsv(bracken_file)
  bracken_abun$Sample <- bracken_file %>% substr(27,38)
  bracken_class <- rbind(bracken_class, bracken_abun)
}

bracken_class <- left_join(bracken_class, sample_names)
bracken_class$Time <-bracken_class$Pond_Time %>% substr(4,4)
bracken_class$Name_Time <- paste(bracken_class$Name, bracken_class$Time, sep = "_")

sample_list <- unique(bracken_class$Name_Time)
bracken_class_simple <- data.frame()

for(i in sample_list){
  sample_abun <- subset(bracken_class, Name_Time == i)
  sample_abun$other <- sum(sample_abun$fraction_total_reads[sample_abun$fraction_total_reads < 0.03])
  sample_abun <- sample_abun %>% add_row(Name_Time = i, fraction_total_reads = unique(sample_abun$other),
                                         name = "Other", Name = unique(sample_abun$Name), Time = unique(sample_abun$Time))
  bracken_class_simple <- rbind(bracken_class_simple, sample_abun)
}

bracken_class_simple <- subset(bracken_class_simple, fraction_total_reads >= 0.03 | name == "Other")
bracken_class_simple$Time <- with(bracken_class_simple, ifelse(Time == "1", "Day 0", "Day 28"))

bracken_class_plot <- ggplot(bracken_class_simple , aes(x = Name, y = fraction_total_reads, fill = reorder(name, fraction_total_reads)))+
  geom_col(position = "stack")+
  labs(x = "Pond", y = "Relative Abundance", fill = "Class")+  
  scale_x_discrete(expand = c(0, 0))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  scale_fill_brewer(palette = "Paired", direction = -1)+
  theme_bw()+
  theme(text = element_text(size = 18), axis.text = element_text(colour = "black"), panel.spacing = unit(0.5,"cm"))+
  facet_wrap(~Time)
ggsave("figures/QC_bracken_class.pdf", bracken_class_plot, limitsize = F, width = 18, height = 6)


bracken_order_list <- list.files(path = "data files/bracken_order/", pattern = "order.bracken", full.names = T)

bracken_order <- data.frame()
for(bracken_file in bracken_order_list){
  bracken_abun <- read_tsv(bracken_file)
  bracken_abun$Sample <- bracken_file %>% substr(27,38)
  bracken_order <- rbind(bracken_order, bracken_abun)
}

bracken_order <- left_join(bracken_order, sample_names)
bracken_order$Time <-bracken_order$Pond_Time %>% substr(4,4)
bracken_order$Name_Time <- paste(bracken_order$Name, bracken_order$Time, sep = "_")

sample_list <- unique(bracken_order$Name_Time)
bracken_order_simple <- data.frame()

for(i in sample_list){
  sample_abun <- subset(bracken_order, Name_Time == i)
  sample_abun$other <- sum(sample_abun$fraction_total_reads[sample_abun$fraction_total_reads < 0.05])
  sample_abun <- sample_abun %>% add_row(Name_Time = i, fraction_total_reads = unique(sample_abun$other),
                                         name = "Other", Name = unique(sample_abun$Name), Time = unique(sample_abun$Time))
  bracken_order_simple <- rbind(bracken_order_simple, sample_abun)
}

bracken_order_simple <- subset(bracken_order_simple, fraction_total_reads > 0.05| name == "Other")

bracken_order_plot <- ggplot(bracken_order_simple , aes(x = Name, y = fraction_total_reads, fill = name))+
  geom_col(position = "stack")+
  labs(x = "Pond", y = "Relative Abundance", fill = "Order")+  
  scale_x_discrete(expand = c(0, 0))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  theme_bw()+
  facet_wrap(~Time)

ggsave("figures/QC_bracken_order.pdf", bracken_order_plot, limitsize = F, width = 16, height = 7)


bracken_genera_list <- list.files(path = "data files/bracken_genera/", pattern = "genus.bracken",full.names = T)

bracken_genera <- data.frame()
for(bracken_file in bracken_genera_list){
  bracken_abun <- read_tsv(bracken_file)
  bracken_abun$Sample <- bracken_file %>% substr(28,39)
  bracken_genera <- rbind(bracken_genera, bracken_abun)
}

bracken_genera <- left_join(bracken_genera, sample_names)
bracken_genera$Time <-bracken_genera$Pond_Time %>% substr(4,4)
bracken_genera$Name_Time <- paste(bracken_genera$Name, bracken_genera$Time, sep = "_")

sample_list <- unique(bracken_genera$Name_Time)
bracken_genera_simple <- data.frame()

for(i in sample_list){
  sample_abun <- subset(bracken_genera, Name_Time == i)
  sample_abun$other <- sum(sample_abun$fraction_total_reads[sample_abun$fraction_total_reads < 0.05])
  sample_abun <- sample_abun %>% add_row(Name_Time = i, fraction_total_reads = unique(sample_abun$other),
                                         name = "Other", Name = unique(sample_abun$Name), Time = unique(sample_abun$Time))
  bracken_genera_simple <- rbind(bracken_genera_simple, sample_abun)
}

bracken_genera_simple <- subset(bracken_genera_simple, fraction_total_reads > 0.01)

bracken_genera_plot <- ggplot(bracken_genera_simple , aes(x = Name, y = fraction_total_reads, fill = name))+
  geom_col(position = "stack")+
  labs(x = "Pond", y = "Relative Abundance", fill = "Genus")+  
  scale_x_discrete(expand = c(0, 0))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  theme_bw()+
  facet_wrap(~Time)

ggsave("figures/QC_bracken_genera.pdf", bracken_genera_plot, limitsize = F, width = 16, height = 7)


### mOTUs

# motus <- read_tsv("data files/merged_sub.motus", skip = 2)
# motus_long <- pivot_longer(motus, cols = c(2:19), names_to = "Sample", values_to = "Abundance")
# motus_long <- left_join(motus_long, sample_names)
# motus_long$Time <- motus_long$Pond_Time %>% substr(4,4)
# motus_long$Name_Time <- paste(motus_long$Name, motus_long$Time, sep = "_")
# colnames(motus_long)[1] <- "Phyla"
# motus_long <- subset(motus_long, Abundance > 0)
# 
# motus_list <- unique(motus_long$Name_Time)
# 
# motus_rel <- data.frame()
# for(i in motus_list){
#   sample_abun <- subset(motus_long, Name_Time == i)
#   sample_abun$unclassified <- sample_abun$Abundance[sample_abun$Phyla == "unassigned"]
#   motus_rel <- rbind(motus_rel, sample_abun)
# }
# 
# motus_rel$classified_abundance <- motus_rel$Abundance / (1 - motus_rel$unclassified)
# 
# 
# motus_phyla <- ggplot(subset(motus_rel, !(Phyla == "unassigned")), aes(x = Name, y = classified_abundance, fill = reorder(Phyla, classified_abundance)))+
#   geom_col(position = "stack")+
#   labs(x = "Pond", y = "Relative Abundance", fill = "Phylum")+
#   facet_wrap(~Time)
# 
# ggsave("figures/mOTUs_phyla.pdf", motus_phyla, limitsize = F, width = 16, height = 7)
  


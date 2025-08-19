library(tidyverse)
library(ggplot2)
library(RColorBrewer)
options(scipen=999)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

sample_names <- read_csv("data files/chapter_1_sample_names.csv")

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
bracken_phyla$Time <- ifelse(grepl("1", bracken_phyla$Name_Time), "Day 0", "Day 28")
bracken_phyla$Treatment <- ifelse(grepl("CTRL", bracken_phyla$Name_Time), "Control", "GBH")
bracken_phyla$Treatment_Time <- paste(bracken_phyla$Treatment, bracken_phyla$Time, sep = " ")

sample_list <- unique(bracken_phyla$Name_Time)

bracken_phyla_simple <- data.frame()
for(i in sample_list){
  sample_abun <- subset(bracken_phyla, Name_Time == i)
  sample_abun$other <- sum(sample_abun$fraction_total_reads[sample_abun$fraction_total_reads < 0.01])
  sample_abun <- sample_abun %>% add_row(Name_Time = i, fraction_total_reads = unique(sample_abun$other),
                                         name = "Other", Name = unique(sample_abun$Name), Time = unique(sample_abun$Time),
                                         Treatment = unique(sample_abun$Treatment), Treatment_Time = unique(sample_abun$Treatment_Time)) 
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
  facet_wrap(~Treatment_Time, scales = "free_x", nrow = 1)

ggsave("figures/QC_bracken_phyla.pdf", bracken_phyla_plot, limitsize = F, width = 16, height = 7)

#class

bracken_class_list <- list.files(path = "data files/class_bracken/", pattern = "class.bracken",full.names = T)

bracken_class <- data.frame()
for(bracken_file in bracken_class_list){
  bracken_abun <- read_tsv(bracken_file)
  bracken_abun$Sample <- bracken_file %>% substr(27,38)
  bracken_class <- rbind(bracken_class, bracken_abun)
}

bracken_class <- left_join(bracken_class, sample_names)
bracken_class$Time <-bracken_class$Pond_Time %>% substr(4,4)
bracken_class$Name_Time <- paste(bracken_class$Name, bracken_class$Time, sep = "_")
bracken_class$Time <- ifelse(grepl("1", bracken_class$Name_Time), "Day 0", "Day 28")
bracken_class$Treatment <- ifelse(grepl("CTRL", bracken_class$Name_Time), "Control", "GBH")
bracken_class$Treatment_Time <- paste(bracken_class$Time, bracken_class$Treatment, sep = " ")

sample_list <- unique(bracken_class$Name_Time)

bracken_class_simple <- data.frame()
for(i in sample_list){
  sample_abun <- subset(bracken_class, Name_Time == i)
  sample_abun$other <- sum(sample_abun$fraction_total_reads[sample_abun$fraction_total_reads < 0.03])
  sample_abun <- sample_abun %>% add_row(Name_Time = i, fraction_total_reads = unique(sample_abun$other),
                                         name = "Other", Name = unique(sample_abun$Name), Time = unique(sample_abun$Time),
                                         Treatment = unique(sample_abun$Treatment), Treatment_Time = unique(sample_abun$Treatment_Time)) 
  bracken_class_simple <- rbind(bracken_class_simple, sample_abun)
}

bracken_class_simple <- subset(bracken_class_simple, fraction_total_reads >= 0.03 | name == "Other")

bracken_class_plot <- ggplot(bracken_class_simple , aes(x = Name, y = fraction_total_reads, fill = reorder(name, fraction_total_reads)))+
  geom_col(position = "stack")+
  labs(x = "Pond", y = "Relative Abundance", fill = "Class")+  
  scale_fill_brewer(palette = "Paired", direction = -1)+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  theme_bw()+
  theme(panel.spacing = unit(0.5, "cm"), text = element_text(size = 16), axis.text = element_text(colour = "black"))+
  facet_wrap(~Treatment_Time, scales = "free_x", nrow = 1)

ggsave("figures/QC_bracken_class.pdf", bracken_class_plot, limitsize = F, width = 18, height = 7)

#order

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

#generna

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
  


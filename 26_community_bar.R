library(tidyverse)
library(RColorBrewer)

options(scipen=999)
setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

sample_names <- read_csv("data files/chapter_1_sample_names.csv")
kraken_files <- list.files("data files/kraken_class", recursive = T, pattern = ".bracken", full.names = T)

all_kraken <- data_frame()
for(i in 1:length(kraken_files)){
  kraken_output <- read_tsv(kraken_files[i])
  kraken_output$Sample <- kraken_files[i] %>% substr(25,36)
  all_kraken <- rbind(all_kraken, kraken_output)
}

write.csv(all_kraken, "data files/all_kraken_class.csv", row.names = F)
all_kraken <- left_join(all_kraken, sample_names)
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
                               "#c0c741","#647d34","#9d303b","#d26471"))+
  labs(x = "Pond", y = "Realtive Abundance", fill = "Class")+
  theme_classic()+
  theme(legend.position = "bottom", legend.key.size = unit(0.4, "cm"),
        axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title = element_text(size = 10, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        legend.title = element_text(size = 10, colour = "black"))+
  guides(fill = guide_legend(reverse=T))+
  facet_wrap(~Time, scales = "fixed")
ggsave("figures/class_community.pdf", kraken_plot, units = "cm", width = 18, height = 13)  
  
  
  
  
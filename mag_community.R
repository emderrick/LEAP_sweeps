library(tidyverse)
library(ggplot2)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

mags <- read_csv("refined data files/MAG_rel_abun.csv")
mags$Name_Time <- paste(mags$Name, mags$Time, sep = " ")
mags$Time <- ifelse(grepl("1", mags$Name_Time), "Day 0", "Day 28")

pond_community <- subset(mags, Relative_Abundance >= 0.01)
pond_community <- pond_community %>% group_by(Name, Time, Treatment) %>% count()
pond_community$Treatment_Time <- paste(pond_community$Treatment, pond_community$Time, sep = " ")


mag_community <- ggplot(pond_community, aes(x = Treatment, group = Treatment, colour = Treatment, fill = Treatment, y = n))+
  geom_boxplot()+
  geom_point()+
  scale_fill_manual(values = c("#98BF64", "#BE93D4"))+
  scale_colour_manual(values = c("darkgreen", "darkmagenta"))+
  theme_classic()+
  labs(y = "MAGs with > 0.01% Relative Abundance")+
  ylim(0,200)+
  facet_wrap(~Time)
ggsave("refined figures/mag_community.pdf", mag_community, limitsize = F, width = 8.5, height = 3.5)

dunn_test(n ~ Treatment_Time, data = pond_community)

sample_names <- read_csv("data files/chapter_1_sample_names.csv")
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
  theme(panel.spacing = unit(0.25, "cm"),
        text = element_text(size = 10),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(size = 6))+
  facet_wrap(~Time, scales = "free_x", nrow = 1)

ggsave("refined figures/community_class.pdf", bracken_class_plot, limitsize = F, width = 8.5, height = 3.5)

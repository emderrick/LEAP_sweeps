library(tidyverse)


total_reads <- read_csv("read_stats.csv")
total_reads <- as.data.frame(apply(total_reads, 2, function(x) gsub('\\s+', '', x)))
total_reads$pond <- total_reads$file %>% substr(1,2)
total_reads$time <- total_reads$file %>% substr(4,4)
total_reads <- total_reads %>% mutate(pulse = case_when(time == 1 ~ 0, time == 2 ~ 1, time == 3 ~ 1, time == 4 ~ 2, time == 5 ~ 2))
total_reads$timepoint <- paste(total_reads$pond, "_pulse", total_reads$pulse, sep = "")
total_reads <- subset(total_reads, pulse == 1)
total_reads_unique <- total_reads[, c(12, 4)] %>% unique()

total_reads_sum <- total_reads_unique %>% group_by(timepoint) %>% summarize(reads = sum(as.numeric(num_seqs)))
total_reads_sum$least <- min(total_reads_sum$reads)
total_reads_sum$percent <- total_reads_sum$least / total_reads_sum$reads


indiv_reads <- total_reads_unique
indiv_reads$num_seqs <- as.numeric(indiv_reads$num_seqs)
indiv_reads$least <- min(indiv_reads$num_seqs)
indiv_reads$percent <- indiv_reads$least / indiv_reads$num_seqs

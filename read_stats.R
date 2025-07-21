library(tidyverse)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

paired_read_stats <- read_tsv("data files/paired_read_stats.tsv")

paired_read_stats$lowest_depth <- min(paired_read_stats$num_seqs)
paired_read_stats$prop_kept <- paired_read_stats$lowest_depth / paired_read_stats$num_seqs


QC_read_stats <- read_tsv("data files/QC_read_stats.tsv")

QC_read_stats$lowest_depth <- min(QC_read_stats$num_seqs)
QC_read_stats$prop_kept <- QC_read_stats$lowest_depth / QC_read_stats$num_seqs

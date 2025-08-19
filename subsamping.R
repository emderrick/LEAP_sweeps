library(tidyverse)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")

mag_contigs <- read_tsv("refined data files/T1_refined.stb", col_names = c("contig", "mag"))
mag_contigs$mag <- mag_contigs$mag %>% str_remove(".fa")
mag_list <- (unique(mag_contigs$mag))

for(i in 1:length(mag_list)){
  contigs <- subset(mag_contigs, mag == mag_list[i])
  contigs <- contigs[, c(1)]
  write_tsv(contigs, paste("refined data files/", mag_list[i], "_contigs.txt", sep = ""), col_names = F)
}

sample_names <- read_csv("refined data files/chapter_1_sample_names.csv")
sample_names$Name_Time <- paste(sample_names$Name, sample_names$Pond_Time %>% str_sub(-1))
mag_SNVs <- read_csv("refined data files/T1_refined_SNV_summary_MAG.csv")

mag_SNVs_5x <- subset(mag_SNVs, mag_coverage >= 5 & mag_breadth >= 0.5)
mag_SNVs_5x <- left_join(mag_SNVs_5x, sample_names[, c(1:4)])
mag_SNVs_5x <- mag_SNVs_5x[, c(1,21,3)]
mag_SNVs_5x$cov_want <- 5
mag_SNVs_5x$prop_keep <- mag_SNVs_5x$cov_want / mag_SNVs_5x$mag_coverage

args <- as.data.frame(paste("cat ", mag_SNVs_5x$mag, "_contigs.txt | xargs samtools view -bh ", mag_SNVs_5x$Sample, "_P_T1_refined.bam", " > ", mag_SNVs_5x$Sample, "_", mag_SNVs_5x$mag, ".bam", sep = ""))
write_tsv(args, "refined data files/samtools_mag_args.tsv", col_names = F) 

sub_args <- as.data.frame(paste("samtools view -s ", mag_SNVs_5x$prop_keep, " -b ", mag_SNVs_5x$Sample, "_", mag_SNVs_5x$mag, ".bam", " > ", mag_SNVs_5x$Sample, "_", mag_SNVs_5x$mag, "_sub.bam", sep = ""))
write_tsv(sub_args, "refined data files/samtools_subsample_args.tsv", col_names = F) 

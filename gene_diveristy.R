library(tidyverse)
library(tidytext)
library(dplyr)
library(ggplot2)
library(cowplot)

setwd("/Users/Emma/Documents/manuscript version/")

mag_files <- list.files(pattern = ".*_gene_matrix_Aug12.csv", full.names = T)

for (mag_file in 1: length(mag_files)){
  
  MAG <- mag_files[mag_file] %>% substr(3,14)
  gene_matrix <- read_csv(mag_files[mag_file])

  Control_Control <- gene_matrix %>% select(matches("gene|Control.*Control.*pass"))
  Control_Control_pass <- Control_Control[rowSums(is.na(Control_Control)) != ncol(Control_Control) - 2, ]
  GBH_GBH <- gene_matrix %>% select(matches("gene|GBH.*GBH.*pass"))
  GBH_GBH_pass <- GBH_GBH[rowSums(is.na(GBH_GBH)) != ncol(GBH_GBH) - 2, ]
  Control_GBH <- gene_matrix %>% select(matches("gene|Control.*GBH.*pass"))
  Control_GBH_pass <- Control_GBH[rowSums(is.na(Control_GBH)) != ncol(Control_GBH) - 2, ]
  Control_GBH_pass <- subset(Control_GBH_pass, !(gene %in% Control_Control_pass$gene | gene %in% GBH_GBH_pass$gene))
  
  Control_GBH_pass_long <- pivot_longer(Control_GBH_pass, cols = -c("gene", "gene_length"), names_to = "comparison", values_to = "direction", values_drop_na = T)
  Control_GBH_pass_long <- subset(Control_GBH_pass_long, select = c("gene", "direction")) %>% unique()
  increase <- subset(Control_GBH_pass_long, direction == "increase")
  decrease <- subset(Control_GBH_pass_long, direction == "decrease")
  
  Control_GBH_pass_parallel <- na.omit(Control_GBH_pass)
  par_increase <- Control_GBH_pass_parallel[rowSums(Control_GBH_pass_parallel == "increase") == ncol(Control_GBH_pass_parallel) - 2, ]
  par_increase$direction <- "increase"
  par_increase <- subset(par_increase, select = c("gene", "direction"))
  par_decrease <- Control_GBH_pass_parallel[rowSums(Control_GBH_pass_parallel == "decrease") == ncol(Control_GBH_pass_parallel) - 2, ]
  par_decrease$direction <- "decrease"
  par_decrease <- subset(par_decrease, select = c("gene", "direction"))
  
  print(MAG)
  print(nrow(Control_GBH_pass_parallel))
  print(nrow(par_increase))
  print(nrow(par_decrease))
  write.csv(increase, paste(MAG, "_increase_pass.csv", sep = ""), row.names = F)
  write.csv(decrease, paste(MAG, "_decrease_pass.csv", sep = ""), row.names = F)
  write.csv(par_increase, paste(MAG, "_increase_par_pass.csv", sep = ""), row.names = F)
  write.csv(par_decrease, paste(MAG, "_decrease_par_pass.csv", sep = ""), row.names = F)

}


inc_gene_files <- list.files(pattern = ".*increase_pass.csv", full.names = T)

inc_sig_genes <- data.frame()
for (gene_file in 1: length(inc_gene_files)){
  
  sig_genes <- read_csv(inc_gene_files[gene_file])
  sig_genes <- subset(sig_genes, select = c(gene, direction))
  sig_genes$mag <- sig_genes$gene %>% substr(1,12)
  inc_sig_genes <- rbind(inc_sig_genes, sig_genes)
  
}

write.csv(inc_sig_genes, "increase_sig_genes.csv", row.names = F)

inc_par_gene_files <- list.files(pattern = ".*increase_par_pass.csv", full.names = T)

inc_par_sig_genes <- data.frame()
for (par_gene_file in 1: length(inc_par_gene_files)){
  
  par_sig_genes <- read_csv(inc_par_gene_files[par_gene_file])
  par_sig_genes <- subset(par_sig_genes, select = c(gene, direction))
  par_sig_genes$mag <- par_sig_genes$gene %>% substr(1,12)
  inc_par_sig_genes <- rbind(inc_par_sig_genes, par_sig_genes)
  
}

write.csv(inc_par_sig_genes, "parallel_increase_sig_genes.csv", row.names = F)


dec_gene_files <- list.files(pattern = ".*decrease_pass.csv", full.names = T)

dec_sig_genes <- data.frame()
for (gene_file in 1: length(dec_gene_files)){
  
  sig_genes <- read_csv(dec_gene_files[gene_file])
  sig_genes <- subset(sig_genes, select = c(gene, direction))
  sig_genes$mag <- sig_genes$gene %>% substr(1,12)
  dec_sig_genes <- rbind(dec_sig_genes, sig_genes)
  
}

write.csv(dec_sig_genes, "decrease_sig_genes.csv", row.names = F)

dec_par_gene_files <- list.files(pattern = ".*decrease_par_pass.csv", full.names = T)

dec_par_sig_genes <- data.frame()
for (par_gene_file in 1: length(dec_par_gene_files)){
  
  par_sig_genes <- read_csv(dec_par_gene_files[par_gene_file])
  par_sig_genes <- subset(par_sig_genes, select = c(gene, direction))
  par_sig_genes$mag <- par_sig_genes$gene %>% substr(1,12)
  dec_par_sig_genes <- rbind(dec_par_sig_genes, par_sig_genes)
  
}

write.csv(dec_par_sig_genes, "parallel_decrease_sig_genes.csv", row.names = F)







# num_cols <- count(as.data.frame(colnames(gene_matrix))) -2
# gene_matrix_long <- pivot_longer(gene_matrix, cols = -c("gene", "gene_length"), names_to = "comparison", values_to = "SNVs_kbp")
# 
# 
# sweepplot <- ggplot(gene_matrix_long, aes(reorder_within(gene, SNVs_kbp, comparison), SNVs_kbp))+
#   geom_point(size = 1)+
#   scale_x_reordered()+
#   theme_classic()+
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
#   facet_wrap(~comparison, scales = "free_x", ncol = 1)
# save_plot(paste(MAG, "_sweep_plot.jpg", sep =""), nrow = num_cols[[1]], sweepplot, base_width = 8, base_height = 4, limitsize = F)
# 


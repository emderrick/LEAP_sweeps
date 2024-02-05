library(tidyverse)
#Fisher's exact test for each COG category (Sig genes in category, sig genes not in category, non-sig genes in category, non-sig genes not in category)

identify_enriched_categories <- function(genes,
                                         background,
                                         gene_to_category_map,
                                         min_category_count = min_category_count,
                                         to_ignore = character()) {
  
  enrichments_out <- data.frame(matrix(NA, nrow = length(gene_to_category_map), ncol = 8))
  rownames(enrichments_out) <- names(gene_to_category_map)
  colnames(enrichments_out) <- c("category", "genes_num_category", "genes_num_other",
                                 "background_num_category", "background_num_other", "OR", "p", "fdr")
  
  enrichments_out[names(gene_to_category_map), "category"] <- names(gene_to_category_map)
  
  for (category in rownames(enrichments_out)) {
    
    if (category %in% to_ignore) { next }
    
    genes_num_category <- length(which(genes %in% gene_to_category_map[[category]]))
    genes_num_other <- length(genes) - genes_num_category
    
    background_num_category <- length(which(background %in% gene_to_category_map[[category]]))
    background_num_other <- length(background) - background_num_category
    
    count_table <- matrix(c(genes_num_category, genes_num_other, background_num_category, background_num_other), nrow = 2, ncol = 2)
    
    if (min(c(genes_num_category + background_num_category, genes_num_other + background_num_other)) < min_category_count) {
      next
    }
    
    fisher_out <- fisher.test(count_table)
    
    enrichments_out[category, c("genes_num_category",
                                "genes_num_other",
                                "background_num_category",
                                "background_num_other", "p")] <- c(genes_num_category,
                                                                   genes_num_other,
                                                                   background_num_category,
                                                                   background_num_other,
                                                                   fisher_out$p.value)
    if (genes_num_other > 0) {
      ratio_numer <- genes_num_category / genes_num_other
    } else {
      ratio_numer <- genes_num_category / 1 
    }
    
    if (background_num_other == 0) {
      ratio_denom <- 1
    } else if(background_num_category == 0) {
      ratio_denom <- 1 / background_num_other
    } else {
      ratio_denom <- background_num_category / background_num_other
    }
    
    enrichments_out[category, "OR"] <- ratio_numer / ratio_denom
  }
  
  if (length(which(rowSums(is.na(enrichments_out)) > 1)) > 0) {
    enrichments_out <- enrichments_out[-which(rowSums(is.na(enrichments_out)) > 1), ]
  }
  
  enrichments_out$fdr <- p.adjust(enrichments_out$p, "fdr")
  
  rownames(enrichments_out) <- NULL
  
  return(enrichments_out)
  
}


COG_categories <- read.table('COG_category_descrip.tsv', header = FALSE, row.names = 1, stringsAsFactors = FALSE, sep = '\t')

# Read in mapping of which COGs are in which COG category (and convert this to list of COG category mappings to COG gene families).
COG_gene_to_category <- read.table("cog-20.to_category.tsv", header = FALSE, sep = "\t")
COG_category_to_COG <- list()
for (category in unique(COG_gene_to_category$V2)) {
  COG_category_to_COG[[category]] <- COG_gene_to_category[which(COG_gene_to_category$V2 == category), "V1"]
}

categories_to_ignore <- c('A', 'B', 'Y', 'Z')

all_background_genes <- read_csv("cog_background_genes.csv")

sig_gene_files <- c("SNV_decrease_genes.csv", "SNV_increase_genes.csv", "parallel_decrease_genes.csv", "parallel_increase_genes.csv", 
                    "threshold_significant_genes_all_subsamp.csv", "gene_cov_sig_increase_all_subsamp.csv", "gene_cov_sig_decrease_all_subsamp.csv")

for(gene_file in sig_gene_files){
  significant_genes_df <- read_csv(gene_file)
  significant_genes <- significant_genes_df[, c("COG_ID")] %>% na.omit()
  significant_genes <- significant_genes[['COG_ID']]
  background_genes_df <- all_background_genes %>% subset(!(gene %in% significant_genes_df$gene))
  background_genes <- background_genes_df[, c("COG_ID")] %>% na.omit()
  background_genes <- background_genes[['COG_ID']] 

  COG_enrichment_output <- identify_enriched_categories(genes = significant_genes,
                                                        background = background_genes,
                                                        gene_to_category_map = COG_category_to_COG,
                                                        min_category_count = 10, #categories without at least 10 COG gene families in the sig set and background are ignored
                                                        to_ignore = categories_to_ignore)

  write.csv(COG_enrichment_output, paste("COG_enrich_", gene_file, sep = ""), row.names = F)
}

EPSPS_class_1 <- list("I4_MAG_00006", "L7_MAG_00028", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")

class_1 <- c("sig_genes_loose_pos_C1_subsamp.csv", "sig_genes_loose_neg_C1_subsamp.csv", "sig_genes_strict_pos_C1_subsamp.csv", "sig_genes_strict_neg_C1_subsamp.csv",
             "sig_genes_threshold_C1_subsamp.csv", "sig_genes_increase_C1_subsamp.csv", "sig_genes_decrease_C1_subsamp.csv")

C1_background_genes <- subset(all_background_genes, mag %in% EPSPS_class_1)

for(gene_file in class_1){
  significant_genes_df <- read_csv(gene_file)
  significant_genes <- significant_genes_df[, c("COG_ID")] %>% na.omit()
  significant_genes <- significant_genes[['COG_ID']]
  background_genes_df <- C1_background_genes %>% subset(!(gene %in% significant_genes_df$gene))
  background_genes <- background_genes_df[, c("COG_ID")] %>% na.omit()
  background_genes <- background_genes[['COG_ID']] 
  
  COG_enrichment_output <- identify_enriched_categories(genes = significant_genes,
                                                        background = background_genes,
                                                        gene_to_category_map = COG_category_to_COG,
                                                        min_category_count = 10, #categories without at least 10 COG gene families in the sig set and background are ignored
                                                        to_ignore = categories_to_ignore)
  
  write.csv(COG_enrichment_output, paste("COG_enrich_", gene_file, sep = ""), row.names = F)
}

EPSPS_class_2 <- list("I4_MAG_00065", "L3_MAG_00058", "L7_MAG_00020", "L7_MAG_00043")

class_2 <- c("sig_genes_loose_pos_C2_subsamp.csv", "sig_genes_loose_neg_C2_subsamp.csv", "sig_genes_strict_pos_C2_subsamp.csv", "sig_genes_strict_neg_C2_subsamp.csv",
             "sig_genes_threshold_C2_subsamp.csv", "sig_genes_increase_C2_subsamp.csv", "sig_genes_decrease_C2_subsamp.csv")

C2_background_genes <- subset(all_background_genes, mag %in% EPSPS_class_2)

for(gene_file in class_2){
  significant_genes_df <- read_csv(gene_file)
  significant_genes <- significant_genes_df[, c("COG_ID")] %>% na.omit()
  significant_genes <- significant_genes[['COG_ID']]
  background_genes_df <- C2_background_genes %>% subset(!(gene %in% significant_genes_df$gene))
  background_genes <- background_genes_df[, c("COG_ID")] %>% na.omit()
  background_genes <- background_genes[['COG_ID']] 
  
  COG_enrichment_output <- identify_enriched_categories(genes = significant_genes,
                                                        background = background_genes,
                                                        gene_to_category_map = COG_category_to_COG,
                                                        min_category_count = 10, #categories without at least 10 COG gene families in the sig set and background are ignored
                                                        to_ignore = categories_to_ignore)
  
  write.csv(COG_enrichment_output, paste("COG_enrich_", gene_file, sep = ""), row.names = F)
}

import itertools
import os
import pickle
import pandas as pd
import csv

strict_pos_MAG_significant_genes = pd.DataFrame()
loose_pos_MAG_significant_genes = pd.DataFrame()

strict_neg_MAG_significant_genes = pd.DataFrame()
loose_neg_MAG_significant_genes = pd.DataFrame()

MAG_gene_files = os.listdir('parevol out subsamp')


for file in MAG_gene_files:
    strict_positive_genes = {}
    strict_negative_genes = {}
    loose_positive_genes = {}
    loose_negative_genes = {}
    MAG = file[0:12]
    input_file = 'parevol out subsamp/' + file
    with open(input_file, 'rb') as handle:
        treatment_combinations_dict = pickle.load(handle)

    combination_list = list(treatment_combinations_dict.keys())
    control_vs_control = []
    control_vs_GBH = []
    control_vs_GBH_delta_n = []
    GBH_vs_GBH = []

    for combination in combination_list:
        genes_list = treatment_combinations_dict[combination]['genes_significant'].tolist()
        delta_n_list = treatment_combinations_dict[combination]['delta_n'].tolist()
        gene_delta_n = list(zip(genes_list, delta_n_list))
    
        if 'Control' in combination[0] and 'Control' in combination[1]:
            control_vs_control.append(genes_list)
        
        if 'Control' in combination[0] and 'GBH' in combination[1]:
            control_vs_GBH.append(genes_list)
            control_vs_GBH_delta_n.append(gene_delta_n)
        
        if 'GBH' in combination[0] and 'GBH' in combination[1]:
            GBH_vs_GBH.append(genes_list)

    all_control_vs_control = list(itertools.chain.from_iterable(control_vs_control))
    unique_control_vs_control = set(all_control_vs_control)

    all_GBH_vs_GBH = list(itertools.chain.from_iterable(GBH_vs_GBH))
    unique_GBH_vs_GBH = set(all_GBH_vs_GBH)

    all_control_vs_GBH = list(itertools.chain.from_iterable(control_vs_GBH))
    unique_control_vs_GBH = set(all_control_vs_GBH)
    
    all_control_vs_GBH_delta_n = list(itertools.chain.from_iterable(control_vs_GBH_delta_n))

    common_genes = set(control_vs_GBH[0]).intersection(*control_vs_GBH[1:])

    strict_divergent_genes = []
    for gene in common_genes:
        if gene not in unique_control_vs_control and gene not in unique_GBH_vs_GBH:
            strict_divergent_genes.append(gene)

    strict_common_gene_delta_n = [[gene, delta_n_count] for gene, delta_n_count in all_control_vs_GBH_delta_n if gene in strict_divergent_genes]    

    strict_sig_genes = {}
    for (gene, count) in strict_common_gene_delta_n:
        if gene in strict_sig_genes:
            strict_sig_genes[gene].append(count)
        else:
            strict_sig_genes[gene] = [count]
    
    for gene, delta_n in strict_sig_genes.items():
        total = len(delta_n)
        positive = len([n for n in delta_n if n > 0])
        negative = len([n for n in delta_n if n < 0])
        if positive == total:
            strict_positive_genes[gene] = delta_n

        elif negative == total:
            strict_negative_genes[gene] = delta_n
        else:
            print(gene, delta_n)

    loose_divergent_genes = []
    for gene in unique_control_vs_GBH:
        if gene not in unique_control_vs_control and gene not in unique_GBH_vs_GBH:
            loose_divergent_genes.append(gene)

    loose_common_gene_delta_n = [[gene, delta_n_count] for gene, delta_n_count in all_control_vs_GBH_delta_n if gene in loose_divergent_genes]    
    
    loose_sig_genes = {}
    for (gene, count) in loose_common_gene_delta_n:
        if gene in loose_sig_genes:
            loose_sig_genes[gene].append(count)
        else:
            loose_sig_genes[gene] = [count]
   
    for gene, delta_n in loose_sig_genes.items():
        
        total = len(delta_n)
        positive = len([n for n in delta_n if n > 0])
        negative = len([n for n in delta_n if n < 0])
        if positive == total:
            loose_positive_genes[gene] = delta_n
        
        elif negative == total:
            loose_negative_genes[gene] = delta_n 
            
        else:
            print(gene, delta_n)     

    strict_pos_significant_genes = pd.DataFrame(strict_positive_genes.items())
    strict_pos_significant_genes['mag'] = MAG
    strict_pos_significant_genes['direction'] = "positive"
    strict_pos_MAG_significant_genes = pd.concat([strict_pos_MAG_significant_genes, strict_pos_significant_genes])
    
    strict_neg_significant_genes = pd.DataFrame(strict_negative_genes.items())
    strict_neg_significant_genes['mag'] = MAG
    strict_neg_significant_genes['direction'] = "negative"
    strict_neg_MAG_significant_genes = pd.concat([strict_neg_MAG_significant_genes, strict_neg_significant_genes])
    
    loose_pos_significant_genes = pd.DataFrame(loose_positive_genes.items())
    loose_pos_significant_genes['mag'] = MAG
    loose_pos_significant_genes['direction'] = "positive"
    loose_pos_MAG_significant_genes = pd.concat([loose_pos_MAG_significant_genes, loose_pos_significant_genes])
    
    loose_neg_significant_genes = pd.DataFrame(loose_negative_genes.items())
    loose_neg_significant_genes['mag'] = MAG
    loose_neg_significant_genes['direction'] = "negative"
    loose_neg_MAG_significant_genes = pd.concat([loose_neg_MAG_significant_genes, loose_neg_significant_genes])


all_strict = pd.concat([strict_pos_MAG_significant_genes, strict_neg_MAG_significant_genes])
all_loose = pd.concat([loose_pos_MAG_significant_genes, loose_neg_MAG_significant_genes])

strict_MAG_list = all_strict.values.tolist()
with open('strict_MAG_significant_genes_subsamp.csv', 'w') as f:
    write = csv.writer(f)
    write.writerow(['mag','direction', 'gene', 'delta_n'])
    write.writerows(strict_MAG_list)

loose_MAG_list = all_loose.values.tolist()
with open('loose_MAG_significant_genes_subsamp.csv', 'w') as f:
    write = csv.writer(f)
    write.writerow(['gene', 'delta_n','mag','direction'])
    write.writerows(loose_MAG_list)

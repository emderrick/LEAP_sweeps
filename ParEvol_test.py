#!/usr/bin/env python
# coding: utf-8

# In[1]:


from __future__ import division

import sys, pickle, math, random, os, itertools, re, collections
from itertools import combinations

import pandas as pd
import numpy as np

#import scipy.stats as stats
from scipy.special import iv

import csv


# In[ ]:


#for running MAGs in a loop
def skellam_pdf(delta_n, lambda_1, lambda_2):

    if delta_n > 0:

        pmf = ((lambda_1/lambda_2)**(delta_n/2)) * iv(delta_n, 2*np.sqrt(lambda_1*lambda_2))
        pmf += ((lambda_2/lambda_1)**(delta_n/2)) * iv(-1*delta_n, 2*np.sqrt(lambda_1*lambda_2))
        pmf *= np.exp((-1*lambda_1) + (-1*lambda_2))

    else:

        pmf = np.exp((-1*lambda_1) + (-1*lambda_2)) * iv(0, 2*np.sqrt(lambda_1*lambda_2))

    return pmf



def calculate_survival(counts_1, counts_2, genes, n_min = 3, alpha = 0.05):

    lambda_1 = sum(counts_1)/len(counts_1)
    lambda_2 = sum(counts_2)/len(counts_2)

    delta_n_original = np.absolute(counts_1-counts_2)
    delta_n = delta_n_original[delta_n_original>n_min]
    
    genes = genes[delta_n_original>n_min]

    delta_n_range = list(range(0,800))
    delta_n_range_array = np.asarray(delta_n_range)



    delta_n_no_absolute = counts_1-counts_2
    delta_n_no_absolute = delta_n_no_absolute[delta_n_original>n_min]


    delta_n_range_array_subset = delta_n_range_array[delta_n_range_array<=max(delta_n)]
    pmf = [skellam_pdf(i, lambda_1, lambda_2) for i in delta_n_range]
    pmf = np.asarray(pmf, dtype=np.float64)

    survival_null = [ 1-sum(pmf[:i]) for i in range(len(pmf)) ]
    survival_null = np.asarray(survival_null)
    survival_null = survival_null[delta_n_range_array<=max(delta_n)]

    survival_obs = [ len(delta_n[delta_n>=i])/len(delta_n) for i in delta_n_range]
    survival_obs = np.asarray(survival_obs)
    survival_obs = survival_obs[delta_n_range_array<=max(delta_n)]

    P_values = [sum(pmf[delta_n_range.index(delta_n_i):]) for delta_n_i in delta_n]
    P_values = np.asarray(P_values)
    P_range = np.linspace(10**-4, 0.05, num=10000)[::-1]


    N_bar_P_star_div_N_P_star = []
    P_stars = []


    for P_range_i in P_range:

        N_P_star = len(P_values[P_values<P_range_i])
        N_bar_P_star = 0

        if N_P_star == 0:
            continue


        for delta_n_j_idx, delta_n_j in enumerate(delta_n):
            if delta_n_j < n_min:
                continue

            P_delta_n_j = sum(pmf[delta_n_range.index(delta_n_j):])

            if P_range_i > P_delta_n_j:
                # no gene specific indices so just multiply the final probability by number of genes
                N_bar_P_star += skellam_pdf(delta_n_j, lambda_1, lambda_2) * len(delta_n_original)

        N_bar_P_star_div_N_P_star.append(N_bar_P_star/N_P_star)
        P_stars.append(P_range_i)


    N_bar_P_star_div_N_P_star = np.asarray(N_bar_P_star_div_N_P_star)
    P_stars = np.asarray(P_stars)
    position_P_star = np.argmax(N_bar_P_star_div_N_P_star<=0.05)
    P_star = P_stars[position_P_star]

    

    return delta_n_range_array_subset, genes, survival_obs, survival_null, delta_n_no_absolute, P_values, P_star


MAG_gene_files = os.listdir('gene matrices')

for file in MAG_gene_files[3:4]:
  
    MAG = file[0:12]
    print(MAG)
    input_file = 'gene matrices/' + file
    gene_count_matrix = pd.read_csv(input_file, sep = ',', index_col =0, header=None)
    treatments = gene_count_matrix.index.tolist()
    all_treatments = treatments[1:]
    gene_count_array = gene_count_matrix.to_numpy()
    my_genes = gene_count_array[0]
    treatment_combinations = list(itertools.combinations(all_treatments, 2))

    all_combinations_dict = {}
    
    for combination in treatment_combinations:
        print(combination)
        counts_1_index = gene_count_matrix.index.get_loc(combination[0])
        counts_1 = np.asarray([int(i) for i in gene_count_array[counts_1_index]])
        counts_2_index = gene_count_matrix.index.get_loc(combination[1])
        counts_2 = np.asarray([int(i) for i in gene_count_array[counts_2_index]])
        n_min = 2
        genes = np.asarray(my_genes)

        delta_n_range_array_subset, genes_keep, survival_obs, survival_null, delta_n, P_values, P_star = calculate_survival(counts_1, counts_2, genes, n_min=n_min)

        # keep significant genes
        P_values_significant = P_values[P_values<P_star]
        genes_significant = genes_keep[P_values<P_star]
        delta_n = delta_n[P_values<P_star]

        gene_dict = {}
        gene_dict['P_values_significant'] = P_values_significant
        gene_dict['genes_significant'] = genes_significant
        gene_dict['delta_n'] = delta_n
        all_combinations_dict[combination] = gene_dict
    
    output_file = MAG + '_300_all_combinations_sig_genes.pickle'
    
    with open(output_file, 'wb') as handle:
        pickle.dump(all_combinations_dict, handle, protocol = pickle.HIGHEST_PROTOCOL)
    
    print(all_combinations_dict)
    print("done " + MAG)


# In[2]:


MAG_gene_files = os.listdir('parevol_300')

MAG_significant_genes = pd.DataFrame()
for file in MAG_gene_files:
    input_file = 'parevol_300/' + file
    with open(input_file, 'rb') as handle:
        treatment_combinations_dict = pickle.load(handle)
    
    MAG = file[0:12]

    combination_list = list(treatment_combinations_dict.keys())
    control_vs_control = []
    control_vs_GBH = []
    GBH_vs_GBH = []
    for combination in combination_list:
        genes_list = treatment_combinations_dict[combination]['genes_significant'].tolist()
        if 'Control' in combination[0] and 'Control' in combination[1]:
            control_vs_control.append(genes_list)
    
        if 'Control' in combination[0] and 'GBH' in combination[1]:
            control_vs_GBH.append(genes_list)
   
        if 'GBH' in combination[0] and 'GBH' in combination[1]:
            GBH_vs_GBH.append(genes_list)
    
    all_control_vs_control = list(itertools.chain.from_iterable(control_vs_control))
    unique_control_vs_control = set(all_control_vs_control)

    all_GBH_vs_GBH = list(itertools.chain.from_iterable(GBH_vs_GBH))
    unique_GBH_vs_GBH = set(all_GBH_vs_GBH)

    #common_genes = list(set(control_vs_GBH[0]).intersection(*control_vs_GBH[1:]))

    all_control_vs_GBH = list(itertools.chain.from_iterable(control_vs_GBH))
    unique_control_vs_GBH = set(all_control_vs_GBH)
    
    divergent_genes = []
    #for gene in common_genes:
        #if gene not in unique_control_vs_control and gene not in unique_GBH_vs_GBH:
            #divergent_genes.append(gene)
    
    for gene in unique_control_vs_GBH:
        if gene not in unique_control_vs_control and gene not in unique_GBH_vs_GBH:
            divergent_genes.append(gene)
    
    divergent_genes = pd.DataFrame(divergent_genes)
    divergent_genes['mag'] = MAG
    MAG_significant_genes = pd.concat([MAG_significant_genes, divergent_genes])


MAG_list = MAG_significant_genes.values.tolist()

with open('not_strict_MAG_significant_genes.csv', 'w') as f:
    write = csv.writer(f)
    write.writerow(['gene','mag'])
    write.writerows(MAG_list)


# In[ ]:


#to check if pmf of 300 gives same results as 1000
file = open("1000_MAG_significant_genes.csv", "r")
genes_1000 = list(csv.reader(file, delimiter=","))
file.close()
file2 = open("MAG_significant_genes.csv", "r")
genes_300 = list(csv.reader(file2, delimiter=","))
file2.close()

def common_member(a, b):
    result = [i for i in a if i in b]
    return result

common_100_300 = common_member(genes_1000, genes_300)
print(len(genes_1000))
print(len(genes_300))
print(len(common_100_300))


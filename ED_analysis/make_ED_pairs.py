# Author: Xiang Zhang

# 12052023: for MC, get a list of all DepMap pairs with a significant +ve and -ve ED score 
# without restricting to gene pairs tested in GIN

import os, sys, re
from math import *
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
from scipy.stats import pearsonr

# Design a function to rename the Genes from the Broad data
def renameGene(x):
    x = x.split()[0]
    return x

# Design a function to rename the Genes from the columns of the Lab data
def renameLabColumn(x):
    x = x.split("_")[0]
    return x

import itertools
def get_all_paris(gene_list, L=2):
    gene1_list = []
    gene2_list = []
    for subset in itertools.combinations(gene_list, L):
        gene1_list.append(subset[0])
        gene2_list.append(subset[1])
    df = pd.DataFrame(data={'gene1': gene1_list, 'gene2': gene2_list})
    return df

def get_fitness_vs_expression(df_fitness, df_expression, depmap_ver='22Q4', output_direc='./'):
    overlap_CLs = np.intersect1d(df_fitness.columns, df_expression.columns)
    overlap_genes = np.intersect1d(df_fitness.index, df_expression.index)
    
    df_ab = get_all_paris(overlap_genes, L=2)
    df_ab.to_csv(output_direc+'all_depmap_'+depmap_ver+'_pairs.tsv', sep='\t', index=None)
    
    pcc_ab_list = []
    pcc_ab_pval_list = []
    pcc_ba_list = []
    pcc_ba_pval_list = []
    for row in df_ab.iterrows():
        a_gene = row[1]['gene1']
        b_gene = row[1]['gene2']
        if (a_gene not in df_fitness.index) or (b_gene not in df_expression.index):
            pcc_ab = np.nan
            pcc_ab_pval = np.nan
        else:
            overlap_CLs_noNA_cl = df_fitness.loc[a_gene, overlap_CLs].dropna().index
            overlap_CLs_noNA_ex = df_expression.loc[b_gene, overlap_CLs].dropna().index
            overlap_CLs_ab = np.intersect1d(overlap_CLs_noNA_cl, overlap_CLs_noNA_ex)
            df_fitness_a = df_fitness.loc[a_gene, overlap_CLs_ab].dropna()
            df_expression_b = df_expression.loc[b_gene, overlap_CLs_ab].dropna()
            pcc_ab_stat = pearsonr(df_fitness_a, df_expression_b)
            pcc_ab = pcc_ab_stat[0]
            pcc_ab_pval = pcc_ab_stat[1]
            
        if (b_gene not in df_fitness.index) or (a_gene not in df_expression.index):
            pcc_ba = np.nan
            pcc_ba_pval = np.nan
        else:
            overlap_CLs_noNA_cl = df_fitness.loc[b_gene, overlap_CLs].dropna().index
            overlap_CLs_noNA_ex = df_expression.loc[a_gene, overlap_CLs].dropna().index
            overlap_CLs_ba = np.intersect1d(overlap_CLs_noNA_cl, overlap_CLs_noNA_ex)
            df_fitness_b = df_fitness.loc[b_gene, overlap_CLs_ba].dropna()
            df_expression_a = df_expression.loc[a_gene, overlap_CLs_ba].dropna()
            pcc_ba_stat = pearsonr(df_fitness_b, df_expression_a)
            pcc_ba = pcc_ba_stat[0]
            pcc_ba_pval = pcc_ba_stat[1]
            
        pcc_ab_list.append(pcc_ab)
        pcc_ab_pval_list.append(pcc_ab_pval)
        pcc_ba_list.append(pcc_ba)
        pcc_ba_pval_list.append(pcc_ba_pval)
    df_ab['PCC_lib_fitness_vs_query_expression'] = pcc_ab_list
    df_ab['AB_pcc_pval'] = pcc_ab_pval_list
    df_ab['PCC_query_fitness_vs_lib_expression'] = pcc_ba_list
    df_ab['BA_pcc_pval'] = pcc_ba_pval_list
    #df_ab['average_PCC'] = (df_ab['PCC_lib_fitness_vs_query_expression'] + df_ab['PCC_query_fitness_vs_lib_expression']) / 2
    df_ab.to_csv(output_direc+'all_depmap_'+depmap_ver+'_ED_scores.tsv', sep='\t', index=None)
    return df_ab



if __name__ == '__main__':
    df_broad_22q4 = pd.read_csv('./22Q4/CRISPRGeneEffect.csv', index_col=0)
    df_broad_22q4_T = df_broad_22q4.T
    broad_gene_22q4 = df_broad_22q4_T.reset_index()["index"].apply(renameGene)
    df_broad_22q4_T.index = broad_gene_22q4

    # Read in CCLE expression file
    df_expression_22q4 = pd.read_csv('./22Q4/OmicsExpressionProteinCodingGenesTPMLogp1.csv', index_col=0)
    df_expression_22q4_T = df_expression_22q4.T
    express_gene_22q4 = df_expression_22q4_T.reset_index()["index"].apply(renameGene)
    df_expression_22q4_T.index = express_gene_22q4

    df_depmap22q4_pairs = get_fitness_vs_expression(df_fitness=df_broad_22q4_T,
                                                df_expression=df_expression_22q4_T,
                                                depmap_ver='22Q4', output_direc='./')

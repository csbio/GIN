# Author: Xiang Zhang
# Date: 11/28/2023
# Description: Co-Annotation enrichment on the levels of clusters

import sys
import os
from os.path import expanduser
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time

import pyreadr
import itertools
from scipy.stats import hypergeom
from statsmodels.stats.multitest import fdrcorrection

# hypergeometric test:
# We can then compute a probability of drawing X red marbles out of N from 
# a jar containing n red marbles out of M in the following way:
def hyper_test(x, M, n, N):
    pval = hypergeom.sf(x-1, M, n, N)
    #print(pval)
    return pval

def get_BH_correct_pval(df, pval_col):
    (reject_list, pval_bh) = fdrcorrection(df[pval_col], alpha=0.05, method='indep', is_sorted=False)
    return pval_bh

def read_coAnno_standard_from_Rdata(input_direc):
    result = pyreadr.read_r(input_direc)
    print("Read in co-annotation standard:", result.keys()) # let's check what objects we got
    df_standard_co = result[list(result.keys())[0]]
    df_standard_co = df_standard_co.set_index('gene1')
    return df_standard_co

def get_all_paris(gene_list, L=2):
    gene1_list = []
    gene2_list = []
    for subset in itertools.combinations(gene_list, L):
        gene1_list.append(subset[0])
        gene2_list.append(subset[1])
    df = pd.DataFrame(data={'gene1': gene1_list, 'gene2': gene2_list})
    return df

'''
# old - abandoned
def get_annotated_pairs(df_pairs, df_standard, output_file=None):
    df_pairs['is_annotated'] = 0
    gene1_idx_Anno1 = np.intersect1d(df_pairs['gene1'], df_standard.index)
    gene2_idx_Anno1 = np.intersect1d(df_pairs['gene2'], df_standard.index)

    # 1. Deal with those whose 'gene1' is in the df_standard.index and 'gene2' is in df_standard['gene2]
    # Get the final gene2 list 
    df_standard_gene1_Anno = df_standard.loc[gene1_idx_Anno1]
    gene1_gene2_Anno = df_standard_gene1_Anno['gene2'].values
    for g1 in gene1_idx_Anno1:
        df_pairs_g1 = df_pairs[df_pairs['gene1']==g1]
        gene1_gene2_Anno = df_standard.loc[[g1], 'gene2'].values.tolist()
        # Get the final gene1 list
        gene1_final_idx = df_pairs_g1[df_pairs_g1['gene2'].isin(gene1_gene2_Anno)].index
        df_pairs.loc[gene1_final_idx, 'is_annotated'] = 1

    # 2. Deal with those whose 'gene2' is in the df_standard.index and 'gene1' is in df_standard['gene2']
    # Get the final gene1 list 
    df_standard_gene2_Anno = df_standard.loc[gene2_idx_Anno1]
    gene2_gene1_Anno = df_standard_gene2_Anno['gene2'].values
    for g2 in gene2_idx_Anno1:
        df_pairs_g2 = df_pairs[df_pairs['gene2']==g2]
        gene2_gene1_Anno = df_standard.loc[[g2], 'gene2'].values.tolist()
        # Get the final gene1 list
        gene2_final_idx = df_pairs_g2[df_pairs_g2['gene1'].isin(gene2_gene1_Anno)].index
        df_pairs.loc[gene2_final_idx, 'is_annotated'] = 1

    df_anno_pairs = df_pairs[df_pairs['is_annotated']==1]

    if output_file != None:
        df_anno_pairs.to_csv(output_file)
        #df_pairs.to_csv(output_file[:-4]+'_all.csv')

    return df_anno_pairs
'''

def get_annotated_pairs(df_pairs, df_standard, gene1_col='gene1', gene2_col='gene2', output_file=None):
    # Idea:
    # 1. Get the background: overlapping pairs with the standard and the input pair list - double index overlap
    # 2. Annotate the pairs: double index overlap between "1"s
    # 3. Annotation in background and in subset
    df_pairs_2idx_ab = df_pairs.reset_index().set_index([gene1_col, gene2_col])
    df_pairs_2idx_ba = df_pairs.reset_index().set_index([gene2_col, gene1_col])
    df_standard_2idx = df_standard.reset_index().set_index(['gene1', 'gene2'])

    ab_overlap = np.intersect1d(df_standard_2idx.index, df_pairs_2idx_ab.index)
    ba_overlap = np.intersect1d(df_standard_2idx.index, df_pairs_2idx_ba.index)
    abba_overlap = ab_overlap.tolist() + ba_overlap.tolist()
    #print(len(ab_overlap))
    #print(len(ba_overlap))
    #print(len(abba_overlap))

    num_bg_pairs = len(abba_overlap)
    #print("Overlap background pairs", num_bg_pairs)

    df_standard_bg = df_standard_2idx.loc[abba_overlap]
    df_standard_bg_anno = df_standard_bg[df_standard_bg['is_annotated']==1]
    num_bg_anno = df_standard_bg_anno.shape[0]
    #print("Overlap background annotations", num_bg_anno)


    df_ab_anno = df_standard_2idx.loc[ab_overlap][df_standard_2idx.loc[ab_overlap, 'is_annotated']==1]
    df_ba_anno = df_standard_2idx.loc[ba_overlap][df_standard_2idx.loc[ba_overlap, 'is_annotated']==1]
    df_ba_anno_reverse_idx = df_ba_anno.reset_index().set_index(['gene2', 'gene1'])

    # Get the overlap bg df for df_pairs (should be the same length as num_bg_pairs)
    df_pairs_bg = pd.concat([df_pairs_2idx_ab.loc[ab_overlap],
                             df_pairs_2idx_ba.loc[ba_overlap].reset_index().set_index([gene1_col, gene2_col])])
    #print(df_pairs_bg.shape)

    df_pairs_bg['is_annotated'] = 0
    df_pairs_bg.loc[df_ab_anno.index, 'is_annotated'] = 1
    df_pairs_bg.loc[df_ba_anno_reverse_idx.index, 'is_annotated'] = 1
    #print(df_pairs_bg['is_annotated'].sum()) # should equal to num_bg_anno

    if output_file != None:
        df_pairs_bg.to_csv(output_file)
        #df_pairs.to_csv(output_file[:-4]+'_all.csv')

    df_anno_pairs = df_pairs_bg[df_pairs_bg['is_annotated']==1]

    return df_anno_pairs


def coAnno_enrichment(df_cluster, df_standard_orig, standard_name, output_direc, shell=True, correction_method=None):
    if not os.path.exists(output_direc):
        os.mkdir(output_direc)
    # Population: all possible pairs from the gene set
    df_allpairs = get_all_paris(df_cluster.index.tolist())
    M = df_allpairs.shape[0]
    print("Population (all possible pairs):", M)

    df_standard = df_standard_orig[df_standard_orig["is_annotated"]==1]

    # Successes of population: all annotated pairs from the gene set
    all_annotated_pairs_direc = output_direc+standard_name+"_annotated_pairs.csv"
    if os.path.exists(all_annotated_pairs_direc):
        print("Loading existing file for all annotated pairs: " + all_annotated_pairs_direc)
        df_all_annotated_pairs = pd.read_csv(all_annotated_pairs_direc) 
    else:
        print("Generating new file for all annotated pairs")
        # Timing
        tic = time.perf_counter()
        df_all_annotated_pairs = get_annotated_pairs(df_allpairs, df_standard, gene1_col='gene1', gene2_col='gene2', output_file=all_annotated_pairs_direc)
        toc = time.perf_counter()
        print(f"Get annotated pair list in {toc - tic:0.4f} seconds")
    n = df_all_annotated_pairs.shape[0]

    df_result = pd.DataFrame(np.zeros([len(df_cluster.columns), 5]), index=df_cluster.columns, columns=[standard_name, "#all_pairs", "#all_annotated_pairs", "#cluster_pairs", "#cluster_annotated_pairs"])
    for col_name in df_cluster.columns:
        cluster_ids = df_cluster[col_name].unique()
        N = 0 # sample_size
        x = 0 # sample_success
        for c_id in cluster_ids:
            df_cid = df_cluster[df_cluster[col_name]==c_id]
            df_allpairs_cid = get_all_paris(df_cid.index.tolist())
            N += df_allpairs_cid.shape[0]
            df_cid_annotated_pairs = get_annotated_pairs(df_allpairs_cid, df_standard)
            x += df_cid_annotated_pairs.shape[0]

        pval = hyper_test(x, M, n, N)
        print(x, M, n, N)
        df_result.loc[col_name, standard_name] = pval
        df_result.loc[col_name, "#all_pairs"] = M
        df_result.loc[col_name, "#all_annotated_pairs"] = n
        df_result.loc[col_name, "#cluster_pairs"] = N
        df_result.loc[col_name, "#cluster_annotated_pairs"] = x

    if correction_method == 'BH':
        df_result[standard_name] = get_BH_correct_pval(df_result, standard_name)

    if shell:
        print("Add shell version")
        df_result["shell_pval"] = np.nan
        df_result["shell_fold"] = np.nan
        for i in range(1, df_result.shape[0]):
            idx = df_result.index[i]
            idx_prev = df_result.index[i-1]
            N = df_result.loc[idx, "#cluster_pairs"] - df_result.loc[idx_prev, "#cluster_pairs"] 
            x = df_result.loc[idx, "#cluster_annotated_pairs"] - df_result.loc[idx_prev, "#cluster_annotated_pairs"] 
            shell_pval = hyper_test(x, M, n, N)
            shell_fold = round(float(x/N / (n/M)), 2)
            df_result.loc[idx, "shell_pval"] = shell_pval
            df_result.loc[idx, "shell_fold"] = shell_fold
        # Distant layer
        N = M - df_result.loc[idx, "#cluster_pairs"]
        x = n - df_result.loc[idx, "#cluster_annotated_pairs"]
        shell_pval = hyper_test(x, M, n, N)
        shell_fold = round(float(x/N / (n/M)), 2)
        df_result.loc['Distant'] = [None, M, n, N, x, shell_pval, shell_fold]

    return df_result

if __name__ == '__main__':
    #df_c1571 = pd.read_csv('./new_exp_4k_clusters_noL2C18_noMitoL4C3_1571genes_c12rm_no507Var2_filtered_addL0_02062023.txt', sep='\t', index_col=0)
    df_c1863 = pd.read_csv('./new_exp_4k_clusters_noL2C18_1863genes_c12rm_no507Var2_filtered_addL0_02062023.txt', sep='\t', index_col=0)

    # Cluster (1571 genes - no mito)
    output_direc1571 = './UPDATED_coAnno_enrichment_1571genes_filtered_noL2C18_noMitoL4C3_withL0/'
    # Cluster (1863 genes - with mito) 
    output_direc1863 = './UPDATED_coAnno_enrichment_1863genes_filtered_noL2C18_withL0/'

    #df_c = df_c1571
    #output_direc = output_direc1571
    df_c = df_c1863
    output_direc = output_direc1863
    if not os.path.exists(output_direc):
        os.makedirs(output_direc)

    # CORUM enrichment
    corum_standard_direc = './corum.Rdata' 
    df_corum_standard_co = read_coAnno_standard_from_Rdata(corum_standard_direc)
    df_CORUM_coAnno_enrich = coAnno_enrichment(df_c, df_corum_standard_co, standard_name='CORUM', output_direc=output_direc, correction_method=None)
    df_CORUM_coAnno_enrich.to_csv(output_direc+'coAnno_enrich_CORUM.csv')

    # GO BP enrichment
    gobp_standard_direc = './GO_BP.Rdata' 
    df_gobp_standard_co = read_coAnno_standard_from_Rdata(gobp_standard_direc)
    df_gobp_coAnno_enrich = coAnno_enrichment(df_c, df_gobp_standard_co, standard_name='GO_BP', output_direc=output_direc, correction_method=None)
    df_gobp_coAnno_enrich.to_csv(output_direc+'coAnno_enrich_GO_BP.csv')

    # Pathway enrichment
    path_standard_direc = './Pathway.Rdata' 
    df_path_standard_co = read_coAnno_standard_from_Rdata(path_standard_direc)
    df_path_coAnno_enrich = coAnno_enrichment(df_c, df_path_standard_co, standard_name='Pathway', output_direc=output_direc, correction_method=None)
    df_path_coAnno_enrich.to_csv(output_direc+'coAnno_enrich_Pathway.csv')

    # localization enrichment
    path_standard_direc = './localization.Rdata' 
    df_path_standard_co = read_coAnno_standard_from_Rdata(path_standard_direc)
    df_path_coAnno_enrich = coAnno_enrichment(df_c, df_path_standard_co, standard_name='Localization', output_direc=output_direc, correction_method=None)
    df_path_coAnno_enrich.to_csv(output_direc+'coAnno_enrich_Localization.csv')

# Author: Xiang Zhang (zhan6668)
# Description: This file starts at 06152024, in order to run the hypergeometric method pipeline on DepMap 22Q4 data
# and qGI_20211111_nRed to check about gene pair set enrichment (whether GIs are enriched in DepMap PCC pairs),
# with 9 unique queries excluded due to their frequency in the intersection > 2%. p-value cutoff is 0.01.

import os, sys, re
from math import *
import numpy as np
#import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats

# hypergeometric test:
# We can then compute a probability of drawing X red marbles out of N from
# a jar containing n red marbles out of M in the following way:
from scipy.stats import hypergeom
def hyper_test(x, M, n, N):
    pval = hypergeom.sf(x-1, M, n, N)
    #print(pval)
    return pval

def get_BH_correct_pval(df, pval_col):
    (reject_list, pval_bh) = fdrcorrection(df[pval_col], alpha=0.05, method='indep', is_sorted=False)
    return pval_bh

# Update on 05152022: Add frequency filter as 2% (qGI qualified pairs / all DepMap significant pairs)
def get_qGI_fc_subset(df_pairs, df_fc, interaction='positive', frequency_cutoff=1.0, topn=5, skipStep1=True, skipStep2=True):
    df_pairs_GI = df_pairs[df_pairs['GI_standard']==interaction]
    df_pairs_GI = df_pairs_GI.set_index('library_gene')
    if interaction == 'positive':
        if not skipStep1:
            df_fc_neg = df_fc[df_fc['rich']<0]
            df_fc_neg_abs = np.abs(df_fc_neg[['rich']])
        else:
            df_fc_neg_abs = df_fc[['rich']]
        df_pairs_GI_fc = df_pairs_GI.loc[np.intersect1d(df_pairs_GI.index, df_fc_neg_abs.index)]
        df_pairs_GI_fc = df_pairs_GI_fc.merge(df_fc_neg_abs, how='left', left_index=True, right_index=True)
        if not skipStep2:
            df_pairs_GI_fc = df_pairs_GI_fc[df_pairs_GI_fc['qGI_score']>df_pairs_GI_fc['rich']]
    elif interaction == 'negative':
        # Update: Don't do any FC filter with negative side
        df_pairs_GI_fc  = df_pairs_GI # [df_pairs_GI['qGI_score']<df_pairs_GI['rich']]
    else:
        print('error')
    df_pairs_GI_fc.index.name = 'library_gene'
    df_pairs_GI_fc = df_pairs_GI_fc.reset_index()
    print("# pairs before frequency filter", len(df_pairs_GI_fc))
    # frequency filter of 2% after the above two filters
    df_pairs_GI_fc_freq = df_pairs_GI_fc.groupby('query_gene').count()['library_gene']/df_pairs_GI_fc.shape[0]
    print("Frequency Top " +str(topn)+":\n",df_pairs_GI_fc_freq.nlargest(topn))
    querys_to_keep = df_pairs_GI_fc_freq[df_pairs_GI_fc_freq<frequency_cutoff].index
    print("Remove pairs with query frequency >= "+str(frequency_cutoff), list(set(df_pairs_GI_fc_freq.index) - set(querys_to_keep)))
    df_pairs_GI_fc = df_pairs_GI_fc[df_pairs_GI_fc['query_gene'].isin(querys_to_keep)]
    print("# pairs after frequency filter", len(df_pairs_GI_fc))
    print("\n\n")
    return df_pairs_GI_fc

def set_enrichment_GI_in_DepMap(df_depmap_bg, df_depmap_subset, df_fc, frequency_filter=1.0, gi_direc='negative'):
    M = df_depmap_bg.shape[0]
    # No frequency filter at this point
    df_demap_bg_qGI = get_qGI_fc_subset(df_depmap_bg, df_fc, interaction=gi_direc, frequency_cutoff=1)
    n = df_demap_bg_qGI.shape[0]
    
    N = df_depmap_subset.shape[0]
    df_depmap_subset_qGI = get_qGI_fc_subset(df_depmap_subset, df_fc, interaction=gi_direc, frequency_cutoff=frequency_filter)
    x = df_depmap_subset_qGI.shape[0]
    
    print(M, n, N, x)
    fc = (x/N) / (n/M)
    print("fold change:", fc)
    pval = hyper_test(x, M, n, N)
    print("p-value:", pval)
    return pval

def exlude_paralogs(df_pairs, df_standard, gene1_col='library_gene', gene2_col='query_gene', output_file=None):
    # Idea:
    # 1. Get the background: overlapping pairs with the standard and the input pair list - double index overlap
    # 2. Annotate the pairs: double index overlap between "1"s
    # 3. Annotation in background and in subset
    df_pairs_2idx = df_pairs.reset_index().set_index([gene1_col, gene2_col])
    df_standard_2idx_ab = df_standard.reset_index().set_index(['gene1', 'gene2'])
    df_standard_2idx_ba = df_standard.reset_index().set_index(['gene2', 'gene1'])
    
    ab_overlap = np.intersect1d(df_standard_2idx_ab.index, df_pairs_2idx.index)
    ba_overlap = np.intersect1d(df_standard_2idx_ba.index, df_pairs_2idx.index)
    abba_overlap = ab_overlap.tolist() + ba_overlap.tolist()
    print(len(ab_overlap))
    print(len(ba_overlap))
    print(len(abba_overlap))
    print(set(abba_overlap))
    df_pairs_no_paralogs = df_pairs_2idx.drop(abba_overlap, axis=0).reset_index()
    return df_pairs_no_paralogs

# Build softer merge function so that
# Either AB or BA should be significant (pass PCC and pval threshold)
# They just could not be both significant but in different directions
def softer_filter_abba(df, pcc_ab_col='PCC_lib_fitness_vs_query_expression',
            pcc_ba_col='PCC_query_fitness_vs_lib_expression', output_direc=None,
                    pcc_cutoff=0.1, pcc_pval=0.01):
    df_abba = df.copy(deep=True)
    print("# pairs before filtering:", df_abba.shape[0])
    
    # Filter significant different directions of AB/BA's PCC
    df_abba['pcc_product'] = df_abba[pcc_ab_col] * df_abba[pcc_ba_col]
    df_abba_oppo = df_abba[df_abba['pcc_product']<0]
    rm_idx = df_abba_oppo[(df_abba_oppo['AB_pcc_pval']<pcc_pval) & (df_abba_oppo['BA_pcc_pval']<pcc_pval)].index
    print("# pairs to remove because of significance in opposite PCC directions:", len(rm_idx))
    df_abba = df_abba.drop(rm_idx)
    
    # Filter based on AB or BA significance
    if pcc_cutoff is not None:
        if pcc_cutoff > 0:
            df_ab_pass = df_abba[(df_abba[pcc_ab_col]>pcc_cutoff) & (df_abba['AB_pcc_pval']<pcc_pval)]
            df_ba_pass = df_abba[(df_abba[pcc_ba_col]>pcc_cutoff) & (df_abba['BA_pcc_pval']<pcc_pval)]
            print("# AB pairs after PCC&pval filter:", df_ab_pass.shape[0])
            print("# BA pairs after PCC&pval filter:", df_ba_pass.shape[0])
            df_abba = pd.concat([df_ab_pass, df_ba_pass]).drop_duplicates()
            print("# pairs after duplicates removed:", df_abba.shape[0])
        else:
            df_ab_pass = df_abba[(df_abba[pcc_ab_col]<pcc_cutoff) & (df_abba['AB_pcc_pval']<pcc_pval)]
            df_ba_pass = df_abba[(df_abba[pcc_ba_col]<pcc_cutoff) & (df_abba['BA_pcc_pval']<pcc_pval)]
            print("# AB pairs after PCC&pval filter:", df_ab_pass.shape[0])
            print("# BA pairs after PCC&pval filter:", df_ba_pass.shape[0])
            df_abba = pd.concat([df_ab_pass, df_ba_pass]).drop_duplicates()
            print("# pairs after duplicates removed:", df_abba.shape[0])
    
    if output_direc != None:
        df_abba.to_csv(output_direc, index=None)
    
    return df_abba

if __name__ == "__main__":
    '''
    df_qGI_nRed_pairs = pd.read_csv('qGI_pairwise_20211111_nRed_cleaned.txt', index_col=0, sep='\t')
    df_qGI_nRed_pairs = df_qGI_nRed_pairs[['Query', 'qGI_score', 'FDR', 'GI_standard']]
    df_qGI_nRed_pairs.columns = ['query_gene', 'qGI_score', 'qGI_FDR', 'GI_standard']
    df_qGI_nRed_pairs.index.name='library_gene'
    '''

    # 08312023: by default we skip Step 1 & Step 2 regarding the FC filter (only do standard GI filter for both pos/neg GIs)
    df_fc = pd.read_csv('../../fc_singlePhenotype_20211111.txt', index_col=0, sep='\t')

    df_qGI_nRed_depmap = pd.read_csv('../DepMap_22Q4_all_qGI_nRed_pairs_noFreq2_9queries.tsv', sep='\t', index_col=None)

    # Merge AB/BA and do the set enrichment 
    df_depmap_abba_softer_merge_posPCC = softer_filter_abba(df_qGI_nRed_depmap, pcc_ab_col='PCC_lib_fitness_vs_query_expression',
            pcc_ba_col='PCC_query_fitness_vs_lib_expression', output_direc='/Users/zhangxiang/Documents/Research_CB/GIN/Code_BROAD_LAB_DATA/suppressor_analysis_04142023/DepMap_22Q4_qGI_enrichment_ABBAsofter/DepMap_positive_ABBAsofter_PCC01_pval001_noFreq2.csv',
                    pcc_cutoff=0.1, pcc_pval=0.01)
    df_depmap_abba_softer_merge_negPCC = softer_filter_abba(df_qGI_nRed_depmap, pcc_ab_col='PCC_lib_fitness_vs_query_expression',
            pcc_ba_col='PCC_query_fitness_vs_lib_expression', output_direc='/Users/zhangxiang/Documents/Research_CB/GIN/Code_BROAD_LAB_DATA/suppressor_analysis_04142023/DepMap_22Q4_qGI_enrichment_ABBAsofter/DepMap_negative_ABBAsofter_PCCn01_pval001_noFreq2.csv',
                    pcc_cutoff=-0.1, pcc_pval=0.01)
    df_depmap_abba_softer_merge_posPCC = pd.read_csv('/Users/zhangxiang/Documents/Research_CB/GIN/Code_BROAD_LAB_DATA/suppressor_analysis_04142023/DepMap_22Q4_qGI_enrichment_ABBAsofter/DepMap_positive_ABBAsofter_PCC01_pval001_noFreq2.csv')
    df_depmap_abba_softer_merge_negPCC = pd.read_csv('/Users/zhangxiang/Documents/Research_CB/GIN/Code_BROAD_LAB_DATA/suppressor_analysis_04142023/DepMap_22Q4_qGI_enrichment_ABBAsofter/DepMap_negative_ABBAsofter_PCCn01_pval001_noFreq2.csv')

    print('Positive ABBA softer (PCC>0.1), pval<0.01, freq 100%): # Pairs:', len(df_depmap_abba_softer_merge_posPCC))
    print('Negative ABBA softer (PCC<-0.1), pval<0.01, freq 100%): # Pairs:', len(df_depmap_abba_softer_merge_negPCC))

    # Optional: get the 4 combos subset
    # Frequency filter comes after the GI filter
    df_posPCC_posGI = get_qGI_fc_subset(df_pairs=df_depmap_abba_softer_merge_posPCC, df_fc=df_fc,
                      interaction='positive', frequency_cutoff=1.0, topn=5)
    df_posPCC_negGI = get_qGI_fc_subset(df_pairs=df_depmap_abba_softer_merge_posPCC, df_fc=df_fc,
                      interaction='negative', frequency_cutoff=1.0, topn=5)
    df_negPCC_posGI = get_qGI_fc_subset(df_pairs=df_depmap_abba_softer_merge_negPCC, df_fc=df_fc,
                      interaction='positive', frequency_cutoff=1.0, topn=5)
    df_negPCC_negGI = get_qGI_fc_subset(df_pairs=df_depmap_abba_softer_merge_negPCC, df_fc=df_fc,
                      interaction='negative', frequency_cutoff=1.0, topn=5)

    print("\n\nSet enrichment:\n\n")
    
    # Only for saving the skip Step12 version final lists
    df_posPCC_posGI.to_csv('./filtered_intersection_pairs_positivePCC_positiveGI_skipStep12_PCC01_pval001_noFreq2.csv', index=None)
    df_posPCC_negGI.to_csv('./filtered_intersection_pairs_positivePCC_negativeGI_skipStep12_PCC01_pval001_noFreq2.csv', index=None)
    df_negPCC_posGI.to_csv('./filtered_intersection_pairs_negativePCC_positiveGI_skipStep12_PCC01_pval001_noFreq2.csv', index=None)
    df_negPCC_negGI.to_csv('./filtered_intersection_pairs_negativePCC_negativeGI_skipStep12_PCC01_pval001_noFreq2.csv', index=None)
    

    # Set enrichment
    pval_posPCC_posGI = set_enrichment_GI_in_DepMap(df_qGI_nRed_depmap, df_depmap_abba_softer_merge_posPCC, df_fc, gi_direc='positive')
    pval_posPCC_negGI = set_enrichment_GI_in_DepMap(df_qGI_nRed_depmap, df_depmap_abba_softer_merge_posPCC, df_fc, gi_direc='negative')
    pval_negPCC_posGI = set_enrichment_GI_in_DepMap(df_qGI_nRed_depmap, df_depmap_abba_softer_merge_negPCC, df_fc, gi_direc='positive')
    pval_negPCC_negGI = set_enrichment_GI_in_DepMap(df_qGI_nRed_depmap, df_depmap_abba_softer_merge_negPCC, df_fc, gi_direc='negative')


    # No paralog version:
    print("\n\nNo paralog version:\n")

    # Ohno (older version)
    '''
    df_paralog = pd.read_csv('/Users/zhangxiang/Documents/Research_CB/GIN/Code_BROAD_LAB_DATA/suppressor_analysis_04142023/hsapiens.Pairs.Relaxed.2R.txt', sep='\t')
    df_paralogs = df_paralog[['Symbol1', 'Symbol2']]
    df_paralogs.columns = ['gene1', 'gene2']
    '''

    # Newer version:
    df_paralog = pd.read_csv('/Users/zhangxiang/Documents/Research_CB/GIN/Code_BROAD_LAB_DATA/suppressor_analysis_04142023/ensembl_pairs_complete_ident_max_g20.tsv', sep='\t')
    df_paralogs = df_paralog[['name1', 'name2']]
    df_paralogs.columns = ['gene1', 'gene2']

    df_qGI_nRed_depmap_noPara = exlude_paralogs(df_qGI_nRed_depmap, df_paralogs)
    df_depmap_abba_softer_merge_posPCC_noPara = exlude_paralogs(df_depmap_abba_softer_merge_posPCC, df_paralogs)
    df_depmap_abba_softer_merge_negPCC_noPara = exlude_paralogs(df_depmap_abba_softer_merge_negPCC, df_paralogs)

    pval_posPCC_posGI_np = set_enrichment_GI_in_DepMap(df_qGI_nRed_depmap_noPara, df_depmap_abba_softer_merge_posPCC_noPara, df_fc, frequency_filter=1.0, gi_direc='positive')
    pval_posPCC_negGI_np = set_enrichment_GI_in_DepMap(df_qGI_nRed_depmap_noPara, df_depmap_abba_softer_merge_posPCC_noPara, df_fc, frequency_filter=1.0, gi_direc='negative')
    pval_negPCC_posGI_np = set_enrichment_GI_in_DepMap(df_qGI_nRed_depmap_noPara, df_depmap_abba_softer_merge_negPCC_noPara, df_fc, frequency_filter=1.0, gi_direc='positive')
    pval_negPCC_negGI_np = set_enrichment_GI_in_DepMap(df_qGI_nRed_depmap_noPara, df_depmap_abba_softer_merge_negPCC_noPara, df_fc, frequency_filter=1.0, gi_direc='negative')

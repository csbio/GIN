# Author: Xiang Zhang (zhan6668)
# Description: This file starts at 07072024, in order to run the hypergeometric method pipeline on DepMap 22Q4 data
# and qGI_20211111 to check about gene pair enrichment (the functional information that the GIN dataset can add to DepMap)
# with 9 unique queries excluded due to their frequency in the intersection > 2%. p-value cutoff is 0.01.

# Update on 07072024: change the input file into "ABBA softer merge no Freq2"
# Will be deployed to lab server for running the co-annotation enrichment test

import os, sys, re
from math import *
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats

import pyreadr

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

def read_coAnno_standard_from_Rdata(input_direc):
    result = pyreadr.read_r(input_direc)
    print("Read in co-annotation standard:", result.keys()) # let's check what objects we got
    df_standard_co = result[list(result.keys())[0]]
    df_standard_co = df_standard_co.set_index('gene1')
    return df_standard_co

def get_annotated_pairs(df_pairs, df_standard, gene1_col='library_gene', gene2_col='query_gene', output_file=None):
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
    print(len(ab_overlap))
    print(len(ba_overlap))
    print(len(abba_overlap))
    
    num_bg_pairs = len(abba_overlap)
    print("Overlap background pairs", num_bg_pairs)
    
    df_standard_bg = df_standard_2idx.loc[abba_overlap]
    df_standard_bg_anno = df_standard_bg[df_standard_bg['is_annotated']==1]
    num_bg_anno = df_standard_bg_anno.shape[0]
    print("Overlap background annotations", num_bg_anno)
    
    
    df_ab_anno = df_standard_2idx.loc[ab_overlap][df_standard_2idx.loc[ab_overlap, 'is_annotated']==1]
    df_ba_anno = df_standard_2idx.loc[ba_overlap][df_standard_2idx.loc[ba_overlap, 'is_annotated']==1]
    df_ba_anno_reverse_idx = df_ba_anno.reset_index().set_index(['gene2', 'gene1'])
    
    # Get the overlap bg df for df_pairs (should be the same length as num_bg_pairs)
    df_pairs_bg = pd.concat([df_pairs_2idx_ab.loc[ab_overlap], 
                             df_pairs_2idx_ba.loc[ba_overlap].reset_index().set_index([gene1_col, gene2_col])])
    print(df_pairs_bg.shape)
    
    df_pairs_bg['is_annotated'] = 0
    df_pairs_bg.loc[df_ab_anno.index, 'is_annotated'] = 1
    df_pairs_bg.loc[df_ba_anno_reverse_idx.index, 'is_annotated'] = 1
    print(df_pairs_bg['is_annotated'].sum()) # should equal to num_bg_anno
    
    if output_file != None:
        df_pairs_bg.to_csv(output_file)
        #df_pairs.to_csv(output_file[:-4]+'_all.csv')
    
    return [num_bg_pairs, num_bg_anno]

def co_annotation_enrichment_pairwise(df, df_subset, df_standard, standard_name='complex', output_direc='./suppressor_analysis_04142023/'):
    
    [num_bg, num_bg_anno] = get_annotated_pairs(df, df_standard=df_standard, 
                                                 gene1_col='library_gene', gene2_col='query_gene', 
                                                  output_file=output_direc+standard_name+'_background_pairs.csv')
    
    [num_sub, num_sub_anno] = get_annotated_pairs(df_subset, df_standard=df_standard, 
                                                 gene1_col='library_gene', gene2_col='query_gene', 
                                                  output_file=output_direc+standard_name+'_subset_pairs.csv')
    
    # We can then compute a probability of drawing X red marbles out of N from 
    # a jar containing n red marbles out of M in the following way:
    pval = hyper_test(x=num_sub_anno, M=num_bg, n=num_bg_anno, N=num_sub)
    print(pval)
    return pval

def group_coanno_test(df_depmap_neg, df_depmap_pos, df_fc, df_coanno_standard, standard_name, output_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    df_depmap_neg_qGI_pos_fc = get_qGI_fc_subset(df_pairs=df_depmap_neg, df_fc=df_fc, interaction='positive')
    df_depmap_neg_qGI_neg_fc = get_qGI_fc_subset(df_pairs=df_depmap_neg, df_fc=df_fc, interaction='negative')
    df_depmap_pos_qGI_pos_fc = get_qGI_fc_subset(df_pairs=df_depmap_pos, df_fc=df_fc, interaction='positive')
    df_depmap_pos_qGI_neg_fc = get_qGI_fc_subset(df_pairs=df_depmap_pos, df_fc=df_fc, interaction='negative')

    print('\nnegPCC posGI')
    co_annotation_enrichment_pairwise(df=df_depmap_neg, df_subset=df_depmap_neg_qGI_pos_fc,
                                      df_standard=df_coanno_standard, standard_name=standard_name,
                                      output_direc=output_dir+'negPCC_posGI_')
    print('\nnegPCC negGI')
    co_annotation_enrichment_pairwise(df=df_depmap_neg, df_subset=df_depmap_neg_qGI_neg_fc,
                                      df_standard=df_coanno_standard, standard_name=standard_name,
                                      output_direc=output_dir+'negPCC_negGI_')
    print('\nposPCC posGI')
    co_annotation_enrichment_pairwise(df=df_depmap_pos, df_subset=df_depmap_pos_qGI_pos_fc,
                                      df_standard=df_coanno_standard, standard_name=standard_name,
                                      output_direc=output_dir+'posPCC_posGI_')
    print('\nposPCC negGI')
    co_annotation_enrichment_pairwise(df=df_depmap_pos, df_subset=df_depmap_pos_qGI_neg_fc,
                                      df_standard=df_coanno_standard, standard_name=standard_name,
                                      output_direc=output_dir+'posPCC_negGI_')
    print("\n\n\n")

def get_qGI_fc_subset(df_pairs, df_fc, interaction='positive', frequency_cutoff=1.0, skipStep1=True, skipStep2=True):
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
        # Skip this suppressor filter Step 2 to be more lenient
        if not skipStep2:
            df_pairs_GI_fc = df_pairs_GI_fc[df_pairs_GI_fc['qGI_score']>df_pairs_GI_fc['rich']]
    elif interaction == 'negative':
        # Update: Don't do any FC filter with negative side
        df_pairs_GI_fc  = df_pairs_GI # [df_pairs_GI['qGI_score']<df_pairs_GI['rich']]
    else:
        print('error')
    df_pairs_GI_fc.index.name = 'library_gene'
    df_pairs_GI_fc = df_pairs_GI_fc.reset_index()
    df_pairs_GI_fc_freq = df_pairs_GI_fc.groupby('query_gene').count()['library_gene']/df_pairs_GI_fc.shape[0]
    querys_to_keep = df_pairs_GI_fc_freq[df_pairs_GI_fc_freq<frequency_cutoff].index
    print("Remove pairs with query frequency >= "+str(frequency_cutoff), list(set(df_pairs_GI_fc_freq.index) - set(querys_to_keep)))
    df_pairs_GI_fc = df_pairs_GI_fc[df_pairs_GI_fc['query_gene'].isin(querys_to_keep)]
    print(len(df_pairs_GI_fc))
    return df_pairs_GI_fc

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
    print("Exclude paralogs:", len(abba_overlap))
    df_pairs_no_paralogs = df_pairs_2idx.drop(abba_overlap, axis=0).reset_index()
    return df_pairs_no_paralogs


if __name__ == "__main__":
    '''
    df_qGI_nRed_pairs = pd.read_csv('qGI_pairwise_20211111_nRed_cleaned.txt', index_col=0, sep='\t')
    df_qGI_nRed_pairs = df_qGI_nRed_pairs[['Query', 'qGI_score', 'FDR', 'GI_standard']]
    df_qGI_nRed_pairs.columns = ['query_gene', 'qGI_score', 'qGI_FDR', 'GI_standard']
    df_qGI_nRed_pairs.index.name='library_gene'
    '''

    # 07072024: softer merge, skip Step 1 & 2 (just do standard GI cutoff for both negative and positive), no 9 Freq2 query genes
    df_depmap_neg_abba_soft = pd.read_csv('DepMap_negative_ABBAsofter_PCCn01_pval001_noFreq2.csv', index_col=None)
    df_depmap_pos_abba_soft = pd.read_csv('DepMap_positive_ABBAsofter_PCC01_pval001_noFreq2.csv', index_col=None)

    print('Negative ABBA softer merge PCC (<-0.1), pval<0.01: # Pairs:', len(df_depmap_neg_abba_soft))
    print('Positive ABBA softer merge PCC (>0.1), pval<0.01: # Pairs:', len(df_depmap_pos_abba_soft))

    df_fc = pd.read_csv('./dataset/fc_singlePhenotype_20211111.txt', index_col=0, sep='\t')

    # 06092023: Exclude paralogs
    # Old version (Ohno)
    #df_paralog = pd.read_csv('./dataset/hsapiens.Pairs.Relaxed.2R.txt', sep='\t')
    #df_paralogs = df_paralog[['Symbol1', 'Symbol2']]
    #df_paralogs.columns = ['gene1', 'gene2']

    # Newer version:
    df_paralog = pd.read_csv('./dataset/ensembl_pairs_complete_ident_max_g20.tsv', sep='\t')
    df_paralogs = df_paralog[['name1', 'name2']]
    df_paralogs.columns = ['gene1', 'gene2']

    # CORUM enrichment
    corum_standard_direc = './dataset/corum.Rdata' 
    df_corum_standard_co = read_coAnno_standard_from_Rdata(corum_standard_direc)

    # GO BP enrichment
    gobp_standard_direc = './dataset/GO_BP.Rdata'
    df_gobp_standard_co = read_coAnno_standard_from_Rdata(gobp_standard_direc)

    # Pathway enrichment
    path_standard_direc = './dataset/Pathway.Rdata'
    df_path_standard_co = read_coAnno_standard_from_Rdata(path_standard_direc)

    # localization enrichment
    loc_standard_direc = './dataset/localization.Rdata'
    df_loc_standard_co = read_coAnno_standard_from_Rdata(loc_standard_direc)

    # PPI enrichment
    df_ppi_standard_co = pd.read_csv('/project/csbio/Xiang/GIN/clustering/ppi_enrich/coannotation_PPI_pairs.txt', sep='\t')
    df_ppi_standard_co = df_ppi_standard_co.set_index('gene1')
    
    # paralog enrichment
    #df_paralog_standard_co = pd.read_csv('./dataset/GIN_screen_bg_pairs_paralogs_annotated.tsv', sep='\t')
    #df_paralog_standard_co = pd.read_csv('./dataset/coannotation_paralog_pairs.txt', sep='\t')
    #df_paralog_standard_co = df_paralog_standard_co.set_index('gene1')

    '''
    df_qGI_nRed_depmap = pd.read_csv('./DepMap_22Q4_all_qGI_nRed_pairs.tsv', sep='\t', index_col=None)

    df_depmap_abba_softer_merge_posPCC = softer_filter_abba(df_qGI_nRed_depmap, pcc_ab_col='PCC_lib_fitness_vs_query_expression',
            pcc_ba_col='PCC_query_fitness_vs_lib_expression', output_direc='DepMap_positive_ABBAsofter_PCC01_pval001.csv',
                    pcc_cutoff=0.1, pcc_pval=0.01)

    df_depmap_abba_softer_merge_negPCC = softer_filter_abba(df_qGI_nRed_depmap, pcc_ab_col='PCC_lib_fitness_vs_query_expression',
            pcc_ba_col='PCC_query_fitness_vs_lib_expression', output_direc='DepMap_negative_ABBAsofter_PCCn01_pval001.csv',
                    pcc_cutoff=-0.1, pcc_pval=0.01)

    ''' 

    # 07282024: Repeat noPara with new paralog standard - temporarily comment out the w/ para version
    '''
    # Set output directory
    output_dir_softer = "./output_softer_skipStep12_noFreq2/"

    # Special paralog enrichment
    group_coanno_test(df_depmap_neg_abba_soft, df_depmap_pos_abba_soft, df_fc, 
            df_coanno_standard=df_paralog_standard_co, standard_name="paralog", output_dir=output_dir_softer)

    group_coanno_test(df_depmap_neg_abba_soft, df_depmap_pos_abba_soft, df_fc, 
            df_coanno_standard=df_corum_standard_co, standard_name="CORUM", output_dir=output_dir_softer)
    group_coanno_test(df_depmap_neg_abba_soft, df_depmap_pos_abba_soft, df_fc, 
            df_coanno_standard=df_gobp_standard_co, standard_name="GOBP", output_dir=output_dir_softer)
    group_coanno_test(df_depmap_neg_abba_soft, df_depmap_pos_abba_soft, df_fc, 
            df_coanno_standard=df_path_standard_co, standard_name="Pathway", output_dir=output_dir_softer)
    group_coanno_test(df_depmap_neg_abba_soft, df_depmap_pos_abba_soft, df_fc, 
            df_coanno_standard=df_loc_standard_co, standard_name="localization", output_dir=output_dir_softer)
    group_coanno_test(df_depmap_neg_abba_soft, df_depmap_pos_abba_soft, df_fc, 
            df_coanno_standard=df_ppi_standard_co, standard_name="PPI", output_dir=output_dir_softer)
    '''

    # Exclude paralogs:
    print("\n\n\n\n\nExclude paralogs:")
    df_depmap_neg_abba_soft = exlude_paralogs(df_depmap_neg_abba_soft, df_paralogs)
    df_depmap_pos_abba_soft = exlude_paralogs(df_depmap_pos_abba_soft, df_paralogs)
    
    output_dir_softer_np = "./output_softer_skipStep12_noFreq2_noPara_new/"

    group_coanno_test(df_depmap_neg_abba_soft, df_depmap_pos_abba_soft, df_fc, 
            df_coanno_standard=df_corum_standard_co, standard_name="CORUM", output_dir=output_dir_softer_np)
    group_coanno_test(df_depmap_neg_abba_soft, df_depmap_pos_abba_soft, df_fc, 
            df_coanno_standard=df_gobp_standard_co, standard_name="GOBP", output_dir=output_dir_softer_np)
    group_coanno_test(df_depmap_neg_abba_soft, df_depmap_pos_abba_soft, df_fc, 
            df_coanno_standard=df_path_standard_co, standard_name="Pathway", output_dir=output_dir_softer_np)
    group_coanno_test(df_depmap_neg_abba_soft, df_depmap_pos_abba_soft, df_fc, 
            df_coanno_standard=df_loc_standard_co, standard_name="localization", output_dir=output_dir_softer_np)
    group_coanno_test(df_depmap_neg_abba_soft, df_depmap_pos_abba_soft, df_fc, 
            df_coanno_standard=df_ppi_standard_co, standard_name="PPI", output_dir=output_dir_softer_np)

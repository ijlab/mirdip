#!/usr/bin/env python

import os
import time
import pandas as pd
import glob
import os
import csv
import re
#from Bio import SeqIO
#from Bio.Seq import Seq
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn.model_selection
import sklearn.datasets
import sklearn.metrics
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import PrecisionRecallDisplay
from sklearn.metrics import roc_curve
from sklearn.metrics import RocCurveDisplay
import math
import types

import seaborn as sns
# for plotting PR, ROC curves

# for figuring out normalization etc.
from pandas_profiling import ProfileReport

def map_genes_with_hgnc(df,
                        left_df_column_name,
                        hgnc_column_name,
                        hgnc=None,
                        check_ensembl=False,
                        check_alias_symbols=False,
                        check_prev_symbols=False,
                        keep_mappings=['symbol','alias_symbol','prev_symbol', 'entrez_id', 
                                        'ensembl_gene_id', 'refseq_accession', 'uniprot_ids', 
                                        'prev_symbols_list', 'alias_symbols_list', 'refseq_accession_list', 'mane_select']):
    """
    Uses the HGNC ID mapping table available at:
    http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt
    """
    if check_alias_symbols:
        if 'symbol_for_merge' not in df.columns: 
            df['symbol_for_merge'] = df[left_df_column_name]
        hgnc_alias = pd.read_csv('/gpfs/lb/mirdip5/hgnc/complete_set_by_ALIAS_symbol_exploded.tsv', sep='\t', header=0)
        df = df.merge(hgnc_alias[[c for c in keep_mappings if c in hgnc_alias.columns]], left_on='symbol_for_merge', right_on='alias_symbols_list', how='left')
        df[left_df_column_name] = df['symbol_for_merge']
        hgnc_alias=None
    if check_prev_symbols:
        if 'symbol_for_merge' not in df.columns: 
            # If isn't already a column called symbol_for_merge, probably due to check_alias=False
            df['symbol_for_merge'] = df[left_df_column_name]
        hgnc_prev = pd.read_csv('/gpfs/lb/mirdip5/hgnc/complete_set_by_PREVIOUS_symbol_exploded.tsv', sep='\t', header=0)
        df = df.merge(hgnc_prev[[c for c in keep_mappings if c in hgnc_prev.columns]], left_on='symbol_for_merge', right_on='prev_symbols_list', how='left')
        df[left_df_column_name] = df['symbol_for_merge']
        hgnc_prev=None

    if not hgnc or hgnc.empty:
        hgnc = pd.read_csv('/gpfs/lb/mirdip5/hgnc/hgnc_complete_set.txt', sep='\t', header=0)
    
    if hgnc_column_name == 'refseq_accession':
        hgnc_refseq = pd.read_csv('/gpfs/lb/mirdip5/hgnc/complete_set_by_REFSEQ_ACCESSION_exploded.tsv', sep='\t', header=0)
        df = df.merge(hgnc_refseq[[c for c in keep_mappings if c in hgnc_refseq.columns]], left_on=left_df_column_name, right_on='refseq_accession_list', how='left')
    else:
        df = df.merge(hgnc[[c for c in keep_mappings if c in hgnc.columns]], left_on=left_df_column_name, right_on=hgnc_column_name, how='left')
    
    # this is for ensembl specifically, using Chiara's mapping downloaded from Biomart.
    if check_ensembl:
        biomart_ids = pd.read_csv('/gpfs/lb/mirdip5/mart_export_ensembl_hgnc.txt', sep='\t', header=0)
        biomart_grouped_map = biomart_ids.groupby('Gene stable ID').agg(set).reset_index()
        biomart_grouped_map = biomart_grouped_map[['Gene stable ID', 'HGNC symbol', 'Transcript stable ID']]
        biomart_grouped_map['symbol'] = biomart_grouped_map['HGNC symbol'].apply(lambda x: list(x)[0] if len(list(x)) == 1 else ','.join(x))
        biomart_grouped_map['ensembl_transcripts'] = biomart_grouped_map['Transcript stable ID'].apply(lambda x: [y for y in list(x)])
        biomart_grouped_map = biomart_grouped_map.rename(columns={'Gene stable ID':'ensembl_gene_id'})
        biomart_grouped_map = biomart_grouped_map[['ensembl_gene_id', 'symbol', 'ensembl_transcripts']]
        if 'ensembl_gene_id' in df.columns:
            df = df.merge(biomart_grouped_map, on='ensembl_gene_id', how='left')
            # if HGNC was null, use ensembl's
            df['symbol'] = df['symbol_x'].combine_first(df['symbol_y'])
            # now check for mane_select transcript in the transcript stable ids
            #df['check_mane_select'] = df['mane_select'].fillna(value='').str.split('|').apply(lambda x: x[0].split('.')[0] if len(x) > 1 else '')
            #df = df
            # Enumerated above it's all of 25 entries. Skip.
    return df

def normalize_scores(df, score_column_name, normalized_score_column_name='', group=1, score_func=None):
    assert group in [1,2], "Group 1 is 1-x, group 2 is x for min/max norm. This is how Tomas did it, and corresponds to config of pipeline."
    df_min = df[score_column_name].astype(float).min()
    df_max = df[score_column_name].astype(float).max()
    if not normalized_score_column_name:
        normalized_score_column_name = score_column_name + '_norm'
    if score_func and isinstance(score_func, types.FunctionType):
        df[normalized_score_column_name] = df[score_column_name].apply(score_func)
    else:
        if group == 2:
            df[normalized_score_column_name] = df[score_column_name].apply(lambda x: (float(x)-df_min)/(df_max-df_min))
        elif group == 1:
            df[normalized_score_column_name] = df[score_column_name].apply(lambda x: 1 - (float(x)-df_min)/(df_max-df_min))
        else:
            # better, just set to the score
            df[normalized_score_column_name] = df[score_column_name]
    return df


def group_mir_gene_pairs_and_take_ranked_product(df,
                                                 score_col,
                                                 gene_col='symbol',
                                                 mir_col='mirdip4_mirbase_id',
                                                 num_scores=3):
    cols = list(df.columns)
    # ['symbol', 'mirdip4_mirbase_id', 'mirmap_score', 'mirmap_score_norm',
    # 'data_source', 'original_gene_symbol', 'original_mirbase_id']
    not_norm_score_col = score_col.replace('_norm', '')
    top3 = df.groupby([gene_col, mir_col]).agg(
        **{
            not_norm_score_col: pd.NamedAgg(column=not_norm_score_col, aggfunc=np.median),
            score_col: pd.NamedAgg(column=score_col, aggfunc=lambda x: np.prod([y for y in sorted(x)][0:num_scores])),
            'data_source': pd.NamedAgg(column='data_source', aggfunc='first'),
            'original_gene_symbol': pd.NamedAgg(column='original_gene_symbol', aggfunc='first'),
            'original_mirbase_id': pd.NamedAgg(column='original_mirbase_id', aggfunc='first')
        }
    )
    #top3 = df.groupby([gene_col, mir_col]).apply(lambda row: np.prod([x for x in sorted(row[score_col])][0:num_scores]))
    top3 = top3.reset_index()
    top3 = top3.drop_duplicates()
    #deduplicated = top3.merge(df, left_on=[gene_col, mir_col], right_on=[gene_col, mir_col], how='left')
    #deduplicated = deduplicated.reset_index(drop=True)
    #deduplicated = deduplicated[cols]
    #deduplicated = deduplicated.drop_duplicates()
    #if 'original_gene_symbol' in cols and 'original_mirbase_id' in cols:
    #    deduplicated = deduplicated.groupby([gene_col, mir_col, score_col]).apply(lambda row: ).reset_index()
    #    for col in ['original_gene_symbol', 'original_mirbase_id']:
    #        deduplicated[col] = deduplicated[col].apply(lambda x: ','.join(x))
    return top3

utrs_to_exclude = [u for u in open('/gpfs/lb/mirdip5/utrs/entries_under_25nt.txt','r')]
def standardize_normalize_integrate_pita_run(predictions,
                                                                                          save_cleaned_results_to,
                                                                                          names=['ensembl_gene_id', 'mirdip4_mirbase_id', 'num_sites', 'pita_score'],
                                                                                          gene_colname='ensembl_gene_id',
                                                                                          mir_colname='mirdip4_mirbase_id',
                                                                                          score_colname='pita_score',
                                                                                          data_source='PITA_May_2021',
                                                                                          max_energy_threshold=-13.0,
                                                                                          utrs_to_exclude=utrs_to_exclude,
                                                                                          group=2):
    pita = pd.read_csv(predictions, sep='\t', header=None, names=names)
    pita = pita.astype({score_colname:'float64'})
    pita = pita[pita['pita_score'] <= max_energy_threshold]
    pita = pita[~(pita[gene_colname].isin(utrs_to_exclude))]
    pita[gene_colname] = pita[gene_colname].str.split('|').apply(lambda x: x[0])
    pita[mir_colname] = pita[mir_colname].str.split(' ').apply(lambda x: x[0])
    pita = normalize_scores(pita, score_colname, group=group)
    pita['original_gene_symbol'] = ''
    pita['original_mirbase_id'] = ''
    pita['data_source'] = data_source
    pita = map_genes_with_hgnc(pita, 'ensembl_gene_id', 'ensembl_gene_id', check_ensembl=True, hgnc=None)
    pita = pita[['symbol', 'mirdip4_mirbase_id', 'pita_score', 'pita_score_norm', 'data_source', 'original_gene_symbol', 'original_mirbase_id']]
    pita.replace(r'^\s*$', np.nan, regex=True)
    pita.dropna(inplace=True)
    pita = pita.drop_duplicates()
    pita = group_mir_gene_pairs_and_take_ranked_product(pita, 'pita_score_norm')    
    pita = normalize_scores(pita, score_colname+'_norm', normalized_score_column_name=score_colname+'_norm', group=group)
    pita.sort_values(by='pita_score_norm', inplace=True, ascending=True)
    pita['original_gene_symbol'] = pita['symbol']
    pita['original_mirbase_id'] = pita['mirdip4_mirbase_id']
    pita.to_csv(save_cleaned_results_to, sep='\t',index=False,header=False, quoting=csv.QUOTE_MINIMAL)
standardize_normalize_integrate_pita_run('/home/waddelld/rnatools/pita/work_known/pita_known_mirdip5.tsv',
                                         '/home/waddelld/rnatools/mirdip5/resources_redo_final/pita_with_neg13_cutoff.txt')
        

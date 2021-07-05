#!/usr/bin/env python

from docopt import docopt
import pandas as pd

import os
import time
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


doc = """Usage: mirdip5_update_data_source_ids.py [--header] [--sep SEP] [--mircol MIRCOL] [--genecol GENECOL] [--deadext DEADEXT] [--changedext CHANGEDEXT] --path PATH --checker CHECKER --outpath OUTPATH

-h --help                     Show this message.
-n --header                   If provided first line will be read as column headers.
-s --sep=<SEP>                Field separator.
-p --path=<PATH>              Path for the dataframe to update IDs of, expected format is mirdip4 final_resources.
-k --checker=<CHECKER>        HGNC symbol checker file.
-m --mircol=<MIRCOL>          Column name or position of micrornas.
-g --genecol=<GENECOL>        Column name or position of genes.
-d --deadext=<DEADEXT>        Extension of "dead" ids from mirbaseconverter, drop these.
-c --changedext=<CHANGEDEXT>  Extension of changed ids from mirbaseconverter, update these.
-o --outpath=<OUTPATH>        Path for the saved dataframe with converted IDs.

"""

def update_mirdip4_gene_symbols(df, hgnc_symbol_checker_results_path, genecolname='mirdip4_gene_symbol'):
    hgnc_mapping = pd.read_csv(hgnc_symbol_checker_results_path, sep=',', header=0, skiprows=1)
    hgnc_mapping = hgnc_mapping[['Input', 'Match type', 'Approved symbol']]
    dead_id_list = [i for i in list(hgnc_mapping[hgnc_mapping['Match type'] == 'Entry withdrawn']['Input'])]
    changed_id_map = {k:v for k, v in zip(list(hgnc_mapping[hgnc_mapping['Match type'] == 'Previous symbol']['Input']),
                                          list(hgnc_mapping[hgnc_mapping['Match type'] == 'Previous symbol']['Approved symbol']))}
    df['original_gene_symbol'] = df[genecolname]
    updated_ids = df[genecolname].map(changed_id_map)
    df[genecolname] = updated_ids.combine_first(df[genecolname])
    df = df[~df[genecolname].isin(dead_id_list)]
    return df

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

def update_mirbase_ids(options):
    df_path = options['--path']
    header = options['--header']
    separator = options['--sep'] if options['--sep'] else '\t'
    mircolname = options['--mircol'] if options['--mircol'] else 'mirdip4_mirbase_id'
    deadext = options['--deadext'] if options['--deadext'] else '.mirbase_v22_IDs.tsv.dead.uniq'
    changedext = options['--changedext'] if options['--changedext'] else '.mirbase_v22_IDs.tsv.changed.uniq'
    df = pd.read_csv(df_path, sep=separator, header=None, quoting=csv.QUOTE_NONE, names=['symbol', 'mirdip4_mirbase_id', 'score', 'score_norm', 'data_source', 'original_gene_symbol', 'original_mirbase_id'])
    # should be mirdip5 now
    dead = pd.read_csv(df_path+deadext,
                        sep=separator, 
                        header=None,
                        quoting=csv.QUOTE_NONE,
                        names=['dead_ids', 'empty_column1', 'empty_column2'])
    # we will use this list to drop any rows found to have one of these
    dead_id_list = [i for i in list(dead['dead_ids'])]
    changed = pd.read_csv(df_path+changedext,
                        sep=separator, 
                        header=0,
                        quoting=csv.QUOTE_NONE,
                        names=['from', 'to', 'accession'])
    changed['to'] = changed['to'].fillna('--')
    changed_id_map = {k:v for k, v in zip(list(changed['from']), list(changed['to'])) if v != '--'}
    # save new "original"
    df['original_mirbase_id'] = df[mircolname]
    updated_ids = df[mircolname].map(changed_id_map)
    df[mircolname] = updated_ids.combine_first(df[mircolname])
    df = df[~df[mircolname].isin(dead_id_list)]
    return df

if __name__ == "__main__":
    arguments = docopt(doc, version='mirdip5_update_mirbase_ids 0.1')
    df = update_mirbase_ids(arguments)
    #df = map_genes_with_hgnc(df, 'symbol', 'symbol', hgnc=None)
    # NOTE: symbol is a different default genecol than mirdip4_gene_symbol, which was used in the juypyter notebook
    genecolname_to_update = arguments['--genecol'] if arguments['--genecol'] else 'symbol'
    df = update_mirdip4_gene_symbols(df, arguments['--checker'], genecolname=genecolname_to_update)
    df.to_csv(arguments['--outpath'], sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE)

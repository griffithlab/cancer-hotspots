#!/usr/bin/env python

import pandas as pd
import os
import sys
import csv

gene = sys.argv[1]

# dataframe of allele frequencies scraped from allelefrequencies.net
freqs = pd.read_csv('/gscmnt/gc2142/griffithlab/abasu/Global_HLA_Freq.csv')
freqs.set_index("Allele", inplace=True)

# directory stores lists of alleles + scores + anchoring filtering info for all 214 genes
directory = '/gscmnt/gc2142/griffithlab/abasu/gene_filtered_neo/'

path = directory + gene + ".tsv"
# dataframe of alleles, binding affinities, seq_length, mutation position, wt score, fold_change for all neoantigens for 1 gene
df = pd.read_csv(path, delimiter='\t', sep='\t',header=None)
df.columns = ["Allele", "Seq_len", "Mut_pos", "MT_score", "WT_score", "Fold_change"]

row_remove = []
for row in range(0,len(df.index)):
    hla = df.iloc[row,0] 
    seq_length = df.iloc[row,1]  
    mutation_position = df.iloc[row,2]
    wt_score = df.iloc[row,4]
    fold_change = df.iloc[row,5]

    # class I - first and last 3 sequence positions are considered anchor positions
    anchor_positions = list(range(1,seq_length+1))[0:3] + list(range(1,seq_length+1))[(seq_length-3):seq_length]
    if 'HLA-A' in hla or 'HLA-B' in hla or 'HLA-C' in hla:
        # mt does not bind better than wt
        if fold_change <= 1:
            # if mutation in an anchor position, remove from df (Group3A)
            if mutation_position in anchor_positions:
                row_remove.append(row)
        # mt does bind better than wt
        elif fold_change > 1:
            # if mutation is in anchor position
            if mutation_position in anchor_positions:
                # and wt score is also a good binder, remove from df (Group3B)
                if wt_score <= 500:
                    row_remove.append(row)


df = df.drop(df.index[row_remove])

freq_data = []
for r in range(0,len(df.index.values)):
    allele = df.iloc[r,0]
    if allele in freqs.index.values:
        freq = freqs.loc[[allele][0]][0]
        freq_data.append(freq)
    else:
        freq_data.append(0)

max_freq = max(freq_data)
neo_count = len(df.index.values)
summary = '/gscmnt/gc2142/griffithlab/abasu/hotspot_graph_summary_data_filtered/' + gene + '_Summary_Fig_Data.tsv'

with open(summary, 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow([gene, max_freq, neo_count])




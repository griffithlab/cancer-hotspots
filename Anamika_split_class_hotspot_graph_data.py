#!/usr/bin/env python

import pandas as pd
import os
import sys
import csv

gene = sys.argv[1]

# directory stores lists of alleles + scores for all 214 genes
directory = '/gscmnt/gc2142/griffithlab/abasu/gene_binding_affinity/'

path = directory + gene + ".tsv"
# dataframe of alleles and binding affinities for 1 gene
allele_score = pd.read_csv(path, delimiter='\t', sep='\t',header=None)
allele_score.columns = ["Allele", "Score"]

# split class 1 and class 2 into different dataframes
class1_allele = []
class2_allele = []
class1_score = []
class2_score = []
for r in range(0, len(allele_score)):
    if "D" in allele_score.iloc[r, 0]:
        class2_allele.append(allele_score.iloc[r, 0])
        class2_score.append(allele_score.iloc[r, 1])
    else: 
        class1_allele.append(allele_score.iloc[r, 0])
        class1_score.append(allele_score.iloc[r, 1])

class1 = pd.DataFrame(class1_allele, columns =['Allele']) 
class1['Score'] = class1_score
class2 = pd.DataFrame(class2_allele, columns =['Allele']) 
class2['Score'] = class2_score

bins = [0, 50, 500]

s1 = pd.cut(class1['Score'], bins=bins).value_counts()
s2 = pd.cut(class2['Score'], bins=bins).value_counts()
s1.index = s1.index.astype(str)
s2.index = s2.index.astype(str)

graph_d_c1 = '/gscmnt/gc2142/griffithlab/abasu/split_class_hotspot_graph_data/' + gene + '_Graph_D_Class1_Data.tsv'
graph_d_c2 = '/gscmnt/gc2142/griffithlab/abasu/split_class_hotspot_graph_data/' + gene + '_Graph_D_Class2_Data.tsv'

with open(graph_d_c1, 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow([gene, s1.loc[['(0, 50]']][0],"0-50 nM"])
    tsv_writer.writerow([gene, s1.loc[['(50, 500]']][0],"50-500 nM"])

with open(graph_d_c2, 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow([gene, s2.loc[['(0, 50]']][0],"0-50 nM"])
    tsv_writer.writerow([gene, s2.loc[['(50, 500]']][0],"50-500 nM"])
    
# dataframe of allele frequencies scraped from allelefrequencies.net
freqs = pd.read_csv('/gscmnt/gc2142/griffithlab/abasu/Global_HLA_Freq.csv')
freqs.set_index("Allele", inplace=True)

allele_freq_c1 = []
for r in range(0,len(class1.index.values)):
    allele = class1.iloc[r,0]
    if allele in freqs.index.values:
        freq = freqs.loc[[allele][0]][0]
        allele_freq_c1.append(freq)
    else:
        allele_freq_c1.append(0)
        
allele_freq_c2 = []
for r in range(0,len(class2.index.values)):
    allele = class2.iloc[r,0]
    if allele in freqs.index.values:
        freq = freqs.loc[[allele][0]][0]
        allele_freq_c2.append(freq)
    else:
        allele_freq_c2.append(0)

class1['Frequency'] = allele_freq_c1
class2['Frequency'] = allele_freq_c2

bins_1 = [-1, 1, 5, 10, 15, 20]
bins_2 = [-1, 1, 5, 10]

s1 = pd.cut(class1['Frequency'], bins=bins_1).value_counts()
s2 = pd.cut(class2['Frequency'], bins=bins_2).value_counts()
s1.index = s1.index.astype(str)
s2.index = s2.index.astype(str)

graph_e_c1 = '/gscmnt/gc2142/griffithlab/abasu/split_class_hotspot_graph_data/' + gene + '_Graph_E_Class1_Data.tsv'
graph_e_c2 = '/gscmnt/gc2142/griffithlab/abasu/split_class_hotspot_graph_data/' + gene + '_Graph_E_Class2_Data.tsv'

with open(graph_e_c1, 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow([gene, s1.loc[['(-1, 1]']][0],"0-1%"])
    tsv_writer.writerow([gene, s1.loc[['(1, 5]']][0],"1-5%"])
    tsv_writer.writerow([gene, s1.loc[['(5, 10]']][0],"5-10%"])
    tsv_writer.writerow([gene, s1.loc[['(10, 15]']][0],"10-15%"])
    tsv_writer.writerow([gene, s1.loc[['(15, 20]']][0],"15-20%"])
    
with open(graph_e_c2, 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow([gene, s2.loc[['(-1, 1]']][0],"0-1%"])
    tsv_writer.writerow([gene, s2.loc[['(1, 5]']][0],"1-5%"])
    tsv_writer.writerow([gene, s2.loc[['(5, 10]']][0],"5-10%"])

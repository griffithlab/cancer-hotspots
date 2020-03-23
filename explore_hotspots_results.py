"""
Usage: explore_hotspots_results.py input_dir output_file 
"""

import os
import sys
import re
import pandas as pd

# input arguments
working_dir=os.getcwd()
input_dir=sys.argv[1]
output_filename=sys.argv[2]
output_file=f'{working_dir}/{output_filename}'
          
# functions
def length_count(df, message):
    unique_seqs = df['MT Epitope Seq'].unique().tolist()
    len_seqs = [len(x) for x in unique_seqs]
    len_dict = {f'{message} Length Count {x}': len_seqs.count(x) for x in set(len_seqs)}
    len_df = pd.DataFrame(len_dict, index=[0])
    return len_df

# list of files to parse
open_files = [f for f in sorted(os.listdir(f'{working_dir}/{input_dir}/'))] # if os.path.isfile(f) and f.endswith("_TUMOR.all_epitopes.tsv")]
for file1 in open_files:
    df_input = pd.read_csv(f'{input_dir}/{file1}', sep='\t')
    # identifying constant and variable portions of the file
    df_constant = df_input[["Chromosome", "Start", "Stop", "Reference", "Variant", "Transcript", "Ensembl Gene ID", "Variant Type", "Mutation", "Protein Position", "Gene Name", "HGVSp"]]
    df_variable = df_input[["HLA Allele", "Peptide Length", "Sub-peptide Position", "Mutation Position", "MT Epitope Seq", "WT Epitope Seq", "Median MT Score", "Median WT Score", "Median Fold Change"]]
    df_variable = df_variable.drop_duplicates()
    # class I and II predictions
    #df_variable_classI = df_variable.loc[df_variable['HLA Allele'].str.contains("HLA-A|HLA-B|HLA-C")]
    #df_variable_classII = df_variable.loc[~df_variable['HLA Allele'].str.contains("HLA-A|HLA-B|HLA-C")]

    # explore df.variable
    variable_dict = {}
    columns = ["HLA Allele", "Peptide Length", "Mutation Position", "MT Epitope Seq"]
    for item in columns:
        item_set = df_variable[item].unique()
        unique_count = len(item_set)
        variable_dict[item] = item_set.tolist()
    
    # counting how many seqs go with each hla allele
    len_dict = {}
    for item in df_variable['HLA Allele'].unique():
        df_subset = df_variable.loc[df_variable['HLA Allele'] == item]
        df_len = len(df_subset.index)
        len_dict[item] = df_len
    
    # how many hla alleles
    alleles_num = len(variable_dict["HLA Allele"])
    # how many neoantigen sequences
    seqs_num = len(variable_dict["MT Epitope Seq"])
    # score category counts from original file
    score_50 = len(df_variable.loc[df_variable["Median MT Score"] < 50])
    score_500 = len(df_variable.loc[(df_variable["Median MT Score"] >= 50) & (df_variable["Median MT Score"] <= 500)])
    score_over_500 = len(df_variable.loc[df_variable["Median MT Score"] > 500])
    # df with only Median MT Score <= 500 nM
    df_neoantigen = df_variable.loc[df_variable["Median MT Score"] <= 500]
    # how many hla alleles in df_neoantigen
    alleles_num_neo = len(df_neoantigen["HLA Allele"].unique().tolist())
    # how many neoantigen sequences in df_neoantigen
    seqs_num_neo = len(df_neoantigen["MT Epitope Seq"].unique().tolist())
    # minimun median score - full line
    min_score_line = df_neoantigen[df_neoantigen["Median MT Score"] == df_neoantigen["Median MT Score"].min()]
    # check out values[0] error
    if len(min_score_line) != 0:
    	min_allele = min_score_line["HLA Allele"].values[0]
    	min_seqs = min_score_line["MT Epitope Seq"].values[0]
    	min_score = min_score_line["Median MT Score"].values[0]
    else:
        print('df_neoantigen length:' + str(len(df_neoantigen)))
        min_allele = "NA"
        min_seqs = "NA"
        min_score = "NA"
    # length counts
    len_df = length_count(df_variable, "Seq")
    neo_df = length_count(df_neoantigen, "Neo")
    # constant variant info
    constant_line = df_constant.drop_duplicates()
    # colnames of file
    final_cols_to_df = ["Total Predictions", "HLA Alleles", "MT Epitope Seqs", "Scores 50", "Scores 50-500", "Scores 500",
                        "Total Neoantigens", "HLA Alleles Neo", "MT Epitope Seqs Neo",
                        "HLA Allele Min", "MT Epitope Seq Min", "Median MT Score Min"]
    # output for each column
    final_cols_output = [len(df_variable), alleles_num, seqs_num, score_50, score_500, score_over_500, len(df_neoantigen), alleles_num_neo, seqs_num_neo, min_allele, min_seqs, min_score]
    # zip the two and create a dictionary
    final_cols_dict = {x:y for x,y in zip(final_cols_to_df, final_cols_output)}
    # create a df from the dictionary
    final_cols_df = pd.DataFrame(final_cols_dict, index=[0])
    # concatenate all dfs
    final_write_df = pd.concat([constant_line, final_cols_df, len_df, neo_df], axis=1, join='inner')
    print(file1)
    # add df to final output file
    if os.path.isfile(output_file):
        final_write_df.to_csv(output_file, mode = 'a', index=False, header=None)
    else:
        final_write_df.to_csv(output_file, mode = 'w', index=False)

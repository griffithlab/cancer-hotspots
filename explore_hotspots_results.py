"""
Usage: python3 explore_hotspots_results.py input_dir output_file

Parses hotspot pVACtools results files for further analysis
"""

import os
import sys
import re
import pandas as pd

def identify_neoantigens(df):
    """Identify class I and II neoantigens from pVACtools total predictions df

    Args:
        df (dataframe): pVACtools total predictions df

    Returns:
        tuple: df for classI neoantigens, df for classII neoantigens, dictionary with score-based output
    """
    # score-based output
    df_low_score = df_variable.loc[df_variable["Median MT Score"] <= 50]
    df_med_score = df_variable.loc[(df_variable["Median MT Score"] > 50) & (df_variable["Median MT Score"] <= 500)]
    df_high_score = df_variable.loc[(df_variable["Median MT Score"] > 500) & (df_variable["Median MT Score"] <= 1000)]
    score = [("low", len(df_low_score)), ("med", len(df_med_score)), ("high", len(df_high_score))]
    score_dict = {x:y for (x,y) in score}
    # anchor neoantigen output
    df_neoantigen = pd.concat([df_low_score, df_med_score])
    class_I = [i for i in df_neoantigen["HLA Allele"] if any(x in i for x in ["HLA-A", "HLA-B", "HLA-C"])]
    class_II = [i for i in df_neoantigen["HLA Allele"] if any(x in i for x in ["DP", "DQ", "DR"])]
    index_list_classI = []
    index_list_classII = []
    for row_index,row in df_neoantigen.iterrows():
        hla, wt_score, fold_change, mutation_position, seq_length = (row["HLA Allele"], row["Median WT Score"], row["Median Fold Change"], row["Mutation Position"], row["Peptide Length"])
        # class I - first and last 3 sequence positions are considered anchor positions
        anchor_positions = list(range(1,seq_length+1))[0:3] + list(range(1,seq_length+1))[(seq_length-3):seq_length]
        if hla in class_I:
            # mt does not bind better than wt
            if fold_change <= 1:
                # if mutation NOT in an anchor position
                if mutation_position not in anchor_positions:
                    #print(f'{row_index} fc <= 1, not in anchor pos')
                    index_list_classI.append(row_index)
            # mt does bind better than wt
            elif fold_change > 1:
                # if mutation NOT in an anchor position
                if mutation_position not in anchor_positions:
                    #print(f'{row_index} fc > 1, not in anchor pos')
                    index_list_classI.append(row_index)
                # if mutation is in anchor position
                elif mutation_position in anchor_positions:
                    # and wt score is at least high than a typical neoantigen
                    if wt_score > 500:
                        #print(f'{row_index} fc > 1, in anchor pos, wt score > 500')
                        index_list_classI.append(row_index)
        # no anchor position parameters - only filter binding affinity < 500nM
        elif hla in class_II:
            index_list_classII.append(row_index)
    df_neo_output_classI = df_neoantigen.loc[index_list_classI, :]
    df_neo_output_classII = df_neoantigen.loc[index_list_classII, :]
    return df_neo_output_classI, df_neo_output_classII, score_dict

def parse_predictions(df, messages):
    """Return general stats of df predictions

    Args:
        df (dataframe): pVACtools df
        messages (list): descriptive name + abbreviation for resulting column ex: [Neoantigens, Neo]

    Returns:
        dataframe: column = name, row = count
    """
    # explore df.variable
    variable_dict = {}
    columns = ["HLA Allele", "MT Epitope Seq"]
    variable_dict[messages[0]] = len(df)
    for item in columns:
        item_set = df[item].unique()
        item_name = f'item {messages[1]}'
        variable_dict[item_name] = [len(item_set)]
        output_df = pd.DataFrame(variable_dict, index=[0])
    return output_df

def find_min_prediction(df):
    """Find the best scoring prediction for class I and II HLA alleles

    Args:
        df (dataframe): pVACtools class I/II neoantigens df

    Returns:
        dataframe: columns = name, row = minimum prediction attributes
    """
    # minimun median mt score - full line
    min_score_line = df[df["Median MT Score"] == df["Median MT Score"].min()]
    columns = ["HLA Allele", "MT Epitope Seq", "Median MT Score"]
    min_dict = {}
    for item in columns:
        min_dict[item] = min_score_line[item].values[0]
    output_df = pd.DataFrame(min_dict, index=[0])
    return output_df if len(df) != 0 else pd.DataFrame({columns[0]: "NA", columns[1]: "NA", columns[2]: "NA"}, index=[0])

def length_count(df, message):
    """Find the number of mut sequences that correspond to each peptide length

    Args:
        df (dataframe): pVACtools df
        message (str): descriptive prefix for resulting column name (ex: Neo for neoantigen)
    
    Return:
        df (dataframe): columns = message + peptide length, row = mut sequence count per peptide length
    """
    unique_seqs = df['MT Epitope Seq'].unique().tolist()
    len_seqs = [len(x) for x in unique_seqs]
    len_dict = {f'{message} Length {x}': len_seqs.count(x) for x in set(len_seqs)}
    len_df = pd.DataFrame(len_dict, index=[0])
    return len_df
        
input_dir=sys.argv[1]
output_filename=sys.argv[2]
working_dir=os.getcwd()
output_file=f'{working_dir}/{output_filename}'

open_files = [f for f in sorted(os.listdir(input_dir)) if f.endswith("_TUMOR.all_epitopes.tsv")]
#open_files = [f for f in sorted(os.listdir(f'{working_dir}/{input_dir}/')) if f.endswith("_TUMOR.all_epitopes.tsv")] 
for file1 in open_files:
    df_input = pd.read_csv(f'{input_dir}/{file1}', sep='\t')
    # constant per file
    df_constant = df_input[["Chromosome", "Start", "Stop", "Reference", "Variant", "Transcript", "Ensembl Gene ID", "Variant Type", "Mutation", "Protein Position", "Gene Name", "HGVSp"]]
    constant_line = df_constant.drop_duplicates()
    # unique for each line
    df_variable = df_input[["HLA Allele", "Peptide Length", "Sub-peptide Position", "Mutation Position", "MT Epitope Seq", "WT Epitope Seq", "Median MT Score", "Median WT Score", "Median Fold Change"]]
    df_variable = df_variable.drop_duplicates()

    file_info = pd.DataFrame({"File Name": file1},index=[0])
    constant_line = df_constant.drop_duplicates()
    
    neo_classI, neo_classII, neo_score_dict = identify_neoantigens(df_variable)
    
    total_pred = parse_predictions(df_variable, ["Predictions", "Total"])
    classI_pred = parse_predictions(neo_classI, ["Neo Class I", "ClassI"])
    classII_pred = parse_predictions(neo_classII, ["Neo Class II", "ClassII"])

    classI_min = find_min_prediction(neo_classI)
    classII_min = find_min_prediction(neo_classII)
    
    total_len = length_count(df_variable, "Seq")
    neo_len = length_count(pd.concat([neo_classI, neo_classII]), "Neo")

    final_df_list = [file_info, constant_line, total_pred, classI_pred, classII_pred, classI_min, classII_min, total_len, neo_len] 
    final_write_df = pd.DataFrame(final_df_list, index=[0])

    if os.path.isfile(output_file):
        final_write_df.to_csv(output_file, mode = 'a', index=False, header=None)
    else:
        final_write_df.to_csv(output_file, mode = 'w', index=False)

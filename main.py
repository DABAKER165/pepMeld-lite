#!/usr/bin/env python

import pandas as pd
import math
import os
from argparse import ArgumentParser, ArgumentTypeError


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise ArgumentTypeError('Boolean value expected.')


def parse_process_arrays_args(parser: ArgumentParser):
    """Parses the python script arguments from bash and makes sure files/inputs are valid"""
    parser.add_argument('--meta_filepath',
                        type=str,
                        help='filepath of the meta data',
                        required=True)
    parser.add_argument('--meta_sep',
                        type=str,
                        default='\t',
                        help='filepath of the meta data',
                        required=False)
    parser.add_argument('--pep_seq_protein_lookup_path',
                        type=str,
                        help='Filepath of the protein lookup file',
                        required=True)
    parser.add_argument('--lookup_sep',
                        type=str,
                        default='\t',
                        help='separator used in the pep_seq_protein_lookup file [ , \t ]',
                        required=False)
    parser.add_argument('--in_dir',
                        type=str,
                        help='Location of the input data files',
                        required=True)
    parser.add_argument('--in_data_ext',
                        type=str,
                        default='Mamu.txt.gz',
                        help='File names  endswith to search for in data files of the in_dir',
                        required=False)
    parser.add_argument('--in_data_sep',
                        type=str,
                        default='\t',
                        help='separator used in the in_data file(s) [ , \t ]')
    parser.add_argument('--sequence_lengths',
                        type=str,
                        default='8,9,10',
                        help='comma delimited list of lengths to include [ 8,9,10 ]',
                        required=False)
    parser.add_argument('--out_dir',
                        type=str,
                        help='output directory',
                        required=True)
    parser.add_argument('--out_prefix',
                        type=str,
                        default='mhc_ppma',
                        help='prefix to use in the output files.',
                        required=False)


def get_process_arrays_args():
    """	Inputs arguments from bash
    Gets the arguments, checks requirements, returns a dictionary of arguments
    Return: args - Arguments as a dictionary
    """
    parser = ArgumentParser()
    parse_process_arrays_args(parser)
    args = parser.parse_args()
    return args


args = get_process_arrays_args()
# same arguments to a local variable by same name as the argument
meta_filepath = args.meta_filepath
meta_sep = args.meta_sep
pep_seq_protein_lookup_path = args.pep_seq_protein_lookup_path
lookup_sep = args.lookup_sep
in_dir = args.in_dir
in_data_ext = args.in_data_ext
in_data_sep = args.in_data_sep
sequence_lengths = args.sequence_lengths
out_dir = args.out_dir
out_prefix = args.out_prefix

# Open Lookup Table
seq_len_list = sequence_lengths.split(',')
seq_len_list = [int(x) for x in seq_len_list]
df_lookup = pd.read_csv(pep_seq_protein_lookup_path,
                        sep=lookup_sep,
                        usecols=['ORIGINAL_SEQ_ID', 'PEPTIDE_SEQUENCE'])
df_lookup.rename(columns={'PEPTIDE_SEQUENCE': 'PROBE_SEQUENCE', 'ORIGINAL_SEQ_ID': 'SEQ_ID'},
                 inplace=True)
# Remove sequences that are in the length list length
df_lookup = df_lookup[df_lookup['PROBE_SEQUENCE'].str.len().isin(seq_len_list)]
df_lookup.drop_duplicates(subset=['PROBE_SEQUENCE'], keep='last', inplace=True)
#############################################################################
# get the data paths in the folder, excluding resource fork files and open #
#############################################################################
data_paths = os.listdir(in_dir)
resrc_fork = '._'
data_paths = [os.path.join(in_dir, x) for x in data_paths if x.endswith(in_data_ext) and not x.startswith(resrc_fork)]
# convert to integer sequence length list from a comma delimited string.

df_data = pd.DataFrame()

for path_i in data_paths:
    df_data_i = pd.read_csv(path_i,
                            sep=in_data_sep)
    # Remove unwanted columns to save space.
    df_data_i.drop(columns=['POSITION', 'PROBE_ID'], inplace=True)
    # Get Peptide Sequence length and filter on it.
    df_data_i['PEP_LEN'] = [len(x) for x in df_data_i['PROBE_SEQUENCE']]
    df_data_i = df_data_i[df_data_i['PEP_LEN'].isin(seq_len_list)]
    # make sure there is data and merge as needed.
    if len(df_data.index) > 0:
        df_data = df_data.merge(df_data_i, on=['PROBE_SEQUENCE', 'PEP_LEN', 'REPL', 'SEQ_ID'], how='inner')
    else:
        df_data = df_data_i.copy()

######################
# open the meta data #
######################
df_meta = pd.read_csv(meta_filepath, sep=meta_sep)
df_meta = df_meta[df_meta['EXCLUDE'] != 'EXCLUDE']
# Create a rename the vendor names to a human readable name
rename_dict = {}
for key_i, val_i in zip(df_meta['VENDOR_NAME'], df_meta['MHC Allele Name']):
    rename_dict[key_i] = val_i
df_data.rename(columns=rename_dict, inplace=True, errors='ignore')
# create a list of actual data columns
df_data_cols = list(set(df_meta['MHC Allele Name']))

# Take the log base 2 of the columns
for col_i in df_data_cols:
    df_data[col_i] = [math.log(x, 2) for x in df_data[col_i]]

# Take the Median by Probe sequences for each data column.
df_data_med = df_data.groupby(["PROBE_SEQUENCE", 'PEP_LEN', 'SEQ_ID'])[df_data_cols].median().reset_index()

# Seperate the _Redundant and Named sequence as they are missing from the correspondence key
df_data_med_no_seq = df_data_med[df_data_med['SEQ_ID'].str.endswith('_REDUNDANT')]
# drop the SEQ_ID column prior to merge so the naming doesn't get messed up.
df_data_med_no_seq.drop(columns=['SEQ_ID'], inplace=True)
df_data_med_seq = df_data_med[~df_data_med['SEQ_ID'].str.endswith('_REDUNDANT')]

# Join only the ones that don't have a sequence key
# outer merge/join on the probe sequence with the Peptide Sequence to Protein Accession Number lookup
df_data_m = df_data_med_no_seq.merge(df_lookup, on='PROBE_SEQUENCE', how='inner')
# Concat the Original SEQ ID with the REDUNDANT now Named ones.
df_data_m = pd.concat([df_data_m, df_data_med_seq], ignore_index=True)

# get column list
df_data_m_cols = list(df_data_m.columns)
# Round the values as you control subtract
for key_i, val_i in zip(df_meta['MHC Allele Name'], df_meta['CONTROL_SUBTRACT']):
    if (key_i in df_data_m_cols) and (val_i in df_data_m_cols):
        df_data_m[key_i] = [round(x - y, 3) for x, y in zip(df_data_m[key_i], df_data_m[val_i])]


df_data_m.to_csv(os.path.join(out_dir, '{0}_control_subtracted.csv.gz'.format(out_prefix)),
                 index=False,
                 compression="gzip")

#################################
# Add IEDB specific information #
#################################

# Drop the duplicate sequences (if they exist as they shouldn't at this point.)
df_data_m.drop_duplicates(subset=['PROBE_SEQUENCE'], keep='last', inplace=True)

df_data_m['SEQ_ID'] = df_data_m['SEQ_ID'].fillna('_')
df_data_m['Epitope Source Protein GenBank ID'] = [x.split('_')[-1] for x in df_data_m['SEQ_ID']]


# Get the control columns and drop them
control_columns = list(df_meta['CONTROL_SUBTRACT'].unique())
for control_column_i in control_columns:
    df_data_m.drop(columns=control_column_i, inplace=True)
    if control_column_i in df_data_cols:
        df_data_cols.remove(control_column_i)

df_data_m.rename(columns={'PROBE_SEQUENCE': 'Epitope Linear Sequence'}, inplace=True)
# get lookup table to protein
df_data_m2 = pd.melt(df_data_m,
                     id_vars=['Epitope Linear Sequence', 'Epitope Source Protein GenBank ID', 'SEQ_ID'],
                     value_vars=df_data_cols,
                     var_name='MHC Allele Name',
                     value_name='Quantitative Measurement')

# Get the thresholds from the metatable for a qualitative measurment
threshold_dict = {}
for key_i, val_i in zip(df_meta['MHC Allele Name'], df_meta['PASSING_THRESHOLD']):
    threshold_dict[key_i] = val_i

df_data_m2['Epitope Name'] = ''
df_data_m2['Epitope Source Organism Taxonomy ID'] = ''
df_data_m2['Peptide Nature'] = 'Peptide from protein'
df_data_m2['Epitope Evidence Code'] = 'Epitope Evidence Code'
df_data_m2['Epitope Structure Defines'] = 'Exact Epitope'
df_data_m2['Epitope Comments'] = ''
df_data_m2['Assay Type'] = 'MHC–peptide binding microarray assay - Fluorescence'
df_data_m2['Assay Response'] = 'Association (or direct binding)'
df_data_m2['Assay Units'] = 'A.U.'
df_data_m2['Qualitative Measurement'] = ['POSITIVE' if x > threshold_dict[y] else 'NEGATIVE' for x, y in zip(df_data_m2['Quantitative Measurement'], df_data_m2['MHC Allele Name'])]
df_data_m2['Measurement Inequality'] = '>'
df_data_m2['Assay Comments'] = 'This assay was performed with 5x replicates with median value shown'


df_data_line1 = pd.DataFrame({'Data Field': ['Required/Optional'],
                              'Epitope Name': ['Optional'],
                              'Epitope Linear Sequence': ['Required'],
                              'Peptide Nature': ['Required'],
                              'Epitope Source Protein GenBank ID': ['Required'],
                              'Epitope Source Organism Taxonomy ID': ['Optional'],
                              'Epitope Evidence Code': ['Required'],
                              'Epitope Structure Defines': ['Required'],
                              'Epitope Comments': ['Optional'],
                              'Assay Type': ['Required'],
                              'Assay Response': ['Required'],
                              'Assay Units': ['Required'],
                              'Qualitative Measurement': ['Required'],
                              'Quantitative Measurement': ['Optional'],
                              'Measurement Inequality': ['Optional'],
                              'MHC Allele Name': ['Required'],
                              'Assay Comments': ['Optional']})
df_data_m2 = pd.concat([df_data_line1, df_data_m2], ignore_index=True)
df_data_m2.drop(columns=['SEQ_ID'], inplace=True)
df_data_m2.to_csv(os.path.join(out_dir, '{0}.txt.gz'.format(out_prefix)),
                  sep='\t',
                  index=False,
                  compression='gzip')
print('Complete')

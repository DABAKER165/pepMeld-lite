# pepMeld-Lite
- Simpler workflow to process Roche-NimbleGen arrays.
- This does not complete the additional processing if the phyical coordinates of the array is present.
  - Use pepMeld for better preprocessing if the physical coordinates of the array are exposed
# Summary
- This python script will:
  - Open the data file, 
    - Take the log base2 of each binding affinity value
    - Subtract the associate control with each value
    - Take the median of the repeated peptide sequences of each array, by MHC type.
  - Open a "meta data" file to 
    - rename sequences to a human readable format. 
    - Apply a qualitative positive/negative "Threshold" for each MHC-Type based on  the distribution.
  - Open the peptide sequence to protein accession number lookup
    - Choose a unique accession number / peptide sequence
    - Join the sequences to the peptide data, for each MHC type.
  - Output a file ready to upload into iedb.org epitope databases.
  
## Running main.py to agregate the data for IEDB.org
### Creates IEDB.org ready files from Roche-NimbleGen Peptide Array datasets Arguments
| argument                    | single_character   | type   | default     | help                                                           | required   |
|:----------------------------|:-------------------|:-------|:------------|:---------------------------------------------------------------|:-----------|
| meta_filepath               | m                  | string | nan         | filepath of the meta data                                      | t          |
| meta_sep                    | e                  | string | \t          | filepath of the meta data                                      | f          |
| pep_seq_protein_lookup_path | p                  | string | nan         | Filepath of the protein lookup file                            | t          |
| lookup_sep                  | l                  | string | \t          | separator used in the pep_seq_protein_lookup file [ , \t ]     | f          |
| in_dir                      | i                  | string | nan         | Location of the input data files                               | t          |
| in_data_ext                 | x                  | string | Mamu.txt.gz | File names  endswith to search for in data files of the in_dir | f          |
| in_data_sep                 | d                  | string | \t          | separator used in the in_data file(s) [ , \t ]                 | nan        |
| sequence_lengths            | s                  | string | 8,9,10      | comma delimited list of lengths to include [ 8,9,10 ]          | f          |
| out_dir                     | o                  | string | nan         | output directory                                               | t          |
| out_prefix                  | f                  | string | mhc_ppma    | prefix to use in the output files.                             | f          |

### Example:
```bash
python3 ~/github/pepMeld-lite/main.py \
--meta_filepath='/Volumes/mhc_dataset/mhc_meta.tsv' \
--meta_sep='\t' \
--pep_seq_protein_lookup_path=/Volumes/mhc_dataset/SIV_SHIV_PEP_QX12_correspondence_key.txt.gz \
--lookup_sep='\t' \
--in_dir='/Volumes/mhc_dataset/input_data' \
--in_data_ext='Mamu.txt.gz' \
--in_data_sep='\t' \
--sequence_lengths='8,9,10' \
--out_dir='/Volumes/mhc_dataset/output_data' \
--out_prefix='mhc_ppma'
```
## Raw data input data format
Required columns:

| ﻿Column Name     |Details                |
|----------------|----------------------|
|PROBE_SEQUENCE| Peptide sequence (all caps)|
|REPL|The corresponding Position of a Replicate (used to match control position)|
|SEQ_ID|Either <Len>MER_REDUNDANT or <Virus Accession Number>_<Protein Accession Number>|
## Meta data input data format
Required columns:

| ﻿Column Name     |Details                |
|----------------|----------------------|
|VENDOR_NAME      |Data column name from the vendor provided data table |
|MHC Allele Name  |This is the MHC Allele Name that will be uploaded to IEDB.org|
|CONTROL_SUBTRACT |This is the column to subtract the values by Peptide_Sequence & Replicate position|
|EXCLUDE| To exclude the column from analysis|
|PASSING_THRESHOLD| Control Subtracted threshold that is considered passing|

## Peptide Sequence to Protein Accession format
Required columns:

| ﻿Column Name |Details|
|----------------|----------------------|
|PEPTIDE_SEQUENCE|Peptide sequence (all caps)|
|ORIGINAL_SEQ_ID| Associated Virus Accession Number _ Protein Accession Number|
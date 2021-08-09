"""
tutorial_proprocess_bulk_files.py 

This tutorial document illustrates pre-processing of bulk TCR beta 
repertoires files. These files were cleaned prior to use with tcrdist3
and tabulation of meta-clonotype conformant TCRs. 
Cleaning steps removes invalid CDR3. 
It also converts from Adaptive gene names to IMGT gene names needed 
to infer CDR1,2,2.5 used to compute pairwise TCRdist.

This data is available via Adaptive Biotechnologies:
https://immunerace.adaptivebiotech.com/data/ 
https://clients.adaptivebiotech.com/pub/covid-2020

The file ImmuneCODE-Repertoires-002.tgz contains bulk repertoires from COVID-19 patients at global sites.
(https://adaptivepublic.blob.core.windows.net/publishedproject-supplements/covid-2020/ImmuneCODE-Repertoires-002.tgz) 

It a large archive of bulk files (over 25 GB compressed).
Thus, one might want to first explore preprocessing on a single file. 
tar --extract --file=ImmuneCODE-Repertoires-002.tgz ImmuneCODE-Review-002/KH20-11642_TCRB.tsv

The mimimum file only contains cdr3_b_aa, v_b_gene, j_b_gene, and productive_frequency fields
The maximum file contains all of the field in a full clone_df.
"""
import os
from tcrdist.adpt_funcs import bulk_adaptive_dataset_to_tcrdist3_clone_df
source_path = 'ImmuneCODE-Review-002/'
dest_path = 'bulk_tcrdist3_ready/'
for file in [f for f in os.listdir(source_path) if f.endswith(".tsv")]:
    print(file)
    df = bulk_adaptive_dataset_to_tcrdist3_clone_df(bulk_filename=os.path.join(source_path, file) ,
        minimum_file=os.path.join(dest_path, f'{file}.tcrdist3.v_min.tsv'),
        maximum_file=os.path.join(dest_path, f'{file}.tcrdist3.v_max.tsv'),
        organism='human',
        chains=['beta'])
    print(df)



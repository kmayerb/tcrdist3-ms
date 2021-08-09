"""
tutorial_tabulate_metaclonotypes.py 

This document illustrates tabulation of meta-clonotype conformant TCRs across 694 
bulk samples. The tutorial use MIRA set 55. Prior to executing this script 
the bulk repertoires have been cleaned and prepared with IMGT names for use 
with tcrdist3 (see tutorial_preprocess_bulk.py).

To reproduce results in the manuscript we will use meta-clontype definition file:
'mira_epitope_55_524_ALRKVPTDNYITTY_KVPTDNYITTY.tcrdist3.csv.ranked_centers_bkgd_ctlr_1E6_manuscript.tsv'

1. To run this example the user must ensure all the bulk files are 
available and specifiy correct 'source_path' variable below. 

2. The outputed counts will be written to folder specified by 'project_path' folder. 
counts by file are written individually.

Two methods are illustrated:

One method uses parmap to run the tabulation in parallel. In this method 
ncpus in the tabulate_metaclonotype() function should be set to 1, 
and the parmap argument 'pm_processes' should be set to number of available 
cpus. This method is ideal when (i) tabulating conformant clones to 
a smaller set of meta-clonotypes (<500) as shown here, and (ii) also 
when a user has access to a computing node with multiple CPUs and >32 GB memory.

The second method is ideal for users running on a laptop 
who wish to perform the tabulation of conformant clones in each 
bulk repertoires in series.  The second method is also appropriate
when tabulating conformant clones against a large number of meta-clonotypes. 

For instance tabulating (>500-5000) meta-clonotypes against 10^6 bulk clones can be 
done efficiently via the second method, where ncpus can be set based on available 
cpu resources.
"""
import multiprocessing
import numpy as np
import os
import pandas as pd
from tcrdist.setup_tests import download_and_extract_zip_file
from tcrdist.repertoire import TCRrep
from tcrdist.breadth import get_safe_chunk
from tcrdist.join import join_by_dist
import re
import time

def tabulate_metaclonotype(
    file,
    metaclonotype_source_path,
    metaclonotype_file, 
    source_path, 
    ncpus =1,
    max_radius = 36,
    write = False,
    project_path = "counts"):
    """
    Tabulate a set of meta-clonotypes in a single bulk repertoires 
    
    Parameters 
    ----------
    metaclonotype_source_path : str
        filepath to metaclonotype file
    metaclonotype_file : str 
        filename containing meta-clonotype definitions
    source_path : str
         filepath to bulk files
    file: str
        filename
    ncpus = 6
        maximum number of cpus to use in meta-clonotype vs. bulk distance computation
    max_radius = 36
        maximum radius to store
    
    Returns
    -------
    df_join : pd.DataFrame
    
    """
    ncpus = min(multiprocessing.cpu_count(), ncpus)

    df_search = pd.read_csv(os.path.join(metaclonotype_source_path, metaclonotype_file), sep = "\t")
    
    df_bulk = pd.read_csv(os.path.join(source_path, file), sep = "\t")
    # When one want to track each clone indivually regardless of identical TRBV-CDR3-TRBJ
    df_bulk = df_bulk.sort_values('count').reset_index(drop = True)
    
    df_bulk['rank'] = df_bulk.index.to_list()
    
    from tcrdist.repertoire import TCRrep
    tr = TCRrep(
        cell_df = df_search,
        organism = "human",
        chains = ['beta'],
        compute_distances= False)
    tr.cpus = ncpus
    
    tr_bulk = TCRrep(
        cell_df = df_bulk,
        organism = "human",
        chains = ['beta'],
        compute_distances= False)
    
    chunk_size = get_safe_chunk(tr.clone_df.shape[0], tr_bulk.clone_df.shape[0])
    
    tr.compute_sparse_rect_distances(
        df = tr.clone_df,
        df2 = tr_bulk.clone_df,
        radius = max_radius ,
        chunk_size = chunk_size)
    
    df_join = join_by_dist(
        how = 'inner',
        csrmat = tr.rw_beta,
        left_df = tr.clone_df,
        right_df = tr_bulk.clone_df,
        left_cols  = tr.clone_df.columns.to_list(),
        right_cols = tr_bulk.clone_df.columns.to_list(),
        left_suffix = '_search',
        right_suffix = '_bulk',
        max_n= 1000,
        radius = max_radius )
    
    # df_join has more results 
    df_join['RADIUS'] = df_join.apply(lambda x: x['dist'] <= x['radius_search'], axis = 1)
    import re
    df_join['MOTIF'] = df_join.apply(lambda x: re.search(string = x['cdr3_b_aa_bulk'],
        pattern = x['regex_search']) is not None, axis = 1)
    df_join['RADIUSANDMOTIF'] =  df_join['RADIUS'] & df_join['MOTIF']
    df_join['EXACT'] = df_join.apply(lambda x: x['dist'] <= 0, axis = 1)
    df_join['unique_clones'] = 1
    
    df_join['feature'] = df_join['v_b_gene_search'] + "+" \
        + df_join['cdr3_b_aa_search'] + "+" \
        + df_join['radius_search'].apply(lambda x : str(x)) + "+"\
        + df_join['regex_search']
    
       # mc_file = 'mira_epitope_55_524_ALRKVPTDNYITTY_KVPTDNYITTY.tcrdist3.csv.ranked_centers_bkgd_ctlr_1E6_manuscript.tsv'
    df_mc = pd.read_csv(os.path.join(metaclonotype_source_path, metaclonotype_file), sep = "\t")
    df_mc['feature'] = df_mc['v_b_gene'] + "+" \
        + df_mc['cdr3_b_aa'] + "+" \
        + df_mc['radius'].apply(lambda x : str(x)) + "+" \
        + df_mc['regex']
        
    radius = df_join.query('RADIUS').groupby('feature')['templates_bulk'].sum()
    motif = df_join.query('RADIUSANDMOTIF').groupby('feature')['templates_bulk'].sum()
    exact = df_join.query('EXACT').groupby('feature')['templates_bulk'].sum()
    
    df_result = pd.concat([pd.DataFrame(index = df_mc['feature'] ), radius, motif, exact], axis = 1)
    
    df_result.columns = ['RADIUS','MOTIF','EXACT']
    df_result['file'] = file
    df_result = df_result.fillna(0)
    if write: 
        outname = os.path.join(project_path, f"{file}.counts.tsv")
        df_result.reset_index(drop = False).to_csv(outname, sep = "\t", index = False)
    
    return (df_join, df_result)
    



if __name__ == "__main__":
    """ < mc_file > contains meta-clonotype definitions"""
    mc_file = 'mira_epitope_55_524_ALRKVPTDNYITTY_KVPTDNYITTY.tcrdist3.csv.ranked_centers_bkgd_ctlr_1E6_manuscript.tsv'
    
    """ < project_path > is where the counts files will be written, 
    YOU WILL NEED TO SPECIFIY THIS FOR YOUR OWN ENVIRONMENT """    
    project_path = 'tutorial_tabulate_metaclonotypes_MIRA_55'

    """ < source_path > points to folder where bulk repertoires are contained 
    YOU WILL NEED TO SPECIFIY THIS FOR YOUR OWN ENVIRONMENT, BASED ON WHERE 
    YOU DEPOSITE THE CLEANED BULK DATAFILES
    """
    source_path = '../ncov_tcrs/bulk2/adpt_bulk_r2clean'
    
    """ < meta_data_file >
    Use a well organized metadata file that contains information about bulk filenames. Note
    that all referenced filenames must be present in folder specified by < source_path > variable above"""
    meta_data_file = "metadata/2020-10-12-metadata.tsv"

    df = pd.read_csv( meta_data_file, sep = "\t")
    """the local path can either be provided in your meta-data file or overwritten """
    df['local_path'] = source_path 
    
    """Check that all files are present"""
    assert np.all([os.path.isfile(os.path.join(source_path,f)) for f in df['filename']]), f"Not all files specified in meta-data file {meta_data_file} are present in {source_path}"
    
    """ Method 1: With more computational resources, you can greatly accelerate this process with parmap"""
    use_method_1 = True 
    if use_method_1:
        import parmap
        results = parmap.map(tabulate_metaclonotype, 
                df['filename'].to_list(),
                metaclonotype_source_path = '.',
                metaclonotype_file =  mc_file ,
                source_path =  source_path  ,
                ncpus = 1, #  <- Note: set this to 1 CPU if using parmap
                write = True,
                project_path = project_path, 
                pm_pbar=True, 
                pm_processes=10) # <- Note: set this to total available CPUs
                
        df_concat = pd.concat([x[1] for x in results]).reset_index()
        df_concat['name'] = df_concat['file'].apply(lambda x : x.replace(".tsv.tcrdist3.v_max.tsv",""))
        df_concat.to_csv(os.path.join(project_path, f"{mc_file}.tabulated_counts_concat.tsv"), sep = "\t", index = False)
            
    """ Method 2: Alternatively, if you are running with minimal resources you can tabulate conformant clones 
        in one bulk repertiore at a time. """
    use_method_2 = False 
    if use_method_2:
        for i,r in df.iterrows():
            df_join, df_result = tabulate_metaclonotype(
                metaclonotype_source_path = '.',
                metaclonotype_file =  mc_file ,
                source_path =  r['local_path'],
                file =  r['filename'],
                ncpus = 2, # <- Note: set this to total available CPUs
                write = False,
                project_path = project_path )
            outname = os.path.join(project_path, f"{r['file']}.counts.tsv")
            df_result.reset_index(drop = False).to_csv(outname, sep = "\t", index = False)



    

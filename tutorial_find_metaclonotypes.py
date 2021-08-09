"""
tutorial_find_metaclonotypes.py 

This functions encapsulates a complete work flow for finding meta-clonotypes in antigen-enriched data.
It uses MIRA Set 55 as an example. It also produced optional visual output in html files showing 
background substracted CDR3 logo motifs for each meta-clonotype.

To ensure that you have necessary tcrsampler files to generate synethetic backgrounds, you may first 
have to run:
python -c "from tcrsampler.setup_db import install_all_next_gen; install_all_next_gen(dry_run = False)"

See more details provided: https://tcrdist3.readthedocs.io/en/latest/metaclonotypes.html

The synthetic backgrounds used to estaimate optimal TCRdist radii 
for each antigen-associated TCR are randomly generated each run, 
so meta-clonotypes derived here may differ slightly from those previously defined. 
If you prefer, to reproduce the analysis for MIRA 55 as it appears in the manuscript 
(Mayer-Blackwell et al. 2020), you can proceed to directly tutorial_tabulate.py using already available file:

mira_epitope_55_524_ALRKVPTDNYITTY_KVPTDNYITTY.tcrdist3.csv.ranked_centers_bkgd_ctlr_1E6_manuscript.tsv.tabulated_counts_concat
"""
import sys
import os
import numpy as np
import pandas as pd
import scipy.sparse
from tcrdist.paths import path_to_base
from tcrdist.repertoire import TCRrep
from tcrsampler.sampler import TCRsampler
from tcrdist.background import get_stratified_gene_usage_frequency
from tcrdist.background import  sample_britanova
from tcrdist.ecdf import distance_ecdf, _plot_manuscript_ecdfs
import matplotlib.pyplot as plt
from tcrdist.neighbors import bkgd_cntl_nn2
from tcrdist.public import _neighbors_variable_radius, _neighbors_sparse_variable_radius
from tcrdist.centers import rank_centers
from tcrdist.background import make_gene_usage_counter, get_gene_frequencies
from progress.bar import IncrementalBar
from tcrdist.public import make_motif_logo
from tcrdist.automate import auto_pgen


def find_metaclonotypes(
    project_path = "tutorial_find_metaclonotypes_MIRA_55",
    source_path = os.path.join(path_to_base,'tcrdist','data','covid19'),
    antigen_enriched_file = 'mira_epitope_55_524_ALRKVPTDNYITTY_KVPTDNYITTY.tcrdist3.csv',
    ncpus = 4, 
    seed = 3434):
    # CONFIGURE OUTPUT PATH
    np.random.seed(seed)
    if not os.path.isdir(project_path):
        os.mkdir(project_path)
    assert os.path.isfile(os.path.join(source_path, antigen_enriched_file))
    # LOAD ANTIGEN ENRICHED DATA
    df = pd.read_csv(os.path.join(source_path, antigen_enriched_file))
    df = df[ (df['v_b_gene'].notna() ) & (df['j_b_gene'].notna()) ]
    # REMOVE V GENES MISSING IN CORD BLOOD DATASET, BECAUSE 
    # WITHOUT ESTIMATE OF PREVALENCE, THESE WOULD LEAD TO POTENTIALLY NON-EPITOPE SPECIFIC 
    # META-CLONOTYPES
    df = df[(df['v_b_gene'] != 'TRBV12-1*01') & (df['v_b_gene'] != 'TRBV6-2*01')].reset_index(drop = True)
    tr = TCRrep(cell_df = df[['subject','cell_type','v_b_gene', 'j_b_gene', 'cdr3_b_aa']], 
                organism = "human", 
                chains = ['beta'], 
                compute_distances = True)
    # COMPUTE PGENs
    auto_pgen(tr)
    tr.cpus = ncpus
    # MAKE SYNTHETIC BACKGROUND 
    ts = TCRsampler(default_background = 'britanova_human_beta_t_cb.tsv.sampler.tsv')
    ts = get_stratified_gene_usage_frequency(ts = ts, replace = True) 
    df_vj_background = tr.synthesize_vj_matched_background(ts = ts, chain = 'beta')
    df_britanova_100K = sample_britanova(size = 100000)
    df_britanova_100K = get_gene_frequencies(ts = ts, df = df_britanova_100K)
    df_britanova_100K['weights'] = 1
    df_britanova_100K['source'] = "stratified_random"
    df_bkgd = pd.concat([df_vj_background.copy(), df_britanova_100K.copy()], axis = 0).\
        reset_index(drop = True)                                              
    assert df_bkgd.shape[0] == 200000
    background_outfile = os.path.join(project_path, f"{antigen_enriched_file}.olga100K_brit100K_bkgd.csv")
    df_bkgd.to_csv(background_outfile, index = False)
    tr_bkgd = TCRrep(
        cell_df = df_bkgd,
        organism = "human", 
        chains = ['beta'], 
        compute_distances = False)
    # COMPUTE DISTANCE BETWEEEN ANTIGEN ENRICHED DATA AND BACKGROUND DATA
    tr.compute_sparse_rect_distances(
        df = tr.clone_df, 
        df2 = tr_bkgd.clone_df,
        radius=50,
        chunk_size = 100)
    scipy.sparse.save_npz(os.path.join(project_path, f"{antigen_enriched_file}.rw_beta.npz"), tr.rw_beta)
    # FIND METACLONOTYPES
    level_tag = '1E6'
    centers_df  = bkgd_cntl_nn2(
        tr               = tr,
        tr_background    = tr_bkgd,
        weights          = tr_bkgd.clone_df.weights,
        ctrl_bkgd        = 10**-6, 
        col              = 'cdr3_b_aa',
        add_cols         = ['v_b_gene', 'j_b_gene'],
        ncpus            = 4,
        include_seq_info = True,
        thresholds       = [x for x in range(0,50,2)],
        generate_regex   = True,
        test_regex       = True,
        forced_max_radius = 36)
    centers_df['neighbors'] = _neighbors_variable_radius(pwmat=tr.pw_beta, radius_list = centers_df['radius'])
    centers_df['K_neighbors'] = centers_df['neighbors'].apply(lambda x : len(x))
    # We determine how many <nsubjects> are in the set of neighbors 
    centers_df['nsubject']  = centers_df['neighbors'].\
            apply(lambda x: tr.clone_df['subject'].iloc[x].nunique())
    centers_outfile = os.path.join(project_path, f'{antigen_enriched_file}.centers_bkgd_ctlr_{level_tag}.tsv')
    centers_df.to_csv(centers_outfile  , sep = "\t" )
    # KEEP ONLY THE NON-REDUNDANT CLONES
    ranked_centers_df = rank_centers(
        centers_df = centers_df, 
        rank_column = 'chi2joint', 
        min_nsubject = 2, 
        min_nr = 1)
    ranked_centers_outfile = os.path.join(project_path, f'{antigen_enriched_file}.ranked_centers_bkgd_ctlr_{level_tag}.tsv')
    ranked_centers_df.to_csv(ranked_centers_outfile  , sep = "\t" )
    
    # OPTIONAL MAKE A VISUAL OUTPUT. HTML OUTPUT
    if ranked_centers_df.shape[0] > 0:
        cdr3_name = 'cdr3_b_aa'
        v_gene_name = 'v_b_gene'
        svgs = list()
        svgs_raw = list()
        bar = IncrementalBar('Processing', max = ranked_centers_df.shape[0])
        for i,r in ranked_centers_df.iterrows():
            bar.next()
            centroid = r[cdr3_name]
            v_gene   = r[v_gene_name]
            svg, svg_raw = make_motif_logo( tcrsampler = ts, 
                                            pwmat = tr.pw_beta,
                                            clone_df = tr.clone_df,
                                            centroid = centroid ,
                                            v_gene = v_gene ,
                                            radius = r['radius'],
                                            pwmat_str = 'pw_beta',
                                            cdr3_name = 'cdr3_b_aa',
                                            v_name = 'v_b_gene',
                                            gene_names = ['v_b_gene','j_b_gene'])
            svgs.append(svg)
            svgs_raw.append(svg_raw)
        bar.next();bar.finish()
        ranked_centers_df['svg']      = svgs
        ranked_centers_df['svg_raw'] = svgs_raw
        # WRITE OUTPUT TO AN HTML FILE
        def shrink(s):
            return s.replace('height="100%"', 'height="20%"').replace('width="100%"', 'width="20%"')
        
        labels =['cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'radius', 
        'regex','pgen','nsubject','K_neighbors', 'bkgd_hits_weighted',
        'chi2dist','chi2re','chi2joint']
        
        output_html_name = os.path.join(project_path, f'{antigen_enriched_file}.ranked_centers_bkgd_ctlr_{level_tag}.html')

        with open(output_html_name, 'w') as output_handle:
            for i,r in ranked_centers_df.iterrows():
                #import pdb; pdb.set_trace()
                svg, svg_raw = r['svg'],r['svg_raw']
                output_handle.write("<br></br>")
                output_handle.write(shrink(svg))
                output_handle.write(shrink(svg_raw))
                output_handle.write("<br></br>")
                output_handle.write(pd.DataFrame(r[labels]).transpose().to_html())
                output_handle.write("<br></br>")
                
                
if __name__ == "__main__":
    find_metaclonotypes(
        project_path = "tutorial_find_metaclonotypes_MIRA_55",
        source_path = os.path.join(path_to_base,'tcrdist','data','covid19'),
        antigen_enriched_file = 'mira_epitope_55_524_ALRKVPTDNYITTY_KVPTDNYITTY.tcrdist3.csv',
        ncpus = 4, 
        seed = 3434)
        
    
    

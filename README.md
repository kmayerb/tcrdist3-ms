# tcrdist3-ms

TCR meta-clonotypes with tcrdist3

### Background 

In Mayer-Blacwell et al. 2020, we tested whether meta-clonotypes carry important antigen-specific signals above and beyond individual clonotypes. Do do so, we searched for meta-clonotype conformant TCRS in COVID-19 patients with repertoires collected 0-30 days after diagnosis (Nolan et al. 2020). 

The MIRA55:ORF1ab set of antigen-associated TCRs was chosen for detailed analysis because, among the MIRA sets, it is comprised of CD8+ TCR β-chains activated by a peptide with the strongest evidence of HLA-restriction, primarily by HLA-A*01. We reasoned that we could compare the abundance of meta-clonotype conformant sequences in independent COVID-19 patients with and without the restricting HLA genotype, and that a significant positive association of abundance with the restricting genotype would provide confirmatory evidence of the meta-clonotype’s SARS-CoV-2 antigen specificity. 

### tcrdist3 documentation

The Python package tcrdist3 has extensive additional documentation:
[https://tcrdist3.readthedocs.io/](https://tcrdist3.readthedocs.io/en/latest/)


### Contents 

Tutorials using tcrdist3 with meta-clonotypes (discovery, tabulation, and beta-binomial counts regression)
These executable tutorials cover the primary computational methods used in Mayer-Blackwell et al. 2020.  

1. `tutorial_find_metaclonotypes.py` - define meta-clonotypes from antigen-associated MIRA TCR β-chain data. 

2. `tutorial_preprocess_bulk_files` - clean bulk TCR β-chains repertoire files from [ImmuneCODE-Repertoires-002.tgz](https://clients.adaptivebiotech.com/pub/covid-2020) for use with tcrdist3.

3. `tutorial_tabulate_metaclonotypes.py` - tabulate metaclonotype conformant TCRs sequences in many bulk repertoires

4. `tutorial_regressions` - beta-binomial regression with corncob (Martin et al. 2020) to find coefficient estimates associated with AGE, SEX, DAYS POST DIAGNOSIS, HLA-A*01


### Data Availability 

[ImmuneRACE data](https://clients.adaptivebiotech.com/pub/covid-2020) is publicly available 

Mayer-Blackwell et al. used data found in the following files: 
* MIRA antigen-associated TCRs : ImmuneCODE-MIRA-Release002.zip - 3.0 MB
* Bulk repertoires: ImmuneCODE-Repertoires-002.tgz - 24 GB

### Acknowledgments

This work was funded by NIH NIAID R01 AI136514-03 (PI Thomas) and ALSAC at St. Jude. The authors thank M. Pogorelyy and A. Minervina for extensive feedback on the manuscript. Scientific Computing Infrastructure at Fred Hutchinson Cancer Research Center was funded by ORIP grant S10OD028685.

### References 

Mayer-Blackwell, K., Schattgen, S., Cohen-Lavi, L., Crawford, J. C., Souquette, A., Gaevert, J. A., ... & Fiore-Gartland, A. (2020). TCR meta-clonotypes for biomarker discovery with tcrdist3: quantification of public, HLA-restricted TCR biomarkers of SARS-CoV-2 infection. bioRxiv.

Martin, B. D., Witten, D., & Willis, A. D. (2020). Modeling microbial abundances and dysbiosis with beta-binomial regression. The annals of applied statistics, 14(1), 94.

Nolan, S., Vignali, M., Klinger, M., Dines, J. N., Kaplan, I. M., Svejnoha, E., ... & Robins, H. S. (2020). A large-scale database of T-cell receptor beta (TCRβ) sequences and binding associations from natural and synthetic exposure to SARS-CoV-2. Research square.


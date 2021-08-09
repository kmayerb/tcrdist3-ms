# tcrdist3-ms

TCR meta-clonotypes with tcrdist3

### Background 

In a manuscript by Mayer-Blackwell et al., 2020, we tested whether meta-clonotypes carry important antigen-specific signals above and beyond individual clonotypes. To do so, we searched for meta-clonotype conformant TCRs in COVID-19 patients with repertoires collected 0-30 days after diagnosis (Nolan et al. 2020). The MIRA55 ORF1ab set of antigen-associated TCRs was chosen for detailed illustrative analysis because, among the MIRA sets, it is comprised of CD8+ TCR β-chains activated by a peptide with the strongest evidence of HLA-restriction, primarily by HLA-A*01. We reasoned that we could compare the abundance of meta-clonotype conformant sequences in an independent cohort of COVID-19 patients with and without the restricting HLA genotype and that a significant positive association of abundance with the restricting genotype would provide confirmatory evidence of the meta-clonotypes' SARS-CoV-2 antigen specificity. 

### tcrdist3 Documentation

The open-source Python package [tcrdist3](https://github.com/kmayerb/tcrdist3) has extensive additional documentation:
[https://tcrdist3.readthedocs.io/](https://tcrdist3.readthedocs.io/en/latest/)


### Contents of this Repository

These executable tutorials cover the primary computational methods used in Mayer-Blackwell et al. 2020, 
using tcrdist3 for meta-clonotypes discovery, tabulation, and testing in beta-binomial counts regressions.

1. `tutorial_find_metaclonotypes.py` - define meta-clonotypes from antigen-associated MIRA TCR β-chain data. 

2. `tutorial_preprocess_bulk_files.py` - clean bulk TCR β-chains repertoire files from [ImmuneCODE-Repertoires-002.tgz](https://clients.adaptivebiotech.com/pub/covid-2020) for use with tcrdist3.

3. `tutorial_tabulate_metaclonotypes.py` - tabulate meta-clonotype conformant TCRs sequences in many bulk repertoires

4. `tutorial_regressions.R` - beta-binomial regression with corncob (Martin et al. 2020) to find coefficient estimates associated with AGE, SEX, DAYS POST DIAGNOSIS, HLA-A*01

#### Files Generated

```
Input: mira_epitope_55_524_ALRKVPTDNYITTY_KVPTDNYITTY.tcrdist3.cs
Outputs: 
├── tutorial_find_metaclonotypes_MIRA_55
    (mira_epitope_55_524_ALRKVPTDNYITTY_KVPTDNYITTY.tcrdist3.csv)
│   ├── (all antigen-associated clonotypes) .centers_bkgd_ctlr_1E6.tsv
│   ├── (synethetic background) .olga100K_brit100K_bkgd.csv.zip
│   ├── (non-redundant meta-clonotypes) .ranked_centers_bkgd_ctlr_1E6.tsv
│   ├── (html logos of meta-clonotypes) .ranked_centers_bkgd_ctlr_1E6.html
│   └── (timing to tabulate in each bulk repertoire) .ranked_centers_bkgd_ctlr_1E6.tsv.benchmark_tabulation.tsv


Input: mira_epitope_55_524_ALRKVPTDNYITTY_KVPTDNYITTY.tcrdist3.csv.ranked_centers_bkgd_ctlr_1E6.tsv
Outputs:
├── tutorial_tabulate_metaclonotypes_MIRA_55
│   ├── 110047437_TCRB.tsv.tcrdist3.v_max.tsv.counts.tsv
│   ├── 110047542_TCRB.tsv.tcrdist3.v_max.tsv.counts.tsv
    ... 
│   ├── KHBR20-00206_TCRB.tsv.tcrdist3.v_max.tsv.counts.tsv
│   └──(all counts concatenated) .ranked_centers_bkgd_ctlr_1E6_manuscript.tsv.tabulated_counts_concat.tsv

Input: .ranked_centers_bkgd_ctlr_1E6_manuscript.tsv.tabulated_counts_concat.tsv
Outputs:
├── tutorial_regression_MIRA_55 
    (mira_epitope_55_524_ALRKVPTDNYITTY_KVPTDNYITTY.tcrdist3.csv.ranked_centers_bkgd_ctlr_1E6_manuscript.tsv)
│   ├── (volcano plots pdf) .tabulated_counts_concat.tsv.beta-binomial-regression.pdf
│   ├── (volcano plots png) .tabulated_counts_concat.tsv.beta-binomial-regression.png
│   └── (tabular results) .tabulated_counts_concat.tsv.beta-binomial-regression.tsv
```

#### Regression Results as Volcano Plots

![mira_epitope_55_524_ALRKVPTDNYITTY_KVPTDNYITTY tcrdist3 csv ranked png_centers_bkgd_ctlr_1E6_manuscript tsv tabulated_counts_concat tsv beta-binomial-regression](https://user-images.githubusercontent.com/46639063/128776168-dc75f23d-e6a9-44d4-845f-97baac478499.png)


### Data Availability 

[ImmuneRACE data](https://clients.adaptivebiotech.com/pub/covid-2020) is publicly available 

Mayer-Blackwell et al. used data found in the following files: 
* MIRA antigen-associated TCRs : ImmuneCODE-MIRA-Release002.zip - 3.0 MB
* Bulk repertoires: ImmuneCODE-Repertoires-002.tgz - 24 GB


### References 

Mayer-Blackwell, K., Schattgen, S., Cohen-Lavi, L., Crawford, J. C., Souquette, A., Gaevert, J. A., ... & Fiore-Gartland, A. (2020). TCR meta-clonotypes for biomarker discovery with tcrdist3: quantification of public, HLA-restricted TCR biomarkers of SARS-CoV-2 infection. bioRxiv.

Martin, B. D., Witten, D., & Willis, A. D. (2020). Modeling microbial abundances and dysbiosis with beta-binomial regression. The annals of applied statistics, 14(1), 94.

Nolan, S., Vignali, M., Klinger, M., Dines, J. N., Kaplan, I. M., Svejnoha, E., ... & Robins, H. S. (2020). A large-scale database of T-cell receptor beta (TCRβ) sequences and binding associations from natural and synthetic exposure to SARS-CoV-2. Research square.



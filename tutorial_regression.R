# tutorial_regression.R 

# This is a tutorial document illustrating beta-binomial regression 
# coefficient estimates associated with AGE, SEX, DAYS POST DIAGNOSIS, HLA-MATCH
# variables using the corncob::bbdml function in R.

# The input is count data of meta-clonotype conformant sequences per sample.
metaclonotype_definition_file = 'mira_epitope_55_524_ALRKVPTDNYITTY_KVPTDNYITTY.tcrdist3.csv.ranked_centers_bkgd_ctlr_1E6_manuscript.tsv'
tabulated_counts_file = 'tutorial_tabulate_metaclonotypes_MIRA_55/mira_epitope_55_524_ALRKVPTDNYITTY_KVPTDNYITTY.tcrdist3.csv.ranked_centers_bkgd_ctlr_1E6_manuscript.tsv.tabulated_counts_concat.tsv'
metadata_file = 'metadata/2020-10-12-metadata.tsv'
hmatch1 = "A\\*01"
project_path = "tutorial_regression_MIRA_55"

# In this example, with MIRA 55, we test the hypothesis that individuals 
# beleived to express HLA-A*01 (see github.com/kmayerb/hla3). 
# have a greater number of meta-clonotype conformant clones 
# than individuals without this HLA-A allele. 

# Meta-clonotypes were derived from MIRA identified TCRS (Nolan et al. 2020)
# acivated by SARS-CoV-2 peptides or sub-peptides within ORF1ab 1316-1330 LRKVPTDNYITTY.
# There was prior evidence that presentation of this peptide is restricted by
# antigen presentation with HLA-A*01.

# The model formulation is 
# 'cbind(W, M - W) ~ AGE+SEX+DAYS+HLA')
# Where AGE is age in years
# Where biological SEX at birth is binary with MALE as reference category
# Where DAYS is 1 if > 2 days post diagnois, 0 otherwise
# Where HLA is 1 if participant was inferred to express HLA-A*01 and 0 otherwise
# Where W is total conformant templates per sample and M is total productive templates.
# where we replace W with RADIUS, MOTIF, or EXACT in regressions below testing 
# association based on each definition. 

require(readr)
require(dplyr)
require(corncob)
require(ggplot2)

########################################################################################################################
### Beta-Binomial Maximum Likelihood Functions
########################################################################################################################
#### FUNCTIONS ####
#' do_corncob 
#' 
#' Define the beta-binomial we are attempting to fit
#' 
#' @param mydata data.frame
do_corncob <- function(mydata, frm = as.formula('cbind(W, M - W) ~ AGE+SEX+DAYS+HLA')){
  cb1 = bbdml(formula = frm,
              phi.formula = ~ 1,
              data = mydata)
  return(cb1)
}
possibly_do_corncob = purrr::possibly(do_corncob, otherwise = NA)
#' assemble bbdml coefficient estimates into a table
#' 
#' @param cb is object result of corncob::bbdml
parse_corncob <- function(cb,i =1){
  y = summary(cb)$coefficients
  rdf = as.data.frame(y) 
  rdf$param = rownames(rdf)
  rdf = rdf %>% mutate(estimate = Estimate,  se = `Std. Error`, tvalue =  `t value`, pvalue = `Pr(>|t|)`, param) %>% 
    mutate(type = ifelse(grepl(param, pattern = "phi"), "phi", "mu")) %>% 
    mutate(type2 = ifelse(grepl(param, pattern = "Intercept"), "intercept", "covariate")) 
  rdf$feature = i
  return(rdf)
}

#' lambda for adding feature column to a dataframe using - applied with purrr::map2(dataframes, names)
f  <- function(df, name){
  df = as.data.frame(df)
  df$variable = rownames(df)
  df['feature'] = name
  return(df)
}


# Load counts of RADIUS, MOTIF, EXACT conformant clones per sample
d = readr::read_tsv(tabulated_counts_file)
meta = readr::read_tsv(metadata_file) %>% select(name, AGE, SEX, DAYS, hla_a1_name, hla_a2_name, M = productive_templates)

# <df_ready> prepare a dataframe with correct HLA and DAYS variables
df_ready = d %>% left_join(meta, by = c("name" = "name")) %>% 
  mutate(hla_a1_name = ifelse(is.na(hla_a1_name), "XXX", hla_a1_name)) %>% 
  mutate(hla_a2_name = ifelse(is.na(hla_a2_name), "XXX", hla_a2_name)) %>%
  mutate(HLA = ifelse(  stringr::str_detect(hla_a1_name, pattern = hmatch1) | 
                              stringr::str_detect(hla_a2_name, pattern = hmatch1) , 1,0)) %>% 
  dplyr::mutate(DAYS  = ifelse( DAYS > 2, 1, 0))

# <dfs> split the data.frame by each feature, so regression is can be computed one feature at a time.
dfs = df_ready %>% split(f = .$index)
# <cb_res_motif> for corncob results, contains the fit BBDML model specified by do_corncob, fit to each feature seperately
cb_res_motif = purrr::map(dfs, ~possibly_do_corncob(.x, frm = as.formula('cbind(MOTIF, M - MOTIF) ~ AGE+SEX+DAYS+HLA') ))
cb_res_table_motif = do.call(rbind, purrr::map2(cb_res_motif, names(cb_res_motif), ~parse_corncob(cb = .x, i = .y))) %>% 
  dplyr::arrange(type,type2, pvalue) %>% 
  dplyr::mutate(`Pr(>|z|)` = `Pr(>|t|)`) %>%
  dplyr::mutate(`z value` = `t value`) %>%
  dplyr::filter(type == "mu", type2 == "covariate") %>% 
  dplyr::mutate(negative_log10_pvalue= -1 * log10(`Pr(>|z|)`)) %>%
  dplyr::mutate(variable = sub(param, pattern = "mu.", replacement = "") ) %>%
  dplyr::mutate(primary = ifelse(variable %in% c("M", "(Intercept)"), 0,1)) %>%
  dplyr::mutate(model = "TCRDIST BETA-BINOMIAL") %>% 
  dplyr::mutate(hmatch1 = hmatch1) %>% 
  dplyr::mutate(type = "MOTIF") %>%
  dplyr::mutate(form = 'cbind(MOTIF, M - MOTIF) ~ AGE+SEX+DAYS+HLA') %>% 
  dplyr::select(c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "variable", 
                  "feature", "primary", "negative_log10_pvalue", "hmatch1", "type", "form"))

# <cb_res_radius> 
cb_res_radius = purrr::map(dfs, ~possibly_do_corncob(.x, frm = as.formula('cbind(RADIUS, M - RADIUS) ~ AGE+SEX+DAYS+HLA') ))
cb_res_table_radius = do.call(rbind, purrr::map2(cb_res_radius, names(cb_res_radius), ~parse_corncob(cb = .x, i = .y))) %>% 
  dplyr::arrange(type,type2, pvalue) %>% 
  dplyr::mutate(`Pr(>|z|)` = `Pr(>|t|)`) %>%
  dplyr::mutate(`z value` = `t value`) %>%
  dplyr::filter(type == "mu", type2 == "covariate") %>% 
  dplyr::mutate(negative_log10_pvalue= -1 * log10(`Pr(>|z|)`)) %>%
  dplyr::mutate(variable = sub(param, pattern = "mu.", replacement = "") ) %>%
  dplyr::mutate(primary = ifelse(variable %in% c("M", "(Intercept)"), 0,1)) %>%
  dplyr::mutate(model = "TCRDIST BETA-BINOMIAL") %>% 
  dplyr::mutate(hmatch1 = hmatch1) %>% 
  dplyr::mutate(type = "RADIUS") %>%
  dplyr::mutate(form = 'cbind(MOTIF, M - MOTIF) ~ AGE+SEX+DAYS+HLA') %>% 
  dplyr::select(c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "variable", 
                  "feature", "primary", "negative_log10_pvalue", "hmatch1", "type", "form"))

# <cb_res_exact> 
cb_res_exact = purrr::map(dfs, ~possibly_do_corncob(.x, frm = as.formula('cbind(EXACT, M - EXACT) ~ AGE+SEX+DAYS+HLA') ))
cb_res_table_exact = do.call(rbind, purrr::map2(cb_res_exact, names(cb_res_exact), ~parse_corncob(cb = .x, i = .y))) %>% 
  dplyr::arrange(type,type2, pvalue) %>% 
  dplyr::mutate(`Pr(>|z|)` = `Pr(>|t|)`) %>%
  dplyr::mutate(`z value` = `t value`) %>%
  dplyr::filter(type == "mu", type2 == "covariate") %>% 
  dplyr::mutate(negative_log10_pvalue= -1 * log10(`Pr(>|z|)`)) %>%
  dplyr::mutate(variable = sub(param, pattern = "mu.", replacement = "") ) %>%
  dplyr::mutate(primary = ifelse(variable %in% c("M", "(Intercept)"), 0,1)) %>%
  dplyr::mutate(model = "TCRDIST BETA-BINOMIAL") %>% 
  dplyr::mutate(hmatch1 = hmatch1) %>% 
  dplyr::mutate(type = "EXACT") %>%
  dplyr::mutate(form = 'cbind(MOTIF, M - MOTIF) ~ AGE+SEX+DAYS+HLA') %>% 
  dplyr::select(c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "variable", 
                  "feature", "primary", "negative_log10_pvalue", "hmatch1", "type", "form"))

# <compare_methods_df> combine regression results from MOTIF, RADIUS, and EXACT conformant clones
compare_methods_df = dplyr::bind_rows(cb_res_table_motif, cb_res_table_radius, cb_res_table_exact) %>% 
  filter(abs(Estimate) < 5 & `Pr(>|z|)` < .99) %>%  
  mutate(type = factor(type, levels = c("MOTIF", "RADIUS", "EXACT"))) %>% 
  mutate(variable = factor(variable, levels = c("AGE","DAYS","SEXMale","HLA")))
# output table as .tsv
compare_methods_df %>% readr::write_tsv(paste0(project_path, "/", tabulated_counts_file,".beta-binomial-regression.tsv"))

# generate a volcano plots by variable and meta-clnootype conformant definition type (i.e, MOTIF, RADIUS, and EXACT)  
gg = ggplot(compare_methods_df, aes(x = Estimate, y = negative_log10_pvalue)) + 
  geom_point(pch = 20, size = .5) + 
  facet_grid(type~ variable, scale= 'free_x') + 
  geom_vline(xintercept = 0, linetype = "dashed", col = 'gray')+
  theme_bw() + 
  ylab(expression(-log[10](p-value)))
# output plots as .pdf
pdf(paste0(project_path, "/", tabulated_counts_file,".beta-binomial-regression.pdf"), width = 11, height = 8)
gg
dev.off()
# output plots as .png

png(paste0(project_path, "/", tabulated_counts_file,".beta-binomial-regression.png"), width = 1100, height = 800)
gg
dev.off()
  
  

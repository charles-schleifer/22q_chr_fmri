# 22q_chr_fmri

## Overview


## Dependencies

## Analysis Steps

#### Preparing individual-level fMRI measures
* see [22q_hoffman/22q_analysis](https://github.com/charles-schleifer/22q_hoffman/tree/main/22q_analysis)

#### Parcellated gene expression
* prep_cabnp_abagen.sh
* abagen_ahba_cabnp.py
* outputs: abagen_vol_methods_report.txt, abagen_surf_methods_report.txt, CAB-NP_subcort_abagen_expression.csv, CAB-NP_surface_abagen_expression.csv

#### Main analysis script
* 22q_chr_fmri.Rmd

#### Permutation comparison of brain maps
* after group difference maps generated in main script, need to export results and run brainSMASH to generate spatial autocorrelation-preserving permutations for null models
* cabnp_prep_brainsmash.sh
* cabnp_distance_mat.py
* cabnp_brainsmash_22q.py
* need to sync repo to hoffman cluster to use job scheduler for computing permuted PLS models. use bash scripts to submit Rscript to scheduler
    * submit_ahba_permute_hoffman.sh --> run_rscript_args.sh --> ahba_permute_hoffman.R
    * collect_permutation_results_hoffman.R
        * output: rsfa_ahba_pls_permuted_yr2.csv

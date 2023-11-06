# 22q_chr_fmri

## Overview
Analysis of parcellated Brain Signal Variability (aka Resting State Fluctuation Amplitude; RSFA), Local Connectivity (aka Network Homogeneity; NetHo), and Global Connectivity (GBC) in 22q11.2 Deletion Syndrome (22qDel), Clinical High Risk for Psychosis (CHR) and typically developing (TD) controls

## Dependencies
* several python packages
  * [abagen](https://abagen.readthedocs.io/en/stable/) for AHBA gene expression parcellation
  * [BrainSMASH](https://brainsmash.readthedocs.io/en/latest/) for generating appropriate spatial null models for brain maps 
* various [R packages](https://github.com/charles-schleifer/22q_chr_fmri/blob/main/package_versions.txt) 
* requires wb_command for ciftiTools functions that read/plot MRI data. 
  * download: https://www.humanconnectome.org/software/get-connectome-workbench
  * script expects the workbench directory to be `/Applications/workbench/bin_macosx64` (either download to this location, or edit the path in the script)
* to read fMRI results from the hoffman2 server, the script expects that the server `hoffman2.idre.ucla.edu:/u/project/cbearden/data` is mounted to your local machine at the path `~/Desktop/hoffman_mount` using an application such as SSHFS (mac download: https://osxfuse.github.io/)
  * requires first-level MRI results to be already computed on server (see https://github.com/charles-schleifer/22q_hoffman)
* paths are all relative to the project directory, which should be detected automatically by the package 'here' if RStudio is opened via the included .Rproj file 


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
        * output: rsfa_ahba_pls_permuted_yr2.csv (full PLS models for each null map are stored on hoffman in git_exclude directory)

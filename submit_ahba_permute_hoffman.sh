#!/bin/bash

# name your job
job="ahba_perm"

# set your study directory
logdir="/u/project/cbearden/data/scripts/charlie/22q_chr_fmri/git_exclude/logs"
mkdir -p ${logdir}

# choose how many jobs in array
length=10000

# path to run script
run="/u/project/cbearden/data/scripts/charlie/22q_chr_fmri/run_rscript_args.sh"

# path to R script
rscript="/u/project/cbearden/data/scripts/charlie/22q_chr_fmri/ahba_permute_hoffman.R"


# qsub command
qsub -cwd -V -N $job -o $logdir -e $logdir -t 1-${length}:1 -l h_data=16G,h_rt=4:00:00 $run $rscript --SGE_TASK_ID=$SGE_TASK_ID --measure="rsfa"

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
qsub -cwd -V -N $job -o $logdir -e $logdir -t 1-${length}:1 -l h_data=16G,h_rt=2:00:00 $run $rscript --measure="rsfa"
qsub -cwd -V -N $job -o $logdir -e $logdir -t 1-${length}:1 -l h_data=16G,h_rt=2:00:00 $run $rscript --measure="netho"
qsub -cwd -V -N $job -o $logdir -e $logdir -t 1-${length}:1 -l h_data=16G,h_rt=2:00:00 $run $rscript --measure="gbc"

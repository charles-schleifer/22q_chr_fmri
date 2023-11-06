#!/bin/sh

source /u/local/Modules/default/init/modules.sh
module load R 

# run R script with all arguments received by this script
Rscript "$@" --SGE_TASK_ID=$SGE_TASK_ID 


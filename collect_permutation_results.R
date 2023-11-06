# C. Schleifer 11/5/23
# script to collect permutation outputs and save as single CSV

# list packages to load
packages <- c("optparse","dplyr","tidyr","pls")

# install packages if not yet installed
all_packages <- rownames(installed.packages())
installed_packages <- packages %in% all_packages
if (any(installed_packages == FALSE)){install.packages(packages[!installed_packages])}

# load packages
invisible(lapply(packages, library, character.only = TRUE))

# set project path
project <- "/u/project/cbearden/data/scripts/charlie/22q_chr_fmri/"


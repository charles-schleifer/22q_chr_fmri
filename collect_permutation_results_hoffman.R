# C. Schleifer 11/5/23
# script to collect permutation outputs and save as single CSV

# list packages to load
packages <- c"dplyr","tidyr")

# install packages if not yet installed
all_packages <- rownames(installed.packages())
installed_packages <- packages %in% all_packages
if (any(installed_packages == FALSE)){install.packages(packages[!installed_packages])}

# load packages
invisible(lapply(packages, library, character.only = TRUE))

# set project path
project <- "/u/project/cbearden/data/scripts/charlie/22q_chr_fmri/"


# measures <- c("rsfa","netho","gbc")
measures <- c("gbc","netho","rsfa")
for(measure in measures){
  # list all files
  fnames <- list.files(file.path(project,"git_exclude/perm_nulls",measure), full.names = TRUE)
  # read CSV and concatenate as df columns
  yr2 <- lapply(fnames, function(f) read.csv(f, header=FALSE, sep=",")) %>% do.call(cbind,.) %>% as.data.frame
  # get column names from file names
  colnames(yr2) <- fnames %>% gsub("^.*\\_","",.) %>% gsub(".csv","",.)
  # save output
  write.table(yr2, file = file.path(project,paste0(measure,"_ahba_pls_permuted_yr2.csv")), row.names = FALSE, col.names = TRUE, sep=",", eol="\r")
}

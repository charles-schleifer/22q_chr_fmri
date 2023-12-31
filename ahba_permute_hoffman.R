# C. Schleifer 11/5/23
# script to be run as array job on hoffman for computing single instance of permuted ahba-mri pls
# saves full pls model for each permutation (output in project/git_exclude/perm_nulls), as well as CSVs with first 5 components Y R2 in project/perm_nulls

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

# get command line options for array index and measure to use
# expects that SGE_TASK_ID is passed as flagged option by bash run script
option_list <- list(
  make_option(c("--SGE_TASK_ID"), type="numeric", default=NULL, 
              help="job array index", metavar="character"),
  make_option(c("--measure"), type="character", default=NULL, 
              help="rsfa, netho, or gbc", metavar="character")) 

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

SGE_TASK_ID=opt$SGE_TASK_ID
measure=opt$measure

# get brainsmash columns to use
cols <- paste0("V",1:10000)
col <- cols[SGE_TASK_ID]

# read CAB-NP network parcellation key
# https://github.com/ColeLab/ColeAnticevicNetPartition
ji_key <- read.table(file.path(project,"CAB-NP/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR_LabelKey.txt"),header=T)

# function to read brainsmash outputs, which are CSVs with one column per parcel and one row per shuffled brain
# needs path to csv, vector of indices (e.g. 1:180 or 181:360), and a key dataframe
# output will be key df merged with n brainsmash columns named V1:Vn, with rows not in indices set to NA
read_brainsmash <- function(path, indices, key=ji_key){
  df <- read.csv(path,header=FALSE)
  # transpose so that rows are indices and columns are replicates
  dft <- as.data.frame(t(df))
  # add index column
  dft$INDEX <- indices
  rownames(dft) <- indices
  # merge with key
  out <- merge(x=dft, all.x=TRUE, y=key, all.y=TRUE, by="INDEX")
  return(out)
}

if(measure=="rsfa"){
  bsmash <- read_brainsmash(path=file.path(project,"CAB-NP/22q_TD_RSFA_L_permuted.csv"), indices=1:180)
}else if(measure=="netho"){
  bsmash <- read_brainsmash(path=file.path(project,"CAB-NP/22q_TD_NetHo_L_permuted.csv"), indices=1:180)
}else if(measure=="gbc"){  
  bsmash <- read_brainsmash(path=file.path(project,"CAB-NP/22q_TD_GBC_L_permuted.csv"), indices=1:180)
}

# load parcellated AHBA data
# this csv is the output of abagen.get_expression_data() python function to extract AHBA expression from CAB-NP surface atlas
# https://abagen.readthedocs.io/en/stable/user_guide/expression.html
ahbaSurfCABNP <- read.csv(file.path(project,"CAB-NP_surface_abagen_expression.csv"), header=T, sep=",") %>% rename("INDEX"="label")
ahbaSurfCABNP_L <- filter(ahbaSurfCABNP, INDEX <= 180)
ahbaSurfCABNP_R <- filter(ahbaSurfCABNP, INDEX > 180 & INDEX <= 360)

# function to compare chosen brain map to AHBA with PLS
ahba_pls <- function(ahba=ahbaSurfCABNP_L, map, val_col, verbose=FALSE){
  if(verbose==TRUE){print(val_col)}
  indices <- ahba$INDEX
  # get list of genes
  genes <- names(ahba)[which(names(ahba) != "INDEX")]
  # filter for only rois in ahba data and merge
  fmap <- filter(map, INDEX %in% indices)[,c("INDEX", val_col)]
  df <- merge(x=ahba, y=fmap, by="INDEX") 
  # make pls formula with val_col on left and all genes on right
  xcols <- paste(genes, collapse=" + ")
  plsformula <- as.formula(paste(val_col,"~", xcols))
  # test model with k-fold cross validation using 6 segments
  pls_test <- plsr(formula=plsformula, data=df, scale=TRUE, center=TRUE, validation="CV", na.action="na.omit", segments=6)
}

# function to return first component Y R2 from pls model
# chose components with ncomp, default is just first component
get_yr2 <- function(input, ncomp=1, estimate="CV"){
  out <- R2(input, estimate = estimate, intercept = FALSE, ncomp=ncomp)$val %>% drop
  return(out)
}

# get pls model
mod <- ahba_pls(map=bsmash, val_col=col, verbose=FALSE)
save(mod, file=file.path(project,"git_exclude/perm_mods/",measure,paste0("null_model_",measure,"_",col,".rda")))

# get explained variance
yr2 <- get_yr2(input=mod, ncomp=1:5, estimate="CV")
write.table(yr2, file=file.path(project,"git_exclude/perm_nulls/",measure,paste0("null_model_yR2_",measure,"_",col,".csv")), row.names=FALSE, col.names=FALSE, sep=",", quote=FALSE, eol="\r")

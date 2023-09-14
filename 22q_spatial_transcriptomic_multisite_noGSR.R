# C Schleifer 9/30/22
# re-run without GSR
# this script reads and plots Allen Human Brain Atlas gene expression data previously mapped to CAB-NP atlas surface with abagen
# also reads RSFA values calculated by (which script?) and tests imaging transcriptomic correlations
# cortical parvalbumin and somatostatin expression are visualized as examples
# https://abagen.readthedocs.io/en/stable/user_guide/parcellations.html
# https://github.com/ColeLab/ColeAnticevicNetPartition

####################################################################################################################################################
### set up workspace
####################################################################################################################################################

# clear workspace
rm(list = ls(all.names = TRUE))

# use SSHFS to mount hoffman2 server (download SSHFS for mac: https://osxfuse.github.io/)
# TODO: set hoffman2 username
uname <- "schleife"
# set local path to mount server
hoffman <- "~/Desktop/hoffman_mount"
# create directory if needed 
if(!file.exists(hoffman)){dir.create(hoffman)}
# make string to run as system command
mntcommand <- paste0("umount -f ", hoffman,"; sshfs ",uname,"@hoffman2.idre.ucla.edu:/u/project/cbearden/data ",hoffman)
# if hoffman directory is empty, use system command and sshfs to mount server, if not empty assume already mounted and skip
if(length(list.files(hoffman)) == 0){system(mntcommand)}else{print(paste(hoffman,"is not empty...skipping SSHFS step"))}


# list packages to load
packages <- c("devtools","conflicted","here","magrittr", "dplyr", "tidyr", "ggplot2","ggpubr","RColorBrewer", "ciftiTools","tableone", "data.table", "reshape2","neuroCombat")

# install packages if not yet installed
# note: ciftiTools install fails if R is started without enough memory on cluster (try 16G)
all_packages <- rownames(installed.packages())
installed_packages <- packages %in% all_packages
if (any(installed_packages == FALSE)){install.packages(packages[!installed_packages])}

# load packages
invisible(lapply(packages, library, character.only = TRUE))

# install neuroComBat from github 
# https://github.com/Jfortin1/neuroCombat_Rpackage
#install_github("jfortin1/neuroCombatData")
#install_github("jfortin1/neuroCombat_Rpackage")

# use the filter function from dplyr, not stats
conflict_prefer("filter", "dplyr")

# get path to project repo directory
project <- here()
print(paste("Project directory:", project))

# set up connectome workbench path for ciftiTools
# https://www.humanconnectome.org/software/get-connectome-workbench
# local wbpath (edit this path if workbench is installed in another location, e.g. on hoffman: /u/project/CCN/apps/hcp/current/workbench/bin_rh_linux64/)
# TODO: edit if necessary
wbpath <- "/Applications/workbench/bin_macosx64/"
ciftiTools.setOption("wb_path", wbpath)

# set local path to SSHFS mount point for hoffman:/u/project/cbearden/data
# TODO: edit if necessary
hoffman <- "~/Desktop/hoffman_mount/"

# load rgl for ciftiTools visualization
# may require XQartz v2.8.1 to be installed locally
if(!require('rgl', quietly=TRUE)){install.packages('rgl')}
rgl::setupKnitr()
rgl::rgl.open(); rgl::rgl.close()

####################################################################################################################################################
### read AHBA and CAB-NP data, define plotting functions, and visualize specific genes
####################################################################################################################################################

# load parcellated AHBA data
# this csv is the output of abagen.get_expression_data() python function to extract AHBA expression from CAB-NP surface atlas
# https://abagen.readthedocs.io/en/stable/user_guide/expression.html
ahbaSurfCABNP <- read.csv(file.path(project,"CAB-NP_surface_abagen_expression.csv"), header=T, sep=",")
ahbaSurfCABNP$label2 <- ahbaSurfCABNP$label
# load AHBA extracted from separated CAB-NP volume in MNI space (no cortical ROIs)
ahbaVolCABNP <- read.csv(file.path(project,"CAB-NP_subcort_abagen_expression.csv"), header=T, sep=",")
ahbaVolCABNP$label2 <- ahbaVolCABNP$label

# load CAB-NP network parcellation
# https://github.com/ColeLab/ColeAnticevicNetPartition
ji_key <- read.table(file.path(project,"CAB-NP/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR_LabelKey.txt"),header=T)
ji_net_keys <- ji_key[,c("NETWORKKEY","NETWORK")] %>% distinct %>% arrange(NETWORKKEY)
# read cifti with subcortical structures labeled 
xii_Ji_parcel <- read_cifti(file.path(project,"CAB-NP/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dscalar.nii"), brainstructures = "all")
xii_Ji_network <- read_cifti(file.path(project,"CAB-NP/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_netassignments_LR.dscalar.nii"), brainstructures = "all")
# read only surface parcels
xii_Ji_parcel_surf <- read_cifti(file.path(project,"CAB-NP/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dscalar.nii"), brainstructures = c("left", "right"))

# function to take xifti atlas (with ROIs denoted by unique values) and return list of xifti matrix indices by brain structure for each ROI
get_roi_atlas_inds <- function(xii){
  # get all unique roi labels
  mxii <- as.matrix(xii)
  vals <- mxii %>% unique %>% sort %>% as.numeric
  # get brain structures to iterate through (cortex_left, cortex_right, subcort)
  xiinames <- names(xii$data)
  # for each ROI, get the indices for each brain structure
  # output is a nested list for each ROI (with "r_" added as prefix) and each brain structure, containing an array of xifti indices corresponding to the ROI in a given brain structure
  out <- lapply(setNames(vals,paste0("r_",vals)), function(v) lapply(setNames(xiinames,xiinames), function(n) which(xii$data[[n]] == v)))
  return(out)
}

# function to create xifti for plotting ROI values on a brain atlas
# input atlas xifti and data frame with at least two cols corresponding to ROI IDs (roi_col) and output values (val_col) e.g. gene expression or functional connectivity
# output modified atlas xifti with ROI IDs replaced with output values (for visualization)
atlas_xifti_new_vals <- function(xii, df, roi_col, val_col){
  # get list of xifti indices for each ROI
  inds <- get_roi_atlas_inds(xii)
  # create blank xii from atlas
  xii_out <- xii
  for (struct in names(xii_out$data)){
    if (!is.null(xii_out$data[[struct]])){
      xii_out$data[[struct]] <- as.matrix(rep(NA,times=nrow(xii_out$data[[struct]])))
    }
  }
  # create new column named roilabel from roi_col
  df$roilabel <- df[,roi_col]
  # for each roi in xii, set all relevant vertices to value from val_col based on roi_col
  for (roi in names(inds)){
    print(roi)
    # get value for roi
    out_val <- as.numeric(filter(df[,c("roilabel",val_col)], roilabel==gsub("r_","",roi))[val_col])
    print(out_val)
    # loop through brain structures, if ROI has any indices in a structure, set those to the output value
    for (struct in names(inds[[roi]])){
      roi_inds <- inds[[roi]][[struct]]
      l <- length(roi_inds)
      if (l > 0){
        xii_out$data[[struct]][roi_inds] <- out_val
      }
    }
  }
  return(xii_out)
}

# get SST/PVALB and PVALB/SST ratios
#ahbaSurfCABNP$SST_PVALB_RATIO <- ahbaSurfCABNP$SST/ahbaSurfCABNP$PVALB
#ahbaSurfCABNP$PVALB_SST_RATIO <- ahbaSurfCABNP$PVALB/ahbaSurfCABNP$SST

# plot surface label values to check that atlas was read and written correctly (should match input atlas)
#plot_label <- atlas_xifti_new_vals(xii=xii_Ji_parcel_surf, df=ahbaSurfCABNP, roi_col="label", val_col="label2")
#view_xifti_surface(plot_label)
#view_xifti_surface(xii_Ji_parcel_surf)

# plot PVALB expression
plot_pv <- atlas_xifti_new_vals(xii=xii_Ji_parcel_surf, df=ahbaSurfCABNP, roi_col="label", val_col="PVALB")
#view_xifti_surface(plot_pv, hemisphere="left", title="PVALB", cex.title=1.3, colors="magma")

# plot SST expression
plot_sst <- atlas_xifti_new_vals(xii=xii_Ji_parcel_surf, df=ahbaSurfCABNP, roi_col="label", val_col="SST")
#view_xifti_surface(plot_sst, hemisphere="left", title="SST", cex.title=1.3, colors="magma")

# plot SST/PVALB ratio
#plot_sst_pvalb_ratio <- atlas_xifti_new_vals(xii=xii_Ji_parcel_surf, df=ahbaSurfCABNP, roi_col="label", val_col="SST_PVALB_RATIO")
#view_xifti_surface(plot_sst_pvalb_ratio,  hemisphere="left", title="SST/PVALB", cex.title=1.3, colors="magma")

#plot_pvalb_sst_ratio <- atlas_xifti_new_vals(xii=xii_Ji_parcel_surf, df=ahbaSurfCABNP, roi_col="label", val_col="PVALB_SST_RATIO")
#view_xifti_surface(plot_pvalb_sst_ratio,  hemisphere="left", title="PVALB/SST", cex.title=1.3, colors="magma")

# combine volume and surface df and plot together
# list of genes that differ between dataframes
#c(setdiff(names(ahbaVolCABNP), names(ahbaSurfCABNP)), setdiff(names(ahbaSurfCABNP), names(ahbaVolCABNP)))
# combine shared columns between vol and surf dataframes (columns may differ slightly if some genes were excluded in only surface/volume)
common_cols <- intersect(names(ahbaVolCABNP), names(ahbaSurfCABNP))
ahbaCombinedCABNP <- rbind(ahbaSurfCABNP[,common_cols], ahbaVolCABNP[,common_cols])
# get SST/PVALB and PVALB/SST
# TODO: when calculating ratio, need to avoid dividing by zero, and figure out why some regions have extremely high vals... 
#ahbaCombinedCABNP$SST_PVALB_RATIO <- ahbaCombinedCABNP$SST/ahbaCombinedCABNP$PVALB
#ahbaCombinedCABNP$PVALB_SST_RATIO <- ahbaCombinedCABNP$PVALB/ahbaCombinedCABNP$SST

# plot volume label values
plot_label_vol <- atlas_xifti_new_vals(xii=xii_Ji_parcel, df=ahbaCombinedCABNP, roi_col="label", val_col="label2")
#view_xifti(plot_label_vol, title="", cex.title=1.3, colors=rainbow(720))
#view_xifti_volume(xii_Ji_parcel, title="orig", cex.title=1.3, colors=rainbow(700))

plot_pv_vol <- atlas_xifti_new_vals(xii=xii_Ji_parcel, df=ahbaCombinedCABNP, roi_col="label", val_col="PVALB")
#view_xifti(plot_pv_vol, title="PVALB", cex.title=1.3, colors="magma")



####################################################################################################################################################
### compute RSFA difference z score in 22qDel vs HCS sample and correlate with PVALB/SST gradient
####################################################################################################################################################

# paths to sessions directories
trio_dir <- file.path(hoffman,"22q/qunex_studyfolder/sessions")
prisma_dir <- file.path(hoffman,"22qPrisma/qunex_studyfolder/sessions")
suny_dir <- file.path(hoffman,"Enigma/SUNY/qunex_studyfolder/sessions")
rome_dir <- file.path(hoffman,"Enigma/Rome/qunex_studyfolder/sessions")
iop_dir <- file.path(hoffman,"Enigma/IoP/qunex_studyfolder/sessions")

# get list of sessions
trio_sessions <- list.files(trio_dir,pattern="Q_[0-9]")
prisma_sessions <- list.files(prisma_dir,pattern="Q_[0-9]")
# exclude Q_0390_09302019 for now due to no AP BOLD
#exclude_sessions <- "Q_0390_09302019"
# for now also exclude several new sessions without rsfa calculated
exclude_sessions <- c("Q_0214_05012017","Q_0390_09302019","Q_0477_01052022","Q_0484_01042022","Q_0508_06232022","Q_0519_05312022","Q_0520_06012022","Q_0521_05202022","Q_0525_06072022","Q_0526_06242022","Q_0527_07112022","Q_0528_07202022","Q_0529_07202022","Q_0541_07182022","Q_0549_10182022","Q_0561_11032022","Q_0568_10252022")
prisma_sessions <- prisma_sessions[! prisma_sessions %in% exclude_sessions]
suny_sessions <- list.files(suny_dir,pattern="X[0-9]")
iop_sessions <- list.files(iop_dir,pattern="GQAIMS[0-9]")
rome_sessions <- c(list.files(rome_dir, pattern="C[0-9]"),list.files(rome_dir, pattern="D[0-9]"))

all_sessions <- c(trio_sessions,prisma_sessions,suny_sessions,rome_sessions,iop_sessions)

# read multisite demo table
demo_multisite <- read.csv("~/Dropbox/PhD/grants/NIMH_f31_2022/scripts/22q_multisite_demo_f31.csv")
# exclude listed sessions if any are in list
demo_multisite <- demo_multisite[! demo_multisite$MRI_S_ID %in% exclude_sessions,]
# subset to baseline
demo_multisite_bl <- filter(demo_multisite, visit_index==1)

# get list of IDs
#use_ids_22q <- filter(ucla_demo_use_hcs_del_xs, SUBJECT_IDENTITY == "PATIENT-DEL")$MRI_S_ID
#use_ids_hcs <- filter(ucla_demo_use_hcs_del_xs, SUBJECT_IDENTITY == "CONTROL")$MRI_S_ID

# read parcel RSFA
# function to read between region results and add columns for roi pair name, site, and ID 
read_csv_results <- function(sdir, fname, sesh, site){
  input <- read.csv(file.path(sdir,sesh,"images/functional",fname))
  session <- rep(sesh, times=nrow(input)) %>% as.data.frame
  site <- rep(site, times=nrow(input)) %>% as.data.frame
  new_cols <- cbind(session,site)
  colnames(new_cols) <- c("MRI_S_ID","site")
  output <- cbind(input,new_cols)
  return(output)
}

# file name to look for
#rsfa_name <- "resting_parcel_RSFA_Atlas_s_hpss_res-mVWMWB1d_lpss_whole_brain_CABNP.csv"
#rsfa_name_prisma <- "restingAP_parcel_RSFA_Atlas_s_hpss_res-mVWMWB1d_lpss_whole_brain_CABNP.csv"
rsfa_name <- "resting_parcel_RSFA_Atlas_s_hpss_res-mVWM1d_lpss_whole_brain_CABNP.csv"
rsfa_name_prisma <- "restingAP_parcel_RSFA_Atlas_s_hpss_res-mVWM1d_lpss_whole_brain_CABNP.csv"

# read rsfa
#rsfa_22q <- lapply(use_ids_22q, function(s) read_csv_results(sesh=s,site="trio",sdir=file.path(hoffman,"22q/qunex_studyfolder/sessions"),fname=rsfa_name)) 
#rsfa_hcs <- lapply(use_ids_hcs, function(s) read_csv_results(sesh=s,site="trio",sdir=file.path(hoffman,"22q/qunex_studyfolder/sessions"),fname=rsfa_name)) 

trio_rsfa <- lapply(filter(demo_multisite_bl, Site=="UCLAtrio")$MRI_S_ID, function(s) read_csv_results(sesh=s,site="UCLAtrio",sdir=trio_dir,fname=rsfa_name)) %>% do.call(rbind,.)
prisma_rsfa <- lapply(filter(demo_multisite_bl, Site=="UCLAprisma")$MRI_S_ID, function(s) read_csv_results(sesh=s,site="UCLAprisma",sdir=prisma_dir,fname=rsfa_name_prisma)) %>% do.call(rbind,.)
suny_rsfa <- lapply(filter(demo_multisite_bl, Site=="SUNY")$MRI_S_ID, function(s) read_csv_results(sesh=s,site="SUNY",sdir=suny_dir,fname=rsfa_name)) %>% do.call(rbind,.)
rome_rsfa <- lapply(filter(demo_multisite_bl, Site=="Rome")$MRI_S_ID, function(s) read_csv_results(sesh=s,site="Rome",sdir=rome_dir,fname=rsfa_name)) %>% do.call(rbind,.)
iop_rsfa <- lapply(filter(demo_multisite_bl, Site=="IoP")$MRI_S_ID, function(s) read_csv_results(sesh=s,site="IoP",sdir=iop_dir,fname=rsfa_name)) %>% do.call(rbind,.)

rsfa_all <- rbind(trio_rsfa,prisma_rsfa,suny_rsfa,rome_rsfa,iop_rsfa)
rsfa_all$INDEX %<>% as.numeric

# harmonize with neuroComBat (use this instead of longComBat for cross-sectional data)
# https://github.com/Jfortin1/neuroCombat_Rpackage
# needs input data matrix with one MRI feature per row, one column per subject

# cast to wide so that rows are parcel indices and columns are MRI_S_IDs
data.table::setDT(rsfa_all)
rsfa_all_subcols <- reshape2::dcast(data=rsfa_all, formula="INDEX ~ MRI_S_ID", value.var="t_sd") 

# get subject ids in order of columns
subcols <- which(names(rsfa_all_subcols)!= "INDEX")
subnames <- names(rsfa_all_subcols)[subcols]

# get covariates from multisite baseline demo df, ordered the same as subject ID columns for input to neurocombat
combat_covars_ordered <- merge(x=data.frame(subnames), y=demo_multisite_bl[,c("MRI_S_ID","AGE","SEX","Site")], by.x="subnames", by.y="MRI_S_ID", all.x=TRUE)

# run neurocombat for parcellated RSFA with site as batch variable and AGE and SEX included in the model
# using default parametric model (this paper shows parametric prior estimates outperform non-parametric in neurocombat: https://www.sciencedirect.com/science/article/pii/S2666956022000605)
rsfa_all_combat <- neuroCombat(dat=rsfa_all_subcols[,subcols], batch=combat_covars_ordered$Site, mod=model.matrix(~combat_covars_ordered$AGE+combat_covars_ordered$SEX), parametric=TRUE)

# get mean of all parcels per subject pre/post combat to look at distributions
parc_means_postcombat <- apply(rsfa_all_combat$dat.combat,  MARGIN=2, FUN="mean")
parc_means_precombat <- apply(rsfa_all_combat$dat.original,  MARGIN=2, FUN="mean")
parc_means <- data.frame(MRI_S_ID=names(parc_means_postcombat), postcombat=as.vector(parc_means_postcombat), precombat=as.vector(parc_means_precombat))
parc_means_demo <- merge(x=parc_means, y=demo_multisite_bl, by="MRI_S_ID", all.x=TRUE, all.y=FALSE)
parc_means_demo$Site <- factor(parc_means_demo$Site, levels=c("UCLAtrio","SUNY","Rome","IoP","UCLAprisma"))


# plot mean RSFA distributions by site pre combat
pl_precombat <- ggplot(parc_means_demo, aes(precombat, fill=Site, y=..count..))+
  geom_density(kernel="gaussian", alpha=0.7)+
  theme_classic()+
  ggtitle("Mean RSFA pre-combat")

# plot mean RSFA distributions by site post combat
pl_postcombat <- ggplot(parc_means_demo, aes(postcombat, fill=Site, y=..count..))+
  geom_density(kernel="gaussian", alpha=0.7)+
  theme_classic()+
  ggtitle("Mean RSFA post-combat")

ggarrange(pl_precombat,pl_postcombat, nrow=1, common.legend=T,legend="right")

# convert combat results  back to long df
rsfa_all_combat_d <- rsfa_all_combat$dat.combat %>% as.data.frame
data.table::setDT(rsfa_all_combat_d)
rsfa_all_combat_d$INDEX <- rownames(rsfa_all_combat_d)
rsfa_all_combat_d$INDEX %<>% as.numeric
rsfa_all_combat_l <- melt.data.table(rsfa_all_combat_d, id.vars="INDEX") %>% rename("MRI_S_ID"="variable")  %>% rename("t_sd_combat"="value")

# merge with original rsfa_all
rsfa_all_combat_l_full <- merge(x=rsfa_all, y=rsfa_all_combat_l, by=c("INDEX","MRI_S_ID"), all.x=TRUE)

## as a check, do the same for data.original from combat output. should be the same as t_sd
#rsfa_all_original_d <- rsfa_all_combat$dat.original %>% as.data.frame
#data.table::setDT(rsfa_all_original_d)
#rsfa_all_original_d$INDEX <- rownames(rsfa_all_original_d)
#rsfa_all_original_d$INDEX %<>% as.numeric
#rsfa_all_original_l <- melt.data.table(rsfa_all_original_d, id.vars="INDEX") %>% rename("MRI_S_ID"="variable")  %>% rename("t_sd_original"="value")
## merge with original rsfa_all
#rsfa_all_combat_l_full_check <- merge(x=rsfa_all_combat_l_full, y=rsfa_all_original_l, by=c("INDEX","MRI_S_ID"), all.x=TRUE)
## check that all match
#all(rsfa_all_combat_l_full_check$t_sd == rsfa_all_combat_l_full_check$t_sd_original)


# get mean and standard dev of hcs for a given parcel to use for normalization
get_parc_combat_hc_mean_sd <- function(df, hcs_ids, parc){
  dff <- filter(df, MRI_S_ID %in% hcs_ids & INDEX == parc)
  out <- data.frame(parcel_hc_mean=mean(as.numeric(dff$t_sd_combat)), parcel_hc_sd=sd(as.numeric(dff$t_sd_combat)))
  return(out)
}

# first get hc stats per parcel
use_ids_hcs <- filter(demo_multisite_bl, SUBJECT_IDENTITY=="CONTROL")$MRI_S_ID
#parc_hc_stats <- lapply(unique(rsfa_all_combat_l_full$INDEX), function(p) get_parc_combat_hc_mean_sd(parc=p, df=rsfa_all_combat_l_full, hcs_ids=use_ids_hcs)) 
parc_hc_stats <- lapply(unique(rsfa_all_combat_l_full$INDEX), function(p) get_parc_combat_hc_mean_sd(parc=p, df=rsfa_all_combat_l_full, hcs_ids=use_ids_hcs)) %>% do.call(rbind,.)
parc_hc_stats$INDEX <- unique(rsfa_all_combat_l_full$INDEX)

# merge control parcel stats with data
rsfa_normed <- merge(x=rsfa_all_combat_l_full, y=parc_hc_stats, by="INDEX", all.x=TRUE)

# normalize data in each parcel based on control mean and sd for that parcel
# t_sd_normed = (t_sd_combat - parcel_hc_mean)/parcel_hc_sd
rsfa_normed$t_sd_normed <- (rsfa_normed$t_sd_combat - rsfa_normed$parcel_hc_mean)/rsfa_normed$parcel_hc_sd

# cast to wide to merge with demographics
data.table::setDT(rsfa_normed)
# add prefix "r_" to parcel indices to use as column names in wide df 
rsfa_normed$index_col <- paste0("r_",rsfa_normed$INDEX)
# get list of indices
parc_cols <- unique(rsfa_normed$index_col)
# make wide df with column form MRI_S_ID and one column per parcel with cells containing t_sd_normed
rsfa_normed_wide <- reshape2::dcast(rsfa_normed, MRI_S_ID ~ index_col, value.var="t_sd_normed") 
# order cols
rsfa_normed_wide <- rsfa_normed_wide[,c("MRI_S_ID", parc_cols)]

# merge with demo table 
df_demo_table_full_mri <- merge(x=demo_multisite_bl, y=rsfa_normed_wide, by="MRI_S_ID")

# function to return beta coefficient for group in a lm predicting MRI from group plus covariates
#lm_parcel_group_covars <- function(df,var){
#  lm(reformulate("SUBJECT_IDENTITY + AGE + SEX", response=var),data=df)$coefficients["SUBJECT_IDENTITYPATIENT-DEL"]
#}
lm_parcel_group_covars <- function(df,var){
  lm(reformulate("SUBJECT_IDENTITY + AGE + SEX", response=var),data=df)
}

# function to apply lm at each parcel
get_parcel_group_lm <- function(df,parc_cols){
  lapply(parc_cols, function(v) lm_parcel_group_covars(var=v, df=df))
}
get_parcel_group_b_lm <- function(df,parc_cols){
  lapply(parc_cols, function(v) lm_parcel_group_covars(var=v, df=df)$coefficients["SUBJECT_IDENTITYPATIENT-DEL"]) 
}

# get group difference effect size for normalized t_sd in each parcel
#rsfa_group_betas <- get_parcel_group_b_lm(df=df_demo_table_full_mri, parc_cols=parc_cols)%>% do.call(rbind,.) %>% as.data.frame
rsfa_group_lms <- get_parcel_group_lm(df=df_demo_table_full_mri, parc_cols=parc_cols)
rsfa_group_betas <- lapply(rsfa_group_lms, function(m)m$coefficients["SUBJECT_IDENTITYPATIENT-DEL"]) %>% do.call(rbind,.) %>% as.vector
rsfa_group_ts <- lapply(rsfa_group_lms, function(m) summary(m)$coefficients["SUBJECT_IDENTITYPATIENT-DEL","t value"]) %>% do.call(rbind,.) %>% as.vector
rsfa_group_ps <- lapply(rsfa_group_lms, function(m) summary(m)$coefficients["SUBJECT_IDENTITYPATIENT-DEL","Pr(>|t|)"]) %>% do.call(rbind,.) %>% as.vector
rsfa_group_fdrp <- p.adjust(rsfa_group_ps, method="fdr")

# add group difference betas to cab-np parcel key (negative number means patients < controls)
ji_key$rsfa_diff_beta <- rsfa_group_betas
# add t- p- and q-values to cab-np parcel key
ji_key$rsfa_diff_t <- rsfa_group_ts
ji_key$rsfa_diff_para_p <- rsfa_group_ps
ji_key$rsfa_diff_para_q <- rsfa_group_fdrp
# create column of betas with non-fdr-significant set to NA
ji_key$rsfa_diff_beta_parasig <- ji_key$rsfa_diff_beta
ji_key$rsfa_diff_beta_parasig[which(ji_key$rsfa_diff_para_q > 0.05)] <- NA


# function to permute group labels to generate null distribution
# TODO: not sure if this is properly accounting for covars (see https://brainder.org/tag/permutation-test/ )
# TODO: first, check parametric p-vals from t-tests above
permute_group <- function(df,group_col){
  # take input df and randomly permute group column
  dfPerm <- df
  dfPerm[,group_col] <- sample(dfPerm[,group_col])
  return(dfPerm)
}

# get permutations of the linear model output for each parcel
nPerm <- 10000
get_parcel_group_b_lm_v <- function(i){
  print(paste0("perm:",i,"/",nPerm))
  get_parcel_group_b_lm(df=permute_group(df=df_demo_table_full_mri,group_col="SUBJECT_IDENTITY"), parc_cols=parc_cols)
}

#perm_betas<- lapply(1:nPerm, get_parcel_group_b_lm_v)
#save(perm_betas, file=file.path(project,"rsfa_permutations.rda"))

# transform list to dataframe with one row per parcel and nperm columns
#perm_beta_df <- lapply(1:length(parc_cols), function(p) lapply(1:nPerm, function(n) as.numeric(perm_betas[[n]][p])) %>% do.call(cbind,.)) %>% do.call(rbind,.) %>% as.data.frame 
#save(perm_beta_df, file=file.path(project,"rsfa_permutations_df.rda"))

# for each parcel, get the percentage of permuted trials that have an effect size with an absolute value greater than the absolute value for the real data
#group_perm_prob <- lapply(1:length(parc_cols), function(p) sum(abs(as.numeric(perm_beta_df[p,])) > abs(as.numeric(rsfa_group_betas[p,])))/nPerm) %>% do.call(rbind,.) %>% as.data.frame

## add p-values and q-values to cab-np parcel key
## TODO: add permuted p-values with proper column names 
#ji_key$rsfa_diff_permp <- group_perm_prob$V1
#ji_key$rsfa_diff_permp_fdr <- p.adjust(ji_key$rsfa_diff_permp, method="fdr")
## add group difference betas (negative number means patients < controls)
#ji_key$rsfa_diff_beta <- rsfa_group_betas[,"SUBJECT_IDENTITYPATIENT-DEL"]
## create column of betas with non-fdr-significant set to NA
#ji_key$rsfa_diff_beta_fdrsig <- ji_key$rsfa_diff_beta
#ji_key$rsfa_diff_beta_fdrsig[which(ji_key$rsfa_diff_permp_fdr > 0.05)] <- NA



# create dataframe for input to atlas_xifti_new_vals() with label column (roi_col) and RSFA outputs (val_col)
#rsfa_bg <- cbind(1:length(rsfa_diff),rsfa_diff,rsfa_delta) %>% as.data.frame
#colnames(rsfa_bg) <- c("label","diff","delta")

# color pals for brain plots
# set up rcolorbrewer palettes
pal_red_yellow_blue <- function(){
  pal <- "RdYlBu"
  ncolors <- 1000
  mycolors <- rev(colorRampPalette(brewer.pal(11,pal))(ncolors))
  return(mycolors)
}

pal_red_blue <- function(){
  pal <- "RdBu"
  ncolors <- 1000
  mycolors <- rev(colorRampPalette(brewer.pal(11,pal))(ncolors))
  return(mycolors)
}

pal_yellow_orange_red <- function(){
  pal <- "YlOrRd"
  ncolors <- 1000
  mycolors <- colorRampPalette(brewer.pal(9,pal))(ncolors)
  return(mycolors)
}

# brain plots
# group difference beta no threshold
plot_rsfa_diff_beta <- atlas_xifti_new_vals(xii=xii_Ji_parcel, df=ji_key, roi_col="INDEX", val_col="rsfa_diff_beta")
view_xifti_surface(plot_rsfa_diff_beta, title="RSFA diff beta (22qDel - HCS)", cex.title=1.3,zlim=c(-0.7,0.7),colors=pal_red_blue())
view_xifti_volume(plot_rsfa_diff_beta, title="RSFA diff beta (22qDel - HCS)", cex.title=1.3,zlim=c(-0.7,0.7),colors=pal_red_blue())


# group difference beta only fdr<0.05
plot_rsfa_diff_beta_fdr <- atlas_xifti_new_vals(xii=xii_Ji_parcel, df=ji_key, roi_col="INDEX", val_col="rsfa_diff_beta_parasig")
view_xifti_surface(plot_rsfa_diff_beta_fdr, title="RSFA diff beta FDR corrected (22qDel - HCS)", cex.title=1.3,zlim=c(-0.7,0.7), colors=pal_red_blue())
view_xifti_volume(plot_rsfa_diff_beta_fdr, title="RSFA diff beta FDR corrected (22qDel - HCS)", cex.title=1.3,zlim=c(-0.7,0.7), colors=pal_red_blue())


## get info for SOBP abstract
# counts and sex
demo_sobp <- CreateTableOne(data=demo_multisite_bl,vars=c("AGE","SEX"),strata="SUBJECT_IDENTITY",addOverall=F)
demo_sobp
# number of significant regions
filter(ji_key, !is.na(rsfa_diff_beta_parasig)) %>% nrow
# decreases network summary cortex
filter(ji_key[1:360,], rsfa_diff_beta_parasig < 0)$NETWORK %>% as.factor %>% summary
filter(ji_key[360:718,], rsfa_diff_beta_parasig < 0)$NETWORK %>% as.factor %>% summary

# increases network summary
filter(ji_key[1:360,], rsfa_diff_beta_parasig > 0)$NETWORK %>% as.factor %>% summary
filter(ji_key[360:718,], rsfa_diff_beta_parasig > 0)$NETWORK %>% as.factor %>% summary

####################################################################################################################################################
### imaging transcriptomics
####################################################################################################################################################


#plot_rsfa_delta <- atlas_xifti_new_vals(xii=xii_Ji_parcel, df=rsfa_bg, roi_col="label", val_col="delta")
#view_xifti_surface(plot_rsfa_delta, title="RSFA delta (22qDel - HCS)", cex.title=1.3, colors="magma")
#view_xifti_surface(plot_rsfa_delta, title="RSFA delta (22qDel - HCS)", cex.title=1.3, colors="magma", hemisphere = "left")


# merge fmri data with ahba by parcel
rsfa_bg_ahba <- merge(x=rsfa_bg, y=ahbaCombinedCABNP, by="label")

# correlate delta with PVALB, SST, and SST/PVALB -- wole brain
cor_wb_delta_pv <- cor.test(rsfa_bg_ahba$delta, rsfa_bg_ahba$PVALB, na.action="omit")
cor_wb_delta_st <- cor.test(rsfa_bg_ahba$delta, rsfa_bg_ahba$SST, na.action="omit")
#cor_wb_delta_st_over_pv <- cor.test(rsfa_bg_ahba$delta, rsfa_bg_ahba$SST_PVALB_RATIO, na.action="omit")

# correlate delta with PVALB, SST, and SST/PVALB -- only LH cortex (RH AHBA has more missing parcels)
rsfa_bg_lh_ahba <- filter(rsfa_bg_ahba, label <= 180)
cor_lh_delta_pv <- cor.test(rsfa_bg_lh_ahba$delta, rsfa_bg_lh_ahba$PVALB, na.action="omit")
cor_lh_delta_st <- cor.test(rsfa_bg_lh_ahba$delta, rsfa_bg_lh_ahba$SST, na.action="omit")
#cor_lh_delta_st_over_pv <- cor.test(rsfa_bg_lh_ahba$delta, rsfa_bg_lh_ahba$SST_PVALB_RATIO, na.action="omit")


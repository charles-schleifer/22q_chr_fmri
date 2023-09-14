# C Schleifer 9/30/22
# this script reads and plots Allen Human Brain Atlas gene expression data previously mapped to CAB-NP atlas surface with abagen
# also reads netho values calculated by (which script?) and tests imaging transcriptomic correlations
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
ji_key_orig <- ji_key
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

# plot surface label values to check that atlas was read and written correctly (should match input atlas)
#plot_label <- atlas_xifti_new_vals(xii=xii_Ji_parcel_surf, df=ahbaSurfCABNP, roi_col="label", val_col="label2")
#view_xifti_surface(plot_label)
#view_xifti_surface(xii_Ji_parcel_surf)

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
#plot_label_vol <- atlas_xifti_new_vals(xii=xii_Ji_parcel, df=ahbaCombinedCABNP, roi_col="label", val_col="label2")
#view_xifti(plot_label_vol, title="", cex.title=1.3, colors=rainbow(720))
#view_xifti_volume(xii_Ji_parcel, title="orig", cex.title=1.3, colors=rainbow(700))

plot_pv <- atlas_xifti_new_vals(xii=xii_Ji_parcel, df=ahbaCombinedCABNP, roi_col="label", val_col="PVALB")
#view_xifti_surface(plot_pv, hemisphere="left", title="PVALB", cex.title=1.3, colors="magma")

plot_sst <- atlas_xifti_new_vals(xii=xii_Ji_parcel, df=ahbaCombinedCABNP, roi_col="label", val_col="SST")
#view_xifti_surface(plot_sst, hemisphere="left", title="SST", cex.title=1.3, colors="magma")

# get PVALB-SST difference
ahbaCombinedCABNP$PVALB_SST_diff <- ahbaCombinedCABNP$PVALB-ahbaCombinedCABNP$SST
ahbaCombinedCABNP$SST_PVALB_diff <- ahbaCombinedCABNP$SST-ahbaCombinedCABNP$PVALB

plot_pv_sst_dif <- atlas_xifti_new_vals(xii=xii_Ji_parcel, df=ahbaCombinedCABNP, roi_col="label", val_col="PVALB_SST_diff")
#view_xifti_surface(plot_pv_sst_dif, hemisphere="left", title="PVALB - SST difference", cex.title=1.3, colors="magma")

plot_sst_pv_dif <- atlas_xifti_new_vals(xii=xii_Ji_parcel, df=ahbaCombinedCABNP, roi_col="label", val_col="SST_PVALB_diff")
#view_xifti_surface(plot_sst_pv_dif, hemisphere="left", title="SST - PVALB difference", cex.title=1.3, colors="magma")


####################################################################################################################################################
### compute netho difference z score in 22qDel vs HCS sample
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
# for now also exclude several new sessions without netho calculated
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
# filter age
demo_multisite_bl <- filter(demo_multisite_bl, AGE >= 8 & AGE <= 45)
  
## get movement info
# function to get mapping between boldn and run name from session_hcp.txt
get_boldn_names <- function(sesh,sessions_dir){
  hcptxt <- read.table(file.path(sessions_dir,sesh,"session_hcp.txt"),sep=":",comment.char="#",fill=T,strip.white=T,col.names=c(1:4)) %>% as.data.frame()
  hcpbolds <- hcptxt %>% filter(grepl("bold[0-9]",X2))
  df_out <- cbind(rep(sesh,times=nrow(hcpbolds)),hcpbolds$X2,hcpbolds$X3)
  colnames(df_out) <- c("sesh","bold_n","bold_name")
  return(df_out)
}

# function to get %udvarsme from images/functional/movement/boldn.scrub
get_percent_udvarsme <- function(sesh,sessions_dir,bold_name_use){
  mov_dir <- file.path(sessions_dir,sesh,"images/functional/movement")
  sesh_bolds <- get_boldn_names(sesh=sesh,sessions_dir=sessions_dir) %>% as.data.frame %>% filter(bold_name == bold_name_use)
  if(nrow(sesh_bolds) > 0){
    boldns_use <- sesh_bolds$bold_n %>% as.vector
    for(i in 1:length(boldns_use)){
      boldn <- boldns_use[i] %>% as.character
      boldn_path <- file.path(mov_dir,paste(boldn,".scrub",sep=""))
      mov_scrub <- read.table(boldn_path, header=T)
      percent_udvarsme <- (sum(mov_scrub$udvarsme == 1)/length(mov_scrub$udvarsme)*100) %>% as.numeric %>% signif(3)
      percent_use <- (sum(mov_scrub$udvarsme == 0)/length(mov_scrub$udvarsme)*100) %>% as.numeric %>% signif(3)
      df_out <- cbind(sesh,boldn,bold_name_use,percent_udvarsme,percent_use)
      colnames(df_out) <- c("sesh","bold_n","bold_name","percent_udvarsme","percent_use")
      return(df_out)
    }
  }
}

# get trio movement
percent_udvarsme_trio <- lapply(trio_sessions,function(s) get_percent_udvarsme(sesh=s,sessions_dir=trio_dir,bold_name_use="resting")) %>% do.call(rbind,.) %>% as.data.frame
# get prisma movement
percent_udvarsme_prisma <- lapply(prisma_sessions,function(s) get_percent_udvarsme(sesh=s,sessions_dir=prisma_dir,bold_name_use="restingAP")) %>% do.call(rbind,.) %>% as.data.frame
# get rome movement
percent_udvarsme_rome <- lapply(rome_sessions,function(s) get_percent_udvarsme(sesh=s,sessions_dir=rome_dir,bold_name_use="resting")) %>% do.call(rbind,.) %>% as.data.frame
# get suny movement
percent_udvarsme_suny <- lapply(suny_sessions,function(s) get_percent_udvarsme(sesh=s,sessions_dir=suny_dir,bold_name_use="resting")) %>% do.call(rbind,.) %>% as.data.frame
# get iop movement
percent_udvarsme_iop <- lapply(iop_sessions,function(s) get_percent_udvarsme(sesh=s,sessions_dir=iop_dir,bold_name_use="resting")) %>% do.call(rbind,.) %>% as.data.frame

# combine
percent_udvarsme_all <- rbind(percent_udvarsme_trio,percent_udvarsme_prisma,percent_udvarsme_rome,percent_udvarsme_suny,percent_udvarsme_iop)
percent_udvarsme_all$percent_udvarsme <- as.numeric(percent_udvarsme_all$percent_udvarsme)
percent_udvarsme_all$percent_use <- as.numeric(percent_udvarsme_all$percent_use)
percent_udvarsme_all$MRI_S_ID <- percent_udvarsme_all$sesh

# add movement to demo bl
demo_multisite_bl <- merge(x=demo_multisite_bl, y=percent_udvarsme_all[,c("MRI_S_ID","percent_udvarsme")], by="MRI_S_ID", all.x=TRUE)

# read parcel netho
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


netho_name <- "resting_ParcelNetHo_Atlas_s_hpss_res-mVWMWB1d_lpss_whole_brain_CABNP.csv"
netho_name_prisma <- "restingAP_ParcelNetHo_Atlas_s_hpss_res-mVWMWB1d_lpss_whole_brain_CABNP.csv"

# read netho
#netho_22q <- lapply(use_ids_22q, function(s) read_csv_results(sesh=s,site="trio",sdir=file.path(hoffman,"22q/qunex_studyfolder/sessions"),fname=netho_name)) 
#netho_hcs <- lapply(use_ids_hcs, function(s) read_csv_results(sesh=s,site="trio",sdir=file.path(hoffman,"22q/qunex_studyfolder/sessions"),fname=netho_name)) 

trio_netho <- lapply(filter(demo_multisite_bl, Site=="UCLAtrio")$MRI_S_ID, function(s) read_csv_results(sesh=s,site="UCLAtrio",sdir=trio_dir,fname=netho_name)) %>% do.call(rbind,.)
prisma_netho <- lapply(filter(demo_multisite_bl, Site=="UCLAprisma")$MRI_S_ID, function(s) read_csv_results(sesh=s,site="UCLAprisma",sdir=prisma_dir,fname=netho_name_prisma)) %>% do.call(rbind,.)
suny_netho <- lapply(filter(demo_multisite_bl, Site=="SUNY")$MRI_S_ID, function(s) read_csv_results(sesh=s,site="SUNY",sdir=suny_dir,fname=netho_name)) %>% do.call(rbind,.)
rome_netho <- lapply(filter(demo_multisite_bl, Site=="Rome")$MRI_S_ID, function(s) read_csv_results(sesh=s,site="Rome",sdir=rome_dir,fname=netho_name)) %>% do.call(rbind,.)
iop_netho <- lapply(filter(demo_multisite_bl, Site=="IoP")$MRI_S_ID, function(s) read_csv_results(sesh=s,site="IoP",sdir=iop_dir,fname=netho_name)) %>% do.call(rbind,.)

netho_all <- rbind(trio_netho,prisma_netho,suny_netho,rome_netho,iop_netho)
netho_all$INDEX %<>% as.numeric

# check for missing netho
#netho_all[which(is.na(netho_all$NetHo)),]$MRI_S_ID %>% as.factor %>% summary(maxsum=500) %>% unique
# remove any NA
netho_all <- filter(netho_all, !is.na(NetHo))

# harmonize with neuroComBat (use this instead of longComBat for cross-sectional data)
# https://github.com/Jfortin1/neuroCombat_Rpackage
# needs input data matrix with one MRI feature per row, one column per subject

# cast to wide so that rows are parcel indices and columns are MRI_S_IDs
data.table::setDT(netho_all)
netho_all_subcols <- reshape2::dcast(data=netho_all, formula="INDEX ~ MRI_S_ID", value.var="NetHo") 

# get subject ids in order of columns
subcols <- which(names(netho_all_subcols)!= "INDEX")
subnames <- names(netho_all_subcols)[subcols]

# get covariates from multisite baseline demo df, ordered the same as subject ID columns for input to neurocombat
combat_covars_ordered <- merge(x=data.frame(subnames), y=demo_multisite_bl[,c("MRI_S_ID","AGE","SEX","Site","percent_udvarsme")], by.x="subnames", by.y="MRI_S_ID", all.x=TRUE)

# run neurocombat for parcellated netho with site as batch variable and AGE and SEX included in the model
# using default parametric model (this paper shows parametric prior estimates outperform non-parametric in neurocombat: https://www.sciencedirect.com/science/article/pii/S2666956022000605)
netho_all_combat <- neuroCombat(dat=netho_all_subcols[,subcols], batch=combat_covars_ordered$Site, mod=model.matrix(~combat_covars_ordered$AGE+combat_covars_ordered$SEX+combat_covars_ordered$percent_udvarsme), parametric=TRUE)

# get mean of all parcels per subject pre/post combat to look at distributions
parc_means_postcombat <- apply(netho_all_combat$dat.combat,  MARGIN=2, FUN="mean")
parc_means_precombat <- apply(netho_all_combat$dat.original,  MARGIN=2, FUN="mean")
parc_means <- data.frame(MRI_S_ID=names(parc_means_postcombat), postcombat=as.vector(parc_means_postcombat), precombat=as.vector(parc_means_precombat))
parc_means_demo <- merge(x=parc_means, y=demo_multisite_bl, by="MRI_S_ID", all.x=TRUE, all.y=FALSE)
parc_means_demo$Site <- factor(parc_means_demo$Site, levels=c("UCLAtrio","SUNY","Rome","IoP","UCLAprisma"))


# plot meean netho distributions by site pre combat
pl_precombat <- ggplot(parc_means_demo, aes(precombat, fill=Site, y=..count..))+
  geom_density(kernel="gaussian", alpha=0.7)+
  #scale_fill_manual(values=c("lightblue","red"))+
  theme_classic()+
  ggtitle("Mean netho pre-combat")

# plot meean netho distributions by site post combat
pl_postcombat <- ggplot(parc_means_demo, aes(postcombat, fill=Site, y=..count..))+
  geom_density(kernel="gaussian", alpha=0.7)+
  #scale_fill_manual(values=c("lightblue","red"))+
  theme_classic()+
  ggtitle("Mean netho post-combat")

ggarrange(pl_precombat,pl_postcombat, nrow=1, common.legend=T,legend="right")

# convert combat results  back to long df
netho_all_combat_d <- netho_all_combat$dat.combat %>% as.data.frame
data.table::setDT(netho_all_combat_d)
netho_all_combat_d$INDEX <- rownames(netho_all_combat_d)
netho_all_combat_d$INDEX %<>% as.numeric
netho_all_combat_l <- melt.data.table(netho_all_combat_d, id.vars="INDEX") %>% rename("MRI_S_ID"="variable")  %>% rename("NetHo.combat"="value")

# merge with original netho_all
netho_all_combat_l_full <- merge(x=netho_all, y=netho_all_combat_l, by=c("INDEX","MRI_S_ID"))

## as a check, do the same for data.original from combat output. should be the same as t_sd
#netho_all_original_d <- netho_all_combat$dat.original %>% as.data.frame
#data.table::setDT(netho_all_original_d)
#netho_all_original_d$INDEX <- rownames(netho_all_original_d)
#netho_all_original_d$INDEX %<>% as.numeric
#netho_all_original_l <- melt.data.table(netho_all_original_d, id.vars="INDEX") %>% rename("MRI_S_ID"="variable")  %>% rename("t_sd_original"="value")
## merge with original netho_all
#netho_all_combat_l_full_check <- merge(x=netho_all_combat_l_full, y=netho_all_original_l, by=c("INDEX","MRI_S_ID"), all.x=TRUE)
## check that all match
#all(netho_all_combat_l_full_check$t_sd == netho_all_combat_l_full_check$t_sd_original)


# get mean and standard dev of hcs for a given parcel to use for normalization
#get_parc_combat_hc_mean_sd <- function(df, hcs_ids, parc){
#  dff <- filter(df, MRI_S_ID %in% hcs_ids & INDEX == parc)
#  out <- data.frame(parcel_hc_mean=mean(as.numeric(dff$t_sd_combat)), parcel_hc_sd=sd(as.numeric(dff$t_sd_combat)))
#  return(out)
#}

get_parc_combat_hc_mean_sd <- function(df, hcs_ids, pt_ids, parc){
  dfh <- filter(df, MRI_S_ID %in% hcs_ids & INDEX == parc)
  dfp <- filter(df, MRI_S_ID %in% pt_ids & INDEX == parc)
  out <- data.frame(parcel_hc_mean=mean(as.numeric(dfh$NetHo.combat)), parcel_hc_sd=sd(as.numeric(dfh$NetHo.combat)),parcel_pt_mean=mean(as.numeric(dfp$NetHo.combat)), parcel_pt_sd=sd(as.numeric(dfp$NetHo.combat)))
  return(out)
}

# first get hc stats per parcel
use_ids_hcs <- filter(demo_multisite_bl, SUBJECT_IDENTITY=="CONTROL")$MRI_S_ID
use_ids_del <- filter(demo_multisite_bl, SUBJECT_IDENTITY=="PATIENT-DEL")$MRI_S_ID
#parc_hc_stats <- lapply(unique(netho_all_combat_l_full$INDEX), function(p) get_parc_combat_hc_mean_sd(parc=p, df=netho_all_combat_l_full, hcs_ids=use_ids_hcs)) 
#parc_hc_stats <- lapply(unique(netho_all_combat_l_full$INDEX), function(p) get_parc_combat_hc_mean_sd(parc=p, df=netho_all_combat_l_full, hcs_ids=use_ids_hcs)) %>% do.call(rbind,.)
parc_hc_stats <- lapply(unique(netho_all_combat_l_full$INDEX), function(p) get_parc_combat_hc_mean_sd(parc=p, df=netho_all_combat_l_full, hcs_ids=use_ids_hcs, pt_ids=use_ids_del)) %>% do.call(rbind,.)
parc_hc_stats$INDEX <- unique(netho_all_combat_l_full$INDEX)

# merge control parcel stats with data
netho_normed <- merge(x=netho_all_combat_l_full, y=parc_hc_stats, by="INDEX", all.x=TRUE)

# normalize data in each parcel based on control mean and sd for that parcel
netho_normed$NetHo.normed <- (netho_normed$NetHo.combat - netho_normed$parcel_hc_mean)/netho_normed$parcel_hc_sd

# cast to wide to merge with demographics
data.table::setDT(netho_normed)
# add prefix "r_" to parcel indices to use as column names in wide df 
netho_normed$index_col <- paste0("r_",netho_normed$INDEX)
# get list of indices
parc_cols <- unique(netho_normed$index_col)
# make wide df with column form MRI_S_ID and one column per parcel with cells containing NetHo.normed
netho_normed_wide <- reshape2::dcast(netho_normed, MRI_S_ID ~ index_col, value.var="NetHo.normed") 
# order cols
netho_normed_wide <- netho_normed_wide[,c("MRI_S_ID", parc_cols)]

## test without normalization
## cast to wide to merge with demographics
#data.table::setDT(netho_all_combat_l_full)
## add prefix "r_" to parcel indices to use as column names in wide df 
#netho_all_combat_l_full$index_col <- paste0("r_",netho_all_combat_l_full$INDEX)
## get list of indices
#parc_cols <- unique(netho_all_combat_l_full$index_col)
## make wide df with column form MRI_S_ID and one column per parcel with cells containing NetHo.normed
#netho_all_combat_wide <- reshape2::dcast(netho_all_combat_l_full, MRI_S_ID ~ index_col, value.var="NetHo.combat") 
## order cols
#netho_all_combat_wide <- netho_all_combat_wide[,c("MRI_S_ID", parc_cols)]

# merge with demo table 
df_demo_table_full_mri <- merge(x=demo_multisite_bl, y=netho_normed_wide, by="MRI_S_ID")
#df_demo_table_full_mri <- merge(x=demo_multisite_bl, y=netho_all_combat_wide, by="MRI_S_ID")

# function to return beta coefficient for group in a lm predicting MRI from group plus covariates
#lm_parcel_group_covars <- function(df,var){
#  lm(reformulate("SUBJECT_IDENTITY + AGE + SEX", response=var),data=df)$coefficients["SUBJECT_IDENTITYPATIENT-DEL"]
#}
lm_parcel_group_covars <- function(df,var){
  lm(reformulate("SUBJECT_IDENTITY + AGE + SEX + Site + percent_udvarsme", response=var),data=df)
}

# function to apply lm at each parcel
get_parcel_group_lm <- function(df,parc_cols){
  lapply(parc_cols, function(v) lm_parcel_group_covars(var=v, df=df))
}
#get_parcel_group_b_lm <- function(df,parc_cols){
#  lapply(parc_cols, function(v) lm_parcel_group_covars(var=v, df=df)$coefficients["SUBJECT_IDENTITYPATIENT-DEL"]) 
#}

# get group difference effect size for normalized t_sd in each parcel
#netho_group_betas <- get_parcel_group_b_lm(df=df_demo_table_full_mri, parc_cols=parc_cols)%>% do.call(rbind,.) %>% as.data.frame
netho_group_lms <- get_parcel_group_lm(df=df_demo_table_full_mri, parc_cols=parc_cols)
netho_group_betas <- lapply(netho_group_lms, function(m)m$coefficients["SUBJECT_IDENTITYPATIENT-DEL"]) %>% do.call(rbind,.) %>% as.vector
netho_group_ts <- lapply(netho_group_lms, function(m) summary(m)$coefficients["SUBJECT_IDENTITYPATIENT-DEL","t value"]) %>% do.call(rbind,.) %>% as.vector
netho_group_ps <- lapply(netho_group_lms, function(m) summary(m)$coefficients["SUBJECT_IDENTITYPATIENT-DEL","Pr(>|t|)"]) %>% do.call(rbind,.) %>% as.vector
netho_group_fdrp <- p.adjust(netho_group_ps, method="fdr")

# set key to only used ROIs
ji_key <- ji_key[as.numeric(gsub("r_","",parc_cols)),]

# add group difference betas to cab-np parcel key (negative number means patients < controls)
ji_key$netho_diff_beta <- netho_group_betas
# add t- p- and q-values to cab-np parcel key
ji_key$netho_diff_t <- netho_group_ts
ji_key$netho_diff_para_p <- netho_group_ps
ji_key$netho_diff_para_q <- netho_group_fdrp
# create column of betas with non-fdr-significant set to NA
ji_key$netho_diff_beta_parasig <- ji_key$netho_diff_beta
ji_key$netho_diff_beta_parasig[which(ji_key$netho_diff_para_q > 0.05)] <- NA

# add group means to key
ji_key <- merge(x=ji_key, y=parc_hc_stats, by="INDEX", all.x=TRUE)

# merge back with original key to have all 718 rows
ji_key_out <- merge(x=ji_key_orig, y=ji_key, by="INDEX", all.x=TRUE)

# write results to csv
#write.csv(ji_key_out, file=file.path(project,"cabnp_key_netho.csv"),row.names = FALSE)

# create dataframe for input to atlas_xifti_new_vals() with label column (roi_col) and netho outputs (val_col)
#netho_bg <- cbind(1:length(netho_diff),netho_diff,netho_delta) %>% as.data.frame
#colnames(netho_bg) <- c("label","diff","delta")

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
plot_netho_diff_beta <- atlas_xifti_new_vals(xii=xii_Ji_parcel, df=ji_key, roi_col="INDEX", val_col="netho_diff_beta")
view_xifti_surface(plot_netho_diff_beta, title="netho diff beta (22qDel - HCS)", cex.title=1.3,zlim=c(-0.7,0.7),colors=pal_red_blue())
view_xifti_volume(plot_netho_diff_beta, title="netho diff beta (22qDel - HCS)", cex.title=1.3,zlim=c(-0.7,0.7),colors=pal_red_blue())


# group difference beta only fdr<0.05
plot_netho_diff_beta_fdr <- atlas_xifti_new_vals(xii=xii_Ji_parcel, df=ji_key, roi_col="INDEX", val_col="netho_diff_beta_parasig")
view_xifti_surface(plot_netho_diff_beta_fdr, title="netho diff beta FDR corrected (22qDel - HCS)", cex.title=1.3,zlim=c(-0.7,0.7), colors=pal_red_blue())
view_xifti_volume(plot_netho_diff_beta_fdr, title="NetHo diff", cex.title=1.3,zlim=c(-0.8,0.8), colors=pal_red_blue(),crop=FALSE, orientation_labels=FALSE, slices=seq(16,42, by=5))

# save fdr corrected maps 
view_xifti_surface(plot_netho_diff_beta_fdr, title="NetHo diff beta FDR corrected (22qDel - HCS)", cex.title=1.3,zlim=c(-0.8,0.8), colors=pal_red_blue(), fname=file.path(project,"figures/netho/netho_fdr_surf.png"), width=1200, height=800)
view_xifti_volume(plot_netho_diff_beta_fdr, title="NetHo diff", cex.title=1.3,zlim=c(-0.8,0.8),slices=seq(16,42, by=5), colors=pal_red_blue(),crop=FALSE, orientation_labels=FALSE, legend_embed=FALSE, fname=file.path(project,"figures/netho/netho_fdr_vol.pdf"), width=1500, height=2000, format="pdf")


# control means
plot_netho_hc <- atlas_xifti_new_vals(xii=xii_Ji_parcel, df=ji_key, roi_col="INDEX", val_col="parcel_hc_mean")
view_xifti_surface(plot_netho_hc, title="netho Control Means", zlim=c(0,25),cex.title=1.3, colors="magma")

# 22qdel means
plot_netho_del <- atlas_xifti_new_vals(xii=xii_Ji_parcel, df=ji_key, roi_col="INDEX", val_col="parcel_pt_mean")
view_xifti_surface(plot_netho_del, title="netho 22qDel Means", zlim=c(0,25),cex.title=1.3, colors="magma")


## get info for SOBP abstract
# counts and sex
demo_sobp <- CreateTableOne(data=demo_multisite_bl,vars=c("AGE","SEX"),strata="SUBJECT_IDENTITY",addOverall=F)
demo_sobp
# number of significant regions
filter(ji_key, !is.na(netho_diff_beta_parasig)) %>% nrow
# decreases network summary cortex
filter(ji_key[1:360,], netho_diff_beta_parasig < 0)$NETWORK %>% as.factor %>% summary
# decreases network summary subcortex
filter(ji_key[360:718,], netho_diff_beta_parasig < 0)$NETWORK %>% as.factor %>% summary

# increases network summary cortex
filter(ji_key[1:360,], netho_diff_beta_parasig > 0)$NETWORK %>% as.factor %>% summary
# increases network summary subcortex
filter(ji_key[360:718,], netho_diff_beta_parasig > 0)$NETWORK %>% as.factor %>% summary

####################################################################################################################################################
### bin by age and test group difference
####################################################################################################################################################

# plot age distributions
ggplot(data=df_demo_table_full_mri, aes(AGE, fill=SUBJECT_IDENTITY, color=SUBJECT_IDENTITY, y=..count..))+
  geom_density(kernel="gaussian", alpha=0.5)+
  scale_fill_manual(values=c("lightblue","red"))+
  scale_color_manual(values=c("lightblue","red"))+
  theme_classic()+
  ggtitle("Age")


#filter(df_demo_table_full_mri, AGE >= 6 & AGE < 14)
#filter(df_demo_table_full_mri, AGE >= 14 & AGE < 22)
#filter(df_demo_table_full_mri, AGE >= 22 & AGE < 30)
#filter(df_demo_table_full_mri, AGE >= 30 & AGE <= 45)

####################################################################################################################################################
### imaging transcriptomics
####################################################################################################################################################


# merge fmri data with ahba by parcel
netho_ahba <- merge(x=ji_key, y=ahbaCombinedCABNP, by.x="INDEX", by.y="label")

# get set of usable parcels from left cortex by excluding rows where a test gene (e.g. PVALB) is NA
#netho_lh_ahba <- filter(netho_ahba, INDEX <= 180 & !is.na(PVALB))
netho_lh_ahba <- filter(netho_ahba, INDEX <= 180)

# correlate delta with PVALB, SST, and PVALB-SST diff
cor_lh_pv <- cor.test(netho_lh_ahba$netho_diff_beta, netho_lh_ahba$PVALB, na.action="omit")
cor_lh_st <- cor.test(netho_lh_ahba$netho_diff_beta, netho_lh_ahba$SST, na.action="omit")
cor_lh_pv_sst_diff <- cor.test(netho_lh_ahba$netho_diff_beta, netho_lh_ahba$PVALB_SST_diff, na.action="omit")

cor_lh_pv
cor_lh_st
cor_lh_pv_sst_diff


# correlate group means with PVALB
cor_hcmean_pv <- cor.test(netho_lh_ahba$parcel_hc_mean, netho_lh_ahba$PVALB, na.action="omit")
cor_delmean_pv <- cor.test(netho_lh_ahba$parcel_pt_mean, netho_lh_ahba$PVALB, na.action="omit")
cor_hcmean_pv
cor_delmean_pv

# correlate group means with SST
cor_hcmean_st <- cor.test(netho_lh_ahba$parcel_hc_mean, netho_lh_ahba$SST, na.action="omit")
cor_delmean_st <- cor.test(netho_lh_ahba$parcel_pt_mean, netho_lh_ahba$SST, na.action="omit")
cor_hcmean_st
cor_delmean_st


# read gene lists from Ziesel et al 2015
zcells <- read.table(file=file.path(project,"Reference_Consistent_Genes_ObtainedBy2StageFiltering.tsv"), header = TRUE, sep="\t", na.strings = "NA", strip.white = TRUE)
zcells_use <- filter(zcells, !is.na(CellType))

# list of genes from AHBA
ahba_names <- names(ahbaCombinedCABNP)[which(!names(ahbaCombinedCABNP) == "label")]

# test if any ziesel cells not in AHBA extracted data
# TODO: ~46 genes in ziesel list not in our AHBA extracted vol+surf data, check why genes excluded by abagen
zcells_use_final <- filter(zcells_use, GeneSymbol %in% ahba_names)
zcells_use_missing <- filter(zcells_use, !GeneSymbol %in% ahba_names)

# list of cell types
zcells_types <- unique(zcells_use_final$CellType)

# function to compute correlation with a given gene
gene_cor <- function(df, mricol, gene){
  out <- cor(df[,mricol], df[,gene], use="complete.obs", method="pearson")
  return(out)
}

# function to compute gene_cor for list of genes
gene_cor_list <- function(df, mricol, genelist){
  out <- lapply(genelist, function(g) gene_cor(df=df, mricol=mricol, gene=g))
  names(out) <- genelist
  return(out)
}

#gene_cor(df=netho_lh_ahba, mricol="netho_diff_beta", gene="PVALB")
#gene_cor_list(df=netho_lh_ahba, mricol="netho_diff_beta", genelist=filter(zcells_use_final, CellType=="Interneuron")$GeneSymbol) %>% do.call(rbind,.) %>% as.vector

zcell_netho_distributions <- lapply(zcells_types, function(t) gene_cor_list(df=netho_lh_ahba, mricol="netho_diff_beta", genelist=filter(zcells_use_final, CellType==t)$GeneSymbol))
names(zcell_netho_distributions) <- zcells_types

# function to plot distribution of correlations
plot_dist <- function(vec,title="title"){
  input <- vec %>% as.data.frame %>% rename("correlation"=".")
  x <- mean(input$correlation)
  ggplot(data=input, aes(correlation, y=..count..))+
    geom_density(kernel="gaussian")+
    geom_vline(xintercept=x, lty="dashed")+
    xlim(-0.5,0.5)+
    theme_classic()+
    ggtitle(title)
}

zcell_netho_plots <- lapply(zcells_types, function(t) plot_dist(vec=unlist(zcell_netho_distributions[[t]]), title=t))

####################################################################################################################################################
### behavior correlations 
####################################################################################################################################################

### First try UCLA only

# function to get psychosis status from SIPS
# to be applied to the row indices of a demographics df, and also given the SIPS df
get_sips <- function(r,demo,sips){
  sub <- demo$SUBJECTID[[r]]
  visit <- demo$CONVERTEDVISITNUM[[r]]
  df_out <- cbind(sub,visit) %>% as.data.frame
  colnames(df_out) <- c("SUBJECTID","CONVERTEDVISITNUM")
  sips_sesh <- sips %>% filter(SUBJECTID == sub & CONVERTEDVISITNUM == visit)
  if(nrow(sips_sesh) < 1){
    # if no match for sub+visit in sips table, set outputs to NA
    df_out[,c("SIPS_p_sum","SIPS_n_sum","SIPS_d_sum","SIPS_g_sum","SIPS_total","SIPS_psychosis_6","SIPS_psspectrum_3")] <- rep(NA,times=7)
  }else if(nrow(sips_sesh) > 1){
    # if more than one match for sub+visit, note error
    df_out[,c("SIPS_p_sum","SIPS_n_sum","SIPS_d_sum","SIPS_g_sum","SIPS_total","SIPS_psychosis_6","SIPS_psspectrum_3")] <- rep("ERROR-duplicates",times=7)
  }else{
    # get SIPS P scores
    sips_p_scores <- sips_sesh[c("P1SEV","P2SEV","P3SEV","P4SEV","P5SEV")]
    # sum SIPS P
    df_out[,"SIPS_p_sum"] <- sum(sips_p_scores)
    # get SIPS N
    sips_n_scores <- sips_sesh[c("N1SEV","N2SEV","N3SEV","N4SEV","N5SEV","N6SEV")]
    df_out[,"SIPS_n_sum"] <- sum(sips_n_scores)
    # get SIPS D
    sips_d_scores <- sips_sesh[c("D1SEV","D2SEV","D3SEV","D4SEV")]
    df_out[,"SIPS_d_sum"] <- sum(sips_d_scores)
    # get SIPS G
    sips_g_scores <- sips_sesh[c("G1SEV","G2SEV","G3SEV","G4SEV")]
    df_out[,"SIPS_g_sum"] <- sum(sips_g_scores)
    # get SIPS total
    df_out["SIPS_total"] <- (df_out["SIPS_p_sum"]+df_out["SIPS_n_sum"]+df_out["SIPS_d_sum"]+df_out["SIPS_g_sum"] )
    # check psychosis criteria of at least one SIPS P score of 6
    count_6 <- length(which(sips_p_scores == 6))
    if(is.na(sum(sips_p_scores))){
      df_out[,"SIPS_psychosis_6"] <- NA
    }else if(count_6 > 0){
      df_out[,"SIPS_psychosis_6"] <- 1
    }else{
      df_out[,"SIPS_psychosis_6"] <- 0
    }
    # check psychosis-spectrum criteria of at least one SIPS P >= 3
    count_3 <- length(which(sips_p_scores >= 3))
    if(is.na(sum(sips_p_scores))){
      df_out[,"SIPS_psspectrum_3"] <- NA
    }else if(count_3 > 0){
      df_out[,"SIPS_psspectrum_3"] <- 1
    }else{
      df_out[,"SIPS_psspectrum_3"] <- 0
    }
  }
  return(df_out)
}

# set location of directory with ucla sistat CSVs
csvdir_ucla <- file.path(project,"demographics/ucla_sistat")
# get list of files_ucla in directory
files_ucla <- list.files(csvdir_ucla)
fpaths <- lapply(files_ucla, function(file) paste(csvdir_ucla,file,sep="/"))
# clean names
fnames <- gsub(".csv","",files_ucla)
fnames <- gsub("Re22Q_","",fnames)
fnames <- gsub("Form_","",fnames)
fnames <- gsub("Qry_","",fnames)
# read all, set to na: "-9999", "-9998","." 
input_all_ucla <- lapply(fpaths, read.csv, header=T, na.strings=c(".","-9999","-9998"), strip.white=T, sep=",")
names(input_all_ucla) <- fnames
df_all_ucla <- lapply(input_all_ucla, function(x) data.frame(x))

# get ucla only subset df_demo_table_full_mri
df_demo_table_ucla_mri <- filter(df_demo_table_full_mri, Site %in% c("UCLAtrio","UCLAprisma"))
df_demo_table_ucla_mri <- merge(x=df_demo_table_ucla_mri, y=df_all_ucla$demo_mri[,c("MRI_S_ID","CONVERTEDVISITNUM")], by="MRI_S_ID", all.x=TRUE)

# get sips 
demo_table_sips <- lapply(1:nrow(df_all_ucla$demo_mri), function(r) get_sips(r=r, demo=df_all_ucla$demo_mri, sips=df_all_ucla$SIPS)) %>% do.call(rbind,.)

# merge sips with demo table
df_demo_table_ucla_mri_beh <- merge(x=df_demo_table_ucla_mri, y=demo_table_sips[,c("SUBJECTID","CONVERTEDVISITNUM","SIPS_total","SIPS_psspectrum_3","SIPS_p_sum","SIPS_n_sum","SIPS_d_sum","SIPS_g_sum")], by=c("SUBJECTID","CONVERTEDVISITNUM"), all.x=T) %>% rename("SIPS_prodromal" = "SIPS_psspectrum_3")

# set sips_total to numeric and sips_prodromal to logical
df_demo_table_ucla_mri_beh %<>% mutate_at(vars("SIPS_prodromal"), ~as.logical(.))
df_demo_table_ucla_mri_beh %<>% mutate_at(vars("SIPS_total"), ~as.numeric(.))


# get indices of regions with significantly decreased netho
ji_key_neg_netho <- filter(ji_key, netho_diff_beta_parasig < 0 )
neg_netho_parcels <- ji_key_neg_netho$INDEX
neg_netho_parces_colnames <- paste0("r_",neg_netho_parcels)

# get column of the mean normed value of the rois with significant negative 22qDel effect for each subject
df_demo_table_ucla_mri_beh$mean_sig_decrease_rois <- apply(df_demo_table_ucla_mri_beh[,neg_netho_parces_colnames],1,mean)

# subset to only 22qDel
df_demo_table_ucla22q_mri_beh <- filter(df_demo_table_ucla_mri_beh, SUBJECT_IDENTITY=="PATIENT-DEL")

# correlate composite negative netho value with SIPS
cor.test(x=df_demo_table_ucla22q_mri_beh$mean_sig_decrease_rois, y=df_demo_table_ucla22q_mri_beh$SIPS_p_sum, method="pearson",  na.action="omit")
cor.test(x=df_demo_table_ucla22q_mri_beh$mean_sig_decrease_rois, y=df_demo_table_ucla22q_mri_beh$SIPS_n_sum, method="pearson",  na.action="omit")
cor.test(x=df_demo_table_ucla22q_mri_beh$mean_sig_decrease_rois, y=df_demo_table_ucla22q_mri_beh$SIPS_g_sum, method="pearson",  na.action="omit")
cor.test(x=df_demo_table_ucla22q_mri_beh$mean_sig_decrease_rois, y=df_demo_table_ucla22q_mri_beh$SIPS_d_sum, method="pearson",  na.action="omit")
cor.test(x=df_demo_table_ucla22q_mri_beh$mean_sig_decrease_rois, y=df_demo_table_ucla22q_mri_beh$SIPS_total, method="pearson",  na.action="omit")

# TODO: in 22q, lm(parcel_netho ~ SIPS_p)







#########
### get 22q11 gene expression for hip and thal regions
#########

# genes from https://academic.oup.com/cercor/article/31/7/3285/6150031 supplement
forsyth_genes <- c("DGCR8", "AIFM3", "SCARF2", "CLDN5", "DGCR2", "P2RX6", "TANGO2", "RANBP1", "HIRA", "UFD1", "ARVCF", "MED15", "COMT", "PRODH", "SLC25A1", "GNB1L", "GP1BB", "MRPL40", "KLHL22", "RIMBP3", "PI4KA", "RTN4R", "C22orf39", "SEPT5", "DGCR6", "SLC7A4", "DGCR6L", "SNAP29")

# n=3 not in AHBA volume CABNP
forsyth_genes_missing <- forsyth_genes[which(!forsyth_genes %in% names(ahbaVolCABNP))]
forsyth_genes_use <- forsyth_genes[which(forsyth_genes %in% names(ahbaVolCABNP))]

# get main structure label from xifti metadata for each subcortical roi (361:718)
#sc_parc_labs <- lapply(361:718, function(r) xii_Ji_parcel$meta$subcort$labels[which(xii_Ji_parcel$data$subcort==r)])
#sc_parc_labs <- lapply(sc_parc_labs, droplevels)
#sc_parc_labs_char <- data.frame(INDEX=361:718,structure=do.call(rbind,lapply(sc_parc_labs, levels)))
#ji_key_sc <- merge(x=filter(ji_key, INDEX >= 361 & INDEX <= 718), y=sc_parc_labs_char, by="INDEX")

# get key for subcortex only
ji_key_sc <- filter(ji_key, INDEX >= 361 & INDEX <= 718)

# get structure name from LABEL column
ji_key_sc$structure <- do.call(rbind,strsplit(ji_key_sc$LABEL, split="_"))[,2]

# get all right hippocampal parcels
ji_key_sc_rhipp <- filter(ji_key_sc, structure=="R-Hippocampus")

# get all right thal frontoparietal parcels
ji_key_sc_rthal_fpn <- filter(ji_key_sc, structure=="R-Thalamus" & NETWORK =="Frontoparietal")
# single parcel for right thal fpn
ind_rthal_fpn <- ji_key_sc_rthal_fpn$INDEX
# 22q gene expression for right thal fpn
exp22q_rthal_fpn <- filter(ahbaVolCABNP, label==ind_rthal_fpn)[,forsyth_genes_use]
exp22q_rthal_fpn <- apply(exp22q_rthal_fpn,2,mean)
# 22q gene expression for right hipp
exp22q_rhipp <- filter(ahbaVolCABNP, label %in% ji_key_sc_rhipp$INDEX, !is.na(DGCR8))[,forsyth_genes_use]
exp22q_rhipp_mean <- apply(exp22q_rhipp,2,mean)

# make data frame for plotting
rhipp <- data.frame(exp=exp22q_rhipp_mean, gene=names(exp22q_rhipp_mean), structure="hipp_mean")
rthal <- data.frame(exp=exp22q_rthal_fpn, gene=names(exp22q_rthal_fpn), structure="thal_fpn")
thal_hip <- rbind(rthal, rhipp)

ggplot(thal_hip, aes(exp, fill=structure))+
  geom_density(kernel="gaussian", alpha=0.3)+
  scale_fill_manual(values=c("lightblue","red"))+
  geom_vline(xintercept = mean(exp22q_rhipp_mean), lty="dashed")+
  geom_vline(xintercept = mean(exp22q_rthal_fpn), lty="dashed")+
  theme_classic()

# use limma package for differential gene expression 
# example paper using method: https://www.nature.com/articles/s41380-022-01489-8
 
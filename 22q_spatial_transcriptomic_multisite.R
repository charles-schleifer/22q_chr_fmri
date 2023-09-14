# C Schleifer 9/30/22
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
packages <- c("devtools","conflicted","here","magrittr", "dplyr", "tidyr", "ggplot2","ggpubr","RColorBrewer", "ciftiTools","tableone", "data.table", "reshape2","neuroCombat","pls")

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
### read sistat data and motion data to get subject lists for analysis 
####################################################################################################################################################

## get motion data for all sessions by reading movement scrubbing files on hoffman
## TODO: confirm that hoffman path is set correctly for local machine
## TODO: currently testing with only trio data, need to add prisma and other sites
#
## function to get mapping between boldn and run name from session_hcp.txt
#get_boldn_names <- function(sesh,sessions_dir){
#  hcptxt <- read.table(file.path(sessions_dir,sesh,"session_hcp.txt"),sep=":",comment.char="#",fill=T,strip.white=T,col.names=c(1:4)) %>% as.data.frame()
#  hcpbolds <- hcptxt %>% filter(grepl("bold[0-9]",X2))
#  df_out <- cbind(rep(sesh,times=nrow(hcpbolds)),hcpbolds$X2,hcpbolds$X3)
#  colnames(df_out) <- c("sesh","bold_n","bold_name")
#  return(df_out)
#}
#
## function to get %udvarsme from images/functional/movement/boldn.scrub (the percentage of frames scrubbed due to either exceeding FD movement threshold or DVARSME signal change threshold) 
#get_percent_udvarsme <- function(sesh,sessions_dir,bold_name_use){
#  mov_dir <- file.path(sessions_dir,sesh,"images/functional/movement")
#  sesh_bolds <- get_boldn_names(sesh=sesh,sessions_dir=sessions_dir) %>% as.data.frame %>% filter(bold_name == bold_name_use)
#  if(nrow(sesh_bolds) > 0){
#    boldns_use <- sesh_bolds$bold_n %>% as.vector
#    for(i in 1:length(boldns_use)){
#      boldn <- boldns_use[i] %>% as.character
#      boldn_path <- file.path(mov_dir,paste(boldn,".scrub",sep=""))
#      mov_scrub <- read.table(boldn_path, header=T)
#      percent_udvarsme <- (sum(mov_scrub$udvarsme == 1)/length(mov_scrub$udvarsme)*100) %>% as.numeric %>% signif(3)
#      percent_use <- (sum(mov_scrub$udvarsme == 0)/length(mov_scrub$udvarsme)*100) %>% as.numeric %>% signif(3)
#      df_out <- cbind(sesh,boldn,bold_name_use,percent_udvarsme,percent_use)
#      colnames(df_out) <- c("sesh","bold_n","bold_name","percent_udvarsme","percent_use")
#      return(df_out)
#    }
#  }
#}
#
## first get all session names on hoffman
#all_22qtrio <- list.files(file.path(hoffman,"22q/qunex_studyfolder/sessions"),pattern="Q_[0-9]")
##all_22qprisma <- list.files(file.path(hoffman,"22qPrisma/qunex_studyfolder/sessions"),pattern="Q_[0-9]")
##all_suny <- list.files(file.path(hoffman,"Enigma/SUNY/qunex_studyfolder/sessions"),pattern="X[0-9]")
##all_iop <- list.files(file.path(hoffman,"Enigma/IoP/qunex_studyfolder/sessions"),pattern="GQAIMS[0-9]")
##all_rome <- list.files(file.path(hoffman,"Enigma/Rome/qunex_studyfolder/sessions"),pattern="[0-9]")
##all_multisite <- c(all_22qtrio,all_22qprisma,all_suny,all_iop,all_rome)
#
## get motion stats (percent udvarsme)
#percent_udvarsme_22qtrio <- lapply(all_22qtrio,function(s) get_percent_udvarsme(sesh=s,sessions_dir=file.path(hoffman,"22q/qunex_studyfolder/sessions"),bold_name_use="resting")) %>% do.call(rbind,.) %>% as.data.frame
##percent_udvarsme_22qprisma <- lapply(all_22qprisma,function(s) get_percent_udvarsme(sesh=s,sessions_dir=file.path(hoffman,"22qPrisma/qunex_studyfolder/sessions"),bold_name_use="restingAP")) %>% do.call(rbind,.) %>% as.data.frame
##percent_udvarsme_suny <- lapply(all_suny,function(s) get_percent_udvarsme(sesh=s,sessions_dir=file.path(hoffman,"Enigma/SUNY/qunex_studyfolder/sessions"),bold_name_use="resting")) %>% do.call(rbind,.) %>% as.data.frame
##percent_udvarsme_iop <- lapply(all_iop,function(s) get_percent_udvarsme(sesh=s,sessions_dir=file.path(hoffman,"Enigma/IoP/qunex_studyfolder/sessions"),bold_name_use="resting")) %>% do.call(rbind,.) %>% as.data.frame
##percent_udvarsme_rome <- lapply(all_rome,function(s) get_percent_udvarsme(sesh=s,sessions_dir=file.path(hoffman,"Enigma/Rome/qunex_studyfolder/sessions"),bold_name_use="resting")) %>% do.call(rbind,.) %>% as.data.frame
#
## TODO: add other sites, for now just test trio
##percent_udvarsme_all <- rbind(percent_udvarsme_22qtrio,percent_udvarsme_22qprisma,percent_udvarsme_suny,percent_udvarsme_iop,percent_udvarsme_rome)
#percent_udvarsme_all <- percent_udvarsme_22qtrio
#
## ensure that values are read as numeric
#percent_udvarsme_all$percent_udvarsme <- as.numeric(percent_udvarsme_all$percent_udvarsme)
#percent_udvarsme_all$percent_use <- as.numeric(percent_udvarsme_all$percent_use)
## get sessions to exclude with over 50% bad frames (this should be empty because all scans with too much movement have been moved to unused_sessions)
#percent_udvarsme_over50 <- filter(percent_udvarsme_all,percent_udvarsme_all$percent_udvarsme > 50)$sesh
#percent_udvarsme_under50 <- filter(percent_udvarsme_all,percent_udvarsme_all$percent_udvarsme < 50)$sesh
#
#
### load sistat data and get lists of scans to use
## all sistat tables should be exported as CSVs into a single directory
## the next several chunks deal with reading, cleaning and annotating the data exported from sistat, and then age matching
## the hcs sample is younger than del due to a large amount of very young hcs subjects. plan is to match samples by using followup timepoints rather than baseline for some younger participants, and dropping several older del subjects, and younger hcs subjects (prioritizing dropping subjects with worse motion stats when possible)
#
## set location of directory with ucla sistat CSVs
#csvdir_ucla <- file.path(project,"demographics/ucla_sistat")
## get list of files_ucla in directory
#files_ucla <- list.files(csvdir_ucla)
#fpaths <- lapply(files_ucla, function(file) paste(csvdir_ucla,file,sep="/"))
## clean names
#fnames <- gsub(".csv","",files_ucla)
#fnames <- gsub("Re22Q_","",fnames)
#fnames <- gsub("Form_","",fnames)
#fnames <- gsub("Qry_","",fnames)
## read all, set to na: "-9999", "-9998","." 
#input_all_ucla <- lapply(fpaths, read.csv, header=T, na.strings=c(".","-9999","-9998"), strip.white=T, sep=",")
#names(input_all_ucla) <- fnames
#df_all_ucla <- lapply(input_all_ucla, function(x) data.frame(x))
#
## filter based on subject lists
## TODO: add prisma 
##ucla_demo <- filter(df_all_ucla$demo_mri, df_all_ucla$demo_mri$MRI_S_ID %in% c(sessions_trio,sessions_prisma))
#ucla_demo <- filter(df_all_ucla$demo_mri, df_all_ucla$demo_mri$MRI_S_ID %in% all_22qtrio)
#
## remove "FAMILY MEMBER" designation from subject identity
#ucla_demo$SUBJECT_IDENTITY <- ucla_demo$SUBJECT_IDENTITY %>% sub("FAMILY MEMBER","",.) %>% sub(",","",.) %>% trimws(which="both") %>% as.factor
## change sex coding from 0/1 to F/M and set to factor
#ucla_demo$SEX <- factor(ucla_demo$SEX,levels=c(0,1),labels=c("F","M"))
## set race=NA to 7 (unknown)
#ucla_demo$RACE[is.na(ucla_demo$RACE)] <- 7
## set race as factor 1=American Indian/Alaska Native; 2=Asian; 3=Native Hawaiian/Pacific Islander; 4=Black or African American; 5=White; 6=Multiple; 7=Unknown
#ucla_demo$RACE <- factor(ucla_demo$RACE,levels=c(1:7),labels=c("1_Native_American","2_Asian","3_Pacific_Island","4_Black","5_White","6_Multiple","7_Unknown"))
## ethnicity as factor with 0=N 1=Y
#ucla_demo$HISPANIC[is.na(ucla_demo$HISPANIC)] <- "Unknown"
#ucla_demo$HISPANIC <- factor(ucla_demo$HISPANIC,levels=c(0,1,"Unknown"),labels=c("N","Y","Unknown"))
## get more accurate age with AGEMONTH/12
#ucla_demo$AGE <- as.numeric(ucla_demo$AGEMONTH)/12 
#
## function to add column to code timepoints relative to sample used (i.e. if visit 1 and 1.12 missing, then 1.24 is baseline)
## trio/prisma coded as T/P-visit_n where T-1 would be the subject's first trio scan and P-1 the first prisma, P-2 the second...
## function should be applied to the indicies of rows (r) in a subset of demo_mri
#gettp <- function(r, df){
#  sub <- df$SUBJECTID[[r]]
#  visit <- df$CONVERTEDVISITNUM[[r]]
#  all_visits <- df$CONVERTEDVISITNUM[which(df$SUBJECTID == sub)] %>% sort
#  n_visits <- length(all_visits)
#  nt_visits <-length(which(all_visits < 2))
#  np_visits <- length(which(all_visits >= 2))
#  visit_index <- which(all_visits == visit)
#  if (visit < 2){
#    label=paste("T-",visit_index,sep="")
#  }else if (visit >= 2){
#    p_visits <- all_visits[which(all_visits >= 2)] %>% sort
#    p_visit_index <- which(p_visits == visit)
#    label=paste("P-",p_visit_index,sep="")
#  }
#  return(c(sub,visit,label,n_visits,nt_visits,np_visits,visit_index))
#}
#
## get timepoints
#timepoints <- sapply(1:nrow(ucla_demo),function(r) gettp(r,ucla_demo)) %>% t %>% as.data.frame
#colnames(timepoints) <- c("SUBJECTID","CONVERTEDVISITNUM","converted_timepoint","n_timepoints","n_trio","n_prisma","visit_index")
#ucla_demo_tp <- cbind(ucla_demo,timepoints[,3:7])
#ucla_demo_tp$visit_index %<>% as.factor
#
## add medication info from summPsych
#ucla_summPsych <- df_all_ucla$summPsych
#ucla_summPsych$PSYTYPE[is.na(ucla_summPsych$PSYTYPE)] <- 5
#ucla_summPsych$medication <- factor(ucla_summPsych$PSYTYPE, levels=c(1,2,3,4,5), labels=c("antipsychotic","antidepressant_or_mood_stabilizer","stimulant","other","none"))
#ucla_summPsych$apd_tf <- ucla_summPsych$medication == "antipsychotic"
#ucla_summPsych$Med_Antipsychotic <- factor(ucla_summPsych$apd_tf, levels=c(T,F), labels=c("Y","N"))
## get psych dx 
#ucla_summPsych$psych_dx <- factor(ucla_summPsych$PSYDIAGNOS, levels=c(1,0,3), labels=c("Y","N","N"))
#
## merge with demo
#ucla_demo_use <- merge(x=ucla_demo_tp, y=ucla_summPsych[,c("SUBJECTID","CONVERTEDVISITNUM","Med_Antipsychotic","psych_dx")], by=c("SUBJECTID","CONVERTEDVISITNUM"), all.x=T)
#
## get IQ
## WASI, WISC-IV, DKEFS and trail making all under df_all_ucla$DKEFS for trio data
## IQSS -- full scale WASI
#ucla_neuro1 <- df_all_ucla$DKEFS[,c("SUBJECTID","CONVERTEDVISITNUM","VOCASS","MATRIXSS","IQSS")] %>% rename("WASI_verbal" = "VOCASS") %>% rename("WASI_matrix" = "MATRIXSS") %>% rename("IQ_full" = "IQSS")
## renewal neuro (prisma) under df_all_ucla$neurocogTest
#ucla_neuro2 <- df_all_ucla$neurocogTest[,c("SUBJECTID","CONVERTEDVISITNUM","VOCA_TSCORE","MATRIX_TSCORE","IQ_SCORE")] %>% rename("WASI_verbal" = "VOCA_TSCORE") %>% rename("WASI_matrix" = "MATRIX_TSCORE") %>% rename("IQ_full" = "IQ_SCORE")
## combine 22q orig and renewal scores before merging with demo table
#ucla_neuro <- rbind(ucla_neuro1, ucla_neuro2)
## merge neuro with demo table
#ucla_demo_use <- merge(x=ucla_demo_use, y=ucla_neuro[,c("SUBJECTID","CONVERTEDVISITNUM","IQ_full","WASI_verbal","WASI_matrix")], by=c("SUBJECTID","CONVERTEDVISITNUM"), all.x=T) 
## record IQ instrument
#ucla_demo_use$IQ_measure <- NA
#ucla_demo_use$IQ_measure[!is.na(ucla_demo_use$IQ_full)] <- "WASI_full_scale"
#
## manually fix missing sex for Q_0381_09102019
#ucla_demo_use[which(ucla_demo_use$MRI_S_ID == "Q_0381_09102019"),"SEX"] <- "F"
#
#
## subset to hcs del
#ucla_demo_use_hcs_del <- ucla_demo_use %>% filter(SUBJECT_IDENTITY=="CONTROL" | SUBJECT_IDENTITY =="PATIENT-DEL")
#
## ugly lines below are to do some convoluted subsetting for age matching for a cross-sectional analysis
## want to take the older timepoints for some of the younger controls
## first get hcs+del scans between age 13 and 30 and use trio timepoint 1
#ucla_demo_use_hcs_del_g13u30_t1 <- ucla_demo_use_hcs_del %>% filter(AGE >= 13 & AGE < 30 & converted_timepoint == "T-1")
## get hcs+del trio scans between age 7 and 13
##ucla_demo_use_hcs_del_g7u13 <- ucla_demo_use_hcs_del %>% filter(AGE >= 7 & AGE < 13 & SCANNER == "trio")
#ucla_demo_use_hcs_del_g7u13 <- ucla_demo_use_hcs_del %>% filter(AGE >= 7 & AGE < 13)
## take timepoint 1 for del 7-13yo scans
#ucla_demo_use_del_g7u13_t1 <- ucla_demo_use_hcs_del_g7u13 %>% filter(converted_timepoint == "T-1",SUBJECT_IDENTITY=="PATIENT-DEL")
## subset hcs
#ucla_demo_use_hcs_g7u13 <- ucla_demo_use_hcs_del_g7u13 %>% filter(SUBJECT_IDENTITY=="CONTROL")
## take timepoint 2 for hcs 7-13yo scans if available
#ucla_demo_use_hcs_g7u13_t2 <- ucla_demo_use_hcs_g7u13 %>% filter(converted_timepoint == "T-2")
## take timepoint 1 for the rest of the 7-13yo hcs
#ucla_demo_use_hcs_g7u13_t1 <- ucla_demo_use_hcs_g7u13 %>% filter(converted_timepoint == "T-1" & !SUBJECTID %in% ucla_demo_use_hcs_g7u13_t2$SUBJECTID)
## combine 13-30 year old both groups T-1 with the 7-13yo deletion T-1 scans and the 7-13yo hcs scans using T-2 when available and T-1 otherwise 
#ucla_demo_use_hcs_del_g7u30_xs <- rbind(ucla_demo_use_hcs_del_g13u30_t1,ucla_demo_use_del_g7u13_t1,ucla_demo_use_hcs_g7u13_t2,ucla_demo_use_hcs_g7u13_t1)
#
## pre-matching demographics summary:
#demo_summary <- CreateTableOne(data=ucla_demo_use_hcs_del_g7u30_xs,vars=c("AGE","SEX"),strata="SUBJECT_IDENTITY",addOverall=F)
#print(demo_summary)
#
#
## need to remove some of the 7-8yo controls to match ages, will order them by %udvarsme and remove worst until ages match
#hcs_7yo <- ucla_demo_use_hcs_del_g7u30_xs %>% filter(AGE < 8 & SUBJECT_IDENTITY == "CONTROL")
#hcs_7yo_mov <- percent_udvarsme_all %>% filter(sesh %in% hcs_7yo$MRI_S_ID)
#hcs_7yo_mov_ordered <- hcs_7yo_mov[order(-hcs_7yo_mov$percent_udvarsme),]
## currently removing the n=4 worst7 year olds
#hcs_7yo_remove <- hcs_7yo_mov_ordered[1:4,1]
#hcs_8yo <- ucla_demo_use_hcs_del_g7u30_xs %>% filter(AGE >= 8 & AGE < 9 & SUBJECT_IDENTITY == "CONTROL")
#hcs_8yo_mov <- percent_udvarsme_all %>% filter(sesh %in% hcs_8yo$MRI_S_ID)
#hcs_8yo_mov_ordered <- hcs_8yo_mov[order(-hcs_8yo_mov$percent_udvarsme),]
## remove the n=1 worst 8 year olds
#hcs_8yo_remove <- hcs_8yo_mov_ordered[1,1]
#hcs_remove <- c(hcs_7yo_remove,hcs_8yo_remove)
#
## filter out subjects in hcs_remove
#ucla_demo_use_hcs_del_xs <- ucla_demo_use_hcs_del_g7u30_xs %>% filter(!MRI_S_ID %in% hcs_remove)
#
## final demo summary
#demo_match_summary <- CreateTableOne(data=ucla_demo_use_hcs_del_xs,vars=c("AGE","SEX"),strata="SUBJECT_IDENTITY",addOverall=F)
#print(demo_match_summary)
#
## final demo table (need to add handedness, medications, psych dx, IQ, SIPS, BOLD movement)
##dir <- "/Users/charlie/Dropbox/PhD/bearden_lab/22q/analyses/striatum_thalamus_fc"
#
## get variables from demo_mri
#df_demo_table_full <- ucla_demo_use_hcs_del_xs[,c("SUBJECTID","CONVERTEDVISITNUM","MRI_S_ID","SUBJECT_IDENTITY","AGE","SEX","EDUDAD","EDUMOM","EDUYEARS","IQ_measure","IQ_full","WASI_verbal","WASI_matrix")]
#
#### hand
## get handedness item scores coded in sistat as 1=L, 2=R, 3=either, 0=no experience
#edin <- df_all_ucla$edin[,c("SUBJECTID","CONVERTEDVISITNUM","EDIN1","EDIN2","EDIN3","EDIN4","EDIN5","EDIN6","EDIN7","EDIN8","EDIN9","EDIN10")]
## function to get total edinburgh score and handedness
## formula is 100*(R-L)/(R+L). score < -40 means left handed, score > 40 right handed
## if more than 2 items NA then score is NA
#get_hand <- function(edin){
#  sub <- edin[c("SUBJECTID")]
#  visit <- edin[c("CONVERTEDVISITNUM")]
#  scores <- edin[c("EDIN1","EDIN2","EDIN3","EDIN4","EDIN5","EDIN6","EDIN7","EDIN8","EDIN9","EDIN10")]
#  l <- sum(scores == 1, na.rm=TRUE)
#  r <- sum(scores == 2, na.rm=TRUE)
#  na <- sum(is.na(scores))
#  if (na < 3){
#    score <- 10*(r-l)
#  }
#  if (na > 2){
#    hand <- NA
#    score <- NA
#  }else if(score > 40){
#    hand <-"R"
#  }else if (score < -40){
#    hand <- "L"
#  }else if (score >= -40 & score <= 40) {
#    hand <- "A"
#  }else{
#    hand <- NA
#  }
#  output <- cbind(sub,visit,score,hand) %>% as.data.frame
#  colnames(output) <- c("SUBJECTID","CONVERTEDVISITNUM","hand_score","hand")
#  return(output)
#}
## get handedness
#edin_result <- lapply(1:nrow(edin), function(r) get_hand(edin[r,])) %>% do.call(rbind,.) %>% as.data.frame
#
## merge handedness with demo table
#df_demo_table_full <- merge(x=df_demo_table_full, y=edin_result[c("SUBJECTID","CONVERTEDVISITNUM","hand")], by=c("SUBJECTID","CONVERTEDVISITNUM"), all.x=T)
#
#### psych dx
## first get SCID columns with Dx (currently only using patient Dx not collateral)
#scid_dx_all <- df_all_ucla$SCID[,c("PATCODE1","PATCODE2","PATCODE3","PATCODE4","PATCODE5","PATCODE6","PATCODE7","PATCODE8")]
#
## get list of unique dx entries
##dx_unique <- scid_dx_all %>% as.matrix %>% as.vector %>% sort %>% unique
#
## create matching key between unique dx and dx groups for demographics table
## first save dx_unique as csv
##write.table(dx_unique, file=file.path(csvdir_ucla,"scid_unique_dx.csv"), row.names=F, col.names=F)
## then manually edit csv so that column 2 contains the dx group for each specific dx. save edited csv as scid_unique_dx_matching.csv
## dx group categories based on DSM-5 https://www.psychiatry.org/File%20Library/Psychiatrists/Practice/DSM/APA_DSM-5-Contents.pdf
## Notes: leave second column blank for non-psych dx (eg Crohn's), code single-episode MDD in full remission as depressive_disorder_past, all other MDD as depressive_disorder
## read matching table back in 
#dx_unique_matching <- read.csv(file.path(csvdir_ucla,"scid_unique_dx_matching.csv"), header=F)
#
## function to take scid patient codes 1-8 for a subject and output binary y/n for each dx in dx_groups based on dx_unique_matching
## should be applied to rows of the scid data frame
#get_general_dx_scid <- function(scid_row,dx_matching){
#  # get subject id and visit columns
#  id_cols <- scid_row[c("SUBJECTID","CONVERTEDVISITNUM")]
#  # get list of all unique dx groups in matching key
#  dx_groups <- dx_matching[,2] %>% sort %>% unique
#  dx_groups <- dx_groups[dx_groups != ""]
#  # get patcodes 1-8
#  patcodes_all <- scid_row[c("PATCODE1","PATCODE2","PATCODE3","PATCODE4","PATCODE5","PATCODE6","PATCODE7","PATCODE8")] %>% as.matrix
#  patcodes <- patcodes_all[patcodes_all != ""]
#  # if subject has data in patcodes, convert to dx groups
#  if(length(patcodes) > 0){
#    # get dx group for each patcode by referencing dx_matching
#    sub_dx <- lapply(patcodes, function(x) filter(dx_matching, V1 == x)$V2) %>% do.call(cbind,.) %>% as.matrix
#    sub_dx <- sub_dx[sub_dx != ""]
#    # check if subject has SCID dx in each dx group, return TRUE for yes, FALSE for no
#    dx_yn <- lapply(dx_groups, function(x) x %in% sub_dx) %>% do.call(cbind,.) %>% as.data.frame
#    # return empty columns if no patcode data (without this will fail for subjects with no data)
#  }else{
#    dx_yn <- matrix(nrow=1,ncol=length(dx_groups)) %>% as.data.frame 
#  }
#  colnames(dx_yn) <- paste("SCID", dx_groups, sep="_")
#  output <- cbind(id_cols, dx_yn)
#  return(output)
#}
#
## get general dx for each scid entry
#scid_general <- lapply(1:nrow(df_all_ucla$SCID), function(r) get_general_dx_scid(scid_row=df_all_ucla$SCID[r,], dx_matching=dx_unique_matching)) %>% do.call(rbind,.) %>% as.data.frame
#
## merge scid general with demo table
#df_demo_table_full <- merge(x=df_demo_table_full, y=scid_general, by=c("SUBJECTID","CONVERTEDVISITNUM"), all.x=T)
#
## count instances of each dx
#dx_counts <- df_demo_table_full %>% dplyr::select(starts_with("SCID_")) %>% colSums(na.rm=T)
#
## get list of dx with more than 2 instances in the data set
#dx_use <- which(dx_counts > 2) %>% names
#
## remove depressive_disorder_past (single episode full remission)
#dx_use <- dx_use[dx_use != "SCID_Depression_Related_Past"]
## remove learning disorder
#dx_use <- dx_use[dx_use != "SCID_Learning_Disorder"]
#
## add info from summPsych
#summpsych <- df_all_ucla$summPsych
#
## meds as factors
#summpsych$PSYTYPE <- factor(summpsych$PSYTYPE, levels=c(1,2,3,4,5), labels=c("antipsychotic","antidepressant_or_mood_stabilizer","stimulant","other","none"))
#
## merge meds with demo table
#df_demo_table_full <- merge(x=df_demo_table_full, y=summpsych[,c("SUBJECTID","CONVERTEDVISITNUM","PSYTYPE")], by=c("SUBJECTID","CONVERTEDVISITNUM"), all.x=T) %>% rename("psych_meds" = "PSYTYPE")
#
## function to mark subject as ASD positive if positive at any visit in summPsych
#get_asd <- function(subject, summ_psych){
#  asd_all <- filter(summ_psych, summ_psych$SUBJECTID==subject)$ASDDIAGNOS
#  # check if any visit coded as 1 (meaning asd=yes)
#  asd_yn <- 1 %in% asd_all
#  return(asd_yn)
#}
#
## add ASD column based on summPsych
#asd_col <- lapply(1:nrow(df_demo_table_full), function(r) get_asd(subject=df_demo_table_full[r,"SUBJECTID"], summ_psych=summpsych)) %>% do.call(rbind,.) %>% as.data.frame
#colnames(asd_col) <- "summPsych_ASD"
#
## merge summPsych ASD with demo table
#df_demo_table_full <- cbind(df_demo_table_full,asd_col)
#
## remove SCID_ASD column, redundant with summPsych
#dx_use <- dx_use[dx_use != "SCID_ASD"]
#
## add percent udvarsme
#df_demo_table_full <- merge(x=df_demo_table_full, y=percent_udvarsme_all[,c("sesh","percent_udvarsme")], by.x="MRI_S_ID", by.y="sesh", all.x=T) %>% rename("percent_BOLD_scrubbed" = "percent_udvarsme")
#
## function to get psychosis status from SIPS
## to be applied to the row indices of a demographics df, and also given the SIPS df
#get_sips <- function(r,demo,sips){
#  sub <- demo$SUBJECTID[[r]]
#  visit <- demo$CONVERTEDVISITNUM[[r]]
#  df_out <- cbind(sub,visit) %>% as.data.frame
#  colnames(df_out) <- c("SUBJECTID","CONVERTEDVISITNUM")
#  sips_sesh <- sips %>% filter(SUBJECTID == sub & CONVERTEDVISITNUM == visit)
#  if(nrow(sips_sesh) < 1){
#    # if no match for sub+visit in sips table, set outputs to NA
#    df_out[,c("SIPS_p_sum","SIPS_n_sum","SIPS_d_sum","SIPS_g_sum","SIPS_psychosis_6","SIPS_psspectrum_3")] <- rep(NA,times=6)
#  }else if(nrow(sips_sesh) > 1){
#    # if more than one match for sub+visit, note error
#    df_out[,c("SIPS_p_sum","SIPS_n_sum","SIPS_d_sum","SIPS_g_sum","SIPS_total", "SIPS_psychosis_6","SIPS_psspectrum_3")] <- rep("ERROR-duplicates",times=7)
#  }else{
#    # get SIPS P scores
#    sips_p_scores <- sips_sesh[c("P1SEV","P2SEV","P3SEV","P4SEV","P5SEV")]
#    # sum SIPS P
#    df_out[,"SIPS_p_sum"] <- sum(sips_p_scores)
#    # get SIPS N
#    sips_n_scores <- sips_sesh[c("N1SEV","N2SEV","N3SEV","N4SEV","N5SEV","N6SEV")]
#    df_out[,"SIPS_n_sum"] <- sum(sips_n_scores)
#    # get SIPS D
#    sips_d_scores <- sips_sesh[c("D1SEV","D2SEV","D3SEV","D4SEV")]
#    df_out[,"SIPS_d_sum"] <- sum(sips_d_scores)
#    # get SIPS G
#    sips_g_scores <- sips_sesh[c("G1SEV","G2SEV","G3SEV","G4SEV")]
#    df_out[,"SIPS_g_sum"] <- sum(sips_g_scores)
#    # get SIPS total
#    df_out["SIPS_total"] <- (df_out["SIPS_p_sum"]+df_out["SIPS_n_sum"]+df_out["SIPS_d_sum"]+df_out["SIPS_g_sum"] )
#    # check psychosis criteria of at least one SIPS P score of 6
#    count_6 <- length(which(sips_p_scores == 6))
#    if(is.na(sum(sips_p_scores))){
#      df_out[,"SIPS_psychosis_6"] <- NA
#    }else if(count_6 > 0){
#      df_out[,"SIPS_psychosis_6"] <- 1
#    }else{
#      df_out[,"SIPS_psychosis_6"] <- 0
#    }
#    # check psychosis-spectrum criteria of at least one SIPS P >= 3
#    count_3 <- length(which(sips_p_scores >= 3))
#    if(is.na(sum(sips_p_scores))){
#      df_out[,"SIPS_psspectrum_3"] <- NA
#    }else if(count_3 > 0){
#      df_out[,"SIPS_psspectrum_3"] <- 1
#    }else{
#      df_out[,"SIPS_psspectrum_3"] <- 0
#    }
#  }
#  return(df_out)
#}
#
## get sips 
#demo_table_sips <- lapply(1:nrow(df_demo_table_full), function(r) get_sips(r=r, demo=df_demo_table_full, sips=df_all_ucla$SIPS)) %>% do.call(rbind,.)
#
## merge sips with demo table
#df_demo_table_full <- merge(x=df_demo_table_full, y=demo_table_sips[,c("SUBJECTID","CONVERTEDVISITNUM","SIPS_total","SIPS_psspectrum_3")], by=c("SUBJECTID","CONVERTEDVISITNUM"), all.x=T) %>% rename("SIPS_prodromal" = "SIPS_psspectrum_3")
#
## set sips_total to numeric and sips_prodromal to logical
#df_demo_table_full %<>% mutate_at(vars("SIPS_prodromal"), ~as.logical(.))
#df_demo_table_full %<>% mutate_at(vars("SIPS_total"), ~as.numeric(.))
#
#
## get IQ
## WASI, WISC-IV, DKEFS and trail making all under df_all_ucla$DKEFS
## VOCASS -- vocab standard score
## MATRIXSS -- matrix standard score
## IQSS -- full scale WASI
## WISC4SS -- WM standard score (missing too many to be useful though)
##neuro <- df_all_ucla$DKEFS[,c("SUBJECTID","CONVERTEDVISITNUM","VOCASS","MATRIXSS","IQSS")] %>% rename("WASI_verbal" = "VOCASS") %>% rename("WASI_matrix" = "MATRIXSS") %>% rename("WASI_IQ_full" = "IQSS")
#
## merge neuro with demo table
##df_demo_table_full <- merge(x=df_demo_table_full, y=neuro, by=c("SUBJECTID","CONVERTEDVISITNUM"), all.x=T) 
#
### make final demographics table to report
## drop unused factor levels
#df_demo_table_full %<>% droplevels
## vector of variables for demo table
#vars_use <- c("AGE","SEX","EDUDAD","EDUMOM","EDUYEARS","hand", "percent_BOLD_scrubbed", "IQ_full","WASI_verbal","WASI_matrix", "SIPS_total","SIPS_prodromal", dx_use, "summPsych_ASD", "psych_meds")
## create table
#demo_match_final <- CreateTableOne(data=df_demo_table_full,vars=vars_use,strata="SUBJECT_IDENTITY",addOverall=F,includeNA=F)
#print(demo_match_final)
## export tableone
##export_demo_table <- print(demo_match_final, quote=F, noSpaces=F, printToggle=T)
##write.csv(export_demo_table, file=file.path(project,"demographics_table1.csv"))
#
#


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
rsfa_name <- "resting_parcellated_voxel_RSFA_Atlas_s_hpss_res-mVWM1d_lpss_whole_brain_CABNP.csv"
rsfa_name_prisma <- "restingAP_parcellated_voxel_RSFA_Atlas_s_hpss_res-mVWM1d_lpss_whole_brain_CABNP.csv"



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
combat_covars_ordered <- merge(x=data.frame(subnames), y=demo_multisite_bl[,c("MRI_S_ID","AGE","SEX","Site","percent_udvarsme")], by.x="subnames", by.y="MRI_S_ID", all.x=TRUE)

# run neurocombat for parcellated RSFA with site as batch variable and AGE and SEX included in the model
# using default parametric model (this paper shows parametric prior estimates outperform non-parametric in neurocombat: https://www.sciencedirect.com/science/article/pii/S2666956022000605)
rsfa_all_combat <- neuroCombat(dat=rsfa_all_subcols[,subcols], batch=combat_covars_ordered$Site, mod=model.matrix(~combat_covars_ordered$AGE+combat_covars_ordered$SEX+combat_covars_ordered$percent_udvarsme), parametric=TRUE)

# get mean of all parcels per subject pre/post combat to look at distributions
parc_means_postcombat <- apply(rsfa_all_combat$dat.combat,  MARGIN=2, FUN="mean")
parc_means_precombat <- apply(rsfa_all_combat$dat.original,  MARGIN=2, FUN="mean")
parc_means <- data.frame(MRI_S_ID=names(parc_means_postcombat), postcombat=as.vector(parc_means_postcombat), precombat=as.vector(parc_means_precombat))
parc_means_demo <- merge(x=parc_means, y=demo_multisite_bl, by="MRI_S_ID", all.x=TRUE, all.y=FALSE)
parc_means_demo$Site <- factor(parc_means_demo$Site, levels=c("UCLAtrio","SUNY","Rome","IoP","UCLAprisma"))


# plot meean RSFA distributions by site pre combat
pl_precombat <- ggplot(parc_means_demo, aes(precombat, fill=Site, y=..count..))+
  geom_density(kernel="gaussian", alpha=0.7)+
  #scale_fill_manual(values=c("lightblue","red"))+
  theme_classic()+
  ggtitle("Mean RSFA pre-combat")

# plot meean RSFA distributions by site post combat
pl_postcombat <- ggplot(parc_means_demo, aes(postcombat, fill=Site, y=..count..))+
  geom_density(kernel="gaussian", alpha=0.7)+
  #scale_fill_manual(values=c("lightblue","red"))+
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
#get_parc_combat_hc_mean_sd <- function(df, hcs_ids, parc){
#  dff <- filter(df, MRI_S_ID %in% hcs_ids & INDEX == parc)
#  out <- data.frame(parcel_hc_mean=mean(as.numeric(dff$t_sd_combat)), parcel_hc_sd=sd(as.numeric(dff$t_sd_combat)))
#  return(out)
#}

get_parc_combat_hc_mean_sd <- function(df, hcs_ids, pt_ids, parc){
  dfh <- filter(df, MRI_S_ID %in% hcs_ids & INDEX == parc)
  dfp <- filter(df, MRI_S_ID %in% pt_ids & INDEX == parc)
  out <- data.frame(parcel_hc_mean=mean(as.numeric(dfh$t_sd_combat)), parcel_hc_sd=sd(as.numeric(dfh$t_sd_combat)),parcel_pt_mean=mean(as.numeric(dfp$t_sd_combat)), parcel_pt_sd=sd(as.numeric(dfp$t_sd_combat)))
  return(out)
}

# first get hc stats per parcel
use_ids_hcs <- filter(demo_multisite_bl, SUBJECT_IDENTITY=="CONTROL")$MRI_S_ID
use_ids_del <- filter(demo_multisite_bl, SUBJECT_IDENTITY=="PATIENT-DEL")$MRI_S_ID
#parc_hc_stats <- lapply(unique(rsfa_all_combat_l_full$INDEX), function(p) get_parc_combat_hc_mean_sd(parc=p, df=rsfa_all_combat_l_full, hcs_ids=use_ids_hcs)) 
#parc_hc_stats <- lapply(unique(rsfa_all_combat_l_full$INDEX), function(p) get_parc_combat_hc_mean_sd(parc=p, df=rsfa_all_combat_l_full, hcs_ids=use_ids_hcs)) %>% do.call(rbind,.)
parc_hc_stats <- lapply(unique(rsfa_all_combat_l_full$INDEX), function(p) get_parc_combat_hc_mean_sd(parc=p, df=rsfa_all_combat_l_full, hcs_ids=use_ids_hcs, pt_ids=use_ids_del)) %>% do.call(rbind,.)
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

# add group means to key
ji_key <- merge(x=ji_key, y=parc_hc_stats, by="INDEX", all.x=TRUE)

# write results to csv
#write.csv(ji_key, file=file.path(project,"cabnp_key_rsfa.csv"),row.names = FALSE)

# function to permute group labels to generate null distribution
# TO DO: not sure if this is properly accounting for covars (see https://brainder.org/tag/permutation-test/ )
# TO DO: first, check parametric p-vals from t-tests above
#permute_group <- function(df,group_col){
#  # take input df and randomly permute group column
#  dfPerm <- df
#  dfPerm[,group_col] <- sample(dfPerm[,group_col])
#  return(dfPerm)
#}
#
## get permutations of the linear model output for each parcel
#nPerm <- 10000
#get_parcel_group_b_lm_v <- function(i){
#  print(paste0("perm:",i,"/",nPerm))
#  get_parcel_group_b_lm(df=permute_group(df=df_demo_table_full_mri,group_col="SUBJECT_IDENTITY"), parc_cols=parc_cols)
#}

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
view_xifti_surface(plot_rsfa_diff_beta, title="RSFA diff beta (22qDel - HCS)", cex.title=1.3,zlim=c(-0.8,0.8),colors=pal_red_blue())
#view_xifti_volume(plot_rsfa_diff_beta, title="RSFA diff", cex.title=1.3,zlim=c(-0.8,0.8),colors=pal_red_blue(),crop=FALSE, orientation_labels=FALSE)


# group difference beta only fdr<0.05
plot_rsfa_diff_beta_fdr <- atlas_xifti_new_vals(xii=xii_Ji_parcel, df=ji_key, roi_col="INDEX", val_col="rsfa_diff_beta_parasig")
view_xifti_surface(plot_rsfa_diff_beta_fdr, title="RSFA diff beta FDR corrected (22qDel - HCS)", cex.title=1.3,zlim=c(-0.8,0.8), colors=pal_red_blue())
#view_xifti_volume(plot_rsfa_diff_beta_fdr, title="RSFA diff", cex.title=1.3,zlim=c(-0.8,0.8), colors=pal_red_blue(),crop=FALSE, orientation_labels=FALSE)

# save fdr corrected maps 
#view_xifti_surface(plot_rsfa_diff_beta_fdr, title="RSFA diff beta FDR corrected (22qDel - HCS)", cex.title=1.3,zlim=c(-0.8,0.8), colors=pal_red_blue(), fname=file.path(project,"figures/rsfa/rsfa_fdr_surf.png"), width=1200, height=800)
#view_xifti_volume(plot_rsfa_diff_beta_fdr, title="RSFA diff", cex.title=1.3,zlim=c(-0.8,0.8),slices=seq(16,42, by=5), colors=pal_red_blue(),crop=FALSE, orientation_labels=FALSE, legend_embed=FALSE, fname=file.path(project,"figures/rsfa/rsfa_fdr_vol.pdf"), format="pdf")

# control means
plot_rsfa_hc <- atlas_xifti_new_vals(xii=xii_Ji_parcel, df=ji_key, roi_col="INDEX", val_col="parcel_hc_mean")
view_xifti_surface(plot_rsfa_hc, title="RSFA Control Means", zlim=c(0,25),cex.title=1.3, colors="magma")

# 22qdel means
plot_rsfa_del <- atlas_xifti_new_vals(xii=xii_Ji_parcel, df=ji_key, roi_col="INDEX", val_col="parcel_pt_mean")
view_xifti_surface(plot_rsfa_del, title="RSFA 22qDel Means", zlim=c(0,25),cex.title=1.3, colors="magma")


## get info for SOBP abstract
# counts and sex
demo_sobp <- CreateTableOne(data=demo_multisite_bl,vars=c("AGE","SEX"),strata="SUBJECT_IDENTITY",addOverall=F)
demo_sobp <- CreateTableOne(data=demo_multisite_bl,vars=c("AGE","SEX","IQ_full","psych_dx","Med_Antipsychotic","percent_udvarsme"),strata="SUBJECT_IDENTITY",addOverall=F)
print(demo_sobp, showAllLevels = TRUE)
# number of significant regions
filter(ji_key, !is.na(rsfa_diff_beta_parasig)) %>% nrow
# decreases network summary cortex
filter(ji_key[1:360,], rsfa_diff_beta_parasig < 0)$NETWORK %>% as.factor %>% summary
# decreases network summary subcortex
filter(ji_key[360:718,], rsfa_diff_beta_parasig < 0)$NETWORK %>% as.factor %>% summary

# increases network summary cortex
filter(ji_key[1:360,], rsfa_diff_beta_parasig > 0)$NETWORK %>% as.factor %>% summary
# increases network summary subcortex
filter(ji_key[360:718,], rsfa_diff_beta_parasig > 0)$NETWORK %>% as.factor %>% summary


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
rsfa_ahba <- merge(x=ji_key, y=ahbaCombinedCABNP, by.x="INDEX", by.y="label")

# get set of usable parcels from left cortex by excluding rows where a test gene (e.g. PVALB) is NA
#rsfa_lh_ahba <- filter(rsfa_ahba, INDEX <= 180 & !is.na(PVALB))
rsfa_lh_ahba <- filter(rsfa_ahba, INDEX <= 180)

# correlate delta with PVALB, SST, and PVALB-SST diff
cor_lh_pv <- cor.test(rsfa_lh_ahba$rsfa_diff_beta, rsfa_lh_ahba$PVALB, na.action="omit")
cor_lh_st <- cor.test(rsfa_lh_ahba$rsfa_diff_beta, rsfa_lh_ahba$SST, na.action="omit")
cor_lh_pv_sst_diff <- cor.test(rsfa_lh_ahba$rsfa_diff_beta, rsfa_lh_ahba$PVALB_SST_diff, na.action="omit")

cor_lh_pv
cor_lh_st
cor_lh_pv_sst_diff


# correlate group means with PVALB
cor_hcmean_pv <- cor.test(rsfa_lh_ahba$parcel_hc_mean, rsfa_lh_ahba$PVALB, na.action="omit")
cor_delmean_pv <- cor.test(rsfa_lh_ahba$parcel_pt_mean, rsfa_lh_ahba$PVALB, na.action="omit")
cor_hcmean_pv
cor_delmean_pv

# correlate group means with SST
cor_hcmean_st <- cor.test(rsfa_lh_ahba$parcel_hc_mean, rsfa_lh_ahba$SST, na.action="omit")
cor_delmean_st <- cor.test(rsfa_lh_ahba$parcel_pt_mean, rsfa_lh_ahba$SST, na.action="omit")
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

#gene_cor(df=rsfa_lh_ahba, mricol="rsfa_diff_beta", gene="PVALB")
#gene_cor_list(df=rsfa_lh_ahba, mricol="rsfa_diff_beta", genelist=filter(zcells_use_final, CellType=="Interneuron")$GeneSymbol) %>% do.call(rbind,.) %>% as.vector

zcell_rsfa_distributions <- lapply(zcells_types, function(t) gene_cor_list(df=rsfa_lh_ahba, mricol="rsfa_diff_beta", genelist=filter(zcells_use_final, CellType==t)$GeneSymbol))
names(zcell_rsfa_distributions) <- zcells_types

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

zcell_rsfa_plots <- lapply(zcells_types, function(t) plot_dist(vec=unlist(zcell_rsfa_distributions[[t]]), title=t))

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


# get indices of regions with significantly decreased RSFA
ji_key_neg_rsfa <- filter(ji_key, rsfa_diff_beta_parasig < 0 )
neg_rsfa_parcels <- ji_key_neg_rsfa$INDEX
neg_rsfa_parces_colnames <- paste0("r_",neg_rsfa_parcels)

# get column of the mean normed value of the rois with significant negative 22qDel effect for each subject
df_demo_table_ucla_mri_beh$mean_sig_decrease_rois <- apply(df_demo_table_ucla_mri_beh[,neg_rsfa_parces_colnames],1,mean)

# subset to only 22qDel
df_demo_table_ucla22q_mri_beh <- filter(df_demo_table_ucla_mri_beh, SUBJECT_IDENTITY=="PATIENT-DEL")

# correlate composite negative RSFA value with SIPS
cor.test(x=df_demo_table_ucla22q_mri_beh$mean_sig_decrease_rois, y=df_demo_table_ucla22q_mri_beh$SIPS_p_sum, method="pearson",  na.action="omit")
cor.test(x=df_demo_table_ucla22q_mri_beh$mean_sig_decrease_rois, y=df_demo_table_ucla22q_mri_beh$SIPS_n_sum, method="pearson",  na.action="omit")
cor.test(x=df_demo_table_ucla22q_mri_beh$mean_sig_decrease_rois, y=df_demo_table_ucla22q_mri_beh$SIPS_g_sum, method="pearson",  na.action="omit")
cor.test(x=df_demo_table_ucla22q_mri_beh$mean_sig_decrease_rois, y=df_demo_table_ucla22q_mri_beh$SIPS_d_sum, method="pearson",  na.action="omit")
cor.test(x=df_demo_table_ucla22q_mri_beh$mean_sig_decrease_rois, y=df_demo_table_ucla22q_mri_beh$SIPS_total, method="pearson",  na.action="omit")

# TODO: in 22q, lm(parcel_rsfa ~ SIPS_p)


# linear models for IQ and psychosis dx from demo table

# function to return beta coefficient for group in a lm predicting MRI from group plus covariates
lm_parcel_beh <- function(df,bhvar,parcel){
  dff <- df
  dff$mri <- dff[,parcel]
  form <- as.formula(paste(bhvar,"~ mri + AGE + SEX + Site + percent_udvarsme"))
  summary(lm(formula=form,data=dff))
}

# function to apply lm at each parcel
get_parcel_beh_lm <- function(df,parc_cols,bhvar){
  lapply(parc_cols, function(v) lm_parcel_beh(df=df, bhvar=bhvar, parcel=v))
}

# get del only
df_demo_table_del_mri <- filter(df_demo_table_full_mri, SUBJECT_IDENTITY=="PATIENT-DEL")
df_demo_table_del_mri$psych_dx <- df_demo_table_del_mri$psych_dx %>% as.factor %>% as.numeric

# get effects of parcel on IQ
del_rsfa_iq_lms <- get_parcel_beh_lm(df=filter(df_demo_table_del_mri, SUBJECT_IDENTITY=="PATIENT-DEL"), parc_cols=parc_cols, bhvar="IQ_full") 

# p-vals for mri effect on IQ
del_rsfa_iq_p <- lapply(del_rsfa_iq_lms, function(d)d$coefficients["mri","Pr(>|t|)"]) %>% do.call(rbind,.) %>% as.numeric
del_rsfa_iq_fdrq <- p.adjust(del_rsfa_iq_p, method = "fdr")

# TODO: partial least squares predicting IQ/psychosis from MRI in 22qDel
# https://cran.r-project.org/web/packages/pls/vignettes/pls-manual.pdf
plsform_iq <- as.formula(paste("IQ_full ~",paste(parc_cols, collapse=" + ")))
plsmod_iq <- plsr(formula=plsform_iq, data=df_demo_table_del_mri, scale=TRUE, center=TRUE, validation="CV", ncomp=10)

validationplot(plsmod_iq,val.type = "RMSE", ncomp=1:7, intercept=FALSE, legendpos="topright")
validationplot(plsmod_iq,val.type = "R2", ncomp=1:7, intercept=FALSE, legendpos="topright")
selectNcomp(plsmod_iq,method = "onesigma", plot=TRUE)

plsform_ps <- as.formula(paste("psych_dx ~",paste(parc_cols, collapse=" + ")))
plsmod_ps <- plsr(formula=plsform_ps, data=df_demo_table_del_mri, scale=TRUE, center=TRUE, validation="CV", ncomp=10)


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
 
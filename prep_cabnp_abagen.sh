#!/bin/sh

# script to convert CIFTI to fsaverage and MNI

# set path to workbench command binary
# TODO: edit this path if workbench is installed somewhere else
wbcommand="/Applications/workbench/bin_macosx64/wb_command"
# update PATH if the environment doesn't include wb_command
checkwb=$(which wb_commandd)
if [ -z $checkwb ]; then 
    echo "wb_command path not set, adding to PATH:"
    export PATH=$PATH:${wbcommand}
    echo $PATH
else 
    echo "wb_command path already set: ${check_wb}" 
fi

# get working directory for githup repo
project=$(git rev-parse --show-toplevel)
echo "Project directory is: ${project}"

# help function for wb_command
#wb_command -gifti-help

# make output directory 
out=${project}/CAB-NP/CAB_NP_converted/
mkdir -p $out

# read the CAB-NP CIFTI dlabel file and convert to LH and RH cortex GIFTIs in fsaverage space, and subcortical NIFTI in MNI volume space
label_in=${project}/CAB-NP/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dlabel.nii
label_out_left=${out}/CAB_NP_surface_orig_L.label.gii
label_out_right=${out}/CAB_NP_surface_orig_R.label.gii
echo "...Separating CIFTI"
wb_command -cifti-separate $label_in COLUMN -label CORTEX_LEFT $label_out_left -label CORTEX_RIGHT $label_out_right

# label resample 
for hemi in L R; do
    ## set up vars for -metric-resample to fsaverage5 
    # metric file to resample
    label_in=${out}/CAB_NP_surface_orig_${hemi}.label.gii
    # path to reference surfaces
    #surf_path=/Users/charlie/Dropbox/PhD/bearden_lab/22q/HCPpipelines-4.1.3/global/templates/standard_mesh_atlases/resample_fsaverage/
    surf_path=${project}/CAB-NP/HCPpipelines-4.1.3_templates_resample_fsaverage
    # appropriate hemisphere and resolution from standard_mesh_atlases/resample_fsaverage (e.g.fs_LR-deformed_to-fsaverage.?.sphere.???.k_fs_LR.surf.gii)
    # confirmed 32k mesh by reading metric_out_left as a text file, Dim0="32492"
    current_sphere=${surf_path}/fs_LR-deformed_to-fsaverage.${hemi}.sphere.32k_fs_LR.surf.gii
    # use appropriate hemi and res (e.g fsaverage5_std_sphere.?.???.k_fsavg_?.surf.gii)
    new_sphere=${surf_path}/fsaverage5_std_sphere.${hemi}.10k_fsavg_${hemi}.surf.gii
    # specify output name *func.gii
    label_out=${out}/CAB_NP_surface_resample_fsaverage5_${hemi}.label.gii
    # use midthickness matching current_sphere
    current_area=${surf_path}/fs_LR.${hemi}.midthickness_va_avg.32k_fs_LR.shape.gii
    # use midthickness matching new_sphere
    new_area=${surf_path}/fsaverage5.${hemi}.midthickness_va_avg.10k_fsavg_${hemi}.shape.gii
    # command to resample surfaces to fsaverage5
    echo "...Resampling"
    wb_command -label-resample $label_in $current_sphere $new_sphere ADAP_BARY_AREA $label_out -area-metrics $current_area $new_area
done


# get volume parcels nii from CIFTI
# use dscalar input rather than dlabel, because dscalar subcort roi values correspond to directly to the parcel IDs but the dlabel ones don't
cifti_in=${project}/CAB-NP/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dscalar.nii
#label_in=${project}/CAB-NP/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dlabel.nii
# volume out name
volume_out=${out}/CAB_NP_vol_separated.nii
echo "...Separating volume"
wb_command -cifti-separate $cifti_in COLUMN -volume-all $volume_out
#wb_command -cifti-separate $label_in COLUMN -volume-all $volume_out

##############################################################################################################################
# Unused code (from previous attempt to use dscalar instead of dlabel inputs, lead to blurring of fsaverage borders)

# create volume nii
#wb_command -label-to-volume-mapping 

# separate CAB-NP CIFTI parcels to single-hemisphere GIFTI files
#cifti_in="/Users/charlie/Dropbox/PhD/bearden_lab/22q/analyses/striatum_thalamus_fc/ColeAnticevicNetPartition-master/data/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dscalar.nii"
#metric_out_left="/Users/charlie/Dropbox/PhD/bearden_lab/AHBA/CAB_NP_converted/CAB_NP_surface_orig_L.func.gii"
#metric_out_right="/Users/charlie/Dropbox/PhD/bearden_lab/AHBA/CAB_NP_converted/CAB_NP_surface_orig_R.func.gii"
#wb_command -cifti-separate $cifti_in COLUMN -metric CORTEX_LEFT $metric_out_left -metric CORTEX_RIGHT $metric_out_right
#
#for hemi in L R; do
#    ## set up vars for -metric-resample to fsaverage5 
#    # metric file to resample
#    metric_in=/Users/charlie/Dropbox/PhD/bearden_lab/AHBA/CAB_NP_converted/CAB_NP_surface_orig_${hemi}.func.gii
#    # path to reference surfaces
#    surf_path=/Users/charlie/Dropbox/PhD/bearden_lab/22q/HCPpipelines-4.1.3/global/templates/standard_mesh_atlases/resample_fsaverage/
#    # appropriate hemisphere and resolution from standard_mesh_atlases/resample_fsaverage (e.g.fs_LR-deformed_to-fsaverage.?.sphere.???.k_fs_LR.surf.gii)
#    # confirmed 32k mesh by reading metric_out_left as a text file, Dim0="32492"
#    current_sphere=${surf_path}/fs_LR-deformed_to-fsaverage.${hemi}.sphere.32k_fs_LR.surf.gii
#    # use appropriate hemi and res (e.g fsaverage5_std_sphere.?.???.k_fsavg_?.surf.gii)
#    new_sphere=${surf_path}/fsaverage5_std_sphere.${hemi}.10k_fsavg_${hemi}.surf.gii
#    # specify output name *func.gii
#    metric_out=/Users/charlie/Dropbox/PhD/bearden_lab/AHBA/CAB_NP_converted/CAB_NP_surface_resample_fsaverage5_${hemi}.func.gii
#    # use midthickness matching current_sphere
#    current_area=${surf_path}/fs_LR.${hemi}.midthickness_va_avg.32k_fs_LR.shape.gii
#    # use midthickness matching new_sphere
#    new_area=${surf_path}/fsaverage5.${hemi}.midthickness_va_avg.10k_fsavg_${hemi}.shape.gii
#    # check files exist
#    #ls $metric_in
#    #ls -d $surf_path
#    #ls $current_sphere
#    #ls $new_sphere
#    #ls $current_area
#    #ls $new_area
#    # command to resample surfaces to fsaverage5
#    wb_command -metric-resample $metric_in $current_sphere $new_sphere ADAP_BARY_AREA $metric_out -area-metrics $current_area $new_area
#done
#
#
#
#ls /Users/charlie/Dropbox/PhD/bearden_lab/AHBA/CAB_NP_converted/
#
#cat /Users/charlie/Dropbox/PhD/bearden_lab/AHBA/CAB_NP_converted/CAB_NP_surface_resample_fsaverage5_L.func.gii| grep "Dim0"


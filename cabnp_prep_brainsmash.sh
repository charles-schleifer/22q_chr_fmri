#!/bin/sh
# separate CAB-NP for input to brainSMASH

path="/Users/charlie/Dropbox/github/22q_chr_fmri/CAB-NP/"
input_dlabel="CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dlabel.nii"
out_name_base="CortexParcelLabels"

wb_command -cifti-separate ${path}/${input_dlabel} COLUMN -label CORTEX_LEFT ${path}/${out_name_base}_L.label.gii
wb_command -cifti-create-label ${path}/${out_name_base}_L.dlabel.nii -left-label ${path}/${out_name_base}_L.label.gii

wb_command -cifti-separate ${path}/${input_dlabel} COLUMN -label CORTEX_RIGHT ${path}/${out_name_base}_R.label.gii
wb_command -cifti-create-label ${path}/${out_name_base}_R.dlabel.nii -right-label ${path}/${out_name_base}_R.label.gii
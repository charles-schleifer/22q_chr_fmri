#!/usr/bin/env python
# coding: utf-8

# import abagen and other libraries
import os
import abagen
from abagen import images
from abagen import reporting
import pandas as pd
pd.set_option('display.max_rows', 1000)
pd.set_option('display.max_columns', 500)

# Download AHBA data
# TODO: uncomment if you need to re-download the data
#files = abagen.fetch_microarray(donors='all', verbose=1)

# get path to parent directory (github repo)
# this works if full script is run, otherwise, need to set project manually to github repo path
# project='/Users/charlie/Dropbox/github/22q_spatial_transcriptomic/'
project=os.path.dirname(__file__)+'/'    
    
# try importing CAB-NP surface in fsLR 32k resolution
# abagen normally expects fsaverage5 but can handle other meshes if provided with the surface file
# see "non-standard parcellations" https://abagen.readthedocs.io/en/stable/user_guide/parcellations.html
#from abagen import images
# for atlas, provide paths to left and right surface GIFTIS (created by prep_cabnp_for_abagen_bash.ipynb)
#atlasCABNPsurf = ('/Users/charlie/Dropbox/PhD/bearden_lab/AHBA/CAB_NP_converted/CAB_NP_surface_orig_L.func.gii','/Users/charlie/Dropbox/PhD/bearden_lab/AHBA/CAB_NP_converted/CAB_NP_surface_orig_R.func.gii')
# for surf, provide path
#surfCABNP = ('/Users/charlie/Dropbox/PhD/bearden_lab/22q/HCPpipelines-4.1.3/global/templates/standard_mesh_atlases/L.sphere.32k_fs_LR.surf.gii','/Users/charlie/Dropbox/PhD/bearden_lab/22q/HCPpipelines-4.1.3/global/templates/standard_mesh_atlases/R.sphere.32k_fs_LR.surf.gii')
#atlasCABNPsurfTree = images.check_atlas(atlas=atlasCABNPsurf, geometry=surfCABNP, space='fslr')
# TODO: figure out why this didn't work

# read CAB-NP parcellation key
ji_parc_key = pd.read_table(project+'/CAB-NP/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR_LabelKey.txt', delimiter="\t")

# edit key to use as atlas info
ji_key_subset = ji_parc_key[['INDEX','HEMISPHERE','LABEL','NETWORK','NETWORKKEY','GLASSERLABELNAME']]
# rename columns that abagen expects
info_ji_parc = ji_key_subset.rename(columns={'INDEX':'id', 'HEMISPHERE':'hemisphere'})
# rename LR to B
info_ji_parc.replace('LR','B',regex=True, inplace=True)
# rename to 'cortex', 'cerebellum', 'subcortex/brainstem', ‘white matter’, or ‘other’
structure = info_ji_parc.LABEL.str.split(pat="-").str[-1]
# rename to 'cortex', 'cerebellum', 'subcortex/brainstem', ‘white matter’, or ‘other’
char_to_replace = {'Ctx':'cortex', 'Cerebellum':'cerebellum', 'Accumbens':'subcortex/brainstem', 'Brainstem':'subcortex/brainstem', 'Caudate':'subcortex/brainstem', 'Diencephalon':'subcortex/brainstem', 'Hippocampus':'subcortex/brainstem', 'Pallidum':'subcortex/brainstem', 'Putamen':'subcortex/brainstem', 'Thalamus':'subcortex/brainstem','Amygdala':'subcortex/brainstem'}
for key, value in char_to_replace.items():
    # Replace key character with value character in string
    structure = structure.replace(key, value)
# add to df
info_ji_parc['structure']=structure
# get only cortex
info_ji_parc_cortex=info_ji_parc[info_ji_parc.structure == 'cortex']
# get subcort
info_ji_parc_subcort=info_ji_parc[info_ji_parc.structure != 'cortex']

# try importing CAB-NP surface converted to fsaverage5 resolution using wb_command
# abagen normally expects fsaverage5 but can handle other meshes if provided with the surface file
# see "non-standard parcellations" https://abagen.readthedocs.io/en/stable/user_guide/parcellations.html
# for atlas, provide paths to left and right surface GIFTIS (created by prep_cabnp_for_abagen_bash.ipynb)
atlasCABNPsurf = (project+'/CAB-NP/CAB_NP_converted/CAB_NP_surface_resample_fsaverage5_L.label.gii',project+'/CAB-NP/CAB_NP_converted/CAB_NP_surface_resample_fsaverage5_r.label.gii')

# surface fsaverage5
atlasCABNPsurfTree = images.check_atlas(atlas=atlasCABNPsurf,atlas_info=info_ji_parc_cortex, space='fsaverage5')

# volume (output of -cifti-separate)
atlasCABNPvol=project+'/CAB-NP/CAB_NP_converted/CAB_NP_vol_separated.nii'
atlasCABNPvolTree = images.check_atlas(atlas=atlasCABNPvol,atlas_info=info_ji_parc_subcort)

# get expression for CABNP surface (fsaverage5 space)
# set norm_structures=True to normalize within brain structures (cortex, cerebellum, subcort/brainstem)
# since vol and surf are computed separately here, normalizing by brain structure is better 
# https://abagen.readthedocs.io/en/stable/user_guide/normalization.html
ahbaSurfCABNP = abagen.get_expression_data(atlasCABNPsurfTree, norm_structures=True, verbose=1)
ahbaSurfCABNP.to_csv(project+'/CAB-NP_surface_abagen_expression.csv')
#return_counts
#return_donors
#return_report
# generate and save methods report
surfReport = reporting.Report(atlasCABNPsurfTree, norm_structures=True)
surfReport_out = surfReport.gen_report()
with open(project+'abagen_surf_methods_report.txt', 'w') as text_file:
    text_file.write(surfReport_out)
    
# volume only
ahbaVolCABNP = abagen.get_expression_data(atlasCABNPvolTree, norm_structures=True, verbose=1)
ahbaVolCABNP.to_csv(project+'/CAB-NP_subcort_abagen_expression.csv')
# generate and save methods report
volReport = reporting.Report(atlasCABNPvolTree, norm_structures=True)
volReport_out = volReport.gen_report()
with open(project+'/abagen_vol_methods_report.txt', 'w') as text_file:
    text_file.write(volReport_out)






















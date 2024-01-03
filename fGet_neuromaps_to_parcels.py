#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 10:17:22 2023

@author: kkumar

Aim: Neuromaps parcellated to different cortical parcellations

Summary:
    Step 1: use "neuromaps" to get various maps in fsaverage5 (FreeSurfer) space
    Step 2: use "ENIGMA_toolbox" to parcellate surface
    *NOTE: For Glasser parcellation use Matlab version of ENIGMA Toolbox
    
    
References:
    1) Markello, Ross D., et al. "Neuromaps: structural and functional 
       interpretation of brain maps." Nature Methods 19.11 (2022): 1472-1479.
    2) neuromaps python package: https://netneurolab.github.io/neuromaps/index.html
    3) ENIGMA Toolbox: https://enigma-toolbox.readthedocs.io/en/latest/#
    4) Larivière, S., Paquola, C., Park, By. et al. The ENIGMA Toolbox: 
       multiscale neural contextualization of multisite neuroimaging datasets. 
       Nat Methods 18, 698–700 (2021). https://doi.org/10.1038/s41592-021-01186-4

"""


#-----------------------------------------------------------------------------
# #--- set directory paths
#-----------------------------------------------------------------------------

# Set the project folder, which has sub-folders for "code", "data_ref", "data_ref/neuromaps_aug2023" etc.
base_folder = "/media/kkumar/SSDN_PostDoc_4TB_1/SSDN_KD/SSDN_Postdoc_Project/Analysis/Project_Folders/CrossCNV_CT_VirtualHist_Aug2023"    # project folder
code_path = base_folder + "/code"
neuromaps_download_path = base_folder + "/data_ref/neuromaps_aug2023"
neuromaps_fsavg5_path = base_folder + "/data_ref/neuromaps_aug2023/neuromaps_fsaverage10k"

#-----------------------------------------------------------------------------
# neuromaps list of 20 brainmaps
#-----------------------------------------------------------------------------

import nibabel as nib
import numpy as np
import os


from neuromaps.datasets import available_annotations
for annotation in available_annotations():
    print(annotation)
   
"""
# 20 brain maps from Markello et al.
('abagen', 'genepc1', 'fsaverage', '10k')
('hcps1200', 'megalpha', 'fsLR', '4k')
('hcps1200', 'megbeta', 'fsLR', '4k')
('hcps1200', 'megdelta', 'fsLR', '4k')
('hcps1200', 'meggamma1', 'fsLR', '4k')
('hcps1200', 'meggamma2', 'fsLR', '4k')
('hcps1200', 'megtheta', 'fsLR', '4k')
#  ('hcps1200', 'megtimescale', 'fsLR', '4k')
('hcps1200', 'myelinmap', 'fsLR', '32k')
('hcps1200', 'thickness', 'fsLR', '32k')
('hill2010', 'devexp', 'fsLR', '164k')
('hill2010', 'evoexp', 'fsLR', '164k')
('margulies2016', 'fcgradient01', 'fsLR', '32k')
('mueller2013', 'intersubjvar', 'fsLR', '164k')
('neurosynth', 'cogpc1', 'MNI152', '2mm')
('raichle', 'cbf', 'fsLR', '164k')
('raichle', 'cbv', 'fsLR', '164k')
('raichle', 'cmr02', 'fsLR', '164k')
('raichle', 'cmruglu', 'fsLR', '164k')
('reardon2018', 'scalingnih', 'civet', '41k')
('reardon2018', 'scalingpnc', 'civet', '41k')

"""

#-----------------------------------------------------------------------------
# Step 1: fetch annotations
#-----------------------------------------------------------------------------

from neuromaps.datasets import fetch_annotation
annotations = fetch_annotation(source=['abagen', 'neurosynth'],data_dir= neuromaps_download_path)
annotations = fetch_annotation(source=['hcps1200', 'hill2010', 'margulies2016', 'mueller2013', 'raichle', 'reardon2018'],data_dir= neuromaps_download_path)


#-----------------------------------------------------------------------------
# Step 2: transform to fsaverage 10k coordinate system
#-----------------------------------------------------------------------------

from neuromaps import transforms
import numpy as np

#-----------------------------------------------------------------------------
# i. Neurosynth MNI space to fsavaerage
neurosynth = fetch_annotation(source='neurosynth',data_dir= neuromaps_download_path)
fsavg = transforms.mni152_to_fsaverage(neurosynth, fsavg_density= '10k')
fsavg_lh, fsavg_rh = fsavg
print(fsavg_lh.agg_data().shape)
print(fsavg_rh.agg_data().shape)

# save the Gifti file
nib.save(fsavg_lh, neuromaps_fsavg5_path + '/neurosynth_fs10k_lh.gii')
nib.save(fsavg_rh, neuromaps_fsavg5_path + '/neurosynth_fs10k_rh.gii')

# Saving the array in a text file
fsavg_array = fsavg_lh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/neurosynth_fs10k_lh.txt', fsavg_array, newline="\n")
fsavg_array = fsavg_rh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/neurosynth_fs10k_rh.txt', fsavg_array, newline="\n")

#-----------------------------------------------------------------------------
# ii. abagen to fsav
abagen = fetch_annotation(source='abagen',data_dir= neuromaps_download_path)
abagen_lh = nib.load(abagen[0])
abagen_rh = nib.load(abagen[1])

# abagen_lh:  'D:\\SSDN_KD\\Temp_KD\\Neuromaps_NatMeth_April2023\\annotations\\abagen\\genepc1\\fsaverage\\source-abagen_desc-genepc1_space-fsaverage_den-10k_hemi-L_feature.func.gii'
# save the Gifti file in current folder (for consistency)
nib.save(abagen_lh, neuromaps_fsavg5_path + '/abagen_fs10k_lh.gii')
nib.save(abagen_rh, neuromaps_fsavg5_path + '/abagen_fs10k_rh.gii')

# Saving the array in a text file
fsavg_array = abagen_lh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/abagen_fs10k_lh.txt', fsavg_array, newline="\n")
fsavg_array = abagen_rh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/abagen_fs10k_rh.txt', fsavg_array, newline="\n")

#-----------------------------------------------------------------------------
# iii. hcps1200

# ('hcps1200', 'megalpha', 'fsLR', '4k')
# ('hcps1200', 'megbeta', 'fsLR', '4k')
# ('hcps1200', 'megdelta', 'fsLR', '4k')
# ('hcps1200', 'meggamma1', 'fsLR', '4k')
# ('hcps1200', 'meggamma2', 'fsLR', '4k')
# ('hcps1200', 'megtheta', 'fsLR', '4k')

#-------------
meg_alpha = fetch_annotation(source='hcps1200',desc='megalpha', space='fsLR', den='4k',data_dir= neuromaps_download_path)
fsavg = transforms.fslr_to_fsaverage(meg_alpha, target_density= '10k')
fsavg_lh, fsavg_rh = fsavg
print(fsavg_lh.agg_data().shape)

# save the Gifti file
nib.save(fsavg_lh, neuromaps_fsavg5_path + '/meg_alpha_fs10k_lh.gii')
nib.save(fsavg_rh, neuromaps_fsavg5_path + '/meg_alpha_fs10k_rh.gii')

# Saving the array in a text file
fsavg_array = fsavg_lh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/meg_alpha_fs10k_lh.txt', fsavg_array, newline="\n")
fsavg_array = fsavg_rh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/meg_alpha_fs10k_rh.txt', fsavg_array, newline="\n")


#-------------
meg_beta = fetch_annotation(source='hcps1200',desc='megbeta', space='fsLR', den='4k',data_dir= neuromaps_download_path)
fsavg = transforms.fslr_to_fsaverage(meg_beta, target_density= '10k')
fsavg_lh, fsavg_rh = fsavg
print(fsavg_lh.agg_data().shape)

# save the Gifti file
nib.save(fsavg_lh, neuromaps_fsavg5_path + '/meg_beta_fs10k_lh.gii')
nib.save(fsavg_rh, neuromaps_fsavg5_path + '/meg_beta_fs10k_rh.gii')

# Saving the array in a text file
fsavg_array = fsavg_lh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/meg_beta_fs10k_lh.txt', fsavg_array, newline="\n")
fsavg_array = fsavg_rh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/meg_beta_fs10k_rh.txt', fsavg_array, newline="\n")


#-------------
meg_delta = fetch_annotation(source='hcps1200',desc='megdelta', space='fsLR', den='4k',data_dir= neuromaps_download_path)
fsavg = transforms.fslr_to_fsaverage(meg_delta, target_density= '10k')
fsavg_lh, fsavg_rh = fsavg
print(fsavg_lh.agg_data().shape)
# save the Gifti file
nib.save(fsavg_lh, neuromaps_fsavg5_path + '/meg_delta_fs10k_lh.gii')
nib.save(fsavg_rh, neuromaps_fsavg5_path + '/meg_delta_fs10k_rh.gii')


# Saving the array in a text file
fsavg_array = fsavg_lh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/meg_delta_fs10k_lh.txt', fsavg_array, newline="\n")
fsavg_array = fsavg_rh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/meg_delta_fs10k_rh.txt', fsavg_array, newline="\n")


#-------------
meg_gamma1 = fetch_annotation(source='hcps1200',desc='meggamma1', space='fsLR', den='4k',data_dir= neuromaps_download_path)
fsavg = transforms.fslr_to_fsaverage(meg_gamma1, target_density= '10k')
fsavg_lh, fsavg_rh = fsavg
print(fsavg_lh.agg_data().shape)
# save the Gifti file
nib.save(fsavg_lh, neuromaps_fsavg5_path + '/meg_gamma1_fs10k_lh.gii')
nib.save(fsavg_rh, neuromaps_fsavg5_path + '/meg_gamma1_fs10k_rh.gii')

# Saving the array in a text file
fsavg_array = fsavg_lh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/meg_gamma1_fs10k_lh.txt', fsavg_array, newline="\n")
fsavg_array = fsavg_rh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/meg_gamma1_fs10k_rh.txt', fsavg_array, newline="\n")


#-------------
meg_gamma2 = fetch_annotation(source='hcps1200',desc='meggamma2', space='fsLR', den='4k',data_dir= neuromaps_download_path)
fsavg = transforms.fslr_to_fsaverage(meg_gamma2, target_density= '10k')
fsavg_lh, fsavg_rh = fsavg
print(fsavg_lh.agg_data().shape)
# save the Gifti file
nib.save(fsavg_lh, neuromaps_fsavg5_path + '/meg_gamma2_fs10k_lh.gii')
nib.save(fsavg_rh, neuromaps_fsavg5_path + '/meg_gamma2_fs10k_rh.gii')

# Saving the array in a text file
fsavg_array = fsavg_lh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/meg_gamma2_fs10k_lh.txt', fsavg_array, newline="\n")
fsavg_array = fsavg_rh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/meg_gamma2_fs10k_rh.txt', fsavg_array, newline="\n")


#-------------
meg_theta = fetch_annotation(source='hcps1200',desc='megtheta', space='fsLR', den='4k',data_dir= neuromaps_download_path)
fsavg = transforms.fslr_to_fsaverage(meg_theta, target_density= '10k')
fsavg_lh, fsavg_rh = fsavg
print(fsavg_lh.agg_data().shape)

# save the Gifti file
nib.save(fsavg_lh, neuromaps_fsavg5_path + '/meg_theta_fs10k_lh.gii')
nib.save(fsavg_rh, neuromaps_fsavg5_path + '/meg_theta_fs10k_rh.gii')

# Saving the array in a text file
fsavg_array = fsavg_lh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/meg_theta_fs10k_lh.txt', fsavg_array, newline="\n")
fsavg_array = fsavg_rh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/meg_theta_fs10k_rh.txt', fsavg_array, newline="\n")



#-----------------------------------------------------------------------------
# ('hcps1200', 'myelinmap', 'fsLR', '32k')
# ('hcps1200', 'thickness', 'fsLR', '32k')
# ('hill2010', 'devexp', 'fsLR', '164k')
# ('hill2010', 'evoexp', 'fsLR', '164k')
# ('margulies2016', 'fcgradient01', 'fsLR', '32k')
# ('mueller2013', 'intersubjvar', 'fsLR', '164k')
# ('neurosynth', 'cogpc1', 'MNI152', '2mm')

myelinmap = fetch_annotation(source='hcps1200',desc='myelinmap', space='fsLR', den='32k',data_dir= neuromaps_download_path)
fsavg = transforms.fslr_to_fsaverage(myelinmap, target_density= '10k')
fsavg_lh, fsavg_rh = fsavg
print(fsavg_lh.agg_data().shape)
# save the Gifti file
nib.save(fsavg_lh, neuromaps_fsavg5_path + '/myelinmap_fs10k_lh.gii')
nib.save(fsavg_rh, neuromaps_fsavg5_path + '/myelinmap_fs10k_rh.gii')

# Saving the array in a text file
fsavg_array = fsavg_lh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/myelinmap_fs10k_lh.txt', fsavg_array, newline="\n")
fsavg_array = fsavg_rh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/myelinmap_fs10k_rh.txt', fsavg_array, newline="\n")



thickness = fetch_annotation(source='hcps1200',desc='thickness', space='fsLR', den='32k',data_dir= neuromaps_download_path)
fsavg = transforms.fslr_to_fsaverage(thickness, target_density= '10k')
fsavg_lh, fsavg_rh = fsavg
print(fsavg_lh.agg_data().shape)
# save the Gifti file
nib.save(fsavg_lh, neuromaps_fsavg5_path + '/thickness_fs10k_lh.gii')
nib.save(fsavg_rh, neuromaps_fsavg5_path + '/thickness_fs10k_rh.gii')

# Saving the array in a text file
fsavg_array = fsavg_lh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/thickness_fs10k_lh.txt', fsavg_array, newline="\n")
fsavg_array = fsavg_rh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/thickness_fs10k_rh.txt', fsavg_array, newline="\n")



#-----------------------------------------------------------------------------
# Dvelopmental and evo expansion
# ---- Right hemisphere only

devexp = fetch_annotation(source='hill2010',desc='devexp', hemi='R',data_dir= neuromaps_download_path)
fsavg = transforms.fslr_to_fsaverage(devexp, target_density= '10k',  hemi='R')
fsavg_rh, = fsavg
print(fsavg_rh.agg_data().shape)
# save the Gifti file
nib.save(fsavg_rh, neuromaps_fsavg5_path + '/devexp_fs10k_rh.gii')
fsavg_array = fsavg_rh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/devexp_fs10k_lh.txt', fsavg_array, newline="\n")
np.savetxt(neuromaps_fsavg5_path + '/devexp_fs10k_rh.txt', fsavg_array, newline="\n")



evoexp = fetch_annotation(source='hill2010',desc='evoexp', hemi='R',data_dir= neuromaps_download_path)
fsavg = transforms.fslr_to_fsaverage(evoexp, target_density= '10k',  hemi='R')
fsavg_rh, = fsavg
print(fsavg_rh.agg_data().shape)
# save the Gifti file
nib.save(fsavg_rh, neuromaps_fsavg5_path + '/evoexp_fs10k_rh.gii')
fsavg_array = fsavg_rh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/evoexp_fs10k_lh.txt', fsavg_array, newline="\n")
np.savetxt(neuromaps_fsavg5_path + '/evoexp_fs10k_rh.txt', fsavg_array, newline="\n")


# evoexp = fetch_annotation(source='hill2010',desc='evoexp', space='fsLR', den='164k',data_dir= neuromaps_download_path)
# fsavg = transforms.fslr_to_fsaverage(evoexp, target_density= '10k')
# fsavg_lh, fsavg_rh = fsavg
# print(fsavg_lh.agg_data().shape)
# # save the Gifti file
# nib.save(fsavg_lh, 'evoexp_fs10k_lh.gii')

#-----------------------------------------------------------------------------
# margulies gradients

margulies = fetch_annotation(source='margulies2016',desc='fcgradient01', space='fsLR', den='32k',data_dir= neuromaps_download_path)
fsavg = transforms.fslr_to_fsaverage(margulies, target_density= '10k')
fsavg_lh, fsavg_rh = fsavg
print(fsavg_lh.agg_data().shape)
# save the Gifti file
nib.save(fsavg_lh, neuromaps_fsavg5_path + '/margulies_fs10k_lh.gii')
nib.save(fsavg_rh, neuromaps_fsavg5_path + '/margulies_fs10k_rh.gii')


# Saving the array in a text file
fsavg_array = fsavg_lh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/margulies_fs10k_lh.txt', fsavg_array, newline="\n")
fsavg_array = fsavg_rh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/margulies_fs10k_rh.txt', fsavg_array, newline="\n")


#-----------------------------------------------------------------------------
# intersubject variability

intersubjvar = fetch_annotation(source='mueller2013',desc='intersubjvar', hemi=('L','R'),data_dir= neuromaps_download_path)
fsavg = transforms.fslr_to_fsaverage(intersubjvar, target_density= '10k',  hemi=('L','R'))
fsavg_lh, fsavg_rh = fsavg
print(fsavg_lh.agg_data().shape)
# save the Gifti file
nib.save(fsavg_lh, neuromaps_fsavg5_path + '/intersubjvar_fs10k_lh.gii')
nib.save(fsavg_rh, neuromaps_fsavg5_path + '/intersubjvar_fs10k_rh.gii')

# Saving the array in a text file
fsavg_array = fsavg_lh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/intersubjvar_fs10k_lh.txt', fsavg_array, newline="\n")
fsavg_array = fsavg_rh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/intersubjvar_fs10k_rh.txt', fsavg_array, newline="\n")


# intersubjvar = fetch_annotation(source='mueller2013',desc='intersubjvar', space='fsLR', den='164k',data_dir= neuromaps_download_path)
# fsavg = transforms.fslr_to_fsaverage(evoexp, target_density= '10k')
# fsavg_lh, fsavg_rh = fsavg
# print(fsavg_lh.agg_data().shape)
# # save the Gifti file
# nib.save(fsavg_lh, 'intersubjvar_fs10k_lh.gii')


#-----------------------------------------------------------------------------
# Metabolic measure maps

# ('raichle', 'cbf', 'fsLR', '164k')
# ('raichle', 'cbv', 'fsLR', '164k')
# ('raichle', 'cmr02', 'fsLR', '164k')
# ('raichle', 'cmruglu', 'fsLR', '164k')


cbf = fetch_annotation(source='raichle',desc='cbf', space='fsLR', den='164k',data_dir= neuromaps_download_path)
fsavg = transforms.fslr_to_fsaverage(cbf, target_density= '10k')
fsavg_lh, fsavg_rh = fsavg
print(fsavg_lh.agg_data().shape)
# save the Gifti file
nib.save(fsavg_lh, neuromaps_fsavg5_path + '/cbf_fs10k_lh.gii')
nib.save(fsavg_rh, neuromaps_fsavg5_path + '/cbf_fs10k_rh.gii')

# Saving the array in a text file
fsavg_array = fsavg_lh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/cbf_fs10k_lh.txt', fsavg_array, newline="\n")
fsavg_array = fsavg_rh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/cbf_fs10k_rh.txt', fsavg_array, newline="\n")


cbv = fetch_annotation(source='raichle',desc='cbv', space='fsLR', den='164k',data_dir= neuromaps_download_path)
fsavg = transforms.fslr_to_fsaverage(cbf, target_density= '10k')
fsavg_lh, fsavg_rh = fsavg
print(fsavg_lh.agg_data().shape)
# save the Gifti file
nib.save(fsavg_lh, neuromaps_fsavg5_path + '/cbv_fs10k_lh.gii')
nib.save(fsavg_rh, neuromaps_fsavg5_path + '/cbv_fs10k_rh.gii')

fsavg_array = fsavg_lh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/cbv_fs10k_lh.txt', fsavg_array, newline="\n")
fsavg_array = fsavg_rh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/cbv_fs10k_rh.txt', fsavg_array, newline="\n")


cmr02 = fetch_annotation(source='raichle',desc='cmr02', space='fsLR', den='164k',data_dir= neuromaps_download_path)
fsavg = transforms.fslr_to_fsaverage(cmr02, target_density= '10k')
fsavg_lh, fsavg_rh = fsavg
print(fsavg_lh.agg_data().shape)
# save the Gifti file
nib.save(fsavg_lh, neuromaps_fsavg5_path + '/cmr02_fs10k_lh.gii')
nib.save(fsavg_rh, neuromaps_fsavg5_path + '/cmr02_fs10k_rh.gii')

fsavg_array = fsavg_lh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/cmr02_fs10k_lh.txt', fsavg_array, newline="\n")
fsavg_array = fsavg_rh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/cmr02_fs10k_rh.txt', fsavg_array, newline="\n")


cmruglu = fetch_annotation(source='raichle',desc='cmruglu', space='fsLR', den='164k',data_dir= neuromaps_download_path)
fsavg = transforms.fslr_to_fsaverage(cmruglu, target_density= '10k')
fsavg_lh, fsavg_rh = fsavg
print(fsavg_lh.agg_data().shape)
# save the Gifti file
nib.save(fsavg_lh, neuromaps_fsavg5_path + '/cmruglu_fs10k_lh.gii')
nib.save(fsavg_rh, neuromaps_fsavg5_path + '/cmruglu_fs10k_rh.gii')

fsavg_array = fsavg_lh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/cmruglu_fs10k_lh.txt', fsavg_array, newline="\n")
fsavg_array = fsavg_rh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/cmruglu_fs10k_rh.txt', fsavg_array, newline="\n")


#-----------------------------------------------------------------------------
# allometric scaling maps
# ('reardon2018', 'scalingnih', 'civet', '41k')
# ('reardon2018', 'scalingpnc', 'civet', '41k')


scalingnih = fetch_annotation(source='reardon2018',desc='scalingnih', space='civet', den='41k',data_dir= neuromaps_download_path)
fsavg = transforms.civet_to_fsaverage(scalingnih, target_density= '10k')
fsavg_lh, fsavg_rh = fsavg
print(fsavg_lh.agg_data().shape)
# save the Gifti file
nib.save(fsavg_lh, neuromaps_fsavg5_path + '/scalingnih_fs10k_lh.gii')
nib.save(fsavg_rh, neuromaps_fsavg5_path + '/scalingnih_fs10k_rh.gii')

fsavg_array = fsavg_lh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/scalingnih_fs10k_lh.txt', fsavg_array, newline="\n")
fsavg_array = fsavg_rh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/scalingnih_fs10k_rh.txt', fsavg_array, newline="\n")


scalingpnc = fetch_annotation(source='reardon2018',desc='scalingpnc', space='civet', den='41k',data_dir= neuromaps_download_path)
fsavg = transforms.civet_to_fsaverage(scalingpnc, target_density= '10k')
fsavg_lh, fsavg_rh = fsavg
print(fsavg_lh.agg_data().shape)
# save the Gifti file
nib.save(fsavg_lh, neuromaps_fsavg5_path + '/scalingpnc_fs10k_lh.gii')
nib.save(fsavg_rh, neuromaps_fsavg5_path + '/scalingpnc_fs10k_rh.gii')

fsavg_array = fsavg_lh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/scalingpnc_fs10k_lh.txt', fsavg_array, newline="\n")
fsavg_array = fsavg_rh.agg_data()
np.savetxt(neuromaps_fsavg5_path + '/scalingpnc_fs10k_rh.txt', fsavg_array, newline="\n")


#-----------------------------------------------------------------------------
# Step 3: parcellate
#-----------------------------------------------------------------------------

# NOTE: needs numpy version <1.24 to avoid error in np.float  # AttributeError: module 'numpy' has no attribute 'float'
from enigmatoolbox.utils.parcellation import surface_to_parcel

# Source code: https://enigma-toolbox.readthedocs.io/en/latest/_modules/enigmatoolbox/utils/parcellation.html#surface_to_parcel

# Parcellate vertexwise data (Example)
#CT_schaefer_200 = surface_to_parcel(CT, 'schaefer_200_conte69')

#------
# Users can use this function to parcellate (fsaverage5 or conte69) vertexwise data
# using different parcellations,
#  including: aparc_fsa5, glasser_fsa5, schaefer_100_fsa5, schaefer_200_fsa5,
#             schaefer_300_fsa5, schaefer_400_fsa5, aparc_conte69, glasser_conte69,
#             schaefer_100_conte69, schaefer_200_conte69, schaefer_300_conte69,
#             schaefer_400_conte69.
#--------


"""
# 20 brain maps from Markello et al.
('abagen', 'genepc1', 'fsaverage', '10k')
('hcps1200', 'megalpha', 'fsLR', '4k')
('hcps1200', 'megbeta', 'fsLR', '4k')
('hcps1200', 'megdelta', 'fsLR', '4k')
('hcps1200', 'meggamma1', 'fsLR', '4k')
('hcps1200', 'meggamma2', 'fsLR', '4k')
('hcps1200', 'megtheta', 'fsLR', '4k')
#  ('hcps1200', 'megtimescale', 'fsLR', '4k')
('hcps1200', 'myelinmap', 'fsLR', '32k')
('hcps1200', 'thickness', 'fsLR', '32k')
('hill2010', 'devexp', 'fsLR', '164k')
('hill2010', 'evoexp', 'fsLR', '164k')
('margulies2016', 'fcgradient01', 'fsLR', '32k')
('mueller2013', 'intersubjvar', 'fsLR', '164k')
('neurosynth', 'cogpc1', 'MNI152', '2mm')
('raichle', 'cbf', 'fsLR', '164k')
('raichle', 'cbv', 'fsLR', '164k')
('raichle', 'cmr02', 'fsLR', '164k')
('raichle', 'cmruglu', 'fsLR', '164k')
('reardon2018', 'scalingnih', 'civet', '41k')
('reardon2018', 'scalingpnc', 'civet', '41k')

"""

array_map_names = ['neurosynth','abagen',
                   'meg_alpha','meg_beta','meg_delta','meg_gamma1','meg_gamma2','meg_theta',
                   'myelinmap','thickness','devexp','evoexp','margulies','intersubjvar',
                   'cbf','cbv','cmr02','cmruglu',
                   'scalingnih','scalingpnc']

#-----------------------------------------------------------------------------
# ##  1. Desikan

parc_name = 'desikan'
neuromaps_parc_path = base_folder + "/data_ref/neuromaps_aug2023/neuromaps_desikan"  # Folder where parcellated maps will be saved

# Check whether the specified path exists or not
isExist = os.path.exists(neuromaps_parc_path)
if not isExist:
   # Create a new directory because it does not exist
   os.makedirs(neuromaps_parc_path)
   
#-- parcellate
for x in range(len(array_map_names)):
    print(array_map_names[x])
    map_name = array_map_names[x]
    data_array_lh_rh = []
    for _, h in enumerate(['lh', 'rh']):
         data_array_lh_rh = np.append(data_array_lh_rh, np.loadtxt(neuromaps_fsavg5_path + '/' + map_name + '_fs10k_{}.txt'.format(h)))
       
    fsavg_to_dk = surface_to_parcel(data_array_lh_rh, 'aparc_fsa5')
    np.savetxt(neuromaps_parc_path + '/' + map_name + '_parc_' + parc_name + '_lhrh.txt', fsavg_to_dk, newline="\n")

#-----------------------------------------------------------------------------
# ##  2 Glasser
# in Matlab only
#-------------------------------
# parc_name = 'glasser'
# neuromaps_parc_path = base_folder + "/data_ref/neuromaps_aug2023/neuromaps_glasser"  # Folder where parcellated maps will be saved

 
# # Check whether the specified path exists or not
# isExist = os.path.exists(neuromaps_parc_path)
# if not isExist:
#    # Create a new directory because it does not exist
#    os.makedirs(neuromaps_parc_path)
   

# #-- parcellate
# for x in range(len(array_map_names)):
#     print(array_map_names[x])
#     map_name = array_map_names[x]
#     data_array_lh_rh = []
#     for _, h in enumerate(['lh', 'rh']):
#          data_array_lh_rh = np.append(data_array_lh_rh, np.loadtxt(neuromaps_fsavg5_path + '/' + map_name + '_fs10k_{}.txt'.format(h)))
       
#     fsavg_to_dk = surface_to_parcel(data_array_lh_rh, 'glasser_fsa5')
#     np.savetxt(neuromaps_parc_path + '/' + map_name + '_parc_' + parc_name + '_lhrh.txt', fsavg_to_dk, newline="\n")


#-----------------------------------------------------------------------------
# ##  3 schaefer_100_fsa5
parc_name = 'schaefer_100'
neuromaps_parc_path = base_folder + "/data_ref/neuromaps_aug2023/neuromaps_schaefer_100" # Folder where parcellated maps will be saved


# Check whether the specified path exists or not
isExist = os.path.exists(neuromaps_parc_path)
if not isExist:
   # Create a new directory because it does not exist
   os.makedirs(neuromaps_parc_path)
   

#-- parcellate
for x in range(len(array_map_names)):
    print(array_map_names[x])
    map_name = array_map_names[x]
    data_array_lh_rh = []
    for _, h in enumerate(['lh', 'rh']):
         data_array_lh_rh = np.append(data_array_lh_rh, np.loadtxt(neuromaps_fsavg5_path + '/' + map_name + '_fs10k_{}.txt'.format(h)))
       
    fsavg_to_dk = surface_to_parcel(data_array_lh_rh, 'schaefer_100_fsa5')
    np.savetxt(neuromaps_parc_path + '/' + map_name + '_parc_' + parc_name + '_lhrh.txt', fsavg_to_dk, newline="\n")


#-----------------------------------------------------------------------------
# ##  4 schaefer_200_fsa5
parc_name = 'schaefer_200'
neuromaps_parc_path = base_folder + "/data_ref/neuromaps_aug2023/neuromaps_schaefer_200" # Folder where parcellated maps will be saved


# Check whether the specified path exists or not
isExist = os.path.exists(neuromaps_parc_path)
if not isExist:
   # Create a new directory because it does not exist
   os.makedirs(neuromaps_parc_path)
   

#-- parcellate
for x in range(len(array_map_names)):
    print(array_map_names[x])
    map_name = array_map_names[x]
    data_array_lh_rh = []
    for _, h in enumerate(['lh', 'rh']):
         data_array_lh_rh = np.append(data_array_lh_rh, np.loadtxt(neuromaps_fsavg5_path + '/' + map_name + '_fs10k_{}.txt'.format(h)))
       
    fsavg_to_dk = surface_to_parcel(data_array_lh_rh, 'schaefer_200_fsa5')
    np.savetxt(neuromaps_parc_path + '/' + map_name + '_parc_' + parc_name + '_lhrh.txt', fsavg_to_dk, newline="\n")

#-----------------------------------------------------------------------------
# ##  5 schaefer_300_fsa5
parc_name = 'schaefer_300'
neuromaps_parc_path = base_folder + "/data_ref/neuromaps_aug2023/neuromaps_schaefer_300" # Folder where parcellated maps will be saved


# Check whether the specified path exists or not
isExist = os.path.exists(neuromaps_parc_path)
if not isExist:
   # Create a new directory because it does not exist
   os.makedirs(neuromaps_parc_path)
   

#-- parcellate
for x in range(len(array_map_names)):
    print(array_map_names[x])
    map_name = array_map_names[x]
    data_array_lh_rh = []
    for _, h in enumerate(['lh', 'rh']):
         data_array_lh_rh = np.append(data_array_lh_rh, np.loadtxt(neuromaps_fsavg5_path + '/' + map_name + '_fs10k_{}.txt'.format(h)))
       
    fsavg_to_dk = surface_to_parcel(data_array_lh_rh, 'schaefer_300_fsa5')
    np.savetxt(neuromaps_parc_path + '/' + map_name + '_parc_' + parc_name + '_lhrh.txt', fsavg_to_dk, newline="\n")




#-----------------------------------------------------------------------------
# ##  6 schaefer_400_fsa5
parc_name = 'schaefer_200'
neuromaps_parc_path = base_folder + "/data_ref/neuromaps_aug2023/neuromaps_schaefer_400"  # Folder where parcellated maps will be saved


# Check whether the specified path exists or not
isExist = os.path.exists(neuromaps_parc_path)
if not isExist:
   # Create a new directory because it does not exist
   os.makedirs(neuromaps_parc_path)
   

#-- parcellate
for x in range(len(array_map_names)):
    print(array_map_names[x])
    map_name = array_map_names[x]
    data_array_lh_rh = []
    for _, h in enumerate(['lh', 'rh']):
         data_array_lh_rh = np.append(data_array_lh_rh, np.loadtxt(neuromaps_fsavg5_path + '/' + map_name + '_fs10k_{}.txt'.format(h)))
       
    fsavg_to_dk = surface_to_parcel(data_array_lh_rh, 'schaefer_400_fsa5')
    np.savetxt(neuromaps_parc_path + '/' + map_name + '_parc_' + parc_name + '_lhrh.txt', fsavg_to_dk, newline="\n")


#-----------------------------------------------------------------------------
# ##  7. economo_koskinas_fsa5
parc_name = 'economo_koskinas'
neuromaps_parc_path = base_folder + "/data_ref/neuromaps_aug2023/neuromaps_economo_koskinas"  # Folder where parcellated maps will be saved


# Check whether the specified path exists or not
isExist = os.path.exists(neuromaps_parc_path)
if not isExist:
   # Create a new directory because it does not exist
   os.makedirs(neuromaps_parc_path)
   

#-- parcellate
for x in range(len(array_map_names)):
    print(array_map_names[x])
    map_name = array_map_names[x]
    data_array_lh_rh = []
    for _, h in enumerate(['lh', 'rh']):
         data_array_lh_rh = np.append(data_array_lh_rh, np.loadtxt(neuromaps_fsavg5_path + '/' + map_name + '_fs10k_{}.txt'.format(h)))
       
    fsavg_to_dk = surface_to_parcel(data_array_lh_rh, 'economo_koskinas_fsa5')
    np.savetxt(neuromaps_parc_path + '/' + map_name + '_parc_' + parc_name + '_lhrh.txt', fsavg_to_dk, newline="\n")




#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------





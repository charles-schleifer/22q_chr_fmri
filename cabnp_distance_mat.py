#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 11:49:14 2023
Script to generate distance matrices for BrainSMASH with CAB-NP parcellation
@author: C. Schleifer
"""

#import pandas as pd
#import numpy as np
#pd.set_option('display.max_rows', 1000)
#pd.set_option('display.max_columns', 500)

# install brainsmash
# pip install brainsmash
# import relevant functions
#from brainsmash.mapgen.base import Base
#from brainsmash.utils.dataio import load
from brainsmash.workbench.geo import cortex
#from brainsmash.workbench.geo import subcortex
#from brainsmash.workbench.geo import parcellate

for hemi in ["L","R"]:
    print(hemi)
    label="CortexParcelLabels_"+hemi+".dlabel.nii"
    surface = "S1200."+hemi+".midthickness_MSMAll.32k_fs_LR.surf.gii"
    cabnp_path = "/Users/charlie/Dropbox/github/22q_chr_fmri/CAB-NP/"
    surface_path = cabnp_path+surface
    label_path = cabnp_path+label
    opath = cabnp_path+"/CABPN_surface_geo_dist_mat_"+hemi+".txt"
    cortex(surface=surface_path, dlabel=label_path, outfile=opath, euclid=False)



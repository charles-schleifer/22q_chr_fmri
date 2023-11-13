#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 18:19:10 2023
Script to generate sets of brain map permutations preserving spatial autocorrelation for null distributions
@author: C. Schleifer
"""
#import numpy as np
import pandas as pd
from brainsmash.mapgen.base import Base
from brainsmash.mapgen.eval import base_fit


# input and output path
path="/Users/charlie/Dropbox/github/22q_chr_fmri/CAB-NP/"

# function to compute n surrogates for a given map
# https://brainsmash.readthedocs.io/en/latest/gettingstarted.html#parcellated-surrogate-maps
def bransmash_generate(n:int, map:str, mat:str):
    print("generating "+str(n)+" permutations")
    print(map)
    print(mat)
    base = Base(x=map, D=mat)
    out = base(n=n)
    print("done")
    return out

# get left hemisphere RSFA permutations
L_22q_TD_RSFA = bransmash_generate(n=10000, map=path+"22q_TD_RSFA_L.txt", mat=path+"CABPN_surface_geo_dist_mat_L.txt")
# save as CSV
pd.DataFrame(L_22q_TD_RSFA).to_csv(path+"22q_TD_RSFA_L_permuted.csv", header=False, index=False)

# get right hemisphere RSFA permutations
R_22q_TD_RSFA = bransmash_generate(n=10000, map=path+"22q_TD_RSFA_R.txt", mat=path+"CABPN_surface_geo_dist_mat_R.txt")
# save as CSV
pd.DataFrame(R_22q_TD_RSFA).to_csv(path+"22q_TD_RSFA_R_permuted.csv", header=False, index=False)

# get left hemisphere NetHo permutations
L_22q_TD_NetHo = bransmash_generate(n=10000, map=path+"22q_TD_NetHo_L.txt", mat=path+"CABPN_surface_geo_dist_mat_L.txt")
# save as CSV
pd.DataFrame(L_22q_TD_NetHo).to_csv(path+"22q_TD_NetHo_L_permuted.csv", header=False, index=False)

# get right hemisphere NetHo permutations
R_22q_TD_NetHo = bransmash_generate(n=10000, map=path+"22q_TD_NetHo_R.txt", mat=path+"CABPN_surface_geo_dist_mat_R.txt")
# save as CSV
pd.DataFrame(R_22q_TD_NetHo).to_csv(path+"22q_TD_NetHo_R_permuted.csv", header=False, index=False)

# get left hemisphere GBC permutations
L_22q_TD_GBC = bransmash_generate(n=10000, map=path+"22q_TD_GBC_L.txt", mat=path+"CABPN_surface_geo_dist_mat_L.txt")
# save as CSV
pd.DataFrame(L_22q_TD_GBC).to_csv(path+"22q_TD_GBC_L_permuted.csv", header=False, index=False)

# get right hemisphere GBC permutations
R_22q_TD_GBC = bransmash_generate(n=10000, map=path+"22q_TD_GBC_R.txt", mat=path+"CABPN_surface_geo_dist_mat_R.txt")
# save as CSV
pd.DataFrame(R_22q_TD_GBC).to_csv(path+"22q_TD_GBC_R_permuted.csv", header=False, index=False)


# function evaluate variogram fits
# https://brainsmash.readthedocs.io/en/latest/gettingstarted.html#evaluating-variogram-fits
def bransmash_eval(n:int, map:str, mat:str):
    print("evaluating "+str(n)+" permutations")
    print(map)
    print(mat)
    out = base_fit(map, mat, nsurr=n)
    return out

#L_22q_TD_RSFA_eval = bransmash_eval(n=100, map=path+"22q_TD_RSFA_L.txt", mat=path+"CABPN_surface_geo_dist_mat_L.txt")
#
#R_22q_TD_RSFA_eval = bransmash_eval(n=100, map=path+"22q_TD_RSFA_R.txt", mat=path+"CABPN_surface_geo_dist_mat_R.txt")
#
#L_22q_TD_NetHo_eval = bransmash_eval(n=100, map=path+"22q_TD_NetHo_L.txt", mat=path+"CABPN_surface_geo_dist_mat_L.txt")
#
#R_22q_TD_NetHo_eval = bransmash_eval(n=100, map=path+"22q_TD_NetHo_R.txt", mat=path+"CABPN_surface_geo_dist_mat_R.txt")


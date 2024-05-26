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


control_22q_rsfa_L = bransmash_generate(n=10000, map=path+"control_22q_rsfa_L.txt", mat=path+"CABPN_surface_geo_dist_mat_L.txt")
pd.DataFrame(control_22q_rsfa_L).to_csv(path+"control_22q_rsfa_L_permuted.csv", header=False, index=False)

control_22q_netho_L = bransmash_generate(n=10000, map=path+"control_22q_netho_L.txt", mat=path+"CABPN_surface_geo_dist_mat_L.txt")
pd.DataFrame(control_22q_netho_L).to_csv(path+"control_22q_netho_L_permuted.csv", header=False, index=False)

control_22q_gbc_L = bransmash_generate(n=10000, map=path+"control_22q_gbc_L.txt", mat=path+"CABPN_surface_geo_dist_mat_L.txt")
pd.DataFrame(control_22q_gbc_L).to_csv(path+"control_22q_gbc_L_permuted.csv", header=False, index=False)

control_22q_rsfa_R = bransmash_generate(n=10000, map=path+"control_22q_rsfa_R.txt", mat=path+"CABPN_surface_geo_dist_mat_L.txt")
pd.DataFrame(control_22q_rsfa_R).to_csv(path+"control_22q_rsfa_R_permuted.csv", header=False, index=False)

control_22q_netho_R = bransmash_generate(n=10000, map=path+"control_22q_netho_R.txt", mat=path+"CABPN_surface_geo_dist_mat_L.txt")
pd.DataFrame(control_22q_netho_R).to_csv(path+"control_22q_netho_R_permuted.csv", header=False, index=False)

control_22q_gbc_R = bransmash_generate(n=10000, map=path+"control_22q_gbc_R.txt", mat=path+"CABPN_surface_geo_dist_mat_L.txt")
pd.DataFrame(control_22q_gbc_R).to_csv(path+"control_22q_gbc_R_permuted.csv", header=False, index=False)

control_n_rsfa_L = bransmash_generate(n=10000, map=path+"control_n_rsfa_L.txt", mat=path+"CABPN_surface_geo_dist_mat_L.txt")
pd.DataFrame(control_n_rsfa_L).to_csv(path+"control_n_rsfa_L_permuted.csv", header=False, index=False)

control_n_netho_L = bransmash_generate(n=10000, map=path+"control_n_netho_L.txt", mat=path+"CABPN_surface_geo_dist_mat_L.txt")
pd.DataFrame(control_n_netho_L).to_csv(path+"control_n_netho_L_permuted.csv", header=False, index=False)

control_n_gbc_L = bransmash_generate(n=10000, map=path+"control_n_gbc_L.txt", mat=path+"CABPN_surface_geo_dist_mat_L.txt")
pd.DataFrame(control_n_gbc_L).to_csv(path+"control_n_gbc_L_permuted.csv", header=False, index=False)

control_n_rsfa_R = bransmash_generate(n=10000, map=path+"control_n_rsfa_R.txt", mat=path+"CABPN_surface_geo_dist_mat_L.txt")
pd.DataFrame(control_n_rsfa_R).to_csv(path+"control_n_rsfa_R_permuted.csv", header=False, index=False)

control_n_netho_R = bransmash_generate(n=10000, map=path+"control_n_netho_R.txt", mat=path+"CABPN_surface_geo_dist_mat_L.txt")
pd.DataFrame(control_n_netho_R).to_csv(path+"control_n_netho_R_permuted.csv", header=False, index=False)

control_n_gbc_R = bransmash_generate(n=10000, map=path+"control_n_gbc_R.txt", mat=path+"CABPN_surface_geo_dist_mat_L.txt")
pd.DataFrame(control_n_gbc_R).to_csv(path+"control_n_gbc_R_permuted.csv", header=False, index=False)




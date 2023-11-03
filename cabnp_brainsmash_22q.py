#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 18:19:10 2023
Script to generate sets of brain map permutations preserving spatial autocorrelation for null distributions
@author: charlie
"""

from brainsmash.mapgen.base import Base

# function to compute n surrogates for a given map
def bransmash_generate(n:int, map:str, mat:str):
    #brain_map_file = "LeftParcelMyelin.txt"  # use absolute paths if necessary!
    #dist_mat_file = "LeftParcelGeodesicDistmat.txt
    print(n)
    print(map)
    print(mat)
    

bransmash_generate(n=1000, map="test", mat="test2")
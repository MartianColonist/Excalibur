#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 12:06:47 2022

@author: arnav
"""

from Cthulhu import HITEMP
import os
import numpy as np
import shutil

def test_create_id_dict():
    assert HITEMP.create_id_dict() == {'H2O': 1, 'ClONO2': 35, 'CO2': 2, 'O3': 3, 'N2O': 4, 'CO': 5, 'CH4': 6, 'O2': 7, 'NO': 8, 'SO2': 9, 'NO2': 10, 'NH3': 11, 'HNO3': 12, 'OH': 13, 'HF': 14, 'HCl': 15, 'HBr': 16, 'HI': 17, 'ClO': 18, 'OCS': 19, 'H2CO': 20, 'HOCl': 21, 'N2': 22, 'HCN': 23, 'CH3Cl': 24, 'H2O2': 25, 'C2H2': 26, 'C2H6': 27, 'PH3': 28, 'COF2': 29, 'SF6': 30, 'H2S': 31, 'HCOOH': 32, 'HO2': 33, 'O': 34, 'NOp': 36, 'HOBr': 37, 'C2H4': 38, 'CH3OH': 39, 'CH3Br': 40, 'CH3CN': 41, 'CF4': 42, 'C4H2': 43, 'HC3N': 44, 'H2': 45, 'CS': 46, 'SO3': 47, 'C2N2': 48, 'COCl2': 49}
    
    
def test_create_pf():
    HITEMP.create_pf(5, 1, './')  # create CO.pf file in 'tests' folder
    file = open('./CO.pf', 'r')
    lines = file.readlines()
    
    Q = [] # partition function
    
    for line in lines[1:]:
        contents = line.split()
        Q.append(float(contents[1]))
        
    file.close()
        
    os.remove('./CO.pf')
        
    Q = np.asarray(Q) 
    
    assert Q.all() > 0
    assert np.isnan(np.sum(Q)) == False
    
  
def test_summon_HITEMP(): # tests that the air broadening file for N2 was created, because that means that the summon function was able to work without errors
    HITEMP.summon_HITEMP(5, 1)  # Summon N2 from HITEMP
    air_broad = './input/CO  ~  (12C-16O)/HITEMP/air.broad'
    broadening_exists = os.path.exists(air_broad)  # True if path exists to air broadening file
    shutil.rmtree('./input')
    assert broadening_exists == True
   
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  5 13:37:11 2022

@author: arnav
"""
from excalibur import HITRAN

import numpy as np
import os
import shutil

def test_create_id_dict():
    assert HITRAN.create_id_dict() == {'H2O': 1, 'ClONO2': 35, 'CO2': 2, 'O3': 3, 'N2O': 4, 'CO': 5, 'CH4': 6, 'O2': 7, 'NO': 8, 'SO2': 9, 'NO2': 10, 'NH3': 11, 'HNO3': 12, 'OH': 13, 'HF': 14, 'HCl': 15, 'HBr': 16, 'HI': 17, 'ClO': 18, 'OCS': 19, 'H2CO': 20, 'HOCl': 21, 'N2': 22, 'HCN': 23, 'CH3Cl': 24, 'H2O2': 25, 'C2H2': 26, 'C2H6': 27, 'PH3': 28, 'COF2': 29, 'SF6': 30, 'H2S': 31, 'HCOOH': 32, 'HO2': 33, 'O': 34, 'NOp': 36, 'HOBr': 37, 'C2H4': 38, 'CH3OH': 39, 'CH3Br': 40, 'CH3CN': 41, 'CF4': 42, 'C4H2': 43, 'HC3N': 44, 'H2': 45, 'CS': 46, 'SO3': 47, 'C2N2': 48, 'COCl2': 49}
    
def test_replace_iso_name():
    assert HITRAN.replace_iso_name('H2(16O)') == '(1H2-16O)'
    assert HITRAN.replace_iso_name('(14N)2') == '(14N2)'
    assert HITRAN.replace_iso_name('(12C)H4') == '(12C-1H4)'
    assert HITRAN.replace_iso_name('H(14N)(16O)3') == '(1H-14N-16O3)'
    
def test_create_pf():  # Tests that the partition function of N2 from HITRAN is non-zero and contains no nan values
    HITRAN.create_pf(22, 1, './')  # create N2.pf file in 'tests' folder
    file = open('./N2.pf', 'r')
    lines = file.readlines()
    
    Q = [] # partition function
    
    for line in lines[1:]:
        contents = line.split()
        Q.append(float(contents[1]))
        
    file.close()
        
    os.remove('./N2.pf')
        
    Q = np.asarray(Q) 
    
    assert Q.all() > 0
    assert np.isnan(np.sum(Q)) == False
    
def test_summon_HITRAN(): # tests that the air broadening file for N2 was created, because that means that the summon function was able to work without errors
    HITRAN.summon_HITRAN(22, 1)  # Summon N2 from HITRAN
    air_broad = './input/N2  ~  (14N2)/HITRAN/air.broad'
    broadening_exists = os.path.exists(air_broad)  # True if path exists to air broadening file
    shutil.rmtree('./input')
    assert broadening_exists == True
    
    
    
    
    
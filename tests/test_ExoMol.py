#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 12:03:57 2022

@author: arnav
"""

from Cthulhu import ExoMol
import os
import numpy as np

    
    
def test_get_default_iso():
    assert ExoMol.get_default_iso('H2_p') == '1H-2H_p'
    assert ExoMol.get_default_iso('FeH') == '56Fe-1H'
    assert ExoMol.get_default_iso('O2') == '16O2'
    assert ExoMol.get_default_iso('HNO3') == '1H-14N-16O3'
    
    
def test_get_default_linelist():
    assert ExoMol.get_default_linelist('CH4', '12C-1H4') == 'MM'
    
def test_summon_ExoMol():
    # test that there is a .pf file at the end of this
    ExoMol.summon_ExoMol('FeH', '56Fe-1H', 'MoLLIST', 'https://www.exomol.com/data/molecules/FeH/56Fe-1H/MoLLIST/')
    pf_exists = os.path.exists('./input/FeH  ~  (56Fe-1H)/ExoMol/MoLLIST/MoLLIST.pf')  # True if path exists to air broadening file
    assert pf_exists == True
    
    
def test_process_files():
    # tests that the partition function file was created properly (since the broadening doesn't always exist)'
    file = open('./input/FeH  ~  (56Fe-1H)/ExoMol/MoLLIST/MoLLIST.pf', 'r')
    lines = file.readlines()
    
    Q = [] # partition function
    
    for line in lines[1:]:
        contents = line.split()
        Q.append(float(contents[1]))
        
    file.close()
        
    Q = np.asarray(Q) 
    
    # Not removing the input directory yet in case I want to test the load_states function below
    
    assert Q.all() > 0
    assert np.isnan(np.sum(Q)) == False
    

'''
    def test_load_states():
'''     
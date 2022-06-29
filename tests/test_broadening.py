#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 21:05:43 2022

@author: arnav
"""

from excalibur import broadening
from excalibur import ExoMol
import os

def test_det_broad():
    ExoMol.summon_ExoMol('HCN', '1H-12C-14N', 'Harris', 'https://www.exomol.com/data/molecules/HCN/1H-12C-14N/Harris/')
    
    broadening_type = broadening.det_broad('./input/HCN  ~  (1H-12C-14N)/ExoMol/Harris')
    
    assert broadening_type == 'H2-He'
    
    H2_He_exists = os.path.exists('./input/HCN  ~  (1H-12C-14N)/ExoMol/Harris/H2.broad') and os.path.exists('./input/HCN  ~  (1H-12C-14N)/ExoMol/Harris/He.broad')
    
    assert H2_He_exists == True
    
    ExoMol.summon_ExoMol('CS', '12C-32S', 'JnK', 'https://www.exomol.com/data/molecules/CS/12C-32S/JnK/')
    
    broadening_type = broadening.det_broad('./input/CS  ~  (12C-32S)/ExoMol/JnK')
    
    assert broadening_type == 'air'
    
    air_exists = os.path.exists('./input/CS  ~  (12C-32S)/ExoMol/JnK/air.broad')
    H2_He_exists = os.path.exists('./input/CS  ~  (12C-32S)/ExoMol/JnK/H2.broad')
    
    assert air_exists == True
    assert H2_He_exists == False
    
    ExoMol.summon_ExoMol('NH', '14N-1H', 'MoLLIST', 'https://www.exomol.com/data/molecules/NH/14N-1H/MoLLIST/')
    
    broadening_type = broadening.det_broad('./input/NH  ~  (14N-1H)/ExoMol/MoLLIST')
    
    assert broadening_type == 'SB07'
    
    air_exists = os.path.exists('./input/NH  ~  (14N-1H)/ExoMol/MoLLIST/air.broad')
    H2_He_exists = os.path.exists('./input/NH  ~  (14N-1H)/ExoMol/MoLLIST/H2.broad')
    
    # Sharp & Burrows broadening not created until the time of computing cross sections. To test, we assert H2-He and air broadening files don't exist
    assert air_exists == False
    assert H2_He_exists == False
    
'''    
def test_create_SB07():
    # ask Ryan what values to test here

def test_read_H2_He():
    # ask Ryan what values to test here
    
    
def test_read_air():
    # ask Ryan what values to test here
    
def test_read_SB07():
    # ask Ryan what values to test here
    
def test_read_custom():
    # ask Ryan what values to test here... would need to add our own broadening file though    '''
    

    
    
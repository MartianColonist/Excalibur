#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  5 13:37:31 2022

@author: arnav
"""

from Excalibur import VALD
import os
import re
import numpy as np
import shutil

def test_summon_VALD():
    # tests that the partition function exists at the end of running summon function
    # this test is intended to be run out of the main folder, not tests. Relevant for the file paths otherwise it will fail.
    
    VALD.summon_VALD('O', 1, './VALD Line Lists/')
    pf_exists = os.path.exists('./input/O  ~  (I)/VALD/O_I.pf')  # True if path exists to air broadening file
    assert pf_exists == True
    #check that the partition function exists at the end of this

    
def test_filter_pf():
    # Tests that the partition function of O from VALD is non-zero and contains no nan values
    
    file = open('./input/O  ~  (I)/VALD/O_I.pf', 'r')
    lines = file.readlines()
    
    Q = [] # partition function
    
    for line in lines[1:]:
        contents = line.split()
        Q.append(float(contents[1]))
        
    file.close()
    
    # Not removing the input directory yet in case I want to test the load_line_list function below
    
        
    Q = np.asarray(Q) 
    
    assert Q.all() > 0
    assert np.isnan(np.sum(Q)) == False
 
'''
def test_load_line_list: gotta ask Ryan what the checks should be for each array that's returned'''
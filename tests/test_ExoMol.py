#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 12:03:57 2022

@author: arnav
"""

from excalibur import ExoMol


def test_check():
    '''
    test that incorrect input parameters raises an HTTPError as expected
    test that correct parameters does not raise an error
    '''
    
    
def test_get_default_iso():
    '''
    test that inputting H2+ gives the correct default isotope
    conditions to test:
         normal molecule (FeH or something)
         more than one of an atom (O2 or something)
         3 or more atoms (HCN, HNO3)
         ions (H3O+)
         H2+ because I have a specific case for this molecule
    '''
    
    
def test_get_default_linelist():
    '''
    conditions to test:
        MgH and CaH make the user input
        Anything from the list returns the correct default
        if not in the list, the system exits (not sure how to test this with an assert statement)
        
    '''
    
    
def test_determine_linelist():
    '''
    This is all user input, don't think it makes sense to test with a test function. 
    Can test manually if desired.'
    '''
    
def test_process_files(input_dir):
    '''
    This function just reformats the broadening and partition function files ExoMol gives. Not sure how to test, 
    outside of adding a separate directory in the github with a .pf and .broad file, and testing to make sure 
    the result matches what we expect (i.e. provide the file it should result in too)
    '''
    
def test_summon_ExoMol(molecule, isotopologue, line_list, URL):
    '''
    Don't see anything to test here, just calls other functions mostly
    '''
    
def test_load_states(input_directory):
    '''
    Just reads in the .states file that ExoMol provides and gets relevant information/arrays from it. Again can't test,
    unless we provide a .states file and the result that is expected'
    '''
    

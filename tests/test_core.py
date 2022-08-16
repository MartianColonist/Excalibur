#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 13:56:11 2022

@author: arnav
"""

from excalibur import core

def test_mass():
    # Tests the mass of a species from VALD, HITRAN, HITEMP, and ExoMol
    
    if core.mass('O', 1, 'vald') > 15.9 and core.mass('O', 1, 'vald') < 16.1: 
        oxygen_mass = True
    else: 
        oxygen_mass = False
        
    assert oxygen_mass == True
    
    if core.mass('CO', '(12C-16O)', 'hitran') > 27.9 and core.mass('CO', '(12C-16O)', 'hitran') < 28.1:
        CO_mass = True
    else: 
        CO_mass = False
        
    assert CO_mass == True
    
    if core.mass('CO', '(13C-16O)', 'hitemp') > 28.9 and core.mass('CO', '(13C-16O)', 'hitemp') < 29.1:
        CO_mass_hitemp = True
    else: 
        CO_mass_hitemp = False
        
    assert CO_mass_hitemp == True
    
    if core.mass('OH', '(16O-1H)', 'MoLLIST') > 16.9 and core.mass('OH', '(16O-1H)', 'MoLLIST') < 17.1:
        OH_mass = True
    else: 
        OH_mass = False
        
    assert OH_mass == True
    
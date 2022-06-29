#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 12:51:11 2022

@author: arnav
"""

from excalibur import misc

def test_check_molecule():
    assert misc.check_molecule('H2O') == True
    assert misc.check_molecule('Hyd') == True  # Looks weird, but should pass based on the details of the check_molecule function
    assert misc.check_molecule('He') == False
    assert misc.check_molecule('oH') == True   # Looks weird, but should pass based on the details of the check_molecule function
    assert misc.check_molecule('C') == False
    assert misc.check_molecule('c') == True    # Looks weird, but should pass based on the details of the check_molecule function
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 22:16:59 2020

@author: arnav
"""

import os
import numpy as np
import pandas as pd
import re
import time
import copy
import requests
import sys
from bs4 import BeautifulSoup
from scipy.interpolate import UnivariateSpline as Interp
from hapi import molecularMass, moleculeName, isotopologueName


from Voigt import Voigt_width, Generate_Voigt_grid_molecules, gamma_L_VALD, gamma_L_impact, analytic_alkali
from calculations import find_index, prior_index, bin_cross_section_atom, bin_cross_section_molecule
from calculations import produce_total_cross_section_EXOMOL, produce_total_cross_section_HITRAN
from calculations import produce_total_cross_section_VALD_atom, produce_total_cross_section_VALD_molecule

from constants import Nu_ref, gamma_0, n_L, P_ref, T_ref
from constants import c, kb, h, m_e, c2, u, pi

from Download_ExoMol import get_default_iso, get_default_linelist


def create_id_dict():
    mol_ID = []
    mol_name = []
    for i in range(1, 50):
        if i == 35: # Skip 35 since moleculeName(35) throws an error from hapi.py
            continue
        else:
            mol_ID.append(i)
            mol_name.append(moleculeName(i))
        
        mol_ID.append(35)
        mol_name.append('ClONO2')

    names_and_IDs = zip(mol_name, mol_ID)

    molecule_dict = dict(names_and_IDs) 
    
    return molecule_dict

def parse_directory(directory):
    """
    Determine which molecule and linelist this directory contains data for (assumes data was downloaded using our script)

    Parameters
    ----------
    directory : TYPE
        DESCRIPTION.

    Returns
    -------
    molecule : TYPE
        DESCRIPTION.
    linelist : TYPE
        DESCRIPTION.

    """
    
    directory_name = os.path.abspath(directory)
    database = os.path.basename(directory_name)
    directory_name = os.path.dirname(directory_name)
    molecule = os.path.basename(directory_name)
    same_molecule = copy.deepcopy(molecule)  # Make a copy of the string because we'll need it for the isotopologue
    molecule = re.sub('[  ~].+', '', molecule)  # Keep molecule part of the folder name
    isotopologue = re.sub('.+[  ~]', '', same_molecule) # Keep isotope part of the folder name
    
    return molecule, isotopologue, database
    
    
def check_molecule(molecule):
    match = re.match('^[A-Z]{1}[a-z]?$', molecule)     # Matches a string containing only 1 capital letter followed by 0 or 1 lower case letters
    
    if match: return False   # If our 'molecule' matches the pattern, it is really an atom
    else: return True        # We did not get a match, therefore must have a molecule


def mass(molecule, isotopologue, linelist):
    if linelist == 'hitran' or linelist == 'hitemp':
        mol_ID = 1
        while moleculeName(mol_ID) != molecule:
            mol_ID += 1
            
        iso_ID = 1
        while True:
            iso_name = isotopologueName(mol_ID, iso_ID) # Need to format the isotopologue name to match ExoMol formatting
    
            # 'H' not followed by lower case letter needs to become '(1H)'
            iso_name = re.sub('H(?![a-z])', '(1H)', iso_name)
    
            # Number of that atom needs to be enclosed by parentheses ... so '(1H)2' becomes '(1H2)'
            matches = re.findall('[)][0-9]{1}', iso_name)
            for match in matches:
                number = re.findall('[0-9]{1}', match)
                iso_name = re.sub('[)][0-9]{1}', number[0] + ')', iso_name)
    
            # replace all ')(' with '-'
            iso_name = iso_name.replace(')(', '-')
            
            if iso_name == isotopologue:
                return molecularMass(mol_ID, iso_ID)
            
            else:
                iso_ID += 1

        
    else:
        isotopologue = isotopologue.replace('(', '')
        isotopologue = isotopologue.replace(')', '')
        url = 'http://exomol.com/data/molecules/' + molecule + '/' + isotopologue + '/' + linelist + '/'
        
        # Parse the webpage to find the .def file and read it
        web_content = requests.get(url).text
        soup = BeautifulSoup(web_content, "lxml")
        def_tag = soup.find('a', href = re.compile("def"))
        new_url = 'http://exomol.com' + def_tag.get('href')
        
        out_file = './def'
        with requests.get(new_url, stream=True) as request:
            with open(out_file, 'wb') as file:
                for chunk in request.iter_content(chunk_size = 1024 * 1024):
                    file.write(chunk)
                    
        data = pd.read_csv(out_file, delimiter = '#', names = ['Value', 'Key'])  # Store the .def file in a pandas DataFrame
        data = data[data['Key'].str.contains('mass')]  # Only use the row that contains the isotopologue mass
        data = data.reset_index(drop = True)  # Reset the index of the DataFrame
        mass = data['Value'][0]
        mass = re.findall('[0-9|.]+', mass)[0]
        
        os.remove(out_file)
        return float(mass)
        
    
def load_ExoMol(input_directory):
    
    # Read in states file (EXOMOL only)
    states_file_name = [filename for filename in os.listdir(input_directory) if filename.endswith('.states')]
    states_file = pd.read_csv(input_directory + states_file_name[0], sep = '\s+', header=None)
    E = np.array(states_file[1])
    g = np.array(states_file[2])
    J = np.array(states_file[3]).astype(np.int64)
    
    del states_file  # Delete file to free up memory    
    
    return E, g, J

def load_VALD():
    # Needs work
    nu_0, gf, E_low, E_up, J_low, l_low, l_up, gamma_nat, gamma_vdw = (np.array([]) for _ in range(9))
    
    for i in range(len(linelist_files)):
        
        # Read in VALD transitions files (also contains broadening parameters)
        trans_file = pd.read_csv(input_directory + linelist_files[i], sep = ' ', header=None, skiprows=1,
                                 dtype={'nu_0': np.float64, 'gf': np.float64, 
                                        'E_low': np.float64, 'E_up': np.float64,
                                        'J_low': np.float64, 'J_up': np.float64,
                                        'l_low': np.int64, 'l_up': np.int64, 
                                        'gamma_nat': np.float64, 
                                        'gamma_vdw': np.float64})
    
        # Merge linelists
        nu_0 = np.append(nu_0, np.array(trans_file[0]))
        gf = np.append(gf, np.array(trans_file[1]))
        E_low = np.append(E_low, np.array(trans_file[2]))
        E_up = np.append(E_up, np.array(trans_file[3]))
        J_low = np.append(J_low, np.array(trans_file[4])).astype(np.int64)
        l_low = np.append(l_low, np.array(trans_file[6]))
        l_up = np.append(l_up, np.array(trans_file[7]))
        gamma_nat = np.append(gamma_nat, np.array(trans_file[8]))
        gamma_vdw = np.append(gamma_vdw, np.array(trans_file[9]))
        
        # If transitions are not in increasing wavenumber order, rearrange
        order = np.argsort(nu_0)  # Indices of nu_0 in increasing order
        nu_0 = nu_0[order]
        gf = gf[order]
        E_low = E_low[order]
        E_up = E_up[order]
        J_low = J_low[order]
        l_low = l_low[order]
        l_up = l_up[order]
        gamma_nat = gamma_nat[order]
        gamma_vdw = gamma_vdw[order]
        
        del trans_file  # Delete file to free up memory
        
    return


def load_pf(input_directory):
    print("Loading partition functions")
    pf_file_name = [filename for filename in os.listdir(input_directory) if filename.endswith('.pf')]
    pf_file = pd.read_csv(input_directory + pf_file_name[0], sep= ' ', header=None, skiprows=1)
    T_pf_raw = np.array(pf_file[0]).astype(np.float64)
    Q_raw = np.array(pf_file[1])

    del pf_file   # Delete file to free up memory
    
    return T_pf_raw, Q_raw


def det_broad(input_directory):
    if 'H2.broad' in os.listdir(input_directory) and 'He.broad' in os.listdir(input_directory):
        broadening = 'H2-He'
        
    elif 'air.broad' in os.listdir(input_directory):
        broadening = 'air'
        
    else:
        broadening = 'Burrows'
        create_Burrows(input_directory)
        # To do: Create a Burrows broadening file and add it to directory
        
    return broadening


def create_Burrows(input_directory):
    burrows_file = input_directory + 'Burrows.broad'
    J = np.arange(31.0)
    gamma_L_0 = np.zeros(31)
    N_L = np.zeros(31)
    
    for i in range(len(J)):
        gamma_L_0[i] = (0.1 - min(J[i], 30) * 0.002) / (1.01325 * 2) # Convert from cm^-1 / atm -> cm^-1 / bar and take width at half-max

    f_out = open(burrows_file, 'w')
    f_out.write('J | gamma_L_0 | n_L \n')
    for i in range(len(J)):
        f_out.write('%.1f %.4f %.3f \n' %(J[i], gamma_L_0[i], N_L[i]))
        
    f_out.close()


def read_H2_He(input_directory):
    
    # Read in H2 broadening file
    broad_file_H2 = pd.read_csv(input_directory + 'H2.broad', sep = ' ', header=None, skiprows=1)
    J_max_H2 = int(np.max(np.array(broad_file_H2[0])))
    gamma_0_H2 = np.array(broad_file_H2[1])
    n_L_H2 = np.array(broad_file_H2[2])
    
    # Read in He broadening file
    broad_file_He = pd.read_csv(input_directory + 'He.broad', sep = ' ', header=None, skiprows=1)
    J_max_He = int(np.max(np.array(broad_file_He[0])))
    gamma_0_He = np.array(broad_file_He[1])
    n_L_He = np.array(broad_file_He[2])

    # Take maximum J'' value for which broadening is a function of J to be lowest for which complete data available
    J_max = np.max(np.array([J_max_H2, J_max_He]))

    # If broadening files not of same length, extend shortest file to same length as longest
    if (J_max_H2 < J_max):
            
        for i in range (J_max_H2, J_max):
                
            gamma_0_H2 = np.append(gamma_0_H2, gamma_0_H2[-1])    # Extended values equal to final value 
            n_L_H2 = np.append(n_L_H2, n_L_H2[-1])                # Extended values equal to final value 
            
    if (J_max_He < J_max):
        
        for i in range (J_max_He, J_max):
                
            gamma_0_He = np.append(gamma_0_He, gamma_0_He[-1])    # Extended values equal to final value 
            n_L_He = np.append(n_L_He, n_L_He[-1])                # Extended values equal to final value 
                
    del broad_file_H2, broad_file_He   # Delete files to free up memory
    
    return J_max, gamma_0_H2, n_L_H2, gamma_0_He, n_L_He
    
    
def read_air(input_directory):
    
    # Read in air broadening file
    broad_file_air = pd.read_csv(input_directory + 'air.broad', sep = ' ', header=None, skiprows = 1)
    J_max = int(np.max(np.array(broad_file_air[0])))
    gamma_0_air = np.array(broad_file_air[1])
    n_L_air = np.array(broad_file_air[2])
    
    del broad_file_air   # Delete file to free up memory  
    
    return J_max, gamma_0_air, n_L_air
    

def read_Burrows(input_directory):
    
    # Read in Burrows broadening file
    broad_file_Burrows = pd.read_csv(input_directory + 'Burrows.broad', sep = ' ', header=None, skiprows=1)
    J_max = int(np.max(np.array(broad_file_Burrows[0])))
    gamma_0_Burrows = np.array(broad_file_Burrows[1])
    #n_L_Burrows = np.array(broad_file_Burrows[2])       # Not really needed, as temperature exponent = 0 for all J''
    
    del broad_file_Burrows   # Delete file to free up memory
    
    return J_max, gamma_0_Burrows


def interpolate_pf(T_pf_raw, Q_raw, T, T_ref):
    
    #***** Interpolate (and extrapolate) partition function to a fine grid *****#
    pf_spline = Interp(T_pf_raw, Q_raw, k=5)
    T_pf_fine = np.linspace(1.0, 10000.0, 9999)       # New temperature grid (extrapolated to 4000K)
    Q_fine = pf_spline(T_pf_fine)                    # Extrapolated partition function
    i_T = np.argmin(np.abs(T_pf_fine - T))           # Index of partition function temperature array closest to desired temperature
    i_T_ref = np.argmin(np.abs(T_pf_fine - T_ref))   # Index of partition function temperature array closest to reference temperature
    Q_T = Q_fine[i_T]                                # Partition function at given temperature, Q(T)
    Q_T_ref = Q_fine[i_T_ref]                        # Partition function at reference temperature, Q(T_ref)
    
    return Q_T, Q_T_ref

def compute_pressure_broadening_atom():
    
    if (species in ['Li', 'Na','K', 'Rb', 'Cs']):  # Special treatments for alkali van der waals widths
                
        gamma_0_H2 = np.zeros(len(nu_0))
        gamma_0_He = np.zeros(len(nu_0))
        n_L_H2 = np.zeros(len(nu_0))
        n_L_He = np.zeros(len(nu_0))
                
        for i in range(len(nu_0)):
                    
            if (gamma_vdw[i] != 0.0):  # For transitions with a VALD broadening value 
                        
                gamma_0_H2[i], n_L_H2[i] = gamma_L_VALD(gamma_vdw[i], (m/u), 'H2')
                gamma_0_He[i], n_L_He[i] = gamma_L_VALD(gamma_vdw[i], (m/u), 'He')
                    
            elif (gamma_vdw[i] == 0.0):  # For transitions without a VALD broadening value 
                        
                gamma_0_H2[i], n_L_H2[i] = gamma_L_impact(E_low[i], E_up[i], l_low[i], l_up[i], species, (m/u), 'H2')
                gamma_0_He[i], n_L_He[i] = gamma_L_impact(E_low[i], E_up[i], l_low[i], l_up[i], species, (m/u), 'He')
                        
                # Remove transitions where Hydrogenic approximation breaks down
          #      accept_condition = np.where(gamma_0_H2 != -1.0)  # Transitions without a VALD broadening value
                        
          #      nu_0 = nu_0[accept_condition]
          #      gf = gf[accept_condition]
          #      E_low = E_low[accept_condition]
          #      E_up = E_up[accept_condition]
          #      J_low = J_low[accept_condition]
          #      l_low = l_low[accept_condition]
          #      l_up = l_up[accept_condition]
          #      n_L_H2 = n_L_H2[accept_condition]
          #      n_L_He = n_L_He[accept_condition]
          #      gamma_0_He = gamma_0_He[accept_condition]
          #      gamma_0_H2 = gamma_0_H2[accept_condition]
                
    else:  # For non-alkali species
                
          #      accept_condition = np.where(log_gamma_vdw != 0.0)  # Transitions without a VALD broadening value 
                
          #      log_gamma_vdw = log_gamma_vdw[accept_condition]
          #      nu_0 = nu_0[accept_condition]
          #      gf = gf[accept_condition]
          #      E_low = E_low[accept_condition]
          #      E_up = E_up[accept_condition]
          #      J_low = J_low[accept_condition]
          #      l_low = l_low[accept_condition]
          #      l_up = l_up[accept_condition]
                
        gamma_0_H2, n_L_H2 = gamma_L_VALD(gamma_vdw, (m/u), 'H2')
        gamma_0_He, n_L_He = gamma_L_VALD(gamma_vdw, (m/u), 'He')
        
    gamma += ((1.0/(4.0*np.pi*(100.0*c))) * gamma_nat)  # Add natural line widths


def compute_H2_He_broadening(gamma_0_H2, T_ref, T, n_L_H2, P, P_ref, X_H2, gamma_0_He, n_L_He, X_He):
    gamma = (gamma_0_H2 * np.power((T_ref/T), n_L_H2) * (P/P_ref) * X_H2 +   # H2+He Lorentzian HWHM for given T, P, and J (ang. mom.)
             gamma_0_He * np.power((T_ref/T), n_L_He) * (P/P_ref) * X_He)    # Note that these are only a function of J''

    return gamma
    
def compute_air_broadening(gamma_0_air, T_ref, T, n_L_air, P, P_ref):
    gamma = (gamma_0_air * np.power((T_ref/T), n_L_air) * (P/P_ref))      # Air-broadened Lorentzian HWHM for given T, P, and J (ang. mom.)

    return gamma
    
def compute_Burrows_broadening(gamma_0_Burrows, P, P_ref):
    gamma = (gamma_0_Burrows * (P/P_ref))      # Equation (15) in Sharp & Burrows (2007)  
    
    return gamma

def create_wavelength_grid_atom():
    # First, we need to find values of gamma_V for reference wavenumber (1000 cm^-1)
    alpha_ref = np.sqrt(2.0*kb*T*np.log(2)/m) * (np.array(nu_ref[2])/c) # Doppler HWHM at reference wavenumber
    gamma_ref = np.min(gamma)                                           # Find minimum value of Lorentzian HWHM
    gamma_V_ref = Voigt_width(gamma_ref, alpha_ref)                     # Reference Voigt width
    
    # Calculate Voigt width for each transition
    alpha = np.sqrt(2.0*kb*T*np.log(2)/m) * (np.array(nu_0)/c)   # Doppler HWHM for each transition
    gamma_V = Voigt_width(gamma, alpha)
    
    # Now compute properties of computational (fine) and output (coarse) wavenumber grid
    
    # Wavenumber spacing of of computational grid (smallest of gamma_V_ref/6 or 0.01cm^-1)
    dnu_fine = np.minimum(gamma_V_ref*Voigt_sub_spacing, dnu_out)      
    
    # Number of points on fine grid (rounded)
    N_points_fine = int((nu_max-nu_min)/dnu_fine + 1)
    
    # Adjust dnu_fine slightly to match an exact integer number of grid spaces
    dnu_fine = (nu_max-nu_min)/(N_points_fine - 1)
    
    cutoffs = np.zeros(len(nu_0))   # Line wing cutoffs for each line
    
    # Line cutoffs at @ min(500 gamma_V, 1000cm^-1)
    for i in range(len(nu_0)):
        cutoffs[i] = dnu_fine * (int((Voigt_cutoff*gamma_V[i])/dnu_fine))

        if (cutoffs[i] >= cut_max): cutoffs[i] = cut_max
                
        # Special cases for alkali resonant lines
        if   ((species == 'Na') and (int(nu_0[i]) in [16978, 16960])): cutoffs[i] = 9000.0   # Cutoff @ +/- 9000 cm^-1
        elif ((species == 'K') and  (int(nu_0[i]) in [13046, 12988])): cutoffs[i] = 9000.0   # Cutoff @ +/- 9000 cm^-1
        #elif ((species == 'Li') and (int(nu_0[i]) in [14908, 14907])): cutoffs[i] = 9000.0   # Cutoff @ +/- 9000 cm^-1
        #elif ((species == 'Rb') and (int(nu_0[i]) in [12820, 12582])): cutoffs[i] = 9000.0   # Cutoff @ +/- 9000 cm^-1
        #elif ((species == 'Cs') and (int(nu_0[i]) in [11735, 9975])):  cutoffs[i] = 9000.0   # Cutoff @ +/- 9000 cm^-1
            
        # Calculate detuning frequencies for Na and K resonance lines
        if (species == 'Na'): nu_detune = 30.0 * np.power((T/500.0), 0.6)
        elif (species == 'K'): nu_detune = 20.0 * np.power((T/500.0), 0.6)
        else: nu_detune = cut_max
        
        # Evaluate number of frequency points for each Voigt function up to cutoff (one tail)
        N_Voigt_points = ((cutoffs/dnu_fine).astype(np.int64)) + 1  
        
        # Define start and end points of fine grid
        #nu_fine = np.linspace(nu_min, nu_max, N_points_fine)
        nu_fine_start = nu_min
        nu_fine_end = nu_max
        
        # Initialise output grid
        N_points_out = int((nu_out_max-nu_out_min)/dnu_out + 1)     # Number of points on coarse grid (uniform)
        nu_out = np.linspace(nu_out_min, nu_out_max, N_points_out)  # Create coarse (output) grid
        
        # Initialise cross section arrays on each grid
        sigma_fine = np.zeros(N_points_fine)    # Computational (fine) grid
        sigma_out = np.zeros(N_points_out)      # Coarse (output) grid

def create_wavelength_grid_molecule(nu_ref, m, T, gamma, Voigt_sub_spacing, dnu_out, cut_max, Voigt_cutoff, nu_out_max, nu_out_min, N_alpha_samples):
    #nu_min = max(1.0, (nu_out_min - cut_max))
    #nu_max = nu_out_max + cut_max
    
    nu_min = 1
    nu_max = 30000

    # First, we need to find values of gamma_V for reference wavenumbers (100, 1000, 10000 cm^-1)
    nu_ref = np.array(nu_ref)   # The reference wavenumbers need to be a numpy array
    alpha_ref = np.sqrt(2.0*kb*T*np.log(2)/m) * (nu_ref/c)    # Doppler HWHM at reference wavenumbers 
    gamma_ref = np.min(gamma)                                           # Find minimum value of Lorentzian HWHM
    gamma_V_ref = Voigt_width(gamma_ref, alpha_ref)                     # Reference Voigt widths
    
    # Now compute properties of computational (fine) and output (coarse) wavenumber grids
            
    # Wavenumber spacing of three regions of computational grid (smallest of gamma_V_ref/6 or 0.01cm^-1)
    dnu_fine = np.minimum(gamma_V_ref*Voigt_sub_spacing, dnu_out*np.ones(3))
    
    # Number of points on fine grid (three regions, rounded)
    N_points_fine_1 = int((nu_ref[1]-nu_min)/dnu_fine[0])
    N_points_fine_2 = int((nu_ref[2]-nu_ref[1])/dnu_fine[1])
    N_points_fine_3 = int((nu_max-nu_ref[2])/dnu_fine[2] + 1)
    N_points_fine = N_points_fine_1 + N_points_fine_2 + N_points_fine_3
    
    # Adjust dnu_fine slightly to match an exact integer number of grid spaces
    dnu_fine[0] = (nu_ref[1]-nu_min)/N_points_fine_1
    dnu_fine[1] = (nu_ref[2]-nu_ref[1])/N_points_fine_2
    dnu_fine[2] = (nu_max-nu_ref[2])/(N_points_fine_3 - 1)
    cutoffs = np.zeros(len(dnu_fine))   # Line wing cutoffs in three regions
    
    # Line cutoffs at @ min(500 gamma_V, 30cm^-1)
    for i in range(len(dnu_fine)):     
        # If grid spacing below maximum dnu (0.01), set to 3000 dnu (~ 500 gamma_V)
        if (dnu_fine[i] < dnu_out):
            cutoffs[i] = Voigt_cutoff*(1.0/Voigt_sub_spacing)*dnu_fine[i] 
            
        # If at maximum dnu (0.01), set to min(500 gamma_V, 30cm^-1)
        else:
            cutoffs[i] = dnu_fine[i] * int(Voigt_cutoff*gamma_V_ref[i]/dnu_fine[i])
            if (cutoffs[i] >= cut_max): cutoffs[i] = cut_max
            
    # Define start and end points of each fine grid
    nu_fine_1_start = nu_min
    nu_fine_1_end = (nu_ref[1] - dnu_fine[0])
    nu_fine_2_start = nu_ref[1]
    nu_fine_2_end = (nu_ref[2] - dnu_fine[1])
    nu_fine_3_start = nu_ref[2]
    nu_fine_3_end = nu_max
    
    #nu_fine_1 = np.linspace(nu_min, nu_ref[1], N_points_fine_1, endpoint=False)
    #nu_fine_2 = np.linspace(nu_ref[1], nu_ref[2], N_points_fine_2, endpoint=False)
    #nu_fine_3 = np.linspace(nu_ref[2], nu_max, N_points_fine_3)
    #nu_fine = np.concatenate([nu_fine_1, nu_fine_2, nu_fine_3])
    
    # Initialise output grid
    N_points_out = int((nu_out_max-nu_out_min)/dnu_out + 1)     # Number of points on coarse grid (uniform)
    nu_out = np.linspace(nu_out_min, nu_out_max, N_points_out)  # Create coarse (output) grid
        
    # Initialise cross section arrays on each grid
    sigma_fine = np.zeros(N_points_fine)    # Computational (fine) grid
    sigma_out = np.zeros(N_points_out)      # Coarse (output) grid
    
    return (sigma_fine, nu_ref, N_points_fine_1, N_points_fine_2, N_points_fine_3, dnu_fine, cutoffs,
            nu_out, sigma_out, nu_fine_1_start, nu_fine_1_end, nu_fine_2_start, nu_fine_2_end, 
            nu_fine_3_start, nu_fine_3_end, N_points_out, nu_min, nu_max, alpha_ref)
    
    
def precompute_Voigt_profiles(nu_ref, nu_max, N_alpha_samples, T, m, cutoffs, dnu_fine, gamma,
                              alpha_ref):
    
    #***** Pre-compute Voigt function array for molecules *****#
    
    print('Pre-computing Voigt profiles...')
    
    t1 = time.perf_counter()
        
    # First, create an array of values of alpha to approximate true values of alpha (N=500 log-spaced => max error of 0.5%)
    nu_sampled = np.logspace(np.log10(nu_ref[0]), np.log10(nu_max), N_alpha_samples)
    alpha_sampled = np.sqrt(2.0*kb*T*np.log(2)/m) * (nu_sampled/c)
    
    # Evaluate number of frequency points for each Voigt function in each spectral region - up to cutoff @ min(500 gamma_V, 30cm^-1)
    N_Voigt_points = ((cutoffs/dnu_fine).astype(np.int64)) + 1  
            
    # Pre-compute and store Voigt functions and first derivatives wrt alpha 
    Voigt_arr = np.zeros(shape=(len(gamma), len(alpha_sampled), np.max(N_Voigt_points)))    # For H2O: V(51,500,3001)
    dV_da_arr = np.zeros(shape=(len(gamma), len(alpha_sampled), np.max(N_Voigt_points)))    # For H2O: V(51,500,3001)
    dV_dnu_arr = np.zeros(shape=(len(gamma), len(alpha_sampled), np.max(N_Voigt_points)))   # For H2O: V(51,500,3001)
    Generate_Voigt_grid_molecules(Voigt_arr, dV_da_arr, dV_dnu_arr, gamma, alpha_sampled, alpha_ref, cutoffs, N_Voigt_points)
        
    t2 = time.perf_counter()
    total1 = t2-t1
            
    print('Voigt profiles computed in ' + str(total1) + ' s')       

    return nu_sampled, alpha_sampled, Voigt_arr, dV_da_arr, dV_dnu_arr, N_Voigt_points
        
        
def write_output_file(cluster_run, output_directory, molecule, T_arr, t, log_P_arr, p, nu_out, sigma_out):
    #***** Now write output files *****#
            
    if (cluster_run == False): f = open(output_directory + str(molecule) + '_T' + str(T_arr[t]) + 'K_log_P' + str(log_P_arr[p]) + '_sigma.txt','w')
    #elif (cluster_run == True): f = open(output_directory + str(molecule) + '_T' + str(T) + 'K_log_P' + str(log_P) + '_sigma.txt','w')
                    
    for i in range(len(nu_out)):
        f.write('%.8f %.8e \n' %(nu_out[i], sigma_out[i]))
                        
    f.close()
    
def replace_iso_name(iso_name):
    # 'H' not followed by lower case letter needs to become '(1H)'
    iso_name = re.sub('H(?![a-z])', '(1H)', iso_name)
    
    # Number of that atom needs to be enclosed by parentheses ... so '(1H)2' becomes '(1H2)'
    matches = re.findall('[)][0-9]{1}', iso_name)
    for match in matches:
        number = re.findall('[0-9]{1}', match)
        iso_name = re.sub('[)][0-9]{1}', number[0] + ')', iso_name)
    
    # replace all ')(' with '-'
    iso_name = iso_name.replace(')(', '-')   
    
    return iso_name
    
    
def find_input_dir(input_dir, database, molecule, isotope, linelist):
    
    if isotope == 'default':
        if database == 'exomol':
            isotope = '(' + get_default_iso(molecule) + ')'
        if database == 'hitran' or database == 'hitemp':
            molecule_dict = create_id_dict()
            mol_id = molecule_dict.get(molecule)
            isotope = isotopologueName(mol_id, 1)
            isotope = replace_iso_name(isotope)
    
    if linelist == 'default':
        if database == 'exomol':
            temp_isotope = re.sub('[(]|[)]', '', isotope)
            linelist = get_default_linelist(molecule, temp_isotope)
        if database == 'hitran':
            linelist = 'HITRAN'
        if database == 'hitemp':
            linelist = 'HITEMP'
            
    input_directory = input_dir + '/' + molecule + '  ~  ' + isotope + '/' + linelist + '/'
    
    if os.path.exists(input_directory):
        return input_directory
    
    else:
        print("You don't seem to have a local folder with the parameters you entered.\n") 
        
        if not os.path.exists(input_dir + '/'):
            print("----- You entered an invalid input directory into the cross_section() function. Please try again. -----")
            sys.exit(0)
        
        elif not os.path.exists(input_dir + '/' + molecule + '  ~  ' + isotope + '/'):
            print("----- There was an error with the molecule + isotopologue you entered. Here are the available options: -----\n")
            for folder in os.listdir(input_dir + '/'):
                if not folder.startswith('.'):
                    print(folder)
            sys.exit(0)
        
        else:
            print("There was an error with the line list. These are the linelists available: \n")
            for folder in os.listdir(input_dir + '/' + molecule + '  ~  ' + isotope + '/'):
                if not folder.startswith('.'):
                    print(folder)
            sys.exit(0)

    
def create_cross_section(input_dir, database, molecule, log_pressure, temperature, isotope = 'default', 
                         linelist = 'default', cluster_run = False, nu_out_min = 200, nu_out_max = 25000, 
                         dnu_out = 0.01, pressure_broadening = 'default', X_H2 = 0.85, X_He = 0.15, 
                         Voigt_cutoff = 500, Voigt_sub_spacing = (1.0/6.0), N_alpha_samples = 500, 
                         S_cut = 1.0e-100, cut_max = 30.0, **kwargs):
    
    print("Beginning cross-section computations...")
    
    database = database.lower()
    
    input_directory = find_input_dir(input_dir, database, molecule, isotope, linelist)
    
    # Use the input directory to define these right at the start
    molecule, isotopologue, database = parse_directory(input_directory)
    if database.lower() != 'hitran' and database.lower() != 'hitemp':
        linelist = database
        database = 'exomol'
    else:
        database = database.lower()
        linelist = database
    
    linelist_files = [filename for filename in os.listdir(input_directory) if filename.endswith('.h5')]
    
    if database == 'exomol':
        print("Loading ExoMol format")
        E, g, J = load_ExoMol(input_directory)
    
    if database == 'hitran':
        print("Loading HITRAN format")
        # Nothing else required
    
    if database == 'hitemp':
        print("Loading HITEMP format")
        # Nothing else required
        
    if database == 'vald':
        print("Loading VALD format")
        # load_VALD or something
    
    T_pf_raw, Q_raw = load_pf(input_directory)
    

    is_molecule = check_molecule(molecule)
    
    if is_molecule and pressure_broadening == 'default':
        broadening = det_broad(input_directory)
        if broadening == 'H2-He':
            J_max, gamma_0_H2, n_L_H2, gamma_0_He, n_L_He = read_H2_He(input_directory)
            
        elif broadening == 'air':
            J_max, gamma_0_air, n_L_air = read_air(input_directory)
            
        elif broadening == 'Burrows':
            J_max, gamma_0_Burrows = read_Burrows(input_directory)
            
    elif is_molecule and pressure_broadening != 'default':
        broadening = pressure_broadening
        if broadening == 'H2-He' and 'H2.broad' in os.listdir(input_directory) and 'He.broad' in os.listdir(input_directory):
            J_max, gamma_0_H2, n_L_H2, gamma_0_He, n_L_He = read_H2_He(input_directory)
        
        elif broadening == 'air' and 'air.broad' in os.listdir(input_directory):
            broadening = 'air'
            J_max, gamma_0_air, n_L_air = read_air(input_directory)
            
        elif broadening == 'Burrows':
            create_Burrows(input_directory)
            J_max, gamma_0_Burrows = read_Burrows(input_directory)
            
        else:
            print("\nYou did not enter a valid type of pressure broadening. Please try again.")
            sys.exit(0)
            
    # Will need an 'else' statement for atomic broadening
    

    # Start clock for timing program
    t_start = time.perf_counter()
    
    # For a single point in P-T space:
    if (cluster_run == False):
        log_P_arr = np.array([log_pressure])        # log_10 (Pressure/bar)
        T_arr = np.array([temperature])         # Temperature (K)
        N_P = len(log_P_arr)
        N_T = len(T_arr)

    #***** Load pressure and temperature for this calculation *****#
    P_arr = np.power(10.0, log_P_arr)   

    for p in range(N_P):
        for t in range(N_T):
            if (cluster_run == False):
            
                P = P_arr[p]   # Atmospheric pressure (bar)
                T = T_arr[t]   # Atmospheric temperature (K)
                
            Q_T, Q_T_ref = interpolate_pf(T_pf_raw, Q_raw, T, T_ref)
            
            m = mass(molecule, isotopologue, linelist) * u
            
            if is_molecule:
                
                if broadening == 'H2-He':
                    gamma = compute_H2_He_broadening(gamma_0_H2, T_ref, T, n_L_H2, P, P_ref, X_H2, gamma_0_He, n_L_He, X_He)
                
                elif broadening == 'air':
                    gamma = compute_air_broadening(gamma_0_air, T_ref, T, n_L_air, P, P_ref)
                    
                elif broadening == 'Burrows':
                    gamma = compute_Burrows_broadening(gamma_0_Burrows, P, P_ref)
                    
                (sigma_fine, nu_ref, N_points_fine_1, N_points_fine_2, 
                 N_points_fine_3, dnu_fine, cutoffs, nu_out, sigma_out, 
                 nu_fine_1_start, nu_fine_1_end, nu_fine_2_start, 
                 nu_fine_2_end, nu_fine_3_start, nu_fine_3_end, 
                 N_points_out, nu_min, nu_max, alpha_ref) = create_wavelength_grid_molecule(Nu_ref, m, T, gamma, Voigt_sub_spacing, 
                                                                              dnu_out, cut_max, Voigt_cutoff, nu_out_max, 
                                                                              nu_out_min, N_alpha_samples)
                                                                                               
                (nu_sampled, alpha_sampled, Voigt_arr, 
                 dV_da_arr, dV_dnu_arr, N_Voigt_points) = precompute_Voigt_profiles(nu_ref, nu_max, N_alpha_samples, T, m, cutoffs, dnu_fine, gamma, alpha_ref)                                                                            
                
            else:
                
                compute_pressure_broadening_atom()
                create_wavelength_grid_atom()
                
            print("Pre-computation complete")
                
            
            print('Generating cross section for ' + molecule + ' at P = ' + str(P) + ' bar, T = ' + str(T) + ' K')
            
            if database == 'exomol':                
                produce_total_cross_section_EXOMOL(linelist_files, input_directory, sigma_fine,
                                           nu_sampled, nu_ref, m, T, Q_T, N_points_fine_1,
                                           N_points_fine_2, N_points_fine_3, dnu_fine,
                                           N_Voigt_points, cutoffs, g, E, J, J_max, 
                                           alpha_sampled, Voigt_arr, dV_da_arr, dV_dnu_arr,
                                           nu_min, nu_max, S_cut)
                
            elif database == 'hitran':
                produce_total_cross_section_HITRAN(linelist_files, input_directory, sigma_fine,
                                           nu_sampled, nu_ref, m, T, Q_T, Q_T_ref,
                                           N_points_fine_1, N_points_fine_2, N_points_fine_3,
                                           dnu_fine, N_Voigt_points, cutoffs, J_max, 
                                           alpha_sampled, Voigt_arr, dV_da_arr, dV_dnu_arr,
                                           nu_min, nu_max, S_cut)
            """    
            elif database == 'vald':
        
                if is_molecule:
                    
                    produce_total_cross_section_VALD_molecule(sigma_fine, nu_sampled, nu_ref, nu_0, E_low, J_low,
                                                      gf, m, T, Q_T, N_points_fine_1, N_points_fine_2,
                                                      N_points_fine_3, dnu_fine, N_Voigt_points, cutoffs, 
                                                      J_max, alpha_sampled, Voigt_arr, dV_da_arr, dV_dnu_arr,
                                                      nu_min, nu_max, S_cut, molecule)
                
                else:
                    
                    produce_total_cross_section_VALD_atom(sigma_fine, nu_0, nu_detune, E_low, gf, m, T, Q_T,
                                                  N_points_fine, N_Voigt_points, alpha, gamma, cutoffs,
                                                  nu_min, nu_max, S_cut, molecule)
                    
            """
            
            # Now bin cross section to output grid    
            print('Binning cross section to output grid...')
        
            if is_molecule:
                bin_cross_section_molecule(sigma_fine, sigma_out, N_points_fine_1, N_points_fine_2,
                                   N_points_fine_3, nu_ref, nu_fine_1_start, nu_fine_1_end,
                                   nu_fine_2_start, nu_fine_2_end, nu_fine_3_start, nu_fine_3_end,
                                   nu_out, N_points_out, 0, nu_min, nu_max)
            
            """
            else:
                bin_cross_section_atom(sigma_fine, sigma_out, nu_fine_start, 
                               nu_fine_end, nu_out, N_points_fine, N_points_out, 0)
            """
            
            #bin_cross_section(sigma_fine, sigma_out_log, nu_fine_1, nu_fine_2, nu_fine_3, nu_out, N_points_out, 1)
        
        t_final = time.perf_counter()
        total_final = t_final-t_start
        
        print('Total runtime: ' + str(total_final) + ' s')
        
        output_directory = re.sub('/input/', '/output/', input_directory)
        
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)
        
        write_output_file(cluster_run, output_directory, molecule, T_arr, t, log_P_arr, p, nu_out, sigma_out)
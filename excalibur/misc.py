import sys
import numpy as np
import pandas as pd
import re
import os
import csv

import excalibur.HITRAN as HITRAN
import excalibur.ExoMol as ExoMol

from .hapi import isotopologueName

def check_molecule(molecule):
    """
    Check if the given string is a molecule

    Parameters
    ----------
    molecule : String
        Molecule name (e.g. H2O).

    Returns
    -------
    True if the given string is a molecule, false otherwise (if it is an atom).

    """
    match = re.match('^[A-Z]{1}[a-z]?$', molecule)     # Matches a string containing only 1 capital letter followed by 0 or 1 lower case letters
    
    if match: return False   # If our 'molecule' matches the pattern, it is really an atom
    else: return True        # We did not get a match, therefore must have a molecule

    
def write_output(output_directory, species, roman_num, T, log_P, broad_type, broadening_file, nu_out, sigma_out):
    """
    Parameters
    ----------
    output_directory : String
        Local directory where the output data is to be stored.
    species : String
        Name of molecule or atomic.
    roman_num : String
        Ionization state, in case the species is an atom. Depicted in Roman numeral form.
    T : int
        Temperature (K) the cross-section was computed at.
    log_P : int
        DESCRIPTION.
    broad_type : String
        The type of broadening used in computing the cross section.
    nu_out : TYPE
        DESCRIPTION.
    sigma_out : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
       
    # Add ionisation state for atoms     
    if (roman_num != ''):
        species = species + '_' + roman_num
        
    if broadening_file == '':
        f = open((output_directory + species + '_T' + str(T) + 'K_log_P' + str(log_P) + '_' + broad_type + '_sigma.txt'),'w')
    else:
        f = open((output_directory + species + '_T' + str(T) + 'K_log_P' + str(log_P) + '_' + broad_type + 
                    '_' + os.path.splitext(broadening_file)[0] + '_sigma.txt'),'w')
                    
    for i in range(len(nu_out)):
        f.write('%.8f %.8e \n' %(nu_out[i], sigma_out[i]))
                        
    f.close()

    
def round_sig_figs(value, sig_figs):
    ''' 
    Round a quantity to a specified number of significant figures.
    '''
    
    return round(value, sig_figs - int(np.floor(np.log10(abs(value)))) - 1)


def int_to_roman(num):
    '''
    Convert an integer to the corresponding Roman numeral.
    '''

    # Define the values and symbols for Roman numerals
    val = [
        1000, 900, 500, 400,
        100, 90, 50, 40,
        10, 9, 5, 4,
        1
    ]
    syb = [
        "M", "CM", "D", "CD",
        "C", "XC", "L", "XL",
        "X", "IX", "V", "IV",
        "I"
    ]
    
    # Initialize an empty string to store the Roman numeral representation
    roman_num = ''
    
    # Start with the largest value and work our way down
    i = 0

    while num > 0:

        # For each value, determine the number of times it can be subtracted from the input number
        for _ in range(num // val[i]):
            roman_num += syb[i]  # Add the corresponding Roman numeral symbol
            num -= val[i]        # Subtract the value from the input number

        i += 1                   # Move to the next value
    
    return roman_num  # Return the Roman numeral representation


def cross_section_collection(new_x, new_y, collection = []):
    '''
    Add new cross section (that is, wavelength and absorption cross section) to the collection

    Parameters
    ----------
    new_x : list
        DESCRIPTION.
    new_y : list
        DESCRIPTION.
    collection : 2d list, optional
        DESCRIPTION. The default is [].

    Returns
    -------
    collection : TYPE
        DESCRIPTION.

    '''

    collection.append([new_x, new_y])

    return collection


def read_cross_section_file(species, database, filename, isotope = 'default',
                            ionization_state = 1, linelist = 'default', output_dir = './output/'):

    """
    Read in a previously computed cross section on a user's machine, and return the 
    wavenumber and absorption cross section columns.

    Parameters
    ----------
    species : String
        Molecule for which the cross section was computed.
    database : String
        Database the line list was downloaded from.
    filename : String
        Full file name of the cross section file(eg. 'this_is_my_file.txt').
    isotope : String
        Isotopologue of the molecule for which the cross-section is to be created. The default is 'default'.
    ionization_state : int
        Ionization state, in case of an atomic species. The default is 1. 
    linelist : String
        Line list that is being used. HITRAN/HITEMP/VALD used as the line list name for these
        databases respectively. ExoMol has its own named line lists. The default is 'default'. 
    output_dir = String
        'Prefix' of the directory containing the cross section files. If the files were downloaded
        using our script with no modifications, output_dir will end in '/output'. The default is './output/'. 

    Returns
    -------
    nu : numpy array
        Wavenumber column of cross section file.
    sigma : numpy array
        Absorption cross section column of cross section file.

    """
    database = database.lower()

    if database in ['hitran', 'hitemp']:
        molecule_dict = HITRAN.create_id_dict()
        mol_id = molecule_dict.get(species)
        if isotope == 'default':
            isotope = isotopologueName(mol_id, 1)
        else:
            isotope = isotopologueName(mol_id, isotope)
        isotope = HITRAN.replace_iso_name(isotope)

    if database == 'exomol':
        if isotope == 'default':
            isotope = ExoMol.get_default_iso(species)

    if database == 'vald':
        ion_roman = int_to_roman(ionization_state)

    if linelist == 'default':
        if database == 'exomol':
            temp_isotope = re.sub('[(]|[)]', '', isotope)
            linelist = ExoMol.get_default_linelist(species, temp_isotope)
        if database == 'hitran':
            linelist = 'HITRAN'
        if database == 'hitemp':
            linelist = 'HITEMP'
        if database == 'vald':
            linelist = 'VALD'

    if database == 'vald':
        tag = '(' + ion_roman + ')'
    else:
        tag = isotope

    if (database == 'exomol'):
        output_directory = (output_dir + species + '  ~  (' + tag + ')/' +
                           'ExoMol' + '/' + linelist + '/')
    else:
        output_directory = (output_dir + species + '  ~  ' + tag + '/' +
                           linelist + '/')

    if os.path.exists(output_directory):
        file  = output_directory + filename
        with open(file, 'r') as f:
            contents = csv.reader(f, delimiter=' ')
            nu = []
            sigma = []
            for cols in contents:
                nu.append(float(cols[0]))
                sigma.append(float(cols[1]))

        nu = np.asarray(nu)
        sigma = np.asarray(sigma)
        return  nu, sigma

    else:
        print("You don't seem to have a local folder with the parameters you entered.\n")

        if not os.path.exists(output_dir + '/'):
            print("----- You entered an invalid output directory into the cross_section() function. Please try again. -----")
            sys.exit(0)

        elif not os.path.exists(output_dir + '/' + species + '  ~  (' + tag + ')/'):
            print("----- There was an error with the molecule + isotope you entered. Here are the available options: -----\n")
            for folder in os.listdir(output_dir + '/'):
                if not folder.startswith('.'):
                    print(folder)
            sys.exit(0)

        else:
            print("There was an error with the line list. These are the linelists available: \n")
            for folder in os.listdir(output_dir + '/' + species + '  ~  (' + tag + ')/'):
                if not folder.startswith('.'):
                    print(folder)
            sys.exit(0)
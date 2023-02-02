import numpy as np
import pandas as pd
import re
import os


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
        Name of molecule or atomc.
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

    
def read_output(output_directory, molecule, T, log_P):
    '''
    Read in wavenumber and cross-section from  file

    Parameters
    ----------
    output_directory : String
        Local directory where the computed cross section is stored.
    molecule : String
        Molecule name.
    T : TYPE
        DESCRIPTION.
    log_P : TYPE
        DESCRIPTION.

    Returns
    -------
    nu : TYPE
        DESCRIPTION.
    sigma : TYPE
        DESCRIPTION.

    '''
    
    file_location = (output_directory + molecule + '_T' + str(T) + 
                     'K_log_P' + str(log_P) + '_sigma_TMP.txt')
    
    file = pd.read_csv(file_location, sep = '\s+', header=None)
    
    nu = np.array(file[0])
    sigma = np.array(file[1])

    return nu, sigma


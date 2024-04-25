import sys
import os
import numpy as np
import h5py
from hapi.hapi import partitionSum, moleculeName, isotopologueName

import Excalibur.downloader as download


def create_id_dict():
    '''
    Recreate the table of molecular names and IDs that HITRAN uses for identification

    Returns
    -------
    molecule_dict : dict
        Dictionary where the keys are HITRAN molecular names and the values are the respective HITRAN molecular IDs,
        as seen on https://hitran.org/lbl/.

    '''
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


def check(molecule, isotope):
    '''
    Checks if the parameters passed into the summon() function by the user are valid to use with the HITEMP database.

    Parameters
    ----------
    molecule : String
        Molecule name (e.g. H2O).
    isotope : int
        Isotopologue number of the molecule (1 = most common isotopologue, 2 = second most common, etc.).

    Returns
    -------
    molecule_id_number : int
        Molecule ID number of the given molecule, as given by the first column in HITEMP here: https://hitran.org/hitemp/.

    '''

    molecule_dict = create_id_dict()
    molecule_id_number = molecule_dict.get(molecule)
    
    if molecule_id_number is None:
        print("\n ----- This molecule/isotopologue ID combination is not valid in HITEMP. Please try calling the summon() function again. Make sure you enter the ID number of the molecule you want. ----- ")
        print("\n ----- A list of supported molecule IDs can be found here: https://hitran.org/hitemp/ -----")
        sys.exit(0)
    
    table = download.HITEMP_table()
    
    # Ensure that the given isotopologue number is less than or equal to the max number of isotopolgues stored by HITRAN
    if molecule_id_number in table['ID'].values:
        row = table.loc[table['ID'] == molecule_id_number]
        isotope_count = row.loc[row.index.values[0], 'Iso Count']
        if isotope <= isotope_count:
            return molecule_id_number
    else:
        print("\n ----- This molecule/isotopologue ID combination is not valid in HITEMP. Please try calling the summon() function again. Make sure you enter the ID number of the molecule you want. ----- ")
        print("\n ----- A list of supported molecule IDs can be found here: https://hitran.org/hitemp/ -----")
        sys.exit(0)
        

def determine_linelist():
    """
    Determines the molecule and isotopologue that would like to be accessed based on user input

    Returns
    -------
    molecule_ID : int
        Molecule ID number of the given molecule, as given by the first column in HITEMP here: https://hitran.org/hitemp/.
    isotopologue_ID : int
        Isotopologue number of the molecule (1 = most common isotopologue, 2 = second most common, etc.).

    """
    
    molecule_dict = create_id_dict()
    table = download.HITEMP_table()
    
    while True:
        molecule = input("What molecule would you like to download the line list for? \n")
        molecule_id_number = molecule_dict.get(molecule)
        if molecule_id_number in table['ID'].values:
            row = table.loc[table['ID'] == molecule_id_number]   # Find the DataFrame row that contains the molecule_ID
            break   
        else:
            print("\n ----- This molecule ID is not valid. Please check to make sure the molecule you want exists on HITEMP. -----")
    
    while True:
        try:
            isotopologue_ID = int(input("What is the isotopologue ID of the isotopologue you would like to download? Enter '1' for the most abundant isotopologue, '2' for next most abundant, etc. More info found here: https://hitran.org/lbl/# \n"))
            isotope_count = row.loc[row.index.values[0], 'Iso Count'] # Find the number of isotopologues of the given molecule ID
            if isotopologue_ID <= isotope_count:
                break
        except ValueError:
            print("\n ----- Please enter an integer for the isotopologue ID number -----")
        else:
            print("\n ----- This molecule/isotopologue ID combination is not valid. Please check to make sure the combo you want exists on HITEMP. ----- ")
        
    return molecule_id_number, isotopologue_ID


def create_pf(mol_ID, iso_ID, folder, T_min = 70, T_max = 3001, step = 1.0):
    """
    Create partition function file using the partitionSum() function already in hapi

    Parameters
    ----------
    mol_ID : int
        Molecule ID number of the given molecule, as given by the first column in HITEMP here: https://hitran.org/hitemp/.
    iso_ID : int
        Isotopologue number of the molecule (1 = most common isotopologue, 2 = second most common, etc.).
    folder : String
        Local directory where the partition function file is to be stored.
    T_min : int, optional
        Minimum temperature (K) for which the partition function is computed by hapy.py. The default is 70.
    T_max : int, optional
        Maximum temperature (K) for which the partition function is computed by hapy.py. The default is 3001.
    step : float, optional
        Increment for temperature between T_min and T_max. The default is 1.0.

    Returns
    -------
    None.

    """
    
    T, Q = partitionSum(mol_ID, iso_ID, [T_min, T_max], step)

    out_file = folder + moleculeName(mol_ID) + '.pf'
    f_out = open(out_file, 'w')
    f_out.write('T | Q \n') 
    for i in range(len(T)):
        f_out.write('%.1f %.4f \n' %(T[i], Q[i]))
        
    
def create_air_broad(input_dir):
    """
    Create an air broadening file using the downloaded line list.

    Parameters
    ----------
    input_dir : String
        Local directory where the broadening file is to be stored.

    Returns
    -------
    None.

    """
    
    # Instantiate arrays which will be needed for creating air broadening file
    J_lower_all, gamma_air, n_air = (np.array([]) for _ in range(3))
    gamma_air_avg, n_air_avg = (np.array([]) for _ in range(2))
    
    for file in os.listdir(input_dir):
        if file.endswith('.h5'):
            with h5py.File(input_dir + file, 'r') as hdf:
                
                # Populate the arrays by reading in each hdf5 file
                J_lower_all = np.append(J_lower_all, np.array(hdf.get('Lower State J')))
                gamma_air = np.append(gamma_air, np.array(hdf.get('Air Broadened Width')))
                n_air = np.append(n_air, np.array(hdf.get('Temperature Dependence of Air Broadening')))
            
    J_sorted = np.sort(np.unique(J_lower_all))
        
    for i in range(len(J_sorted)):
        
        gamma_air_i = np.mean(gamma_air[np.where(J_lower_all == J_sorted[i])])
        n_air_i = np.mean(n_air[np.where(J_lower_all == J_sorted[i])])
        gamma_air_avg = np.append(gamma_air_avg, gamma_air_i)
        n_air_avg = np.append(n_air_avg, n_air_i)
        
    # Write air broadening file
    out_file = input_dir + 'air.broad'
    f_out = open(out_file,'w')
    
    f_out.write('J | gamma_L_0 | n_L \n')
    
    for i in range(len(J_sorted)):
        f_out.write('%.1f %.4f %.4f \n' %(J_sorted[i], gamma_air_avg[i], n_air_avg[i]))
    
    f_out.close()
   
    
def summon_HITEMP(molecule, isotopologue):
    """
    Main function, uses calls to other functions to perform the download

    Parameters
    ----------
    molecule : int
        Molecule ID number of the given molecule, as given by the first column in HITEMP here: https://hitran.org/hitemp/.
    isotopologue : int
        Isotopologue number of the molecule (1 = most common isotopologue, 2 = second most common, etc.).

    Returns
    -------
    None.

    """
    print("\nFetching data from HITEMP...\nMolecule:", moleculeName(molecule), "\nIsotopologue", isotopologueName(molecule, isotopologue), "\n")
    
    output_folder = download.create_directories(mol_ID = molecule, iso_ID = isotopologue,
                                                database = 'HITEMP')
    create_pf(molecule, isotopologue, output_folder)
    download.download_HITEMP_line_list(molecule, isotopologue, output_folder)
    print("\n\nNow creating air broadening file ...")
    create_air_broad(output_folder)
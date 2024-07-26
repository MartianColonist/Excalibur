''' Functions to repackage calculated cross section grids as HDF5 files '''

import os
import numpy as np
import pandas as pd
import h5py
import re
from tqdm import tqdm
from datetime import date

from Excalibur.misc import read_cross_section_file, int_to_roman


def read_cross_section_grid(species, database, log_P, T, ionization_state,
                            isotope, linelist, broadening,
                            output_dir = './output/'):
    '''
    Read pre-computed cross section .txt files into a single large numpy array.
    '''

    N_PT = len(log_P) * len(T)

    if (database.lower() == 'vald'):
        ion_append = '_' + int_to_roman(ionization_state)
    else:
        ion_append = ''

    # Create progress bar
    with tqdm(total=N_PT, desc = "Progress") as pbar:

        # Loop over pressures and temperatures
        for i in range(len(log_P)):
            for j in range(len(T)):
        
                # Extract cross sections for this P,T
                nu, sigma = read_cross_section_file(species, database, 
                                                    filename = (species + 
                                                                ion_append + 
                                                                '_T' + str(T[j]) + 'K_log_P' + 
                                                                str(log_P[i]) + '_' + 
                                                                broadening + '_sigma.txt'),
                                                    ionization_state = ionization_state,
                                                    isotope = isotope, linelist = linelist, 
                                                    output_dir = output_dir)
                
                # Initialise cross section array
                if ((i == 0) and (j == 0)):
                    log_sigma = np.zeros(shape=(len(log_P),len(T),len(nu)))
                
                # Convert cross section into m^2
                log_sigma_PT = np.log10(np.abs(sigma)*1.0e-4 + 1.0e-250)  # Absolute for (v. rare) numerical -ve values
                
                del sigma
                
                # Reverse direction of array (increasing nu-> increasing wl)
                #sigma_PT = sigma_PT[::-1]
                
                log_sigma[i,j,:] = log_sigma_PT

                del log_sigma_PT

                # Update progress bar
                pbar.update(1)
                 
    return nu, log_sigma


def make_single_species_HDF5(species, database, linelist, log_P, T,
                             isotope = 'default', broadening = 'H2-He', 
                             output_dir = './output/HDF5/', version = 1):
    '''
    Save cross section grid for a single chemical species as a HDF5 file.
    '''

    if (isotope == 'default'):
        print("Now preparing HDF5 file for: " + species)
    else:
        print("Now preparing HDF5 file for: " + species + " (" + isotope + ")")
    
    # Make output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Find ionisation state
    if ('+' in species):
        if (database.lower() == 'vald'):
            ionization_state = 2
            species_no_charge = species[:-1]   # Remove charge from string
        elif (database.lower() == 'exomol'):
            ionization_state = 2
            species_no_charge = re.sub('[+]', '_p', species)
    else:
        ionization_state = 1
        species_no_charge = species

    # Get current date (for file metadata)
    today = date.today()

    if (isotope == 'default'):

        # Create new empty HDF5 file
        HDF5_database = h5py.File(output_dir + species + '.hdf5', 'w')

        # Create new group
        g = HDF5_database.create_group(species)

    else:

        # Create new empty HDF5 file
        HDF5_database = h5py.File(output_dir + isotope + '.hdf5', 'w')    

        # Create new group
        g = HDF5_database.create_group(isotope)

    # Store opacity characteristics
    g.attrs["Source"] = database
    g.attrs["Linelist"] = linelist
    g.attrs["Broadening"] = broadening

    # Denote version and modification time
    g.attrs["Version"] = "v" + str(version)
    g.attrs["Modified"] = str(today)

    # Read in pre-computed cross section files into numpy array
    nu, log_sigma = read_cross_section_grid(species_no_charge, database, 
                                            log_P, T, ionization_state, isotope, 
                                            linelist, broadening)
    
    # Make subgroups to hold data
    g1 = g.create_dataset('nu', data = nu, dtype = 'float64')
    g2 = g.create_dataset('log(P)', data = log_P, dtype = 'float64')
    g3 = g.create_dataset('T', data = T, dtype = 'float64')
    g4 = g.create_dataset('log(sigma)', data = log_sigma, compression = "gzip", 
                          dtype = 'float32', shuffle = True)
    
    # Define meaning of each array
    g1.attrs["Variable"] = "Wavenumber"
    g2.attrs["Variable"] = "Pressure (log10 scale)"
    g3.attrs["Variable"] = "Temperature"
    g4.attrs["Variable"] = "Cross section (log10 scale)"
    
    # Define units of arrays
    g1.attrs["Units"] = "cm^-1"
    g2.attrs["Units"] = "log10(P/bar)"
    g3.attrs["Units"] = "K"
    g4.attrs["Units"] = "log10(sigma/m^2)"
    
    del log_sigma
    
    # Finally, close file to write database to disk
    HDF5_database.close()

    print(species + " done")


def extend_HDF5_database(database_name, species):
    '''
    Add a chemical species from an individual HDF5 file to the combined database.
    '''
    
    # If database doesn't already exist, create file, otherwise load existing database
    if (os.path.isfile('./' + database_name)):
        HDF5_database_all = h5py.File('./' + database_name, 'r+')
    else:
        HDF5_database_all = h5py.File('./' + database_name, 'w')

    # Load HDF5 file for species to add to composite database
    HDF5_database_species = h5py.File('./output/HDF5/' + species + '.hdf5', 'r')

    # Check if this species already exists in database (e.g. an old linelist)
    if ('/' + species in HDF5_database_all):

        # Delete old dataset
        del HDF5_database_all[species]

    # Copy the individual species HDF5 file contents in the composite database
    HDF5_database_species.copy(HDF5_database_species[species], HDF5_database_all)

    # Close both HDF5 files
    HDF5_database_all.close()
    HDF5_database_species.close()

    print(species + " done")


#***** Example script to run EXCALIBUR *****#

import numpy as np

from excalibur.core import summon
from excalibur.core import compute_cross_section
from excalibur.plot import plot_cross_section, cross_section_collection

# Parameters

species = 'CO'
isotope = 'default'
database = 'exomol'

P = 1       # Pressure (bar)
T = 1000       # Temperature (K)


# Download line list
summon(species = species, database = database, 
       isotope = isotope, VALD_data_dir = './VALD Line Lists/')


# Create cross section
nu, sigma = compute_cross_section(input_dir = './input/', database = database, 
                                  species = species, log_pressure = np.log10(P), 
                                  temperature = T, isotope = isotope,
                                  N_cores = 1, nu_out_min = 5000, nu_out_max = 50000, dnu_out = 0.01)

'''
# Download line list
summon(species = 'NaCl', database = database, 
       isotope = isotope, VALD_data_dir = './VALD Line Lists/')

# Create cross section
nu2, sigma2 = compute_cross_section(input_dir = './input/', database = database, 
                                  species = 'NaCl', log_pressure = np.log10(P), 
                                  temperature = T, isotope = isotope,
                                  N_cores = 1, nu_out_min = 200, nu_out_max = 5000, dnu_out = 0.01)
'''

cross_sections = cross_section_collection(nu, sigma, collection = [])


#spectra = spectra_collection(nu2, sigma2, collection = spectra)



# Plot cross section
spectrum = plot_cross_section(collection = cross_sections, labels = ['CO'], filename = 'CO_Li2015_T1000_P1')


'''
spectra = spectra_collection(nu, sigma, collection = [])

# Download line list
summon(species = 'NaCl', database = database, 
       isotope = isotope, VALD_data_dir = './VALD Line Lists/')


# Create cross section
nu, sigma2 = compute_cross_section(input_dir = './input/', database = database, 
                                  species = 'NaCl', log_pressure = np.log10(P), 
                                  temperature = T, isotope = isotope,
                                  N_cores = 1, nu_out_min = 200, nu_out_max = 50000, dnu_out = 0.01)

spectra = spectra_collection(nu, sigma2)
'''

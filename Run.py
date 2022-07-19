#***** Example script to run EXCALIBUR *****#

import numpy as np

from excalibur.core import summon
from excalibur.core import compute_cross_section
from excalibur.plot import plot_sigma_wl

# Parameters

species = 'FeH'
isotope = 'default'
database = 'exomol'

P = 1e-3       # Pressure (bar)
T = 2000       # Temperature (K)

# Download line list
summon(species = species, database = database, 
       isotope = isotope, VALD_data_dir = './VALD Line Lists/')


# Create cross section
nu, sigma = compute_cross_section(input_dir = './input/', database = database, 
                                  species = species, log_pressure = np.log10(P), 
                                  temperature = T, isotope = isotope,
                                  N_cores = 1, nu_out_min = 200, nu_out_max = 50000, dnu_out = 0.01)


# Plot cross section
plot_sigma_wl(nu_arr = nu, sigma_arr = sigma, species = species, temperature = T, 
              log_pressure = np.log10(P), database = database, plot_dir = './plots/')

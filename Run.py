#***** Example script to run EXCALIBUR *****#

from excalibur.core import summon, compute_cross_section

species = 'CO'
database = 'ExoMol'

# Download line list
summon(database=database, species = species)

P = [1, 2]  # Pressure in bars
T = [1000, 2000]  # Temperature in Kelvin
input_directory = './input/' # Top level directory containing line lists

compute_cross_section(species = species, database = database, temperature = T, input_dir = input_directory, 
                      pressure = P)


'''

from excalibur.core import summon
from excalibur.core import compute_cross_section
from excalibur.plot import plot_sigma_wl
# Parameters

species = 'FeH'
isotope = 'default'
database = 'exomol'

P = 1       # Pressure (bar)
T = 2000       # Temperature (K)

# Download line list
summon(species = species, database = database, 
       isotope = isotope, VALD_data_dir = './VALD Line Lists/')


# Create cross section
nu, sigma = compute_cross_section(input_dir = './input/', database = database, 
                                  species = species, pressure = P, 
                                  temperature = T, isotope = isotope,
                                  N_cores = 1, nu_out_min = 200, nu_out_max = 50000, dnu_out = 0.01)


# Plot cross section
plot_sigma_wl(nu_arr = nu, sigma_arr = sigma, species = species, temperature = T, 
              log_pressure = np.log10(P), database = database, plot_dir = './plots/')
'''

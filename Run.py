#***** Example script to run EXCALIBUR *****#

'''Download main isotopologue of NO from ExoMol'''

from excalibur.core import summon, compute_cross_section
from excalibur.plot import read_cross_section_file, cross_section_collection, plot_cross_section

species = 'NO'
#species2 = 'K'
#species3 = 'Ca'
#species4 = 'Ca'

database = 'hitran'
database2 = 'exomol'

'''
summon(database=database, species = species)  # Download line list
summon(database=database2, species = species)  # Download line list
#summon(database=database, species = species2, VALD_data_dir = './VALD Line Lists/') # Download line list
#summon(database=database, species = species3, VALD_data_dir = './VALD Line Lists/')  # Download line list
#summon(database=database, species = species4, VALD_data_dir = './VALD Line Lists/', ionization_state=2)  # Download line list


P = 1  # Pressure in bars
T = 300  # Temperature in Kelvin
input_directory = './input/' # Top level directory containing line lists

# Make cross section for 12C-16O
compute_cross_section(species = species, database = database,
                                  temperature = T, pressure = P, input_dir = input_directory)

# Make cross section for 12C-16O
compute_cross_section(species = species, database = database2,
                                  temperature = T, pressure = P, input_dir = input_directory)
'''
# Make cross section for 13C-16O
#compute_cross_section(species = species, database = database, isotope = '13C-16O', linelist = 'Li2015', temperature = T, pressure = P, input_dir = input_directory)

#compute_cross_section(database=database, species = species2, pressure = P, temperature = T, input_dir = input_directory)

#compute_cross_section(database=database, species = species3, pressure = P, temperature = T, input_dir = input_directory)

#compute_cross_section(database=database, species = species4, pressure = P, temperature = T, input_dir = input_directory, ionization_state=2)

nu, sigma = read_cross_section_file(species = species, database = database,
                                    filename = 'NO_T300K_log_P0.0_air_sigma.txt')

nu2, sigma2 = read_cross_section_file(species = species, database = database2,
                                    filename = 'NO_T300K_log_P0.0_air_sigma.txt')
'''
nu2, sigma2 = read_cross_section_file(species = species2, database = database, 
                                      filename = 'K_I_T1000K_log_P0.0_H2-He_sigma.txt')

nu3, sigma3 = read_cross_section_file(species = species3, database = database,  
                                      filename = 'Ca_I_T1000K_log_P0.0_H2-He_sigma.txt')

nu4, sigma4 = read_cross_section_file(species = species4, database = database, 
                                      ionization_state= 2,
                                      filename = 'Ca_II_T1000K_log_P0.0_H2-He_sigma.txt')
'''
cross_sections = []

# Add first cross section to collection
cross_sections = cross_section_collection(new_x = nu, new_y = sigma, collection = cross_sections)
cross_sections = cross_section_collection(new_x = nu2, new_y = sigma2, collection = cross_sections)
#cross_sections = cross_section_collection(new_x = nu3, new_y = sigma3, collection = cross_sections)
#cross_sections = cross_section_collection(new_x = nu4, new_y = sigma4, collection = cross_sections)

plot_cross_section(collection = cross_sections, labels = ['NO HITRAN' ,'NO ExoMol'], filename = 'NO_HITRAN_ExoMol_test', 
                   linewidth = 0.7, smooth_data = False, x_min=7.200, x_max=7.201)

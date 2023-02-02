#***** Example script to run EXCALIBUR *****#

'''Download main isotopologue of NO from ExoMol'''

from excalibur.core import summon, compute_cross_section
from excalibur.plot import read_cross_section_file, cross_section_collection, plot_cross_section

species = 'CO'
species2 = 'K'
species3 = 'Ca'
species4 = 'Ca'

database = 'exomol'


#summon(database=database, species = species)  # Download line list
#summon(database=database, species = species, isotope = '13C-16O', linelist='Li2015')  # Download line list
#summon(database=database, species = species3, VALD_data_dir = './VALD Line Lists/')  # Download line list
#summon(database=database, species = species4, VALD_data_dir = './VALD Line Lists/', ionization_state=2)  # Download line list


P = 1  # Pressure in bars
T = 1000  # Temperature in Kelvin
input_directory = './input/' # Top level directory containing line lists

#compute_cross_section(database=database, species = species, pressure = P, temperature = T, input_dir = input_directory)

#compute_cross_section(database=database, species = species, pressure = P, temperature = T, input_dir = input_directory, linelist='Li2015', isotope = '13C-16O')

#compute_cross_section(database=database, species = species, pressure = P, temperature = T, input_dir = input_directory)

#compute_cross_section(database=database, species = species4, pressure = P, temperature = T, input_dir = input_directory, ionization_state=2)

nu, sigma = read_cross_section_file(species = species, database = database,
                                    filename = 'CO_T1000K_log_P0.0_H2-He_sigma.txt')

nu2, sigma2 = read_cross_section_file(species = species, database = database, 
                                      isotope = '13C-16O', linelist = 'Li2015', 
                                      filename = 'CO_T1000K_log_P0.0_H2-He_sigma.txt')

cross_sections = []

# Add first cross section to collection
cross_sections = cross_section_collection(new_x = nu, new_y = sigma, collection = cross_sections)

# Add second cross section to collection, making sure to specify the previous collection as a parameter
cross_sections = cross_section_collection(new_x = nu2, new_y = sigma2, collection = cross_sections)

plot_cross_section(collection = cross_sections, labels = ['12C-16O', '13C-16O'], filename = 'Isotopologue_cross_sections_for_joss_figure', smooth_data = True)



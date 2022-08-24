#***** Example script to run EXCALIBUR *****#


'''Download main isotopologue of NO from ExoMol'''

from excalibur.core import summon
from excalibur.core import compute_cross_section

species = 'CO'

database = 'ExoMol'

'''
summon(database=database, species = species)  # Download line list

summon(database=database, species = species, isotope = '13C-16O', linelist = 'Li2015')  # Download line list
'''
P = 1  # Pressure in bars
T = 1000  # Temperature in Kelvin
input_directory = './input/' # Top level directory containing line lists

compute_cross_section(database=database, species = species, isotope = '12C-16O', linelist = 'Li2015', pressure = P,
                      temperature = T, input_dir = input_directory)
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

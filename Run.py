#***** Example script to run EXCALIBUR *****#


'''Download main isotopologue of NO from ExoMol'''

from excalibur.core import summon
from excalibur.core import compute_cross_section

species = 'Na'

database = 'Vald'

summon(database=database, species = species, ionization_state = 2, VALD_data_dir='./VALD Line Lists/')  # Download line list

P = 1  # Pressure in bars
T = 1000  # Temperature in Kelvin
input_directory = './input/' # Top level directory containing line lists

compute_cross_section(database=database, species = species, pressure = P, ionization_state = 2,
                      temperature = T, input_dir = input_directory)
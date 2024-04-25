import os
import numpy as np
import pandas as pd
import re
import time
import requests
import sys
import numba
from bs4 import BeautifulSoup
from scipy.interpolate import UnivariateSpline as Interp
from hapi.hapi import molecularMass, moleculeName, isotopologueName

from .calculate import produce_total_cross_section_VALD_atom

from Excalibur.constants import c, kb, u, P_ref, T_ref

import Excalibur.ExoMol as ExoMol
import Excalibur.HITRAN as HITRAN
import Excalibur.HITEMP as HITEMP
import Excalibur.VALD as VALD
import Excalibur.downloader as download
import Excalibur.broadening as broadening
import Excalibur.Voigt as Voigt
import Excalibur.calculate as calculate

from Excalibur.misc import write_output, check_molecule

def mass(species, isotopologue, linelist):
    """
    Determine the mass of a given chemical species-isotopologue combination in atomic mass units (amu).

    Parameters
    ----------
    species : String
        Molecule/atom we are calculating the mass for.
    isotopologue : String
        Isotopologue of the species we are calculating the mass for.
    linelist : String
        The line list this species' cross-section will later be calculated for. Used to 
        identify between ExoMol, HITRAN/HITEMP, and VALD

    Returns
    -------
    float
        Mass in amu of the given species-isotopologue combination.

    """
    
    # Determine the mass of HITRAN and HITEMP species by using built-in molecularMass() function in hapi.py
    if linelist == 'hitran' or linelist == 'hitemp':
        mol_ID = 1
        while moleculeName(mol_ID) != species:
            mol_ID += 1
            
        iso_ID = 1
        
        while True:
            iso_name = isotopologueName(mol_ID, iso_ID) # Need to format the isotopologue name to match ExoMol formatting
            
            iso_name = HITRAN.replace_iso_name(iso_name)
            
            if iso_name == isotopologue:
                return molecularMass(mol_ID, iso_ID)
            
            else:
                iso_ID += 1
      
    # For VALD atomic line lists, just use the weighted average atomic mass
    elif linelist == 'vald':   
        
        # Atomic masses - Weighted average based on isotopic natural abundances found here: 
        # https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
        mass_dict = {'H': 1.00794072, 'He': 4.00260165, 'Li': 6.94003706, 'Be': 9.012182, 
                     'B': 10.81102777, 'C': 12.0107359, 'N': 14.00674309, 'O': 15.9994053, 
                     'F': 18.998403, 'Ne': 20.1800463, 'Na': 22.989770, 'Mg': 24.30505187, 
                     'Al': 26.981538, 'Si': 28.0853852, 'P': 30.973762, 'S': 32.06608499,
                     'Cl': 35.45653261, 'Ar': 39.94767659, 'K': 39.09830144, 'Ca': 40.07802266, 
                     'Sc': 44.955910, 'Ti': 47.86674971, 'V': 50.941472, 'Cr': 51.99613764,
                     'Mn': 54.938050, 'Fe': 55.84515013, 'Co': 58.933200, 'Ni': 58.69335646, 
                     'Cu': 63.5456439, 'Zn': 65.3955669, 'Ga': 69.72307155, 'Ge': 72.61275896, 
                     'As': 74.921596, 'Se': 78.95938897, 'Br': 79.90352862, 'Kr': 83.79932508, 
                     'Rb': 85.46766375, 'Sr': 87.61664598, 'Y': 88.905848, 'Zr': 91.22364739, 
                     'Nb': 92.906378, 'Mo': 95.93129084, 'Tc': 97.907216, 'Ru': 101.06494511, 
                     'Rh': 102.905504, 'Pd': 106.41532721, 'Ag': 107.8681507, 'Cd': 112.41155267, 
                     'In': 114.81808585, 'Sn': 118.71011064, 'Sb': 121.7597883, 'Te': 127.60312538, 
                     'I': 126.904468, 'Xe': 131.29248065, 'Cs': 132.905447, 'Ba': 137.32688569, 
                     'La': 138.90544868, 'Ce': 140.11572155, 'Pr': 140.907648, 'Nd': 144.23612698, 
                     'Pm': 144.912744, 'Sm': 149.46629229, 'Eu': 151.96436622, 'Gd': 157.25211925, 
                     'Tb': 158.925343, 'Dy': 162.49703004, 'Ho': 164.930319, 'Er': 167.25630107, 
                     'Tm': 168.934211, 'Yb': 173.0376918, 'Lu': 174.96671757, 'Hf': 178.48497094, 
                     'Ta': 180.94787594, 'W': 183.84177868, 'Re': 186.20670567, 'Os': 190.22755215, 
                     'Ir': 192.21605379, 'Pt': 194.73875746, 'Au': 196.966552, 'Hg': 200.59914936, 
                     'Tl': 204.38490867, 'Pb': 207.21689158, 'Bi': 208.980383, 'Po': 208.980383,
                     'At': 209.987131, 'Rn': 222.017570, 'Fr': 223.019731, 'Ra': 226.025403,
                     'Ac': 227.027747, 'Th': 232.038050, 'Pa': 231.035879, 'U': 238.02891307,
                     'Np': 237.048167, 'Pu': 244.064198, 'Am': 243.061373, 'Cm': 247.070347,
                     'Bk': 247.070299, 'Cf': 251.079580, 'Es': 252.082972, 'Fm': 257.095099,
                     'Md': 258.098425, 'No': 259.101024, 'Lr': 262.109692, 'Rf': 263.118313,
                     'Db': 262.011437, 'Sg': 266.012238, 'Bh': 264.012496, 'Hs': 269.001341,
                     'Mt': 268.001388, 
                     }
        
        return mass_dict.get(species)

    # For ExoMol line lists, get the mass contained within the molecule's .def file
    else:
        
        isotopologue = isotopologue.replace('(', '')
        isotopologue = isotopologue.replace(')', '')

        url = 'http://exomol.com/data/molecules/' + species + '/' + isotopologue + '/' + linelist + '/'

        # Parse the webpage to find the .def file and read it
        web_content = requests.get(url).text
        soup = BeautifulSoup(web_content, "lxml")
        def_tag = soup.find('a', href = re.compile("def"))
        new_url = url + def_tag.get('href')
        
        ### If the code above breaks again for some reason, one of these lines may fix the problem:
                #new_url = 'http://exomol.com' + def_tag.get('href') //// this code was prior to the href tag for the .def files changing randomly
                #new_url = 'http://www.exomol.com/db/' + species + '/' + isotopologue + '/' + linelist + '/' + def_tag.get_text()
                #new_url = url + def_tag.get('href')  //// code that should work but I have no idea why it doesn't
        
        
        out_file = './def'
        with requests.get(new_url, stream=True) as request:
            with open(out_file, 'wb') as file:
                for chunk in request.iter_content(chunk_size = 1024 * 1024):
                    file.write(chunk)
                    
        data = pd.read_csv(out_file, delimiter = '#', names = ['Value', 'Key'])  # Store the .def file in a pandas DataFrame
        data = data[data['Key'].str.contains('mass')]  # Only use the row that contains the isotopologue mass
        data = data.reset_index(drop = True)  # Reset the index of the DataFrame
        mass = data['Value'][0]
        mass = re.findall('[0-9|.]+', mass)[0]
        
        os.remove(out_file)
        
        return float(mass)
        

def load_pf(input_directory):
    '''
    Read in a pre-downloaded partition function.

    Parameters
    ----------
    input_directory : String
        Local directory where the .pf file is stored.

    Returns
    -------
    T_pf_raw : numpy array
        RYAN.
    Q_raw : numpy array
        RYAN.

    '''
    
    print("Loading partition function")
    
    # Look for files in input directory ending in '.pf'
    pf_file_name = [filename for filename in os.listdir(input_directory) if filename.endswith('.pf')]
    
    # Read partition function
    pf_file = pd.read_csv(input_directory + pf_file_name[0], sep= ' ', header=None, skiprows=1)
    
    # First column in standard format is temperature, second is the partition function
    T_pf_raw = np.array(pf_file[0]).astype(np.float64)
    Q_raw = np.array(pf_file[1])
    
    return T_pf_raw, Q_raw


def interpolate_pf(T_pf_raw, Q_raw, T, T_ref):
    '''
    Interpolate partition function to the temperature of the cross section computation.
    RYAN

    Parameters
    ----------
    T_pf_raw : TYPE
        RYAN.
    Q_raw : TYPE
        RYAN.
    T : TYPE
        RYAN.
    T_ref : TYPE
        RYAN.

    Returns
    -------
    Q_T : TYPE
        RYAN.
    Q_T_ref : TYPE
        RYAN.

    '''
    
    # Interpolate partition function onto a fine grid using a 5th order spline
    pf_spline = Interp(T_pf_raw, Q_raw, k=5)
    
    # Define a new temperature grid (extrapolated to 10,000K)
    T_pf_fine = np.linspace(1.0, 10000.0, 9999)      
    
    # Using spline, interpolate and extrapolate the partition function to the new T grid
    Q_fine = pf_spline(T_pf_fine)                    
    
    # Find the indices in the fine temperature grid closest to the user specified and reference temperatures
    idx_T = np.argmin(np.abs(T_pf_fine - T))
    idx_T_ref = np.argmin(np.abs(T_pf_fine - T_ref))   
    
    # Find partition function at the user specified and reference temperatures
    Q_T = Q_fine[idx_T]                               
    Q_T_ref = Q_fine[idx_T_ref]                       
    
    return Q_T, Q_T_ref


def create_nu_grid_atom_OLD(atom, T, m, gamma, nu_0, Voigt_sub_spacing, 
                            dnu_out, nu_out_min, nu_out_max, Voigt_cutoff, cut_max):
    '''
    Create the computational (fine) and output (coarse) wavenumber grids for
    an atomic cross section calculation.
    
    Note: for atoms a single grid is used over the entire wavenumber range.

    Parameters
    ----------
    atom : TYPE
        RYAN.
    T : TYPE
        RYAN.
    m : TYPE
        RYAN.
    gamma : TYPE
        RYAN.
    nu_0 : TYPE
        RYAN.
    Voigt_sub_spacing : TYPE
        RYAN.
    dnu_out : TYPE
        RYAN.
    nu_out_min : TYPE
        RYAN.
    nu_out_max : TYPE
        RYAN.
    Voigt_cutoff : TYPE
        RYAN.
    cut_max : TYPE
        RYAN.

    Returns
    -------
    sigma_fine : TYPE
        RYAN.
    sigma_out : TYPE
        RYAN.
    nu_detune : TYPE
        RYAN.
    N_points_fine : TYPE
        RYAN.
    N_Voigt_points : TYPE
        RYAN.
    alpha : TYPE
        RYAN.
    cutoffs : TYPE
        RYAN.
    nu_min : TYPE
        RYAN.
    nu_max : TYPE
        RYAN.
    nu_fine_start : TYPE
        RYAN.
    nu_fine_end : TYPE
        RYAN.
    nu_out : TYPE
        RYAN.
    N_points_out : TYPE
        RYAN.

    '''
    
    # Define the minimum and maximum wavenumber on grid to go slightly beyond user's output limits
    nu_min = 1
    nu_max = nu_out_max + 1000
    
    # First, we need to find values of gamma_V for reference wavenumber (1000 cm^-1)
    alpha_ref = np.sqrt(2.0*kb*T*np.log(2)/m) * (np.array(1000.0)/c)  # Doppler HWHM at reference wavenumber
    gamma_ref = np.min(gamma)                                            # Find minimum value of Lorentzian HWHM
    gamma_V_ref = Voigt.HWHM(gamma_ref, alpha_ref)                       # Reference Voigt width
    
    # Calculate Voigt width for each transition
    alpha = np.sqrt(2.0*kb*T*np.log(2)/m) * (np.array(nu_0)/c)   # Doppler HWHM for each transition
    gamma_V = Voigt.HWHM(gamma, alpha)   # Voigt HWHM
    
    #**** Now compute properties of computational (fine) and output (coarse) wavenumber grid *****
    
    # Wavenumber spacing of of computational grid (smallest of gamma_V_ref/6 or 0.01cm^-1)
    dnu_fine = np.minimum(gamma_V_ref*Voigt_sub_spacing, dnu_out)      
    
    # Number of points on fine grid (rounded)
    N_points_fine = int((nu_max-nu_min)/dnu_fine + 1)
    
    # Adjust dnu_fine slightly to match an exact integer number of grid spaces
    dnu_fine = (nu_max-nu_min)/(N_points_fine - 1)
    
    cutoffs = np.zeros(len(nu_0))   # Line wing cutoffs for each line
    
    # Line cutoffs at min(500 gamma_V, 1000cm^-1)
    for i in range(len(nu_0)):
        
        cutoffs[i] = dnu_fine * (int((Voigt_cutoff*gamma_V[i])/dnu_fine))

        if (cutoffs[i] >= cut_max): cutoffs[i] = cut_max
                
        # Special cases for alkali resonant lines
  #      if ((atom == 'Na') and (int(nu_0[i]) in [16978, 16960])):
  #          cutoffs[i] = 9000.0   # Cutoff @ +/- 9000 cm^-1
  #      elif ((atom == 'K') and  (int(nu_0[i]) in [13046, 12988])): 
  #          cutoffs[i] = 9000.0   # Cutoff @ +/- 9000 cm^-1
    if (atom in ['Na', 'K']):
        cutoffs[i] = 9000.0

      
    # Calculate detuning frequencies for Na and K resonance lines (after Baudino+2015)
    if (atom == 'Na'): 
        nu_detune = 30.0 * np.power((T/500.0), 0.6)
    elif (atom == 'K'): 
        nu_detune = 20.0 * np.power((T/500.0), 0.6)
    else: 
        nu_detune = cut_max
    
    # Evaluate number of frequency points for each Voigt function up to cutoff (one tail)
    N_Voigt_points = ((cutoffs/dnu_fine).astype(np.int64)) + 1  
    
    # Define start and end points of fine grid
    nu_fine_start = nu_min
    nu_fine_end = nu_max
    
    # Initialise output grid
    N_points_out = int((nu_out_max-nu_out_min)/dnu_out + 1)     # Number of points on coarse grid (uniform)
    nu_out = np.linspace(nu_out_min, nu_out_max, N_points_out)  # Create coarse (output) grid
    
    # Initialise cross section arrays on each grid
    sigma_fine = np.zeros(N_points_fine)    # Computational (fine) grid
    sigma_out = np.zeros(N_points_out)      # Coarse (output) grid
    
    return (sigma_fine, sigma_out, nu_detune, N_points_fine, N_Voigt_points, alpha, 
            cutoffs, nu_min, nu_max, nu_fine_start, nu_fine_end, nu_out, N_points_out)

def create_nu_grid(nu_out_min, nu_out_max, dnu_out):
    '''
    RYAN

    Parameters
    ----------
    nu_out_min : TYPE
        RYAN.
    nu_out_max : TYPE
        RYAN.
    dnu_out : TYPE
        RYAN.

    Returns
    -------
    nu_compute : TYPE
        RYAN.

    '''

    # Define the minimum and maximum wavenumber on grid to go slightly beyond user's output limits
    nu_min = min(1, nu_out_min)
    nu_max = nu_out_max + 1000
        
    # Initialise computational grid
    N_compute = int((nu_max - nu_min)/dnu_out + 1)       # Number of points on computational grid (uniform)
    nu_compute = np.linspace(nu_min, nu_max, N_compute)  # Create computational (output) grid
      
    return nu_compute
    

def summon(database = '', species = '', isotope = 'default', VALD_data_dir = '',
           linelist = 'default', ionization_state = 1):
    '''
    Makes calls to other downloader files to retrieve data from the desired database

    Parameters
    ----------
    database : String, optional
        The database from which to download the line list (ExoMol, HITRAN, HITEMP, VALD). The default is ''.
    species : String, optional
        Molecular or atomic species to . The default is ''.
    isotope : String, optional
        Isotopologue of the molecular species. For ExoMol, this is specified on the webpage. For HITRAN/HITEMP use 
        1 for the most naturally abundant isotopologue, 2 for the next most abundant, etc. The default is 'default'.
    VALD_data_dir : String, optional
        Local directory where the VALD line list is stored. The default is ''.
    linelist : String, optional
        Species line list to download. Used mainly for ExoMol, where there are multiple line lists for the same species.
        For HITRAN/HITEMP/VALD, the line list is just the name of the database. The default is 'default'.
    ionization_state : int, optional
        Ionization state to be used, in case of an atomic species. A neutral atom corresponds to an ionization state
        of 1. The default is 1.

    Returns
    -------
    None.

    '''

    
    # If the user has specified a chemical species and line list database, they do not want to be guided via the terminal
    # and we set user_prompt = False
    if database != '' and species != '': 
        user_prompt = False
    else: 
        user_prompt = True
        
    # The user does want to be guided via terminal prompts
    if user_prompt: 
        
        while True:
            database = input('Which line list database do you wish to download from (ExoMol, HITRAN, HITEMP, or VALD)?\n')
            database = database.lower()
            if database == 'exomol' or database == 'hitran' or database == 'hitemp' or database == 'vald' :
                break
            else:
                print("\n ----- This is not a supported database, please try again ----- ")
        
        if database == 'exomol': 
            mol, iso, lin, URL = ExoMol.determine_linelist()
            ExoMol.summon_ExoMol(mol, iso, lin, URL)
            
        if database == 'hitran':
            mol, iso = HITRAN.determine_linelist()
            HITRAN.summon_HITRAN(mol, iso)
            
        if database == 'hitemp':
            mol, iso = HITEMP.determine_linelist()
            HITEMP.summon_HITEMP(mol, iso)
            
        if database == 'vald':
            mol, ion = VALD.determine_linelist(VALD_data_dir)
            VALD.summon_VALD(mol, ion, VALD_data_dir)
            
    # The user has already entered parameters
    if not user_prompt: 
        
        db = database.lower()
        spe = species
        
        if isinstance(isotope, str):
            try:
                isotope = int(isotope)
            except ValueError:
                pass
            
        iso = isotope
        lin = linelist
        ion = ionization_state
        
        if db == 'exomol':
            
            spe = re.sub('[+]', '_p', spe)  # Handle ions
            iso = re.sub('[+]', '_p', iso)  # Handle ions
            
            if isotope == 'default':
                ExoMol.check(spe)
                iso = ExoMol.get_default_iso(spe)
            if linelist == 'default':
                ExoMol.check(spe, iso)
                lin = ExoMol.get_default_linelist(spe, iso)

            ExoMol.check(spe, iso, lin)
            URL = "http://exomol.com/data/molecules/" + spe + '/' + iso + '/' + lin + '/'
            ExoMol.summon_ExoMol(spe, iso, lin, URL)
            
        elif db == 'hitran':
            
            if isotope == 'default':
                iso = 1
            else:
                if isinstance(isotope, str):
                    print("----- Error: for HITRAN / HITEMP isotopes are specified by integers.\n")
                    print("----- 1 = most common, 2 = 2nd most common etc. ----- ")
                    sys.exit(0)
            
            spe = HITRAN.check(spe, iso)
            HITRAN.summon_HITRAN(spe, iso)
            
        elif db == 'hitemp':
            
            if isotope == 'default':
                iso = 1
            else:
                if isinstance(isotope, str):
                    print("----- Error: for HITRAN / HITEMP isotopes are specified by integers.\n")
                    print("----- 1 = most common, 2 = 2nd most common etc. ----- ")
                    sys.exit(0)

            spe = HITEMP.check(spe, iso)
            HITEMP.summon_HITEMP(spe, iso)
            
        elif db == 'vald':
            
            VALD.check(spe, ion, VALD_data_dir)
            VALD.summon_VALD(spe, ion, VALD_data_dir)
        
        else:
            print("\n ----- You have not passed in a valid database. Please try calling the summon() function again. ----- ")
            sys.exit(0)
        
    print("\nLine list ready.\n")
    
    
def compute_cross_section(input_dir, database, species, temperature, pressure = None, 
                          log_pressure = None, isotope = 'default', ionization_state = 1, 
                          linelist = 'default', cluster_run = False, set_mass = None,
                          nu_out_min = 200, nu_out_max = 25000, dnu_out = 0.01, 
                          broad_type = 'default', broadening_file = '',
                          gamma_0_fixed = 0.07, n_L_fixed = 0.50, X_H2 = 0.85, 
                          X_He = 0.15, Voigt_cutoff = 500, Voigt_sub_spacing = (1.0/6.0), 
                          N_alpha_samples = 500, S_cut = 1.0e-100, cut_max = 30.0, 
                          N_cores = 1, verbose = True):
    '''
    
    Main function to compute cross section, given that the requisite line list has already been downloaded.

    Parameters
    ----------
    input_dir : String
        Local directory where the line list is stored.
    database : String
        Database from which the line list was downloaded.
    species : String
        Molecule or atom for which the cross section is to be computed.
    temperature : TYPE
        DESCRIPTION.
    pressure : TYPE
        DESCRIPTION.
    log_pressure : TYPE
        DESCRIPTION.
    isotope : String, optional
        Isotopologue of the species (if species is a molecule). For ExoMol, this is specified on the webpage. For HITRAN/HITEMP use 
        1 for the most naturally abundant isotopologue, 2 for the next most abundant, etc. The default is 'default'.
    ionization_state : int, optional
        Ionization state to be used, in case of an atomic species. A neutral atom corresponds to an ionization state
        of 1. The default is 1.
    linelist : TYPE, optional
        Species line list to download. For ExoMol there are often multiple line lists for the same species.
        For HITRAN/HITEMP/VALD, the line list is just the name of the database. The default is 'default'.
    cluster_run : bool, optional
        DESCRIPTION. The default is False.
    set_mass : double, optional
        Mass of the species in atomic mass units (amu). The default is None.
    nu_out_min : double, optional
        Minimum wavenumber for which the cross section is computed. The default is 200.
    nu_out_max : double, optional
        Maximum wavenumber for which the cross section is computed.. The default is 25000.
    dnu_out : double, optional
        RYAN. The default is 0.01.
    broad_type : String, optional
        Type of pressure broadening to be used in the computation. The default is 'default'.
    broadening_file : String, optional
        Broadening file to be used for the calculation. Usually used for custom broadening files. The default is ''.
    gamma_0_fixed : double, optional
        RYAN. The default is 0.07. 
    n_L_fixed : double, optional
        RYAN. The default is 0.50.
    X_H2 : double, optional
        RYAN. The default is 0.85.
    X_He : double, optional
        RYAN. The default is 0.15.
    Voigt_cutoff : double, optional
        RYAN. The default is 500.
    Voigt_sub_spacing : double, optional
        RYAN. The default is (1.0/6.0).
    N_alpha_samples : RYAN, optional
        RYAN. The default is 500.
    S_cut : double, optional
        RYAN. The default is 1.0e-100.
    cut_max : double, optional
        RYAN. The default is 30.0.
    N_cores : int, optional
        Number of parallel cores to be used in the computation. The default is 1.
    verbose : bool, optional
        Controls whether or not intermittent print statements are displayed relaying the status of the cross section computation.
        The default is True.

    Returns
    -------
    nu_out : double
        Wavenumber .
    sigma_out : double
        DESCRIPTION.

    '''
    
    print("Beginning cross-section computations...")
    
    # Start clock for timing program
    t_start = time.perf_counter()
    
    # Configure numba to parallelise with the user specified number of cores
    numba.set_num_threads(N_cores)
    
    if pressure is None and log_pressure is None:
        print("\nYou need to specify a pressure to run the cross section computation. Please try again.")
        sys.exit(0)
    
    # Cast pressure to np arrays
    if log_pressure != None:
        if not isinstance(log_pressure, list):  
            log_pressure = [log_pressure]
        log_pressure = np.array(log_pressure)
        pressure = np.power(10, log_pressure)
    else:
        if not isinstance(pressure, np.ndarray) and not isinstance(pressure, list):  
            pressure = [pressure]
        pressure = np.array(pressure)
        log_pressure = np.log10(pressure)
    
    if not isinstance(temperature, np.ndarray) and not isinstance(temperature, list):  
        temperature = np.array([temperature])
        
    # Cast all temperatures and pressures to floats
    for i in range(len(pressure) - 1):
        pressure[i] = float(pressure[i])
    
    for i in range(len(temperature) - 1):
        temperature[i] = float(temperature[i])
    
    database = database.lower()
    
    # Locate the input_directory where the line list is stored
    input_directory = download.find_input_dir(input_dir, database, species, isotope, ionization_state, linelist)
    
    # Use the input directory to define parameters right at the start
    linelist, isotopologue = download.parse_directory(input_directory, database)
    
    # HITRAN, HITEMP, and VALD do not have seperate line list names
    if database != 'exomol':
        linelist = database
    
    # Load full set of downloaded line list files
    linelist_files = [filename for filename in os.listdir(input_directory) if filename.endswith('.h5')]
    
    if database == 'exomol':
        print("Loading ExoMol format")
        E, g, J = ExoMol.load_states(input_directory)  # Load from .states file
    
    elif database == 'hitran':
        print("Loading HITRAN format")
        # Nothing else required at this stage
    
    elif database == 'hitemp':
        print("Loading HITEMP format")
        # Nothing else required at this stage
        
    elif database == 'vald':
        print("Loading VALD format")
        nu_0, gf, E_low, E_up, J_low, l_low, l_up, \
        Gamma_nat, Gamma_vdw, alkali = VALD.load_line_list(input_directory, species)
        
    # Load partition function
    T_pf_raw, Q_raw = load_pf(input_directory)
    
    # Find mass of the species
    if (set_mass is None):
        m = mass(species, isotopologue, linelist) * u
    else:
        m = set_mass * u
    
    # Check if we have a molecule or an atom
    is_molecule = check_molecule(species)
    
    # Store ionization state as roman numeral for later (atoms only)
    roman_num = ''
    if is_molecule == False:
        for i in range(ionization_state):
            roman_num += 'I'
        
    # If user didn't specify a type of pressure broadening, determine based on available broadening data
    if is_molecule and broad_type == 'default':
        
        broad_type = broadening.det_broad(input_directory)
        
        if broad_type == 'H2-He':
            J_max, gamma_0_H2, n_L_H2, gamma_0_He, n_L_He = broadening.read_H2_He(input_directory)
            
        elif broad_type == 'air':
            J_max, gamma_0_air, n_L_air = broadening.read_air(input_directory)
            
        elif broad_type == 'SB07':
            J_max, gamma_0_SB07 = broadening.read_SB07(input_directory)
            
    # If user specifed a pressure broadening prescription, proceed to load the relevant broadening file
    elif is_molecule and broad_type != 'default':
        
        if (broad_type == 'H2-He' and 'H2.broad' in os.listdir(input_directory) 
                                  and 'He.broad' in os.listdir(input_directory)):
            J_max, gamma_0_H2, n_L_H2, gamma_0_He, n_L_He = broadening.read_H2_He(input_directory)
        
        elif broad_type == 'air' and 'air.broad' in os.listdir(input_directory):
            J_max, gamma_0_air, n_L_air = broadening.read_air(input_directory)
            
        elif broad_type == 'SB07':

            # Create Sharp & Burrows broadening file if it doesn't already exist
            if os.path.isfile(input_directory + 'SB07.broad') == False:
                broadening.create_SB07(input_directory)
            
            J_max, gamma_0_SB07 = broadening.read_SB07(input_directory)
            
        elif broad_type == 'custom' and broadening_file in os.listdir(input_directory):
            J_max, gamma_0_air, n_L_air = broadening.read_custom(input_directory, broadening_file)
        
        elif broad_type == 'fixed':
            J_max = 0
        
        else:
            print("\nYou did not enter a valid type of pressure broadening. Please try again.")
            sys.exit(0)
            
    # For atoms, only H2-He pressure broadening is currently supported
    elif is_molecule == False:
        
        if broad_type != 'default' and broad_type != 'H2-He':
            print("You did not specify a valid choice of pressure broadening.\n" 
                  "For atoms the only supported option is 'H2-He', so we will continue by using that." )
        
        broad_type = 'H2-He'
        
        gamma_0_H2, gamma_0_He, \
        n_L_H2, n_L_He = broadening.read_atom(species, nu_0, gf, E_low, E_up, 
                                              J_low, l_low, l_up, Gamma_nat, 
                                              Gamma_vdw, alkali, m)
    
    #***** Load pressure and temperature for this calculation *****#
    P_arr = pressure     # Pressure array (bar)
    T_arr = temperature  # Temperature array (K)
    
    # If conducting a batch run on a cluster
    if (cluster_run == True):
            
        try:
            idx_PT = int(sys.argv[1])
            
        except IndexError:
            print("\n----- You need to enter a command line argument if cluster_run is set to True. ----- ")
            sys.exit(0)
            
        except ValueError:
            print("\n----- The command line argument needs to be an int. -----")
            sys.exit(0)
            
        if idx_PT >= len(P_arr) * len(T_arr):
            print("\n----- You have provided a command line argument that is out of range for the specified pressure and temperature arrays. -----")
            sys.exit(0)
        
        P = P_arr[idx_PT//len(T_arr)]   # Atmospheric pressure (bar)
        T = T_arr[idx_PT%len(T_arr)]    # Atmospheric temperature (K)
        
        # For a cluster run, each core separately handles a single (P,T) combination
        N_P = 1
        N_T = 1
        
    # If running on a single machine, compute a cross section for each (P,T) pair sequentially
    else:
        
        N_P = len(P_arr)
        N_T = len(T_arr)
    
    # Bool which is used to see if more than one (P,T) pair is 
    grid = True if N_P * N_T > 1 else False

    # Compute cross section for each pressure and temperature point
    for p in range(N_P):
        for t in range(N_T):
            
            # When not running on a cluster, select the next (P,T) pair
            if (cluster_run == False):
                
                P = P_arr[p]   # Atmospheric pressure (bar)
                T = T_arr[t]   # Atmospheric temperature (K)
            
            # Interpolate the tabulated partition function to the desired temperature and reference temperature
            Q_T, Q_T_ref = interpolate_pf(T_pf_raw, Q_raw, T, T_ref)
            
            # Handle pressure broadening, wavenumber grid creation and Voigt profile pre-computation for molecules
            if is_molecule:
                
                # Compute Lorentzian HWHM as a function of J_low
                if broad_type == 'H2-He':
                    gamma = broadening.compute_H2_He(gamma_0_H2, T_ref, T, 
                                                     n_L_H2, P, P_ref, X_H2, 
                                                     gamma_0_He, n_L_He, X_He)
                elif broad_type == 'air':
                    gamma = broadening.compute_air(gamma_0_air, T_ref, T, 
                                                   n_L_air, P, P_ref)
                elif broad_type == 'SB07':
                    gamma = broadening.compute_SB07(gamma_0_SB07, P, P_ref)
                elif broad_type == 'custom': 
                    gamma = broadening.compute_air(gamma_0_air, T_ref, T,    # Computation step is the same as for air broadening
                                                   n_L_air, P, P_ref)
                elif broad_type == 'fixed':
                    gamma = np.array([(gamma_0_fixed * np.power((T_ref/T), n_L_fixed) * (P/P_ref))])  # Fixed Lorentizian HWHM (1 element array)
                    
                # Create wavenumber grid for cross section compuation
                nu_compute = create_nu_grid(nu_out_min, nu_out_max, dnu_out)
                                                                                                                            
                # Initialise cross section arrays for computations
                sigma_compute = np.zeros(len(nu_compute))    # Computational grid
                
                #***** Pre-compute Voigt function array for molecules *****#
    
                print('Pre-computing Voigt profiles...')
    
                t1 = time.perf_counter()    
                
                # Pre-compute template Voigt profiles
                (nu_sampled, alpha_sampled, 
                 cutoffs, N_Voigt, Voigt_arr, 
                 dV_da_arr, dV_dnu_arr, 
                 dnu_Voigt) = Voigt.precompute_molecules(nu_compute, dnu_out, m, T, 
                                                         Voigt_sub_spacing, Voigt_cutoff, 
                                                         N_alpha_samples, gamma, cut_max)
                
                t2 = time.perf_counter()
                time_precompute = t2-t1
            
                print('Voigt profiles computed in ' + str(time_precompute) + ' s')  
                
            # Handle pressure broadening and wavenumber grid creation for atoms
            elif is_molecule == False:
                
                # Compute Lorentzian HWHM line-by-line for atoms
                gamma = broadening.compute_H2_He(gamma_0_H2, T_ref, T, n_L_H2, 
                                                 P, P_ref, X_H2, gamma_0_He, 
                                                 n_L_He, X_He)
                
                # Add natural broadening for each line
                gamma += ((1.0/(4.0*np.pi*(100.0*c))) * Gamma_nat)  
            
                # Create wavenumber grid properties for cross section calculation      
                nu_compute = create_nu_grid(nu_out_min, nu_out_max, dnu_out)
                
                # Initialise cross section arrays for computations
                sigma_compute = np.zeros(len(nu_compute))    # Computational grid
                
                #***** Pre-compute Voigt function array for molecules *****#
    
                print('Pre-computing Voigt profiles...')
    
                t1 = time.perf_counter()    
                
                # Pre-compute Voigt profiles for each line on computational grid
                cutoffs, N_Voigt, Voigt_arr = Voigt.precompute_atoms(species, nu_compute, m, T, gamma, 
                                                                     nu_0, Voigt_cutoff, cut_max)
                
                t2 = time.perf_counter()
                time_precompute = t2-t1
            
                print('Voigt profiles computed in ' + str(time_precompute) + ' s')  
                

                                                                                                                                        
            print("Pre-computation steps complete")
            
            if is_molecule:
                if grid:
                    index = p * N_T + (t + 1)
                    
                    print('Generating cross section for ' + species + ' at P = ' + str(P) + ' bar, T = ' + str(T) + 
                          ' K' + '   [' + str(index) + ' of ' + str(N_P * N_T) + ']')
                else:
                    print('Generating cross section for ' + species + ' at P = ' + str(P) + ' bar, T = ' + str(T) + ' K')
            else:
                print('Generating cross section for ' + species + ' ' + roman_num + ' at P = ' + str(P) + ' bar, T = ' + str(T) + ' K')

            # Call relevant cross section computation function for given line list
            if database == 'exomol':    
                calculate.cross_section_EXOMOL(linelist_files, input_directory, 
                                               nu_compute, sigma_compute, alpha_sampled, 
                                               m, T, Q_T, g, E, J, J_max, N_Voigt, cutoffs,
                                               Voigt_arr, dV_da_arr, dV_dnu_arr, dnu_Voigt, S_cut,
                                               verbose)
                
            elif database in ['hitran', 'hitemp']:
                calculate.cross_section_HITRAN(linelist_files, input_directory, 
                                               nu_compute, sigma_compute, alpha_sampled, 
                                               m, T, Q_T, Q_T_ref, J_max, N_Voigt, cutoffs,
                                               Voigt_arr, dV_da_arr, dV_dnu_arr, dnu_Voigt, S_cut,
                                               verbose)
                
            elif database == 'vald':
                produce_total_cross_section_VALD_atom(nu_compute, sigma_compute, nu_0, 
                                                      E_low, gf, m, T, Q_T, N_Voigt, 
                                                      cutoffs, Voigt_arr, S_cut)
                                    
            # Clip ends from computational grid to leave output wavenumber and cross section grids            
            nu_out = nu_compute[(nu_compute >= nu_out_min) & (nu_compute <= nu_out_max)]
            sigma_out = sigma_compute[(nu_compute >= nu_out_min) & (nu_compute <= nu_out_max)]
        
            # Create output directory (if not already present)
            output_directory = re.sub('/input/', '/output/', input_directory)
    
            if not os.path.exists(output_directory):
                os.makedirs(output_directory)
    
            # Write cross section to file
            write_output(output_directory, species, roman_num, 
                         T, np.log10(P), broad_type, broadening_file, nu_out, sigma_out)
    
    # Print final runtime
    t_final = time.perf_counter()
    total_final = t_final-t_start
    
    print('Total runtime: ' + str(total_final) + ' s')

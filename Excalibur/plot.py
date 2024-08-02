import os
import numpy as np
import matplotlib
import matplotlib.style
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, \
                              ScalarFormatter, NullFormatter, LogLocator
from scipy.ndimage import gaussian_filter1d

from Cthulhu.misc import round_sig_figs

import contextlib
with contextlib.redirect_stdout(None): #suppress HITRAN automatic print statement
    from hapi.hapi import isotopologueName
import Cthulhu.HITRAN as HITRAN
import Cthulhu.ExoMol as ExoMol
import sys
import csv
import re

plt.style.use('classic')
plt.rc('font', family = 'serif')
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['figure.facecolor'] = 'white'


def plot_cross_section(collection, labels, filename, plot_dir = './plots/',
                       x_min = None, x_max = None, y_min = None, y_max = None,
                       color_list = [], smooth_data = False, std = 1000,
                       x_unit = 'micron', x_axis_scale = 'log', 
                       save_fig = True, ax = None, **kwargs):
    """
    Generate a plot of cross section file[s], in both wavelength and wavenumber

    Parameters
    ----------
    molecule : TYPE
        DESCRIPTION.
    temperature : TYPE
        DESCRIPTION.
    log_pressure : TYPE
        DESCRIPTION.
    nu_arr : TYPE, optional
        DESCRIPTION. The default is [].
    sigma_arr : TYPE, optional
        DESCRIPTION. The default is [].
    file : TYPE, optional
        DESCRIPTION. The default is ''.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """


    if not os.path.exists('./plots/'):
        os.mkdir('./plots')


    # Add small value to each computed cross-section to avoid log(0) on log plots
    for xs in collection:
        xs[1] = xs[1] + 1.0e-250

    nu_min, nu_max, sigma_min, sigma_max = find_min_max_nu_sigma(collection)

    # Set xlims so that user is able to view a windowed region of the cross section
    if x_min != None:
        if (x_unit == 'micron'):
            nu_max = 1.0e4/x_min
        else:
            nu_min = x_min
    else:
        x_min = 1e4/nu_max

    if x_max != None:
        if (x_unit == 'micron'):
            nu_min = 1.0e4/x_max
        else:
            nu_max = x_max
    else:
        x_max = 1e4/nu_min

    #***** Format x ticks *****#

    x_range = x_max - x_min

    # Create x formatting objects
    if (x_range >= 0.2):                      
        if (x_max < 1.0):                         # Plotting the optical range
            xmajorLocator = MultipleLocator(0.1)
            xminorLocator = MultipleLocator(0.02)
        else:                                      # Plot extends into the infrared
            xmajorLocator = MultipleLocator(1.0)
            xminorLocator = MultipleLocator(0.1)
    elif ((x_range < 0.2) and (x_range >= 0.02)):   # High-resolution zoomed plots
        xmajorLocator = MultipleLocator(0.01)
        xminorLocator = MultipleLocator(0.002)
    else:                                             # Super high-resolution
        xmajorLocator = MultipleLocator(0.001)
        xminorLocator = MultipleLocator(0.0002)

    xmajorFormatter = FormatStrFormatter('%g')
    xminorFormatter = NullFormatter()

    fig = plt.figure()

    if (ax == None):
        ax = plt.gca()
    else:
        ax = ax

    # Set x axis to be linear or logarithmic
    ax.set_yscale('log')
    ax.set_xscale(x_axis_scale)

    # Assign formatter objects to axes
    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_major_formatter(xmajorFormatter)
    ax.xaxis.set_minor_locator(xminorLocator)
    ax.xaxis.set_minor_formatter(xminorFormatter)

    # Define colors for plotted spectra (default or user choice)
    if color_list == []:   # If user did not specify a custom colour list
        colors = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']
    else:
        colors = color_list

    # plot each spectra in the array
    for i in range(len(collection)):
        (nu, spec) = collection[i]
        wl = 1.0e4/nu # Convert wavenumber (in cm^-1) to wavelength (in um)

        if smooth_data == True:
            spec = gaussian_filter1d(spec, std)

        ax.plot(wl, spec, lw=0.5, alpha = 0.5, color= colors[i], label = labels[i])

    # Set y range
    if (y_min == None):
        y_min = 1.0e-50
    if (y_max == None):
        y_max = sigma_max

    ax.set_xlim([1e4/nu_max, 1e4/nu_min]) # Convert wavenumber (in cm^-1) to wavelength (in um)
    ax.set_ylim([y_min, y_max * 10])

    ax.set_xlabel(r'Wavelength (μm)', fontsize = 14)
    ax.set_ylabel(r'Cross Section (cm$^2$)', fontsize = 14)

    # Decide at which wavelengths to place major tick labels
    wl_ticks = set_wl_ticks((1e4/nu_max), (1e4/nu_min), x_axis_scale)
    ax.set_xticks(wl_ticks)

    legend = ax.legend(loc='upper right', shadow=False, frameon=False, prop={'size':10})

    for legline in legend.legendHandles:
        legline.set_linewidth(1.0)

    plt.tight_layout()

    if (save_fig == True):
        plt.savefig('./plots/' + filename, dpi = 300)

    # print("\nPlotting complete.")


def find_min_max_nu_sigma(spectra):
    '''
    Find min and max wavenumber for which the cross sections were computed, and
    determine the min and max values that the cross section takes on in that range.

    Parameters
    ----------
    spectra : TYPE
        DESCRIPTION.

    Returns
    -------
    wl_min : TYPE
        DESCRIPTION.
    wl_max : TYPE
        DESCRIPTION.
    sigma_min : TYPE
        DESCRIPTION.
    sigma_max : TYPE
        DESCRIPTION.

    '''

    nu_min = 1e10   # Dummy value

    # Loop over each model, finding the minimum wavenumber value
    for i in range(len(spectra)):

        nu_min_i = np.min(spectra[i][0])
        nu_min = min(nu_min, nu_min_i)


    nu_max = 1e-10  # Dummy value

    # Loop over each model, finding the maximum wavenumber value
    for i in range(len(spectra)):

        nu_max_i = np.max(spectra[i][0])
        nu_max = max(nu_max, nu_max_i)


    sigma_min = 1e10   # Dummy value

    # Loop over each model, finding the minimum cross section value
    for i in range(len(spectra)):

        sigma_min_i = np.min(spectra[i][1])
        sigma_min = min(sigma_min, sigma_min_i)


    sigma_max = 1e-200  # Dummy value

    # Loop over each model, finding the maximum cross section value
    for i in range(len(spectra)):

        sigma_max_i = np.max(spectra[i][1])
        sigma_max = max(sigma_max, sigma_max_i)

    return nu_min, nu_max, sigma_min, sigma_max


def set_wl_ticks(wl_min, wl_max, wl_axis = 'log'):
    '''
    Calculates default x axis tick spacing for cross section plots.
    
    Args:
        wl_min (float):
            The minimum wavelength to plot.
        wl_max (float):
            The maximum wavelength to plot.
        wl_axis (str, optional):
            The type of x-axis to use ('log' or 'linear').
            
    Returns:
        np.array:
            The x axis tick values for the given wavelength range.

    '''

    wl_range = wl_max - wl_min

    if (wl_max < wl_min):
        raise Exception("Error: max wavelength must be greater than min wavelength.")

    # For plots over a wide wavelength range
    if (wl_range > 0.2):
        if (wl_max <= 1.0):
            wl_ticks_1 = np.arange(round_sig_figs(wl_min, 1), round_sig_figs(wl_max, 2)+0.01, 0.1)
            wl_ticks_2 = np.array([])
            wl_ticks_3 = np.array([])
            wl_ticks_4 = np.array([])
        elif (wl_max <= 2.0):
            if (wl_min < 1.0):
                wl_ticks_1 = np.arange(round_sig_figs(wl_min, 1), 1.0, 0.2)
                wl_ticks_2 = np.arange(1.0, round_sig_figs(wl_max, 2)+0.01, 0.2)
            else:
                wl_ticks_1 = np.array([])
                wl_ticks_2 = np.concatenate(([round_sig_figs(wl_min, 2)], np.arange(np.ceil(wl_min*5)/5.0, round_sig_figs(wl_max, 2)+0.01, 0.2)))
            wl_ticks_3 = np.array([])
            wl_ticks_4 = np.array([])
        elif (wl_max <= 3.0):
            if (wl_min < 1.0):
                wl_ticks_1 = np.arange(round_sig_figs(wl_min, 1), 1.0, 0.2)
                wl_ticks_2 = np.arange(1.0, round_sig_figs(wl_max, 2)+0.01, 0.5)
            else:
                wl_ticks_1 = np.array([])
                wl_ticks_2 = np.concatenate(([round_sig_figs(wl_min, 2)], np.arange(np.ceil(wl_min*5)/5.0, round_sig_figs(wl_max, 2)+0.01, 0.2)))
            wl_ticks_3 = np.array([])
            wl_ticks_4 = np.array([])
        elif (wl_max <= 5.0):
            if (wl_min < 1.0):
                wl_ticks_1 = np.arange(round_sig_figs(wl_min, 1), 1.0, 0.2)
                wl_ticks_2 = np.arange(1.0, 3.0, 0.5)
                wl_ticks_3 = np.arange(3.0, round_sig_figs(wl_max, 2)+0.01, 1.0)
            elif (wl_min < 3.0):
                wl_ticks_1 = np.array([])
                wl_ticks_2 = np.concatenate(([round_sig_figs(wl_min, 2)], np.arange(np.ceil(wl_min*5)/5.0, 3.0, 0.2)))
                if (wl_axis == 'log'):
                    wl_ticks_3 = np.arange(3.0, round_sig_figs(wl_max, 2)+0.01, 0.5)
                elif (wl_axis == 'linear'):
                    wl_ticks_3 = np.arange(3.0, round_sig_figs(wl_max, 2)+0.01, 0.2)
            else:
                wl_ticks_1 = np.array([])
                wl_ticks_2 = np.array([])
                wl_ticks_3 = np.concatenate(([round_sig_figs(wl_min, 2)], np.arange(np.ceil(wl_min*5)/5.0, round_sig_figs(wl_max, 2)+0.01, 0.2)))
            wl_ticks_4 = np.array([])
        elif (wl_max <= 10.0):
            if (wl_min < 1.0):
                if (wl_axis == 'log'):
                    wl_ticks_1 = np.arange(round_sig_figs(wl_min, 1), 1.0, 0.2)
                    wl_ticks_2 = np.arange(1.0, 3.0, 0.5)
                    wl_ticks_3 = np.arange(3.0, round_sig_figs(wl_max, 2)+0.01, 1.0)
                elif (wl_axis == 'linear'):
                    wl_ticks_1 = np.arange(round_sig_figs(wl_min, 1), 1.0, 1.0)
                    wl_ticks_2 = np.arange(1.0, 3.0, 1.0)            
                    wl_ticks_3 = np.arange(3.0, round_sig_figs(wl_max, 2)+0.01, 1.0)
            elif (wl_min < 3.0):
                wl_ticks_1 = np.array([])
                wl_ticks_2 = np.concatenate(([round_sig_figs(wl_min, 2)], np.arange(np.ceil(wl_min*2)/2.0, 3.0, 0.5)))
                if (wl_axis == 'log'):
                    wl_ticks_3 = np.arange(3.0, round_sig_figs(wl_max, 2)+0.01, 1.0)
                elif (wl_axis == 'linear'):
                    wl_ticks_3 = np.arange(3.0, round_sig_figs(wl_max, 2)+0.01, 0.5)
            else:
                wl_ticks_1 = np.array([])
                wl_ticks_2 = np.array([])
                if (wl_axis == 'log'):
                    wl_ticks_3 = np.concatenate(([round_sig_figs(wl_min, 2)], np.arange(np.ceil(wl_min*2)/2.0, round_sig_figs(wl_max, 2)+0.01, 0.5)))
                elif (wl_axis == 'linear'):
                    wl_ticks_3 = np.concatenate(([round_sig_figs(wl_min, 2)], np.arange(np.ceil(wl_min*2)/2.0, round_sig_figs(wl_max, 2)+0.01, 0.5)))
            wl_ticks_4 = np.array([])
        else:
            if (wl_min < 1.0):
                if (wl_axis == 'log'):
                    wl_ticks_1 = np.arange(round_sig_figs(wl_min, 1), 1.0, 0.2)
                    wl_ticks_2 = np.arange(1.0, 4.0, 1.0)
                    wl_ticks_3 = np.arange(4.0, 10.0, 2.0)
                    wl_ticks_4 = np.arange(10.0, round_sig_figs(wl_max, 2)+0.01, 10.0)
                elif (wl_axis == 'linear'):
                    wl_ticks_1 = np.arange(round_sig_figs(wl_min, 1), 1.0, 1.0)
                    wl_ticks_2 = np.arange(1.0, 4.0, 1.0)            
                    wl_ticks_3 = np.arange(4.0, 10.0, 2.0)
                    wl_ticks_4 = np.arange(10.0, round_sig_figs(wl_max, 3)+0.01, 10.0)
            elif (wl_min < 3.0):
                wl_ticks_1 = np.array([])
                if (wl_axis == 'log'): 
                    wl_ticks_2 = np.concatenate(([round_sig_figs(wl_min, 2)], np.arange(np.ceil(wl_min*2)/2.0, 4.0, 0.5)))
                    wl_ticks_3 = np.arange(4.0, 10.0, 2.0)
                    wl_ticks_4 = np.arange(10.0, round_sig_figs(wl_max, 2)+0.01, 10.0)
                elif (wl_axis == 'linear'):
                    wl_ticks_2 = np.concatenate(([round_sig_figs(wl_min, 2)], np.arange(np.ceil(wl_min), 4.0, 1.0)))
                    wl_ticks_3 = np.arange(4.0, 10.0, 2.0)
                    wl_ticks_4 = np.arange(10.0, round_sig_figs(wl_max, 3)+0.01, 10.0)
            elif (wl_min < 10.0):
                wl_ticks_1 = np.array([])
                wl_ticks_2 = np.array([])
                if (wl_axis == 'log'):
                    wl_ticks_3 = np.concatenate(([round_sig_figs(wl_min, 2)], np.arange(np.ceil(wl_min), 10.0, 1.0)))
                    wl_ticks_4 = np.arange(10.0, round_sig_figs(wl_max, 3)+0.01, 10.0)
                elif (wl_axis == 'linear'):
                    wl_ticks_3 = np.concatenate(([round_sig_figs(wl_min, 2)], np.arange(np.ceil(wl_min), 10.0, 1.0)))
                    wl_ticks_4 = np.arange(10.0, round_sig_figs(wl_max, 3)+0.01, 10.0)
            else:
                wl_ticks_1 = np.array([])
                wl_ticks_2 = np.array([])
                wl_ticks_3 = np.array([])
                if (wl_axis == 'log'):
                    wl_ticks_4 = np.concatenate(([round_sig_figs(wl_min, 3)], np.arange(np.ceil(wl_min), round_sig_figs(wl_max, 3)+0.01, 10.0)))
                elif (wl_axis == 'linear'):
                    wl_ticks_4 = np.concatenate(([round_sig_figs(wl_min, 3)], np.arange(np.ceil(wl_min*2)/2.0, round_sig_figs(wl_max, 3)+0.01, 0.5)))

        wl_ticks = np.concatenate((wl_ticks_1, wl_ticks_2, wl_ticks_3, wl_ticks_4))

    # For high-resolution (zoomed in) spectra
    else:

        # Aim for 10 x-axis labels
        wl_spacing = round_sig_figs((wl_max - wl_min), 1)/10
        
        major_exponent = round_sig_figs(np.floor(np.log10(np.abs(wl_spacing))), 1)
        
        # If last digit of x labels would be 3,6,7,8,or 9, bump up to 10
        if (wl_spacing > 5*np.power(10, major_exponent)):
            wl_spacing = 1*np.power(10, major_exponent+1)
        elif (wl_spacing == 3*np.power(10, major_exponent)):
            wl_spacing = 2*np.power(10, major_exponent)

        wl_ticks = np.arange(round_sig_figs(wl_min, 3), round_sig_figs(wl_max, 3)+0.0001, wl_spacing)

    return wl_ticks


def compare_cross_sections(molecule, label_1, label_2, nu_arr_1 = [], nu_arr_2 = [],
                           sigma_arr_1 = [], sigma_arr_2 = [], **kwargs):

    wl_1 = 1.0e4/nu_arr_1
    wl_2 = 1.0e4/nu_arr_2

    if not os.path.exists('./plots/'):
        os.mkdir('./plots')

    print("\nComparing cross-sections of", molecule)

    #***** Make wavenumber plot *****#
    fig = plt.figure()
    ax = plt.gca()

    ax.set_yscale("log")
 #   ax.set_xscale("log")

    ax.plot(nu_arr_1, sigma_arr_1, lw=0.3, alpha = 0.5, color= 'crimson', label = (molecule + r' Cross Section ' + label_1))
    ax.plot(nu_arr_2, sigma_arr_2, lw=0.3, alpha = 0.5, color= 'royalblue', label = (molecule + r' Cross Section ' + label_2))

 #   ax.set_xlim([1800.2, 1804.1])
 #   ax.set_ylim([3.0e-22, 3.0e-18])
    ax.set_xlim([995.0, 1005.0])
    ax.set_ylim([1.0e-23, 1.0e-20])


    ax.set_ylabel(r'Cross Section (cm$^2$)', size = 14)
    ax.set_xlabel(r'Wavenumber (cm$^{-1}$)', size = 14)

    legend = plt.legend(loc='upper right', shadow=False, frameon=False, prop={'size':10})

    plt.tight_layout()

    plt.savefig('./plots/' + molecule + '_compare_' + label_1 + '_' + label_2 + '_Region_2.pdf')

    if (1 == 2):

        #***** Make wavelength plot *****#
        fig = plt.figure()
        ax = plt.gca()

        ax.set_yscale("log")
        ax.set_xscale("log")

        min_sigma = 1.0e-30
        max_sigma = 10.0**(np.max(np.ceil(np.log10(sigma_arr_1 + 1.0e-250) / 2.0) * 2.0))

        ax.plot(wl_1, sigma_arr_1, lw=0.3, alpha = 0.5, color= 'crimson', label = (molecule + r' Cross Section ' + label_1))
        ax.plot(wl_2, sigma_arr_2, lw=0.3, alpha = 0.5, color= 'royalblue', label = (molecule + r' Cross Section ' + label_2))

        ax.set_xticks([0.4, 0.6, 0.8, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0])
        ax.set_xticklabels(['0.4', '0.6', '0.8', '1', '2', '4', '6', '8', '10'])

     #   ax.set_ylim([min_sigma, max_sigma])
     #   ax.set_xlim([0.4, 10])
     #   ax.set_ylim([1.0e-25, 1.0e-18])
     #   ax.set_xlim([1.1, 1.7])

        ax.set_ylim([8.0e-25, 1.0e-22])
        ax.set_xlim([1.2, 1.201])


        ax.set_ylabel(r'Cross Section (cm$^2$)', size = 14)
        ax.set_xlabel(r'Wavelength (μm)', size = 14)

    #    ax.text(0.5, 5.0e-12, (r'T = ' + str(temperature) + r' K, P = ' + str(pressure) + r' bar'), fontsize = 10)

        legend = plt.legend(loc='upper right', shadow=False, frameon=False, prop={'size':10})

        plt.tight_layout()

        plt.savefig('./plots/' + molecule + '_comparison_' + label_1 + '_' + label_2 + '_zoom_a0.pdf')

    print("\nPlotting complete.")


def read_cross_section_file(species, database, filename, isotope = 'default',
                            ionization_state = 1, linelist = 'default', output_dir = './output/'):

    """
    Read in a previously computed cross section on a user's machine, and return the 
    wavenumber and absorption cross section columns.

    Parameters
    ----------
    species : String
        Molecule for which the cross section was computed.
    database : String
        Database the line list was downloaded from.
    filename : String
        Full file name of the cross section file(eg. 'this_is_my_file.txt').
    isotope : String
        Isotopologue of the molecule for which the cross-section is to be created. The default is 'default'.
    ionization_state : int
        Ionization state, in case of an atomic species. The default is 1. 
    linelist : String
        Line list that is being used. HITRAN/HITEMP/VALD used as the line list name for these
        databases respectively. ExoMol has its own named line lists. The default is 'default'. 
    output_dir = String
        'Prefix' of the directory containing the cross section files. If the files were downloaded
        using our script with no modifications, output_dir will end in '/output'. The default is './output/'. 

    Returns
    -------
    nu : numpy array
        Wavenumber column of cross section file.
    sigma : numpy array
        Absorption cross section column of cross section file.

    """
    database = database.lower()

    if database in ['hitran', 'hitemp']:
        molecule_dict = HITRAN.create_id_dict()
        mol_id = molecule_dict.get(species)
        if isotope == 'default':
            isotope = isotopologueName(mol_id, 1)
        else:
            isotope = isotopologueName(mol_id, isotope)
        isotope = HITRAN.replace_iso_name(isotope)

    if database == 'exomol':
        if isotope == 'default':
            isotope = ExoMol.get_default_iso(species)

    if database == 'vald':
        ion_roman = ''
        for i in range(ionization_state):
            ion_roman += 'I'


    if linelist == 'default':
        if database == 'exomol':
            temp_isotope = re.sub('[(]|[)]', '', isotope)
            linelist = ExoMol.get_default_linelist(species, temp_isotope)
        if database == 'hitran':
            linelist = 'HITRAN'
        if database == 'hitemp':
            linelist = 'HITEMP'
        if database == 'vald':
            linelist = 'VALD'

    if database == 'vald':
        tag = '(' + ion_roman + ')'
    else:
        tag = isotope

    if (database == 'exomol'):
        output_directory = (output_dir + species + '  ~  (' + tag + ')/' +
                           'ExoMol' + '/' + linelist + '/')
    else:
        output_directory = (output_dir + species + '  ~  ' + tag + '/' +
                           linelist + '/')

    if os.path.exists(output_directory):
        file  = output_directory + filename
        with open(file, 'r') as f:
            contents = csv.reader(f, delimiter=' ')
            nu = []
            sigma = []
            for cols in contents:
                nu.append(float(cols[0]))
                sigma.append(float(cols[1]))

        nu = np.asarray(nu)
        sigma = np.asarray(sigma)
        return  nu, sigma

    else:
        print("You don't seem to have a local folder with the parameters you entered.\n")

        if not os.path.exists(output_dir + '/'):
            print("----- You entered an invalid output directory into the cross_section() function. Please try again. -----")
            sys.exit(0)

        elif not os.path.exists(output_dir + '/' + species + '  ~  (' + tag + ')/'):
            print("----- There was an error with the molecule + isotope you entered. Here are the available options: -----\n")
            for folder in os.listdir(output_dir + '/'):
                if not folder.startswith('.'):
                    print(folder)
            sys.exit(0)

        else:
            print("There was an error with the line list. These are the linelists available: \n")
            for folder in os.listdir(output_dir + '/' + species + '  ~  (' + tag + ')/'):
                if not folder.startswith('.'):
                    print(folder)
            sys.exit(0)

def cross_section_collection(new_x, new_y, collection = []):
    '''
    Add new cross section (that is, wavelength and absorption cross section) to the collection

    Parameters
    ----------
    new_x : list
        DESCRIPTION.
    new_y : list
        DESCRIPTION.
    collection : 2d list, optional
        DESCRIPTION. The default is [].

    Returns
    -------
    collection : TYPE
        DESCRIPTION.

    '''

    collection.append([new_x, new_y])

    return collection
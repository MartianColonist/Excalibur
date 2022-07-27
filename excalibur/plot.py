import pandas as pd
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoLocator, FormatStrFormatter, \
                              FuncFormatter, ScalarFormatter, NullFormatter
from scipy.ndimage import gaussian_filter1d


def plot_spectra(spectra, spectra_labels, filename, plot_dir = './plots/', 
                 x_min = None, x_max = None, y_min = None, y_max = None,
                 color_list = [], smooth_data = True, std = 1000, **kwargs):
    """
    Generate a plot of a cross_section file, in both wavelength and wavenumber

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
    for spectrum in spectra:
        spectrum[1] = spectrum[1] + 1.0e-250
    
    nu_min, nu_max, sigma_min, sigma_max = find_min_max_nu_sigma(spectra)
    
    fig, ax = plt.subplots()
    
    # Define colors for plotted spectra (default or user choice)
    if (color_list == []):   # If user did not specify a custom colour list
        colors = ['green', 'crimson', 'black', 'darkgrey', 'navy', 'brown']
    else:
        colors = color_list
    
    # plot each spectra in the array
    for i in range(len(spectra)):
        (nu, spec) = spectra[i]
        wl = 1.0e4/nu # Convert wavenumber (in cm^-1) to wavelength (in um)
        
        if smooth_data == True:
            spec = gaussian_filter1d(spec, std)
        
        ax.loglog(wl, spec, lw=0.3, alpha = 0.5, color= colors[i], label = spectra_labels[i]) 
    
        # Set y range
    if (y_min == None):
        y_min = 1.0e-30
    if (y_max == None):
        y_max = sigma_max

    ax.set_ylim(y_min, y_max * 10)
    ax.set_xlim(1e4/nu_max, 1e4/nu_min) # Convert wavenumber (in cm^-1) to wavelength (in um)
    
    ax.set_ylabel(r'Cross Section (cm$^2$)', size = 14)
    ax.set_xlabel(r'Wavelength (μm)', size = 14)
    
    legend = plt.legend(loc='upper right', shadow=False, frameon=False, prop={'size':10})
    
    for legline in legend.legendHandles:
        legline.set_linewidth(1.0)
        
    plt.tight_layout()
    plt.savefig('./plots/' + filename)
    
    print("\nPlotting complete.")


def spectra_collection(new_x, new_y, collection = []):
    '''
    Add new spectrum (that is, wavelength and cross section) to the collection

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
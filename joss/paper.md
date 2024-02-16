---
title: '`Excalibur`: An Open Source Molecular and Atomic Cross Section Computation Code for Substellar Atmospheres'
tags:
  - Python
  - astronomy
  - exoplanets
  - spectroscopy
  - atmospheric retrieval
  - atmospheric models
  - JWST
authors:
  - name: Arnav Agrawal
    orcid:
    affiliation: "3"
  - name: Ryan J. MacDonald
    orcid: 0000-0003-4816-3469
    affiliation: "1, 2, 3"
affiliations:
  - name: Department of Astronomy, University of Michigan, 1085 S. University Ave., Ann Arbor, MI 48109, USA
    index: 1
  - name: NHFP Sagan Fellow
    index: 2
  - name: Department of Astronomy and Carl Sagan Institute, Cornell University, 122 Sciences Drive, Ithaca, NY 14853, USA
    index: 3
date:
bibliography: paper.bib

--- 

# Summary

Atmospheric studies of exoplanets and brown dwarfs are a cutting-edge and rapidly-evolving area of astrophysics research. Powerful new telescopes, such as the James Webb Space Telescope (JWST) and the upcoming Extremely Large Telescopes (ELTs), are able to capture in detail spectra of planets and brown dwarfs and thereby probe their chemical composition and physical properties. Calculating models of exoplanet or brown dwarf spectra requires knowledge of the wavelength-dependent absorption of light (cross sections) by the molecules and atoms in the atmosphere. Without accurate cross sections, one cannot reliably measure the chemical composition of substellar atmospheres. 

Cross sections are typically pre-computed on a grid of pressure and temperatures from large databases of quantum mechanical transitions (line lists), such as ExoMol [@Tennyson:2020], HITRAN [@Gordon:2022], HITEMP [@Rothman:2010], and VALD [@Pakhomov:2017]. However, the process of calculating cross sections from line lists is often computationally demanding and has required complex and specialised tools. We aim here to lower the access barrier for users to learn how to calculate molecular and atomic cross sections. 

`Excalibur` is a fully Python package that rapidly calculates cross sections from atomic and molecular line lists. `Excalibur` includes modules to automatically download molecular line lists from online databases and compute cross sections on a user-specified temperature, pressure, and wavenumber grid. `Excalibur` requires only CPUs and can run on a user's laptop (for smaller line lists) or on a large cluster in parallel (for billions of lines). `Excalibur` includes in-depth Jupyter tutorials in the online documentation. Finally, `Excalibur` is intended not only for research purposes, but as an educational tool to demystify the process of making cross sections for atmospheric models.

# Computing Molecular and Atomic Cross Sections with `Excalibur`

The purpose of the `Excalibur` package is schematically represented in \autoref{fig:EXCALIBUR_architecture}. Here we walk through this flowchart, highlighting 3 major use cases of, and explain `Excalibur`'s role in the broader process of modelling exoplanetary atmospheres.

![The role and applications of the `Excalibur` Python package. `Excalibur` can download molecular and atomic line lists and calculate the corresponding absorption cross sections as a function of temperature, pressure, and wavenumber. Cross sections made by `Excalibur` can be used in radiative transfer codes to calculate model spectra of exoplanet and brown dwarf atmospheres. \label{fig:EXCALIBUR_architecture}](figures/EXCALIBUR_JOSS_Figure.png){width=100%}

The first use of `Excalibur` is to download existing molecular line lists from online databases. `Excalibur`'s `summon` function can automatically download lines lists from ExoMol and HITRAN/HITEMP (for the latter a user must make an account on https://hitran.org/) and reformat the line lists into space-efficient HDF5 files. Ancillary input files required to calculate cross sections, such as partition functions and pressure broadening files, are also downloaded automatically. Alternatively, the user may manually download a line list from their respective websites and point `Excalibur` to the directory hosting the files. VALD line lists must be downloaded manually by a user with an account on http://vald.astro.uu.se/ (given the terms of use for VALD3), but we provide instructions on how to do this in the `Excalibur` documentation. `Excalibur` currently supports ExoMol, HITRAN, HITEMP, and VALD line lists, though we welcome user requests for additional line list databases support. Once a line list has been downloaded, the user can move onto the next major use case of `Excalibur`, computing cross sections.

The foremost feature of `Excalibur` is its ability to straightforwardly compute atomic and molecular cross sections at high speeds (typically > 100,000 lines per second on one CPU). `Excalibur` is widely accessible, as it does not require GPUs, can run on a standard laptop, and as a fully Python code it is easy for beginners to install and use. To compute a cross section, user's simply call `Excalibur`'s `compute_cross_section` function, specifying the location of the line list, the temperature and pressure, and the wavenumber range. More advanced users can specify custom settings via optional arguments (e.g. Voigt wing cutoffs, intensity cutoffs, or a user-provided pressure broadening file). The documentation and function doc strings explains the various arguments users can provide to `compute_cross_section`.  The computed cross section is output by default as a .txt file in the `output` folder on the user's machine, but `Excalibur` also offers utility functions to combine multiple cross sections (e.g. a grid of temperature and pressures for one or more chemical species) into a HDF5 cross section database. These HDF5 files can then be readily plugged into the user's favourite exoplanet or brown dwarf modelling or retrieval code.

Ryan up to here.

To bring the cross sections to life visually and to complete the flow from 'Inputs' to 'Outputs' in \autoref{fig:EXCALIBUR_architecture}, Excalibur lends the capability to create publication-quality cross section plots. This can be as simple as reading in a single file, and letting Excalibur plot the cross section with default settings. For users that choose to dive in deeper, Excalibur provides the functionality to overplot multiple cross sections, apply a Gaussian filter to smooth the data, change the axes between linear and logarithmic, and more. Depicted in \autoref{fig:EXCALIBUR_architecture} are 3 examples of the beautiful plots Excalibur can produce. From left to right, they show cross section plots produced for molecules, atoms and ions, and molecular isotopologues.

The cross sections Excalibur produces are directly inputted into radiative transfer codes, as shown in \autoref{fig:EXCALIBUR_architecture}. These codes take all the cross sections we have computed for various species and determine which combination of atoms and molecules "best" (in a probabilistic sense) models the exoplanetary atmosphere. That is, by using radiative transfer algorithms to sample different factor spaces of atoms/molecules, we can determine the statistical confidence for a molecule being present in the exoplanetary atmosphere (as compared to the spectra observed by a telescope). This application is depicted in the bottom of \autoref{fig:EXCALIBUR_architecture}.

# Statement of Need
The recent launch of JWST has with it the goal of discovering hundreds of new exoplanets. Naturally, scientists will want to analyze transmission spectra of these newly discovered exoplanets to learn more about their composition and atmosphere. By allowing scientists to easily create molecular and atomic cross sections (the building block for models of transmission spectra), Excalibur bridges the gap between observation and theory.

There are other codes that compute cross sections (CITE), but there are a couple of key features that are unique to Excalibur. Primarily, Excalibur is intended to be used as an educational stepping stone for scientists to begin to use cross sections in their own research. While the focus of cross sections in this paper has been around their use in studying exoplanets, there are applications in LIST 2-3 OTHER FIELDS. Excalibur is written completely in Python, making it easier for users to run the code, and even change certain algorithms or scripts if they choose. Thorough documentation and in-depth tutorials are provided and maintained on the [website](https://excalibur-alpha.readthedocs.io/en/latest/), showing users how to customise parameters and other variables to their preference.

# Future Developments

Connects with POSEIDON retrieval code?

# Documentation

Documentation for `Excalibur`, with step-by-step tutorials illustrating research applications, is available at [https://excalibur-alpha.readthedocs.io/en/latest/](https://excalibur-alpha.readthedocs.io/en/latest/). 

# Similar Tools

PyExoCross, HeliosK, RADIS
- install and use
- compare speeds (lines/sec)

# Acknowledgements

# References
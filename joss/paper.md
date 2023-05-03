---
title: '`Excalibur`: Description'
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
    affiliation: "4"
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
  - name: Cornell University, Ithaca, NY 14853, USA
    index: 4
date:
bibliography: paper.bib

aas-doi: 
aas-journal: Astrophysical Journal

--- 

# Summary

The field of exoplanet research is quickly growing and requires software that can keep up with the demands of new scientific data. Powerful new telescopes, such as the James Webb Space Telescope (JWST), are able to capture in detail the transmission spectra of exoplanets by measuring the fraction of light the planet blocks as it crosses in front of its host star. By comparing observed transmission spectra to known models, scientists are able to learn about the planet's atmosphere, chemical composition, and more.  However, in order for the spectra to be useful, scientists must have accurate models to compare to. The process for creating these models has historically been convoluted and inaccessible for researchers. The core feature of these models is the effective size of an atom or molecule at different wavelengths, which we term the cross section. Atomic and molecular cross sections are derived from line lists. 

LINE LIST IS NOT SIMPLY 2 COLUMNS. Simply put, a line list is composed of two columns, one being the wavelength at which the rotational, vibrational, or electronic transition takes place, and the other being the probability of the transition occurring. DOES THE END OF THIS PARAGRAPH MAKE SENSE TO A RANDOM PERSON?

Excalibur is a Python package that opens up the world of cross sections for the average researcher. It allows scientists to download atomic and molecular line lists from online databases, compute cross sections with various factor spaces, and create beautiful plots of the cross sections. With in-depth tutorials and instructive scientific explanations, Excalibur is intended not just as a research package, but as an educational tool. 

# Computing and Plotting Cross Sections with Excalibur
IN GENERAL, IN THIS PARAGRAPH I WILL WANT TO TIE BACK INTO THE FIGURE A BIT MORE.
As mentioned, Excalibur provides 3 major functions that together fully encompass the process of creating and displaying a cross section. We briefly describe what each of these functions entails, referring to the Excalibur workflow shown in \autoref{fig:EXCALIBUR_architecture} throughout.

The first application of Excalibur is to download molecular and atomic line lists from online databases. Excalibur currently supports downloads from three databases, ExoMol, HITRAN (including it's high-temperature database HITEMP), and VALD. A user is able to specify various parameters, like the name of the molecule or database, to indicate which line list is desired. The downloaded line list is placed in an 'input' folder on the user's computer and converted to an HDF5 file format to save storage space. This initial step is shown in the top row of \autoref{fig:EXCALIBUR_architecture}.

With a line list now available, the user can move onto the next major use case of Excalibur, computing cross sections at high speeds (~ 10,000 transitions per second). Excalibur provides the user with as much freedom as desired in the calculation of the cross section. A user can choose to specify just a few basic parameters (molecule, temperature, pressure, etc.) or define every aspect of how the cross section is computed, from the type of broadening used to (DESCRIBE ANOTHER PARAMETER THAT IS RELEVANT). The Excalibur package, and this step in particular, was designed to be as flexible as possible to fit the wants of any researcher. Once computed, the cross section is placed as a text file in an 'output' folder on the user's machine. This use case of Excalibur represents the flow from line list to cross section in \autoref{fig:EXCALIBUR_architecture}.

Finally, Excalibur lends the capability to create quality plots to be used in research publications, as illustrated by the various colorful plots in \autoref{fig:EXCALIBUR_architecture}, one step down from the cross sections. There are many ways the user can customize plots, from overplotting multiple cross sections to using a Gaussian filter to smooth data.

The overall architecture and use of the Excalibur package is diagrammed in \autoref{fig:EXCALIBUR_architecture}. The top half represents the steps needed to produce a complex scientific model of multiple molecular or atomic cross sections. The bottom half highlights the real spectrum a telescope might observe from watching a transiting exoplanet, and how this spectrum is defined by the molecules and atoms within the exoplanet's atmosphere. Comparison of these two models is necessary for a complete understanding of the exoplanet.

![WRITE DESCRIPTION OF FIGURE HERE. \label{fig:EXCALIBUR_architecture}](figures/EXCALIBUR_JOSS_Figure.png){width=100%}

# Statement of Need
The recent launch of JWST has with it the goal of discovering hundreds of new exoplanets. Naturally, scientists will want to analyze transmission spectra of all of these newly discovered exoplanets to learn more about their composition and atmosphere. By allowing scientists to easily create molecular and atomic cross sections (the building block for models of transmission spectra), Excalibur bridges the gap between observation and theory.

There are other codes that compute cross sections (CITE), but there are a couple of key features that are unique to Excalibur. Primarily, Excalibur is intended to be used as an educational stepping stone for scientists to begin to use cross sections in their own research. While the focus of cross sections in this paper has been around their use in studying exoplanets, there are applications in LIST 2-3 OTHER FIELDS. Excalibur is written completely in Python, making it easier for users to run the code, and even analyze certain algorithms or scripts if they choose. Thorough documentation and in-depth tutorials are provided on the website, showing users how to customize parameters and other variables to their preference.

# Future Developments

ANY FUTURE DEVELOPMENTS TO MENTION?

# Documentation

Documentation for `Excalibur`, with step-by-step tutorials illustrating research applications, is available at [https://excalibur-alpha.readthedocs.io/en/latest/](https://excalibur-alpha.readthedocs.io/en/latest/). 

# Similar Tools

ASK RYAN.

# Acknowledgements

# References
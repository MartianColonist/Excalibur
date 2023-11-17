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

The field of exoplanet research is quickly growing and requires software that can keep up with the demands of new scientific data. Powerful new telescopes, such as the James Webb Space Telescope (JWST), are able to capture in detail the transmission spectra of exoplanets by measuring the fraction of light the planet blocks as it crosses in front of its host star. By comparing observed transmission spectra to known models, scientists are able to learn about the planet's atmosphere, chemical composition, and more.  However, in order for the spectra to be useful, scientists must have accurate models to compare to. The process for creating these models has historically been convoluted and inaccessible for researchers. The core input of these models is the effective size of an atom or molecule at different wavelengths, which we term the cross section. Atomic and molecular cross sections are derived from line lists- complex data sets that contain multitudes of quantitative information on the way an atom or molecule interacts with different wavelengths of light. Line lists, and the cross sections they generate, are essential in devising models to study exoplanetary atmospheres. 

Excalibur is a Python package that opens up the world of cross sections for the average researcher. It allows scientists to download atomic and molecular line lists from online databases, compute cross sections with various factor spaces, and create beautiful plots of the cross sections. With in-depth tutorials and instructive scientific explanations, Excalibur is intended not just as a research package, but as an educational tool. 

Critically, a line list details the wavelength at which a rotational, vibrational, or electronic transition takes place within a species, and the probability of the transition occurring. (don't know if this last sentence is necessary) In short, line lists and the cross sections they generate are essential in devising models to study exoplanetary atmospheres. 

# Computing and Plotting Cross Sections with Excalibur
IN GENERAL, IN THIS PARAGRAPH I WILL WANT TO TIE BACK INTO THE FIGURE A BIT MORE.
As mentioned, Excalibur provides 3 major functions that together encompass the process of creating and displaying a cross section. We briefly describe what each of these functions entails, referring frequently to the Excalibur workflow shown in \autoref{fig:EXCALIBUR_architecture} throughout. In general, \autoref{fig:EXCALIBUR_architecture} can be visualized as a 

The overall architecture and workflow of the Excalibur package is diagrammed in \autoref{fig:EXCALIBUR_architecture}. We walk through this flowchart, describe each of the major functions that Excalibur performs, and explain Excalibur's role in the broader process of modelling exoplanetary atmospheres.

Starting from the top of \autoref{fig:EXCALIBUR_architecture}, the first application of Excalibur is to download existing molecular and atomic line lists from online databases. Excalibur currently supports downloads from three databases, ExoMol, HITRAN (including it's high-temperature extension HITEMP), and VALD. A user specifies various parameters, like the name of the molecule or database, to indicate which line list is desired. The downloaded line list is placed in an 'input' folder on the user's computer and converted to an HDF5 file format to save storage space. With a line list now available, the user can move onto the next major use case of Excalibur, computing cross sections.

The foremost feature of Excalibur is its ability to compute atomic and molecular cross sections at high speeds (~ 10,000 transitions per second) using the algorithm defined in (IS THERE A PAPER/THESIS TO REFERENCE HERE?). The entire cross section calculation process is represented by the Excalibur logo in \autoref{fig:EXCALIBUR_architecture}. As depicted in the figure, a line list goes in, Excalibur performs its computations, and a cross section comes out. There are many subsidiary functions involved in the computation step, and \autoref{fig:EXCALIBUR_architecture}'s depiction of all this into one logo is meant for simplification and clarity, not as a hidden 'black box'. Excalibur actually provides the user with significant freedom in the calculation of the cross section. A user can choose to specify just a few basic parameters (molecule, temperature, pressure, etc.) or define finer aspects of how the cross section is computed, from the type of broadening used to the Voigt line wing cutoff. The computed cross section is placed as a text file in an 'output' folder on the user's machine. This output text file brings us to the next use case of Excalibur, plotting cross sections.

To bring the cross sections to life visually, and to complete the flow from 'Inputs' to 'Outputs' in \autoref{fig:EXCALIBUR_architecture}, Excalibur lends the capability to create publication-quality cross section plots. This can be as simple as reading in a single file, and letting Excalibur plot the cross section with default settings. For users that choose to dive in deeper, Excalibur provides the functionality to overplot multiple cross sections, apply a Gaussian filter to smooth the data, change the axes between linear and logarithmic, and more. Depicted in \autoref{fig:EXCALIBUR_architecture} are 3 examples of the beautiful plots Excalibur can produce. From left to right, they show cross section plots produced for molecules, atoms and ions, and molecular isotopologues.

The cross sections Excalibur produces are absolutely essential to understanding exoplanetary atmospheres. 
FINAL SENTENCE CONNECTING THE OUTPUTS TO RADIATIVE TRANSFER AND RADIATIVE TRANSFER TO TRANSMISSION SPECTRA.

The top half of the figure, above the dashed line, represents the steps Excalibur follows to produce a model of multiple molecular or atomic cross sections. The bottom half, below the dashed line, highlights the real spectrum a telescope might observe from watching a transiting exoplanet, and how this spectrum is defined by the molecules and atoms within the exoplanet's atmosphere. Excalibur serves the top half of \autoref{fig:EXCALIBUR_architecture} by providing researchers with the tools needed to quickly and accurately create cross section models. SENTENCE CONNECTING THE TWO HALVES OF THE DIAGRAM. For the remainder of this section, we walk through \autoref{fig:EXCALIBUR_architecture} from the top-down, while describing in more detail the major functions Excalibur can perform. 

![WRITE DESCRIPTION OF FIGURE HERE. \label{fig:EXCALIBUR_architecture}](figures/EXCALIBUR_JOSS_Figure.png){width=100%}

# Statement of Need
The recent launch of JWST has with it the goal of discovering hundreds of new exoplanets. Naturally, scientists will want to analyze transmission spectra of these newly discovered exoplanets to learn more about their composition and atmosphere. By allowing scientists to easily create molecular and atomic cross sections (the building block for models of transmission spectra), Excalibur bridges the gap between observation and theory.

There are other codes that compute cross sections (CITE), but there are a couple of key features that are unique to Excalibur. Primarily, Excalibur is intended to be used as an educational stepping stone for scientists to begin to use cross sections in their own research. While the focus of cross sections in this paper has been around their use in studying exoplanets, there are applications in LIST 2-3 OTHER FIELDS. Excalibur is written completely in Python, making it easier for users to run the code, and even change certain algorithms or scripts if they choose. Thorough documentation and in-depth tutorials are provided and maintained on the website, showing users how to customize parameters and other variables to their preference.

WORTH MENTIONING POSEIDON HERE?

# Future Developments

ANY FUTURE DEVELOPMENTS TO MENTION?

# Documentation

Documentation for `Excalibur`, with step-by-step tutorials illustrating research applications, is available at [https://excalibur-alpha.readthedocs.io/en/latest/](https://excalibur-alpha.readthedocs.io/en/latest/). 

# Similar Tools

ASK RYAN.

# Acknowledgements

# References
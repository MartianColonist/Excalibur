Introduction
------------

This page is split into two sections. The first part assumes that the user has
little, if any, prior knowledge of molecular and atomic cross sections.
It explains what a cross section is, how it is calculated, and why it is useful.
The second part bullet points the many functions that the Excalibur package
can perform. Even for more knowledgeable users, we recommend taking
the extra few minutes to read both sections when working with Excalibur for the first time.

Part 1: A Beginner's Guide to Cross Sections
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To introduce the science of cross sections, we first need to define some terms. For simplicity,
we'll focus on atoms, but the same concepts can be applied directly to the study of molecules as well.
A 'line list' is a list of all electron transitions within an atom. These include not only the familiar
transitions of electrons between atomic energy levels, but also transitions between different rotational and
vibrational modes of the atom. A line list can be as simple as a two column table, one being the wavelength
at which the transition takes place and the other the probability for the transition to take place.

Thus, a line list can be used to identify the wavelengths at which the atom absorbs more or less light. A higher
transition probability leads to more light absorbed at that wavelength, which leads to the atom appearing more opaque
and thus having a larger effective cross sectional area. A plot of the effective cross sectional area as a function of wavelength
is what we refer to as the 'cross section' of an atom. See the image below for an example of how different cross sections may appear on a plot.
Notice how the Na cross section, which is the only atom plotted, has fewer intricacies and detail to the plot. This comes from the fact 
that atomic line lists are generally much smaller than molecular line lists.

.. image:: ./images/Example_cross_section.pdf
  :width: 600
  :alt: Image overplotting 3 example cross sections for TiO, Na, and CH:subscript:`4`

There are a couple of other parameters that are important to the creation of cross sections. First of all, line lists can be
very temperature dependent, so a line list that was calculated in the laboratory at 300 Kelvin (K) may not be suitable to calculate
a cross-section for an atom at 1000K [`Hargreaves & Gordon<https://meetingorganizer.copernicus.org/EPSC-DPS2019/EPSC-DPS2019-919-1.pdf>`_]. 
Secondly, it is vital to keep in mind that the wavelengths at which electron transitions take place are not quite 
singular values as may appear in a line list, but a small range centered around those values.
There are thus many factors that can lead to broadening of the atomic spectra; the two biggest causes of spectral broadening are
pressure and temperature, which is why Excalibur computes cross sections at a certain pressure and temperature 
[`NIST<https://www.nist.gov/pml/atomic-spectroscopy-compendium-basic-ideas-notation-data-and-formulas/atomic-spectroscopy-6>`_].

Cross sections are used extensively, especially in atmospheric and climate science. Their most prominent purpose in astronomy
is to interpret spectra obtained from observations of exoplanets [`Gordon<https://hitran.org/media/refs/HITRAN-2020.pdf>`_]. As an exoplanet
crosses in front of its host star, part of the star's light is blocked by the exoplanet's atmosphere. By measuring exactly
how much light is blocked at different wavelengths and comparing that to molecular cross section models, scientists can determine the
molecular composition of the exoplanet's atmosphere.

Part 2: Excalibur's Features
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Now that we understand the importance of cross sections to exoplanet research, we elaborate on what Excalibur is and what
it does. We highlight the many unique features that allow Excalibur to help us achieve our science goals.

Functionality
"""""""""""""

* With Excalibur, users can download a line list from an online database, compute an atomic, molecular, or ionic cross section at a certain pressure and temperature, and then plot the cross section.
* Line list databases currently supported are ExoMol, HITRAN, HITEMP, and VALD.
* Excalibur plotting scripts offer the ability to create beautiful, organized plots to include in research papers and presentations.

User-Friendly
"""""""""""""

* Requires minimal coding knowledge
* Detailed, step-by-step tutorials provided for every feature of Excalibur
* Written in Python, making it easy for users to understand the underlying code
* Thorough, well-maintained documentation

High-performance
""""""""""""""""

* Cross section computations are performed at ≈10,000 transitions per second (COMPARE TO OTHER CODES OUT THERE?)
* File compression makes it possible to download and store multiple line lists and cross sections on just a personal laptop

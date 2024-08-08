Cthulhu's documentation
=============================================

Welcome to Cthulhu's lair!

Cthulhu is a Python package that allows users to effortlessly create molecular
and atomic absorption cross sections from online line list databases for 
exoplanet or brown dwarf radiative transfer application.

Cross sections encode the interaction strength between light and gases in an 
atmosphere, making their accurate computation a core requirement for 
calculating spectra of substellar atmospheres.

Cthulhu's features currently include:

* Fully Python, lowering the barrier to learn how cross sections are made.
* Automatically download line lists from online databases, including ExoMol, HITRAN, HITEMP, and VALD.
* Calculate cross sections at high speeds (≈ 100,000 transition per second).
* Generate HDF5 opacity files to plug into your favourite radiative transfer code.
* Produce visually-appealing plots out of the box. 
* Runs on your laptop (for smaller line lists) or on a supercomputer.

If this is your first time using Cthulhu, please proceed to the Introduction page on the sidebar.

Cthulhu is available under the BSD 3-Clause License. If you use Cthulhu in your own work, 
please cite Agrawal & MacDonald (2024).

.. toctree::
   :maxdepth: 1
   :hidden:

   content/introduction
   content/installation

.. toctree::
  :maxdepth: 1
  :caption: Guide:

  content/getting_started
  content/tutorials

.. toctree::
 :maxdepth: 2
 :caption: Code Documentation:

 content/contributing
 autoapi/index

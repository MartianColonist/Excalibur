.. Excalibur-alpha documentation master file, created by
   sphinx-quickstart on Thu Aug 27 17:19:57 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Excalibur-alpha's documentation.
=============================================

Welcome to Excalibur! 

Excalibur is a Python package that allows users to effortlessly create molecular
and atomic cross sections to aid in their own research, including for the creation of
accurate atmospheric retrieval models. 

Excalibur's features currently include:

* Capability to download line lists from online databases, including ExoMol, HITRAN, HITEMP, and VALD
* Cross section computations performed at high speeds (≈ 10,000 transition per second)
* Plotting scripts to create beautiful, organized graphs to include in research papers and presentations
* Completely written in Python, making it easy for researchers to use and adapt the code
* Designed to be used from a personal laptop, no supercomputers needed

If this is your first time using Excalibur, please proceed to the Introduction page on the sidebar.

Excalibur is available under the BSD 3-Clause License. If you use Excalibur in your own work, please cite [INSERT CITATION].

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

 autoapi/index

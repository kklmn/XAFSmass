XAFSmass
========

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.1219124.svg)](http://dx.doi.org/10.5281/zenodo.1219124)

<p align="center">
  <img src="doc/_images/1powder_150.png" width=320 />
</p>

XAFSmass is a python program for calculating the mass of XAFS [X-ray Absorption Fine Structure] samples.
The chemical formula parser understands parentheses and weight percentage, also in nested form.
XAFSmass reports the quantity (weight, thickness or pressure) together with the expected height
of the absorption edge. The GUI is provided by Qt.

Copyright (c) 2015 Konstantin Klementiev and Roman Chernikov under the MIT License terms

Dependencies
------------

numpy, pyparsing and matplotlib are required. Qt must be provided by either
PyQt5, PySide2, PyQt6 or PySide6 by means of qtpy.

Installation
------------

Unzip the .zip file into a suitable directory and run ``python XAFSmassQt.py``.
On Windows, run ``pythonw XAFSmassQt.py`` or give it a .pyw extension to
suppress the console window.

You may want to run ``python setup.py install`` in order to put the XAFSmass
package to the standard location.

Documentation
-------------

See it on http://xafsmass.readthedocs.io
[![Documentation Status](https://readthedocs.org/projects/xafsmass/badge/?version=latest)](https://xafsmass.readthedocs.io)

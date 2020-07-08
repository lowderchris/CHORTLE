# CHORTLE
Coronal Hole Observer and Regional Tracker for Long-term Examination

## Overiew
Manifesting as regions of decreased emission in extreme ultraviolet (EUV) and x-ray wavelengths, coronal holes are the observational signatures of open solar magnetic field. Coronal plasma within these regions is free to flow outward along open magnetic field lines, resulting in reduced density and emission. Studying these coronal hole regions provides useful insights into their connection with open magnetic field and the evolution of the solar cycle.

A previous incarnation of this code, GACHD ([Lowder et al 2014](http://adsabs.harvard.edu/abs/2014ApJ...783..142L), [Lowder et al 2017](http://adsabs.harvard.edu/abs/2017SoPh..292...18L)), was developed and tested. This is a rewritten and revamped version of this code, written to utilize Python and SunPy.

## Goal
- Automated and adaptive coronal hole detection
- Input from multiple data sources, including SOHO/EIT, SDO/AIA, and STEREO/EUVI A&B
- Simple input - provide a CR number, and chortle will do the rest!
- Merging and bridging these datasets to cover solar cycles 23 and 24
- Tracking of coronal hole features to provide statistics on lifetimes and evolution

## Installation

Currently I've uploaded this code as a VERY ROUGH DRAFT of porting the basic methodology of the original code over from IDL to Python. It exists as a script, with all kinds of terrible work-in-progress qualities. Use at your own peril. Notably some... *improvements* to the SunPy library of code have broken critical aspects of the reprojection portion of this code, which needs to be corrected. Notably a separate version of the code, designed to work with EIT data, is in development and is not fully functional at the moment. Eventually this will be merged and modularized into one set of routines.

TL;DL - Don't use this code quite yet, but watch this space as I continue to develop it towards a more user-friendly and accurate version.

## Notes


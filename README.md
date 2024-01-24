# planetMagFields
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5140421.svg)](https://doi.org/10.5281/zenodo.5140421)

Routines to easily access information about magnetic fields of planets in our solar system and visualize them in both 2D and 3D. These require [NumPy](https://numpy.org/), [Matplotlib](https://matplotlib.org/) and [SciPy](https://www.scipy.org/) (pronounced "Sigh Pie"). Other than that, the following external libraries are used for a few different functions:

 - 2D plotting for map projections other than Hammer : [Cartopy](https://scitools.org.uk/cartopy/docs/latest/) library ( see more under [Projections](#projections) )
 - Potential extrapolation: [SHTns](https://bitbucket.org/nschaeff/shtns) library
 - Writing vts files for 3D visualisation: [PyEVTK](https://github.com/paulo-herrera/PyEVTK) library

# Documentation

The documentation is available here: [](https://ankitbarik.github.io/planetMagFields/)


# Acknowledgements

I would like to thank [Regupathi Angappan](https://github.com/reguang) for motivating me to create this package. I thank [Jon Aurnou](https://epss.ucla.edu/people/faculty/543/) for testing it out and pointing out runtime errors. Lastly, I would like to thank [Thomas Gastine](https://github.com/tgastine) for comparing the plots with real data and pointing out a normalization error which has been fixed.

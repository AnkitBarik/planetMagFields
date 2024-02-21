# planetMagFields
![Build workflow](https://github.com/AnkitBarik/planetMagFields/actions/workflows/main.yml/badge.svg)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/AnkitBarik/planetMagFields/HEAD?labpath=%2FExploreFieldsInteractively.ipynb)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5140421.svg)](https://doi.org/10.5281/zenodo.5140421)

Software to easily access and analyze information about magnetic fields of planets in our solar system and visualize them in both 2D and 3D.

[Prerequisites](#prerequisites)

[Installation](#installation)

[Features and examples](#features-and-examples)

[Jupyter notebook](#jupyter-notebook)

[Documentation](#documentation)

[Acknowledgements](#acknowledgements)

# Prerequisites

`planetMagFields` requires [NumPy](https://numpy.org/), [Matplotlib](https://matplotlib.org/) and [SciPy](https://www.scipy.org/). Other than that, the following external libraries are used for a few different functions:

 - 2D plotting for map projections other than Hammer : [Cartopy](https://scitools.org.uk/cartopy/docs/latest/) library
 - Potential extrapolation: [SHTns](https://bitbucket.org/nschaeff/shtns) library
 - Writing vts files for 3D visualisation: [PyEVTK](https://github.com/paulo-herrera/PyEVTK) library

# Installation

`planetMagFields` can be installed in a few different ways:

## Using `PYTHONPATH`

Download the package from the GitHub repository and add it
to `PYTHONPATH`:

```bash

$ git clone https://github.com/AnkitBarik/planetMagFields
$ export PYTHONPATH=$PYTHONPATH:/path/to/planetMagFields
```

## Using `setup.py`

You can also use `setup.py` to install `planetMagFields`:

```bash

$ git clone https://github.com/AnkitBarik/planetMagFields
$ cd planetMagFields
$ pytho3 setup.py install --user
```

Or using `pip`:

```bash

$ git clone https://github.com/AnkitBarik/planetMagFields
$ cd planetMagFields
$ python3 -m pip install . --user
```

# Features and examples

## The `Planet` class

This gives access to all the relevant properties of a planet and has methods to plot
the field and write a `vts` file for 3D visualization. Usage:

```python
from planetmagfields import Planet
p = Planet(name='earth',datDir='planetmagfields/data/')
```

This displays the some information about the planet

```bash
Planet: Earth
l_max = 13
Dipole tilt (degrees) = -9.410531
```

and gives access to
variables associated with the planet such as:

  * ``p.lmax`` : maximum spherical harmonic degree till which data is available
  * ``p.glm``, ``p.hlm``: the Gauss coefficients
  * ``p.Br`` : computed radial magnetic field at surface
  * ``p.dipTheta`` : dipole tilt with respect to the rotation axis
  * ``p.dipPhi`` : dipole longitude ( in case zero longitude is known, applicable to Earth )
  * ``p.idx`` : indices to get values of Gauss coefficients
  * ``p.model`` : the magnetic field model used. Available models can be obtained using the `get_models` function. Selects the latest available model when unspecified.

Example using ``IPython``:

```python
In [1]: from planetmagfields import Planet

In [2]: p = Planet(name='jupiter',model='jrm09')
Planet: Jupiter
Model: jrm09
l_max = 10
Dipole tilt (degrees) = 10.307870

In [3]: p.glm[p.idx[2,0]]      # g20
Out[3]: 11670.4

In [4]: p.hlm[p.idx[4,2]]      # h42
Out[4]: 27811.2
```

## 2D and 3D visualizations

`planetMagFields` can be used to produce both 2D and 3D visualizations of a planetary field.
While doing so, it also provides the option of potential extrapolation using the [SHTns](https://bitbucket.org/nschaeff/shtns/)
library. Two examples are shown below: Earth's surface field and a 3D visualization of Jupiter's field using
[Paraview](https://www.paraview.org/). This is done by writing a `.vts` file using the [PyEVTK](https://github.com/paulo-herrera/PyEVTK)
library.

<p align="center" width="100%">
<img src="https://raw.githubusercontent.com/AnkitBarik/planetMagFields/main/doc/_static/images/2d/earth2d.png" width="500">
</p>

<p align="center" width="100%">
<img src="https://raw.githubusercontent.com/AnkitBarik/planetMagFields/main/doc/_static/images/3d/jupiter3d.png" width=500>
</p>

# Jupyter notebook

For quick and easy visualization we include a [Jupyter](https://jupyter.org/) notebook with a binder link (see badge at the top).
This makes use of [Jupyter widgets](https://ipywidgets.readthedocs.io/) to provide dropdown lists of planets and available magnetic
field models for each as well as a slider for radial level, as shown below

<p align="center" width="100%">
<img src="https://raw.githubusercontent.com/AnkitBarik/planetMagFields/main/doc/_static/images/jupyter_screenshot2.png" width="500">
</p>

This plots the radial magnetic field at the chosen radial level and the corresponding magnetic field spectrum,

<p align="center" width="100%">
<img src="https://raw.githubusercontent.com/AnkitBarik/planetMagFields/main/doc/_static/images/jupyter_screenshot3.png" width="500">
</p>

# Documentation

Full list of features with examples as well as the magnetic field models used are described in detail in the documentation,
available here: https://ankitbarik.github.io/planetMagFields/


# Acknowledgements

I would like to thank [Regupathi Angappan](https://github.com/reguang) for motivating me to create this package, testing it and for writing the Jupyter notebook. I thank [Jon Aurnou](https://epss.ucla.edu/people/faculty/543/) for testing it out and pointing out runtime errors. I would like to thank [Thomas Gastine](https://github.com/tgastine) for comparing the plots with real data and pointing out a normalization error which has been fixed. Thanks a lot to [Arthus](https://github.com/arthus701) for adding the `setup.py`.

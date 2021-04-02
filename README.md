# planetMagFields
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Routines to visualize magnetic fields of planets in our solar system, in both 2D and 3D. For 2D plots, the routines use a Mollweide projection by default using the [cartopy](https://scitools.org.uk/cartopy/docs/latest/) library, but in the absence of it, it falls back to the internal Hammer projection.

The following external libraries are used for a few different functions:

 - Potential extrapolation: [SHTns](https://bitbucket.org/nschaeff/shtns) library
 - Writing vts files for 3D visualisation: [PyEVTK](https://github.com/paulo-herrera/PyEVTK) library

# The `planet` class

This gives access to all the relevant properties of a planet and has methods to plot
the field and write a `vts` file for 3D visualization. Usage:

```python
from planetmagfields import *
p = planet(name='earth',datDir='planetmagfields/data/')
```

This displays the some information about the planet

```
Planet: Earth
Depth (fraction of surface radius) = 1.00
l_max = 13
Dipole tilt (degrees) = -9.410531
```

and gives access to
variables associated with the planet such as:

 - `p.lmax` : maximum spherical harmonic degree till which data is available
 - `p.glm`, `p.hlm`: the Gauss coefficients
 - `p.Br` : computed radial magnetic field at surface
 - `p.dipTheta` : dipole tilt with respect to the rotation axis
 - `p.dipPhi` : dipole longitude ( in case zero longitude is known, applicable to Earth )

as well as the functions:

## `planet.plot()`

This function plots a 2D surface plot of the radial magnetic field at radius `r` given in terms of the surface radius.
For example,

```python
from planetmagfields import *
p = planet(name='earth',datDir='planetmagfields/data/')
p.plot(r=1,proj='moll')
```

produces the info mentioned above first and then the following plot of Earth's magnetic field using a Mollweide projection

<img src="planetmagfields/images/2d/earth2d.png" width="500">

while

```python
from planetmagfields import *
p = planet(name='jupiter',datDir='planetmagfields/data/')
p.plot(r=0.85,proj='moll')
```
produces the following info about Jupiter and then plot that follows

```
Planet: Jupiter
Depth (fraction of surface radius) = 1.00
l_max = 10
Dipole tilt (degrees) = 10.307870
```

<img src="planetmagfields/images/jupiter_r085.png" width="500">

This can be compared with Fig. 1 g from [Moore et al. 2018](https://doi.org/10.1038/s41586-018-0468-5)

## `planet.writeVtsFile`

This function writes a vts file that can be used to produce 3D visualizations of field lines with Paraview/VisIt. Usage:

```python
p.writeVtsFile(potExtra=True, ratio_out=2, nrout=32)
```
where,

  - `potExtra` : bool, whether to use potential extrapolation
  - `ratio_out`: float, radius till which the field would be extrapolated in terms of the surface radius
  - `nrout`: radial resolution for extrapolation

Example of a 3D image produced using Paraview for Neptune's field, extrapolated till 5 times the surface radius.

<img src="planetmagfields/images/3d/neptune3d.png" width="500">

# Field filtering using `planet.plot_filt`

The `planet` class also provides a function for producing a filtered view of the radial magnetic field using the function `plot_filt`.
This function can take in either an arbitrary array of spherical harmonic degrees and orders or cut-off values. This is illustrated
below with examples, assuming the user is in the repository directory.

## Saturn's small-scale magnetic field at a depth of 0.75 planetary radius ( degree > 3 )

```python
from planetmagfields import *
p = planet(name='saturn')
p.plot_filt(r=0.75,lCutMin=4,proj='moll')
```

<img src="planetmagfields/images/saturn_lgeq4_2d.png" width="500">

Compare this with Fig. 20 B from [Cao et al. 2020](https://doi.org/10.1016/j.icarus.2019.113541)

## Jupiter's surface field restricted to degrees 1,2,3 and order 3

```python
from planetmagfields import *
p = planet(name='jupiter')
p.plot_filt(r=1,larr=[1,2,3],marr=[3],proj='moll')
```

<img src="planetmagfields/images/jupiter_l123m3_2d.png" width="500">

## Earth's smaller scale surface field for degree > 4 and order > 3

```python
from planetmagfields import *
p = planet(name='earth')
p.plot_filt(r=1,lCutMin=5,mmin=4,proj='moll')
```

<img src="planetmagfields/images/earth_lgeq5mgeq4_2d.png" width="500">

# Potential extrapolation

The `potextra` module provides a method for potential extrapolation of a planet's magnetic field.
This uses the [SHTns](https://bitbucket.org/nschaeff/shtns) library for spherical harmonic transforms.
Usage example:

```python
import numpy as np
from planetmagfields import *
p = planet('saturn')
ratio_out = 5 # Ratio (> 1) in terms of surface radius to which to extrapolate
nrout = 32 # Number of grid points in radius between 1 and ratio_out
rout = np.linspace(1,ratio_out,nrout)
brout, btout, bpout = potextra.extrapot(p.lmax,1.,p.Br,rout)
```

# Quickplot using the `magField` script:

```
$ ./magField <planet> <radius> <projection>
```

This will plot the radial magnetic field of a planet (any of the names from the list
below, case insensitive) at a radius given in terms of the surface radius. The default
is the surface field. For example,

```
$ ./magField earth moll
```

displays the same information as above about Earth's field and produces the surface field of Earth while

```
$ ./magField jupiter 0.85 moll
```

produces the same plot of Jupiter's field as shown before.

```
$ ./magField all <radius> <projection>
```

would produce a table of information about dipole tilt for each planet and magnetic field maps of all different planets at the given radius in a single figure.

For example:

```
$ ./magField all 0.9 moll
```

would give

```
|=========|======|=======|
|Planet   | Theta| Phi   |
|=========|======|=======|
|Mercury  | 0.0  | 0.0   |
|Earth    | -9.4 | -72.7 |
|Jupiter  | 10.3 | -16.6 |
|Saturn   | 0.0  | 0.0   |
|Uranus   | 58.6 | -53.6 |
|Neptune  | 46.9 | -72.0 |
|Ganymede | -4.2 | 25.5  |
|---------|------|-------|
```

followed by the following plot

![All fields r=0.9](planetmagfields/images/magField_all_09.png)

# Spherical harmonic normalization and Condon-Shortley phase

All the Gauss coefficients in the collected data are Schmidt semi-normalized.
Only the data for Earth uses a Condon-Shortley phase, the others do not.

# Data sources

Mercury : [Anderson et al. 2012](https://doi.org/10.1029/2012JE004159)

Earth   : [IGRF 13](https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html)

Jupiter : [JRM09, Connerney et al. 2018](https://doi.org/10.1002/2018GL077312)

Saturn  : [Cassini 11+, Cao et al. 2020](https://doi.org/10.1016/j.icarus.2019.113541)

Uranus  : [Connerny et al. 1987](https://doi.org/10.1029/JA092iA13p15329)

Neptune : [Connerny et al. 1991](https://doi.org/10.1029/91JA01165)

Ganymede: [Kivelson et al. 2002](https://doi.org/10.1006/icar.2002.6834)

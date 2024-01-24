.. planetMagFields documentation master file, created by
   sphinx-quickstart on Mon Jan 22 08:32:05 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _secExamples:

Features and examples
=======================

The `Planet` class
+++++++++++++++++++

This gives access to all the relevant properties of a planet and has methods to plot
the field and write a `vts` file for 3D visualization. Usage:

.. code-block:: python

   from planetmagfields import *
   p = Planet(name='earth',datDir='planetmagfields/data/')


This displays the some information about the planet

.. code-block:: bash

   Planet: Earth
   l_max = 13
   Dipole tilt (degrees) = -9.410531


and gives access to
variables associated with the planet such as:

 - `p.lmax` : maximum spherical harmonic degree till which data is available
 - `p.glm`, `p.hlm`: the Gauss coefficients
 - `p.Br` : computed radial magnetic field at surface
 - `p.dipTheta` : dipole tilt with respect to the rotation axis
 - `p.dipPhi` : dipole longitude ( in case zero longitude is known, applicable to Earth )
 - `p.idx` : indices to get values of Gauss coefficients
 - `p.model` : the magnetic field model used. Available models can be obtained using the `get_models()` function. Selects the latest available model when unspecified.

Example using `IPython`:

.. code-block:: python

   In [1]: from planetmagfields import *

   In [2]: p = Planet(name='jupiter',model='jrm09')
   Planet: Jupiter
   Model: jrm09
   l_max = 10
   Dipole tilt (degrees) = 10.307870

   In [3]: p.glm[p.idx[2,0]]      # g20
   Out[3]: 11670.4

   In [4]: p.hlm[p.idx[4,2]]      # h42
   Out[4]: 27811.2

as well as the functions:


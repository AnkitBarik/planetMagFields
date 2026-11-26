.. _seccoords:

######################
Coordinate conventions
######################

The same conventions are used for every planet and every part of
`planetMagFields`: field maps, :py:meth:`Planet.orbit_path`, the
``--trajectory`` option of the :ref:`command line tool <seccli>`, potential
extrapolation, VTK output and the 3D web app.

Positions
=========

Positions are given in planet-centred spherical coordinates
:math:`(r, \theta, \phi)` in the body-fixed frame of the field model:

.. list-table::
   :header-rows: 1
   :widths: 12 88

   * - Symbol
     - Meaning
   * - :math:`r`
     - Distance from the planet's centre, in units of the reference radius of
       the field model (:math:`r = 1` is the reference surface). For example,
       6371.2 km for Earth (IGRF) and 71,492 km for Jupiter (JRM09, JRM33).
   * - :math:`\theta`
     - Co-latitude in radians: :math:`0` at the north pole, :math:`\pi/2` at the
       equator, :math:`\pi` at the south pole.
   * - :math:`\phi`
     - **East** longitude in radians, increasing eastward, in the body-fixed
       frame of the model. Any value is accepted; it is periodic in
       :math:`2\pi`.

Converting from latitude and longitude in degrees:

.. math::

   \theta = \frac{\pi}{180}\,(90^\circ - \mathrm{latitude}), \qquad
   \phi = \frac{\pi}{180}\,\lambda_E = \frac{\pi}{180}\,(360^\circ - \lambda_W),

where :math:`\lambda_E` is east longitude and :math:`\lambda_W` is west
longitude.

Field components
================

Field components are spherical components at each point:

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Component
     - Positive direction
   * - :math:`B_r`
     - Outward, away from the planet's centre
   * - :math:`B_\theta`
     - Towards increasing :math:`\theta`, i.e. **southward**
   * - :math:`B_\phi`
     - **Eastward**

Geomagnetic data are often given as north, east and down components
:math:`(X, Y, Z)`. These relate to the spherical components by

.. math::

   X = -B_\theta, \qquad Y = B_\phi, \qquad Z = -B_r .

For example, at Earth's surface at 0° N, 0° E in 2020,
:py:meth:`Planet.orbit_path` gives :math:`B_\theta \approx -27{,}600` nT (the
horizontal field points north) and :math:`B_\phi \approx -2{,}250` nT (a
westward declination of about 5°).

Cartesian coordinates
=====================

3D output (VTK files and 3D plots) uses a right-handed Cartesian frame fixed to
the planet:

* :math:`z` along the rotation axis, towards the north pole,
* :math:`x` in the equatorial plane towards :math:`\phi = 0`,
* :math:`y` in the equatorial plane towards :math:`\phi = 90^\circ` E,

so that

.. math::

   x = r\sin\theta\cos\phi, \qquad y = r\sin\theta\sin\phi, \qquad z = r\cos\theta .

Gauss coefficients
==================

For every planet the Gauss coefficients :math:`g_l^m, h_l^m` are for Schmidt
semi-normalised associated Legendre functions :math:`P_l^m` **without** the
Condon-Shortley phase, so that

.. math::

   B_r = \sum_{l=1}^{l_{max}} (l+1)\, r^{-(l+2)} \sum_{m=0}^{l}
   \left(g_l^m \cos m\phi + h_l^m \sin m\phi\right) P_l^m(\cos\theta),

with :math:`\phi` the east longitude above. See :ref:`secMath` for details.

Maps
====

On 2D maps, longitude is :math:`\phi` wrapped into
:math:`[-180^\circ, 180^\circ]`, with east to the right and the prime meridian
of the model's frame at the centre.

Body-fixed frames
=================

What "longitude" refers to depends on the planet:

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Planet
     - Frame
   * - Earth
     - Geographic coordinates. :math:`\phi` is geographic east longitude.
   * - Jupiter
     - System III (1965). It is conventionally quoted as west longitude
       :math:`\lambda_W`, so :math:`\phi = 360^\circ - \lambda_W`. For JRM33 this
       puts the dipole tilt towards 196.38° W and the Great Blue Spot near
       90° W, as in `Connerney et al. (2022)
       <https://doi.org/10.1029/2021JE007055>`_.
   * - Ganymede
     - Body-fixed, tidally locked frame.
   * - Uranus, Neptune
     - Longitude systems based on Voyager 2 rotation periods. Their
       uncertainty makes longitudes drift over decades, so they are best
       suited to comparisons with Voyager-era data.
   * - Saturn
     - The models are axisymmetric, so longitude does not matter.
   * - Mercury
     - IAU body-fixed frame. The models are nearly axisymmetric, so longitude
       matters little.

Example: field along a trajectory
=================================

Computing the field at points given as latitude and west longitude in System
III, 1.5 Jupiter radii from the centre:

.. code-block:: python

   import numpy as np
   from planetmagfields import Planet

   lat   = np.array([ 0.0, 30.0, -45.0])     # degrees
   lon_w = np.array([90.0, 196.4, 300.0])    # System III west longitude, degrees
   r     = np.full(3, 1.5)                   # Jupiter radii

   theta = np.deg2rad(90.0 - lat)
   phi   = np.deg2rad(360.0 - lon_w)         # east longitude

   p = Planet(name='jupiter', model='jrm33', units='nT', info=False)
   p.orbit_path(r, theta, phi)
   print(p.br_orb, p.btheta_orb, p.bphi_orb)  # outward, southward, eastward

For Earth, use geographic latitude and east longitude directly
(:math:`\phi = \lambda_E`).

.. note::

   Before version 1.9, maps labelled longitude as :math:`\phi - 180^\circ`.
   For planets other than Earth this rotated maps by 180°: Jupiter's Great
   Blue Spot, for example, appeared near 96° E instead of near 90° W. Field
   values along trajectories were correct for these planets.

   For Earth, the field was evaluated with the Condon-Shortley phase, which
   rotates it by 180° and so cancelled the map offset: maps were correct, but
   :py:meth:`Planet.orbit_path` and ``--trajectory`` needed
   :math:`\phi` = geographic longitude + 180° to give the right field.

   Now :math:`\phi` is east longitude everywhere, for maps and trajectories.

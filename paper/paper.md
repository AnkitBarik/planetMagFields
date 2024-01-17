---
title: >-
    planetMagFields: A Python package for analyzing and plotting planetary magnetic field data
authors:
  - name: Ankit Barik
    email: ankit.barik@gmail.com
    affiliation: [1]
    orcid: 0000-0001-5747-669X
    corresponding: true
  - name: Regupathi Angappan
    email: rangapp1@jhu.edu
    affiliation: [1]
    orcid: 0000-0002-6258-0659
affiliations:
  - index: 1
    name: Johns Hopkins University
date: 2024-01-04
bibliography: paper.bib
tags:
  - Python
  - planetary science
  - magnetic fields
  - visualization
---

# Summary
Long term observations and space missions have generated a wealth of data on the magnetic fields of the Earth and other solar system planets [@IGRF13;@Connerney1987;@Connerney1991;@Kivelson2002;@Anderson2012;@Cao2020;@Connerney2022]. `planetMagfields` is a Python package designed to have all the planetary magnetic field data currently available in one place and to provide an easy interface to access the data. `planetMagfields` focuses on planetary bodies that generate their own magnetic field, namely Mercury, Earth, Jupiter, Saturn, Uranus, Neptune and Ganymede. `planetMagfields` provides functions to compute as well as plot the magnetic field on the planetary surface or at a distance above or under the surface. It also provides functions to filter out the field to large or small scales as well as to produce `.vts` files to visualize the field in 3D using Paraview [@Ahrens2005;@Ayachit2015], VisIt [@visit] or similar rendering software.

# Statement of need



# Mathematics

Magnetic fields in planets are generated by electric currents in a fluid region inside them [@Jones2011;@Schubert2011;@Stanley2014]. Outside this region, in the absence of current sources, the magnetic field $\vec{B}$ can be written as the gradient of a scalar potential $V$,

$$\vec{B}(r,\theta,\phi) = -\nabla V = R_p \sum_{l,m} \left(\dfrac{R_p}{r}\right)^{l+1} [g_l^m \cos(m\phi) + h_l^m \sin(m\phi)] P_l^m (\cos\theta) \label{eq:gaussCoeff},$$

where, $g_l^m$ and $h_l^m$ are called the Gauss coefficients. Here, $(r,\theta,\phi)$ are spherical coordinates representing radial distance from the center of a planet, co-latitude and longitude, respectively. $R_p$ represents the radius of the planet and $P_l^m$ are associated Legendre functions of order $l$ and degree $m$, where $l$ and $m$ are integers. The Gauss coefficients represent the multipole modes of a planet's magnetic field. For example, $g_1^0$ represents the axial dipole (along the rotation axis) while $g_1^1$ and $h_1^1$ represent orthogonal components of the equatorial dipole. Thus, the dipole tilt of a planet, or the angle between the dipole and the rotation axis is given by:

$$\theta_{dip} = \tan^{-1}\dfrac{\sqrt{\left(g_1^1\right)^2 + \left(h_1^1\right)^2}}{g_1^0}$$

while the longitude of the dipole is given by:

$$\phi_{dip} = \tan^{-1}h_1^1/g_1^1$$

The raw data obtained from satellites or space missions are usually inverted to obtain these Gauss coefficients. These coefficients are the key to describing the surface magnetic field of a planet as well as how that magnetic field looks like at a certain altitude from the surface. The magnetic energy content on the surface in a certain degree $l$ is given by the Lowes spectrum:

$$R_{l} = (l + 1) \sum_{m}\left( \left(g_l^m\right)^2 + \left(h_l^m\right)^2\right),$$

$l$ plays the role of a wavenumber. Low degrees represent large spatial features in the field while high degrees represent small scale features. The maximum available degree $l_{max}$ of data for a particular planet depends on the quality of observations. For example, for Earth $l_{max} = 13$ because beyond that the magnetic field of magnetized rocks on the crust obscures any signal coming from the self generated field. Similarly, Jupiter's field was known only well constrained till $l_{max} = 4$ [@Connerney1998] before the Juno mission provided excellent observations of finer scale structure to extend the well constrained $l_{max}$ to 13 [@Connerney2022].

# Description of the software

`planetMagfields` has datafiles containing Gauss coefficients from various studies of planetary magnetic models for different planets. These coefficients are then used to obtain the magnetic field on a grid of latitude and longitude

# References


.. planetMagFields documentation master file, created by
   sphinx-quickstart on Mon Jan 22 08:32:05 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to planetMagFields' documentation!
===========================================

`planetMagFields` is a package that provides an easy interface to plot and analyze planetary magnetic field data. The raw data from space missions are typically inverted to obtain 'Gauss coefficients', given by:

.. math::
   :label: eqGauss

   V = R_p \sum \left(\frac{R_p}{r}\right)^{l+1} \left[g_{lm}\cos(m\phi) + h_{lm}\sin(m\phi)\right] ,

where :math:`V` is the magnetic scalar potential, such that the magnetic field is given by,

.. math::
   :label: eqBgradv

   \vec{B} = -\nabla V .

Here, :math:`g_{lm}` and :math:`h_{lm}` are called the Gauss coefficients, usually defined on a planetary surface. They are the key to extrpolating the magnetic field of a planet anywhere in regions that are free of electrical currents. All the Gauss coefficients in the collected data are `Schmidt semi-normalized <https://adsabs.harvard.edu/full/2005GeoJI.160..487W>`_. Only the data for Earth uses a `Condon-Shortley phase <https://mathworld.wolfram.com/Condon-ShortleyPhase.html>`_, the others do not. The prerequisites, installation and features are described in the pages below. However, for easy and quick visualization, skip ahead to the :ref:`Jupyter notebook <secjup>` section.

.. toctree::
   :maxdepth: 2
   :numbered:

   prereq
   installation
   examples
   models
   projections
   jupyter
   api



.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`

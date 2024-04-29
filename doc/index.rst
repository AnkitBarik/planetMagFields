.. planetMagFields documentation master file, created by
   sphinx-quickstart on Mon Jan 22 08:32:05 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to planetMagFields' documentation!
===========================================

`planetMagFields` is a package that provides an easy interface to plot and analyze planetary magnetic field data. `planetMagFields` is free and open source and is available at https://github.com/AnkitBarik/planetMagFields .

The raw data from space missions are typically inverted to obtain 'Gauss coefficients', given by:

.. math::
   :label: eqGauss

   V = R_p \sum \left(\frac{R_p}{r}\right)^{l+1} \left[g_{lm}\cos(m\phi) + h_{lm}\sin(m\phi)\right] ,

where :math:`V` is the magnetic scalar potential, such that the magnetic field is given by,

.. math::
   :label: eqBgradv

   \vec{B} = -\nabla V .

Here, :math:`g_{lm}` and :math:`h_{lm}` are called the Gauss coefficients, usually defined on a planetary surface. They are the key to extrpolating the magnetic field of a planet anywhere in regions that are free of electrical currents. All the Gauss coefficients in the collected data are `Schmidt semi-normalized <https://adsabs.harvard.edu/full/2005GeoJI.160..487W>`_. Only the data for Earth uses a `Condon-Shortley phase <https://mathworld.wolfram.com/Condon-ShortleyPhase.html>`_, the others do not. The prerequisites, installation and features are described in the pages below. However, for easy and quick visualization, skip ahead to the :ref:`Jupyter notebook <secjup>` section.

Code contribution and reporting issues
--------------------------------------

`planetMagFields` is an open source project and anyone is welcome to contribute to it. If you wish to contribute to this project, please follow the guidelines below:

 - Please make sure the tests pass before opening a pull request.
 - Please follow Python `PEP-8 <https://peps.python.org/pep-0008/>`_ guidelines when it comes to code style.
 - If you wish to add new data file, please name it following the convention `<planet>_<model>.dat`, where `planet` and `model` denote the names of the planet and the magnetic field model being used. See the `data` directory for examples.
 - If you implement a new feature or data source, please update the documentation accordingly.

Please report any bugs or other issues through `GitHub Issues <https://github.com/AnkitBarik/planetMagFields/issues>`_.

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

.. planetMagFields documentation master file, created by
   sphinx-quickstart on Mon Jan 22 08:32:05 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _secInstall:

#########################
Installation
#########################


Prerequisites
****************

Most basic functionalities of `planetMagFields` require only `NumPy <https://numpy.org/>`_, `Matplotlib <https://matplotlib.org/>`_ and `SciPy <https://www.scipy.org/>`_. Other than that, the following external libraries are used for a few different functions:

- 2D plotting for map projections other than Hammer : `Cartopy <https://scitools.org.uk/cartopy/docs/latest/>`_ library (see more under :ref:`Projections <secproj>`)
- Potential extrapolation (optional): `SHTns <https://bitbucket.org/nschaeff/shtns>`_ library for fast spherical harmonic transforms, else it falls back to SciPy.
- 3D visualisation: `PyVista <https://docs.pyvista.org/>`_ library
- Writing vts files for 3D visualisation: `PyEVTK <https://github.com/paulo-herrera/PyEVTK>`_ library

`planetMagFields` can be installed in a few different ways:

Using `pip`
***********

This is probably the easiest way. `planetMagFields` is available on `PyPI <https://pypi.org/project/planetMagFields/>`_ and can be installed with

.. code-block:: bash

   $ python3 -m pip install planetMagFields

Using `setup.py`
*****************

You can also use `setup.py` to install `planetMagFields`:

.. code-block:: bash

   $ git clone https://github.com/AnkitBarik/planetMagFields
   $ cd planetMagFields
   $ python3 setup.py install --user

Or using `pip`:

.. code-block:: bash

   $ git clone https://github.com/AnkitBarik/planetMagFields
   $ cd planetMagFields
   $ python3 -m pip install . --user

Using `PYTHONPATH`
******************

Download the package from the `GitHub repository <https://github.com/AnkitBarik/planetMagFields>`_ and it
to `PYTHONPATH`:

.. code-block:: bash

   $ git clone https://github.com/AnkitBarik/planetMagFields
   $ export PYTHONPATH=$PYTHONPATH:/path/to/planetMagFields



.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`

.. planetMagFields documentation master file, created by
   sphinx-quickstart on Mon Jan 22 08:32:05 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to planetMagFields' documentation!
===========================================

`planetMagFields` is a package that provides an easy interface to plot and analyze planetary magnetic field data. `planetMagFields` is free and open source and is available at https://github.com/AnkitBarik/planetMagFields .

Documentation
-------------

The prerequisites, installation and features are described in the pages below. However, for easy and quick visualization, skip ahead to the :ref:`Jupyter notebook <secjup>` section. The details of the math are explained in the :ref:`Mathematics <secMath>` section.

.. toctree::
   :maxdepth: 2
   :numbered:

   prereq
   installation
   examples
   math
   models
   projections
   jupyter
   api


Citing `planetMagFields`
------------------------

If you're using `planetMagFields` for your work, please cite the `JOSS paper <https://joss.theoj.org/papers/10.21105/joss.06677#>`_:

   Barik et al., (2024). planetMagFields: A Python package for analyzing and plotting planetary magnetic field data. Journal of Open Source Software, 9(97), 6677, https://doi.org/10.21105/joss.06677

.. code-block:: bibtex

   @article{Barik2024,
       doi = {10.21105/joss.06677},
       url = {https://doi.org/10.21105/joss.06677},
       year = {2024},
       publisher = {The Open Journal},
       volume = {9},
       number = {97},
       pages = {6677},
       author = {Barik, Ankit and Angappan, Regupathi},
       title = {planetMagFields: A Python package for analyzing and plotting planetary magnetic field data},
       journal = {Journal of Open Source Software}
   }

Code contribution and reporting issues
--------------------------------------

`planetMagFields` is an open source project and anyone is welcome to contribute to it. If you wish to contribute to this project, please follow the guidelines below:

 - Please make sure the tests pass before opening a pull request.
 - Please follow Python `PEP-8 <https://peps.python.org/pep-0008/>`_ guidelines when it comes to code style.
 - If you wish to add new data file, please name it following the convention `<planet>_<model>.dat`, where `planet` and `model` denote the names of the planet and the magnetic field model being used. See the `data` directory for examples.
 - If you implement a new feature or data source, please update the documentation accordingly.

Please report any bugs or other issues through `GitHub Issues <https://github.com/AnkitBarik/planetMagFields/issues>`_.


.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`

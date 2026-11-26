.. planetMagFields documentation master file, created by
   sphinx-quickstart on Mon Jan 22 08:32:05 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to planetMagFields' documentation!
===========================================

`planetMagFields` is a package that provides an easy interface to plot and analyze planetary magnetic field data. `planetMagFields` is free and open source and is available at https://github.com/AnkitBarik/planetMagFields .

.. raw:: html

   <div class="pmf-showcase">
     <a class="pmf-shot" href="https://ankitbarik.github.io/planetMagFields/app/" title="Open the 3D explorer">
       <img src="_static/images/webapp_explorer.jpg" width="2400" height="1350"
            alt="The planetMagFields 3D explorer showing Earth's magnetic field in 2025: the radial field on the surface with coastlines, 100 field lines, and the Lowes spectrum.">
     </a>
     <div class="pmf-cta">
       <a class="pmf-button" href="https://ankitbarik.github.io/planetMagFields/app/" target="_blank">Launch the 3D explorer</a>
       <span class="pmf-cta-note">Runs in your browser. Nothing to install.</span>
     </div>
   </div>

Explore planetary magnetic fields in 3D
---------------------------------------

The interactive `3D explorer <https://ankitbarik.github.io/planetMagFields/app/>`_ computes everything in your browser from the same Gauss coefficients as the package:

* **Every planet and model:** Mercury, Earth, Jupiter, Saturn, Uranus, Neptune and Ganymede, with all :ref:`supported field models <secmodels>`.
* **Field lines** traced live from the surface out to 10 planetary radii, with the surface field at any radius from 0.5 to 5.
* **Earth through time:** animate the IGRF field from 1900 to 2030.
* **Live filtering** by spherical harmonic degree and order, alongside the Lowes spectrum.
* **Spacecraft trajectories:** load a CSV of positions to see the path on the globe and plot :math:`B_r`, :math:`B_\theta`, :math:`B_\phi` and :math:`|B|` along it, then download the values. Your data stays on your computer. Built-in examples compare the models with real magnetometer data from Voyager 2 at Uranus, Juno at Jupiter and Cassini's Grand Finale at Saturn.
* **Shareable views:** the page address records the planet, model, radius and settings.

Documentation
-------------

The prerequisites, installation and features are described in the pages below. The details of the math are explained in the :ref:`Mathematics <secMath>` section.

.. toctree::
   :maxdepth: 2
   :numbered:

   installation
   examples
   cli
   math
   models
   coordinates
   projections
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
 - If you wish to add new data file, please name it following the convention `<planet>_<model>.dat`, where `planet` and `model` denote the names of the planet and the magnetic field model being used. See the `data` directory for examples. Then add an entry for the model under that planet's `"models"` in `planetmagfields/data/models.json` with its name and reference (`"name"`, `"url"`), plus `"format"` if the file layout differs from the planet's default.
 - If you implement a new feature or data source, please update the documentation accordingly.

Please report any bugs or other issues through `GitHub Issues <https://github.com/AnkitBarik/planetMagFields/issues>`_.


.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`

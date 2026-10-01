.. _seccli:

######################
Command Line Interface
######################

Once ``planetMagFields`` is installed, the ``magfield`` command is available directly
from the terminal without writing any Python code.

.. code-block:: bash

   magfield --help

Basic usage
-----------

Plot the surface magnetic field of Earth (default):

.. code-block:: bash

   magfield

Plot Jupiter's field at 85 % of its radius using a different colormap:

.. code-block:: bash

   magfield -p jupiter -r 0.85 -c viridis

Plot all planets at once:

.. code-block:: bash

   magfield -p all

Save a plot to a file instead of opening an interactive window:

.. code-block:: bash

   magfield -p saturn -s saturn_field.png

Compute field components along a spacecraft trajectory (outputs CSV to stdout):

.. code-block:: bash

   magfield -p saturn --trajectory orbit.csv

Same but save the result to a file:

.. code-block:: bash

   magfield -p saturn --trajectory orbit.csv -s orbit_field.csv

Full option reference
---------------------

Planet Selection
^^^^^^^^^^^^^^^^

.. option:: -p <name>, --planet <name>

   Planet to plot. Use ``all`` to plot every planet on one figure.
   (default: ``earth``)

.. option:: -o <model>, --model <model>

   Magnetic field model to use. Uses the latest available model when omitted.
   Available models can be listed with :py:func:`~planetmagfields.get_models`.

Data Options
^^^^^^^^^^^^

.. option:: -r <value>, --radius <value>

   Radial level at which to evaluate the field, scaled to the planetary radius.
   Must be positive. (default: ``1.0``)

.. option:: -u <unit>, --unit <unit>

   Unit for the magnetic field values. Accepted values: ``muT``, ``nT``, ``Gauss``.
   (default: ``muT``)

.. option:: --trajectory <file>

   Path to a CSV file with columns ``r``, ``theta``, ``phi`` (co-latitude and
   longitude both in **radians**). When supplied, the CLI switches to
   *trajectory mode*: it computes Br, Btheta and Bphi along the path using
   :py:meth:`Planet.orbit_path <planetmagfields.Planet.orbit_path>` instead of
   producing a map plot.

   Output is a CSV written to stdout, or to the file given by ``--save`` if
   that flag is also provided.

   .. note::

      ``orbit_path`` uses the `SHTns <https://bitbucket.org/nschaeff/shtns>`_
      library when available and falls back to SciPy otherwise. Install the
      optional ``potextra`` extra for best performance:
      ``pip install planetMagFields[potextra]``

   Example trajectory file (``orbit.csv``)::

      r,theta,phi
      17.829,1.1865,4.7087
      17.821,1.1863,4.7088
      17.813,1.1861,4.7090

Plot Options
^^^^^^^^^^^^

.. option:: --dark

   Use a dark background style for the figure.

.. option:: -c <name>, --cmap <name>

   Matplotlib colormap name. (default: ``RdBu_r``)

.. option:: -l <n>, --levels <n>

   Number of contour levels. (default: ``20``)

.. option:: -m <projection>, --mapproj <projection>

   Map projection. See :ref:`Projections <secproj>` for supported values.
   (default: ``Mollweide``)

.. option:: --vmin <value>

   Fix the lower bound of the colorscale.

.. option:: --vmax <value>

   Fix the upper bound of the colorscale.

Output Options
^^^^^^^^^^^^^^

.. option:: -s <file>, --save <file>

   Save the figure to *file* (e.g. ``field.png``) instead of opening an
   interactive window. The parent directory is created automatically if needed.

.. option:: -v, --verbose

   Print debug-level information during execution.

.. option:: --config <file>

   Path to a JSON configuration file (see :ref:`Config file <cliconfigfile>`
   below). Defaults to the ``config.json`` bundled with the package.

.. option:: --version

   Print the version number and exit.

.. _cliconfigfile:

Configuration file
------------------

All options can be set in a JSON file so you do not have to repeat long
command lines. Values in the config file act as defaults; any flag supplied on
the command line still takes precedence.

Pass a custom config file with ``--config``:

.. code-block:: bash

   magfield --config /path/to/myconfig.json

An example configuration file (matching the bundled ``planetmagfields/config.json``):

.. code-block:: json

   {
       "_comment": "Copy and edit this file, then pass it with --config",

       "planet": "earth",
       "model": null,

       "r": 1.0,
       "unit": "muT",

       "cmap": "RdBu_r",
       "levels": 20,
       "proj": "Mollweide",
       "vmin": null,
       "vmax": null,

       "dark": false,
       "save_path": null,
       "verbose": false,
       "trajectory": null
   }

Keys beginning with ``_`` are treated as comments and ignored. Unknown keys
produce a warning. ``null`` values leave the corresponding option at its
built-in default.


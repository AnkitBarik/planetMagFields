#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import importlib
import importlib.util
import numpy as np
from .models import planetlist

stdDatDir = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data/')

UNITS = {
    'nt':    (1.,   'nT'),
    'mut':   (1e-3, r'$\mu$T'),
    'gauss': (1e-5, 'Gauss'),
}

UNIT_TEXT = {'nt': 'nT', 'mut': 'μT', 'gauss': 'G'}

INSTALL_HINTS = {
    'shtns':   'https://bitbucket.org/nschaeff/shtns',
    'pyevtk':  'pip install pyevtk',
    'pyvista': 'pip install pyvista',
    'cartopy': 'pip install cartopy',
}


def get_unit(units):
    """Returns (factor, label) converting nT to the requested units."""
    try:
        return UNITS[units.lower()]
    except KeyError:
        raise ValueError("Unknown units '%s', must be one of 'nT', 'muT' or 'Gauss'"
                         % units) from None


def unit_text(units):
    """Plain text unit label, e.g. for plotly figures."""
    get_unit(units)
    return UNIT_TEXT[units.lower()]


def has_module(name):
    return importlib.util.find_spec(name) is not None


def require(name):
    """Imports an optional dependency, raising an ImportError with an
    installation hint if it is missing."""
    try:
        return importlib.import_module(name)
    except ImportError as exc:
        pkg = name.split('.')[0]
        raise ImportError("This requires the %s library: %s"
                          % (pkg, INSTALL_HINTS.get(pkg, 'pip install ' + pkg))) from exc


def get_models(planetname,datDir=stdDatDir):
    """Prints available models for a planet.

    Parameters
    ----------
    datDir : str
        Directory where the data file is present. Files are assumed to be named
        as <planetname>_<modelname>.dat,
        e.g.: earth_igrf13.dat, jupiter_jrm09.dat etc.
    planetname : str
        Name of the planet

    Returns
    -------
    models : str array
        Array of available model names
    """

    from glob import glob
    planetname = planetname.lower()
    dataFiles = glob(os.path.join(datDir, planetname+"_*.dat"))
    models = [os.path.basename(f)[len(planetname)+1:-len('.dat')] for f in dataFiles]
    return np.sort(models)

def is_dark_color(color):
    """
    Determine if a color is dark using relative luminance.

    Parameters
    ----------
    color : str or tuple
        Color as hex string ('#RRGGBB'), RGB tuple (r, g, b) with values 0-255,
        or normalized RGB tuple (r, g, b) with values 0-1.

    Returns
    -------
    bool
        True if color is dark, False if light.
    """

    # Parse color to RGB values (0-255)
    if isinstance(color, str):
        color = color.lstrip('#')
        if len(color) == 6:
            r, g, b = tuple(int(color[i:i+2], 16) for i in (0, 2, 4))
        elif len(color) == 3:
            r, g, b = tuple(int(c*2, 16) for c in color)
        else:
            # Named colors
            import matplotlib.colors as mcolors
            rgb = mcolors.to_rgb(color)
            r, g, b = int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255)
    elif isinstance(color, (tuple, list)):
        if all(0 <= c <= 1 for c in color):
            # Normalized RGB
            r, g, b = int(color[0]*255), int(color[1]*255), int(color[2]*255)
        else:
            # 0-255 RGB
            r, g, b = color[0], color[1], color[2]
    else:
        raise ValueError(f"Unknown color format: {color}")

    # Calculate relative luminance (ITU-R BT.709)
    # Human eye is most sensitive to green, least to blue
    luminance = (0.2126 * r + 0.7152 * g + 0.0722 * b) / 255

    # Threshold of 0.5 is standard; can adjust based on preference
    return luminance < 0.5
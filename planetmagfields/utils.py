#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np

stdDatDir = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data/')
planetlist = ["mercury", "earth", "jupiter", "saturn", "uranus", "neptune",
              "ganymede"]

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
    dataFiles = glob(datDir+'/'+planetname+"*.dat")
    models = []
    for k,filename in enumerate(dataFiles):
        modelname = filename.split('_')[1].split('.dat')[0]
        models.append(modelname)
    models = np.sort(models)
    return models

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
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
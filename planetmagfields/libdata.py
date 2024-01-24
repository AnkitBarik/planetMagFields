#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import os
from .libgauss import gen_idx
from .utils import stdDatDir

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
    dataFiles = glob(datDir+'/'+planetname+"*.dat")
    models = []
    for k,filename in enumerate(dataFiles):
        modelname = filename.split('_')[1].split('.dat')[0]
        models.append(modelname)
    models = np.sort(models)
    return models

def get_data(datDir,planetname="earth",model=None,year=2020):
    """
    Reads data file for a planet and rearranges to create arrays of Gauss coefficients,
    glm and hlm.

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
    glm : array_like
        Coefficients of real part of spherical harmonics (often called glm in
        literature)
    hlm : array_like
        Coefficients of imaginary part of spherical harmonics (often called hlm in
        literature)
    lmax : int
        Maximum spherical harmonic degree
    idx : int array
        Array of indices that correspond to an (l,m) combination. For example,
        g(0,0) -> 0, g(1,0) -> 1, g(1,1) -> 2 etc.
    """

    models_avail = get_models(datDir,planetname)
    datfile = datDir + planetname + '_' + model.lower() + '.dat'
    try:
        tmpdat = np.loadtxt(datfile,dtype=object)
        del tmpdat
    except FileNotFoundError:
        print("Could not read datafile, please double check path, planet and model!")
        print("For selected planet %s, the following models are available:" %planetname)
        print(models_avail)
        print("Use get_models to get a list of models!")

    m0file = ( (planetname == "mercury" and model == "anderson2012")
            or (planetname == "saturn") )

    if planetname == "jupiter" and model in ['jrm09','jrm33']:
        dat = np.loadtxt(datfile,usecols=[1,-2])
        l_dat = dat[:,1]
        ghlm = dat[:,0]

        lmax = np.int32(l_dat.max())
        if model == 'jrm33':
            lmax = 18

        g = []
        h = []

        ########################
        # Separate glm and hlm
        ########################

        for i in range(1,lmax+1):
            mask = l_dat == i
            n = len(l_dat[mask])
            half = int(n/2)

            g.append(ghlm[mask][:half+1])
            h.append(np.concatenate([[0.],ghlm[mask][half+1:]]))

        g  = np.concatenate(g)
        h  = np.concatenate(h)

    elif m0file:

        dat = np.loadtxt(datfile,usecols=[3])
        g   = dat.flatten()
        lmax = len(g)
        h = np.zeros_like(g)

    else:
        dat = np.loadtxt(datfile,dtype=object)
        gh  = dat[:,0]
        dat = np.float32(dat[:,1:])
        lmax = np.int32(dat[:,0]).max()

        if planetname != "earth":

            mask = gh == 'g'
            gDat = dat[mask,:]
            g   = gDat[:,-1]

            mask = gh == 'h'
            hDat = dat[mask,:]
            h   = hDat[:,-1]

        else: # Earth, accounting for linear SV

            if year > 2025:
                print("IGRF is only defined till 2025,please be careful while selecting year!")

            years = 1900 + 5*np.arange(dat.shape[1])
            idx = np.argmin(np.abs(years - year))

            if years[idx] < year:
                selected_idx = idx
            else:
                selected_idx = idx-1

            if years[selected_idx] == 2020:
                sv = dat[:,-1]
            else:
                sv = (dat[:,selected_idx+1] - dat[:,selected_idx])/5
            selected_dat = ( dat[:,selected_idx] +
                            sv * (year - years[selected_idx]) )

            mask = gh == 'g'
            g = selected_dat[mask]
            mask = gh == 'h'
            h = selected_dat[mask]

            m = dat[mask,1]
            m1Idx = np.where(m == 1.)[0]
            h   = np.insert(h,m1Idx,0.)



    # Insert (0,0) -> 0 for less confusion

    glm = np.insert(g,0,0.)
    hlm = np.insert(h,0,0.)

    if m0file:
        idx = np.arange(lmax+1)
    else:
        idx = gen_idx(lmax)

    return glm,hlm,lmax,idx

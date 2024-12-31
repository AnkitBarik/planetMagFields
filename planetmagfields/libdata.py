#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import os
from .libgauss import gen_idx
from .utils import get_models

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
    mmax : int
        Maximum spherical harmonic order. This is required to distinguish cases
        with maximum order of zero.
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

    if ( (planetname == "mercury" and model == "anderson2012")
        or (planetname == "saturn") ):
        mmax = 0
    else:
        mmax = None

    if planetname == "jupiter" and model in ['jrm09','jrm33']:
        dat = np.loadtxt(datfile,dtype=object)
        gh = dat[:,3]
        l_dat = np.int32(dat[:,-2])
        ghlm = np.float32(dat[:,1])

        lmax = l_dat.max()
        if model == 'jrm33':
            lmax = 18

        gmask = gh == 'g'
        hmask = gh == 'h'
        g = ghlm[gmask]
        h = ghlm[hmask]

        m = np.int32(dat[hmask,-1])
        m1Idx = np.where(m == 1.)[0]
        h   = np.insert(h,m1Idx,0.)

    elif mmax == 0:

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

            m = dat[mask,1]
            m1Idx = np.where(m == 1.)[0]
            h   = np.insert(h,m1Idx,0.)

        else: # Earth, accounting for linear SV

            if year > 2030:
                print("IGRF-14 is only defined till 2030,please be careful while selecting year!")

            years = 1900 + 5*np.arange(dat.shape[1])

            from scipy import interpolate
            f = interpolate.interp1d(years,dat,fill_value='extrapolate')
            selected_dat = f(year)

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

    # Ensure type int
    lmax = int(lmax)

    if mmax == 0:
        idx = np.arange(lmax+1)
    else:
        idx = gen_idx(lmax)
        mmax = lmax

    return glm,hlm,lmax,idx,mmax

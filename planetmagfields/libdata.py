#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import warnings
import numpy as np
from .libgauss import gen_idx
from .models import model_format, model_info
from .utils import get_models


def _insert_hm0(h, m):
    """Inserts h(l,0) = 0 before every h(l,1) entry."""
    return np.insert(h, np.where(m == 1.)[0], 0.)


def _read_jrm(datfile):
    dat = np.loadtxt(datfile,dtype=object)
    gh = dat[:,3]
    ghlm = np.float32(dat[:,1])
    gmask = gh == 'g'
    hmask = gh == 'h'
    g = ghlm[gmask]
    h = _insert_hm0(ghlm[hmask], np.int32(dat[hmask,-1]))
    lmax = np.int32(dat[:,-2]).max()
    return g, h, lmax


def _read_axisymmetric(datfile):
    g = np.loadtxt(datfile,usecols=[3]).flatten()
    return g, np.zeros_like(g), len(g)


def _split_gh(gh, dat, vals):
    gmask = gh == 'g'
    hmask = gh == 'h'
    return vals[gmask], _insert_hm0(vals[hmask], dat[hmask,1])


def _read_generic(datfile):
    dat = np.loadtxt(datfile,dtype=object)
    gh  = dat[:,0]
    dat = np.float32(dat[:,1:])
    g, h = _split_gh(gh, dat, dat[:,-1])
    return g, h, np.int32(dat[:,0]).max()


def _read_igrf(datfile, year):
    if year < 1900:
        warnings.warn("IGRF-14 is only defined from 1900, please be careful while selecting year!")
    elif year > 2030:
        warnings.warn("IGRF-14 is only defined till 2030, please be careful while selecting year!")

    dat = np.loadtxt(datfile,dtype=object)
    gh  = dat[:,0]
    dat = np.float32(dat[:,1:])

    # Columns 0,1 are l,m; columns 2:-1 are year data; last column is secular variation
    year_data = dat[:, 2:-1]
    years = 1900 + 5 * np.arange(year_data.shape[1])

    from scipy import interpolate
    f = interpolate.interp1d(years, year_data, fill_value='extrapolate')
    g, h = _split_gh(gh, dat, f(year))
    return g, h, np.int32(dat[:,0]).max()


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
    model : str
        Name of the model
    year : float
        Year for time dependent models (Earth)

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
    mmax : int
        Maximum spherical harmonic order. This is required to distinguish cases
        with maximum order of zero.
    """

    model = model.lower()
    datfile = os.path.join(datDir, planetname + '_' + model + '.dat')
    if not os.path.exists(datfile):
        raise FileNotFoundError(
            "Could not find %s. Models available for %s: %s"
            % (datfile, planetname, [str(m) for m in get_models(planetname, datDir)]))

    fmt = model_format(planetname, model)
    if fmt == 'jrm':
        g, h, lmax = _read_jrm(datfile)
    elif fmt == 'axisymmetric':
        g, h, lmax = _read_axisymmetric(datfile)
    elif fmt == 'igrf':
        g, h, lmax = _read_igrf(datfile, year)
    else:
        g, h, lmax = _read_generic(datfile)

    lmax = int(model_info(model).get('lmax', lmax))

    # Insert (0,0) -> 0 for less confusion
    glm = np.insert(g,0,0.)
    hlm = np.insert(h,0,0.)

    idx = gen_idx(lmax)

    if fmt == 'axisymmetric':
        # Expand to the full triangular size so idx[l,0] is always in bounds.
        mmax = 0
        ncoeff = (lmax+1)*(lmax+2)//2
        glm_full = np.zeros(ncoeff)
        hlm_full = np.zeros(ncoeff)
        glm_full[idx[:,0]] = glm[:lmax+1]
        hlm_full[idx[:,0]] = hlm[:lmax+1]
        glm, hlm = glm_full, hlm_full
    else:
        mmax = lmax

    return glm,hlm,lmax,idx,mmax

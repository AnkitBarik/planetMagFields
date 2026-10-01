#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.special import sph_harm_y
from scipy.special import roots_legendre
from .models import phase_factor

def gen_idx(lmax):

    """
    Generate index to convert from (l,m) to [idx]

    Parameters
    ----------
    lmax : int
        This defines the maximum spherical harmonic degree

    Returns
    -------
    idx : int
        Array index of the spherical harmonic (l,m)
    """

    idx = np.zeros([lmax+1,lmax+1])
    count = 0

    for l in range(lmax+1):
        for m in range(l+1):

            idx[l,m] = count
            count += 1

    return np.int32(idx)

def get_grid(nphi=256,ntheta=128):
    """
    Generates 2D grid of longitude and co-latitude. The longitude grid is equally
    spaced, the co-latitude grid uses the zeros of Legendre polynomial with degree
    ntheta.

    Parameters
    ----------
    nphi : int, optional
        Number of points in longitude, by default 256
    ntheta : int, optional
        Number of poitns in co-latitude, by default 128

    Returns
    -------
    p2D : ndarray(float, ndim=2)
        Longitude at every point on a (longitude,co-latitude) grid
    th2D : ndarray(float, ndim=2)
        Co-latitude at every point on a (longitude,co-latitude) grid
    phi : ndarray(float, ndim=1)
        Array of longitude values
    theta : ndarray(float, ndim=1)
        Array of co-latitude values
    """

    phi    = np.linspace(0.,2*np.pi,nphi)
    x,w    = roots_legendre(ntheta)
    theta  = np.sort(np.arccos(x))

    p2D, th2D = np.meshgrid(phi,theta,indexing='ij')

    return p2D, th2D, phi, theta

def getB(lmax,mmax,glm,hlm,idx,r,p2D,th2D,planetname="earth"):
    """
    This function computes the radial magnetic field from arrays of Gauss
    coefficients glm and hlm. It uses scipy's sph_harm function for spherical
    harmonics.

    Parameters
    ----------
    lmax : int
        Maximum degree of spherical harmonic
    lmax : int
        Maximum order of spherical harmonic
    glm : array_like
        Gauss coefficients of cos(m*phi)
    hlm : array_like
        Gauss coefficients of sin(m*phi)
    idx : int array
        Array of indices that correspond to an (l,m) combination. For example,
        g(0,0) -> 0, g(1,0) -> 1, g(1,1) -> 2 etc.
    r : float
        Radial level scaled to planetary surface
    p2D : ndarray(float, ndim=2)
        2D array defining longitude (phi).
        This ranges from 0 to 2*pi and has a shape (nphi,ntheta)
    th2D : ndarray(float, ndim=2)
        2D array defining co-latitude (theta).
        This ranges from 0 to pi and has a shape (nphi,ntheta)
    planetname : str, optional
        Name of the planet, by default "earth"

    Returns
    -------
    Br : ndarray(float, ndim=2)
        Radial magnetic field on the grid defined by p2D and th2D
    """

    Br = np.zeros_like(p2D)

    if mmax > 0:
        for l in range(1,lmax+1):
            for m in range(l+1):
                ylm = sph_harm_y(l, m, th2D, p2D)

                # Scipy includes the Condon-Shortley phase, which only Earth models use
                fac_m = phase_factor(planetname, m)

                if m != 0:
                    fac_m *= np.sqrt(2)

                fac = fac_m * (l+1) * r**(-l-2) * np.sqrt((4.*np.pi)/(2*l+1))

                G = glm[idx[l,m]] * np.real(ylm)
                H = hlm[idx[l,m]] * np.imag(ylm)

                Br +=   np.real(fac * (G + H))
    else:
        for l in range(1,lmax+1):
            ylm = sph_harm_y(l, 0, th2D, p2D)
            fac = (l + 1) * r**(-l-2) * np.sqrt((4.*np.pi)/(2*l+1))

            Br += fac * glm[idx[l,0]] * np.real(ylm)

    return Br

def get_spec(glm,hlm,idx,lmax,mmax,r=1.0):
    """
    Computes Lowes spectrum of a planet with Gauss coefficients glm and hlm at
    a radial level r, scaled to the planetary radius.

    Parameters
    ----------
    glm : array_like
        Gauss coefficients of cos(m*phi)
    hlm : array_like
        Gauss coefficients of sin(m*phi)
    idx : int array
        Array of indices that correspond to an (l,m) combination. For example,
        g(0,0) -> 0, g(1,0) -> 1, g(1,1) -> 2 etc.
    lmax : int
        Maximum degree of spherical harmonic
    mmax : int
        Maximum order of spherical harmonic
    r : float, optional
        Radial level scaled to planetary surface, by default 1

    Returns
    -------
    E : array_like
        Magnetic energy in spherical harmonic degrees
    emag_10 : float
        Magnetic energy in the axial dipole
    E_symm : float
        Equatorially symmetric magnetic energy
    E_antisymm : float
        Equatorially anti-symmetric magnetic energy
    E_axisymm : float
        Axisymmetric magnetic energy
    """

    E = np.zeros(lmax+1)
    E_symm = 0.
    E_antisymm = 0.
    E_axisymm = 0.

    if mmax == 0:
        for l in range(1,lmax+1):
            E[l] = (l+1) * r**(-2*l-4) *(np.abs(glm[idx[l,0]])**2 + np.abs(hlm[idx[l,0]])**2)
        emag_10 = E[1]
        E_symm = np.sum(E[::2])
        E_antisymm = np.sum(E[1::2])
        E_axisymm = np.sum(E)
    else:
        for l in range(1,lmax+1):
            for m in range(l+1):
                emag_lm = (l+1) * r**(-2*l-4) *(np.abs(glm[idx[l,m]])**2 + np.abs(hlm[idx[l,m]])**2)
                E[l] += emag_lm
                if ( (l-m)%2 == 0 ):
                    E_symm += emag_lm
                else:
                    E_antisymm += emag_lm
            E_axisymm += (l+1) * r**(-2*l-4) *(np.abs(glm[idx[l,0]])**2 + np.abs(hlm[idx[l,0]])**2)

        emag_10 = 2 * r**(-2*1-4)* np.abs(glm[idx[1,0]])**2
    return E, emag_10, E_symm, E_antisymm, E_axisymm

def filt_Gauss(glm,hlm,lmax,model_mmax,idx,larr=None,marr=None,
               lCutMin=0,lCutMax=None,mmin=0,mmax=None):
    """
    Filters Gauss coefficients by using either a fixed array of spherical harmonic
    degrees or orders or a minimum or maximum degree or order. Coefficients are either
    restricted to degree and order values defined by larr and marr or range defined
    by lCutMin, lCutMax and mmin, mmax, respectively.

    Parameters
    ----------
    glm : array_like
        Gauss coefficients of cos(m*phi)
    hlm : array_like
        Gauss coefficients of sin(m*phi)
    idx : int array
        Array of indices that correspond to an (l,m) combination. For example,
        g(0,0) -> 0, g(1,0) -> 1, g(1,1) -> 2 etc.
    lmax : int
        Maximum degree of spherical harmonic
    model_mmax : int
        Maximum order of spherical harmonic of the field model
    larr : int array, optional
        Array of desired spherical harmonic degrees, by default None
    marr : int array, optional
        Array of desired spherical harmonic orders, by default None
    lCutMin : int, optional
        Minimum spherical harmonic degree to retain, by default 0
    lCutMax : int, optional
        Maximum spherical harmonic degree to retain, by default None
        If None, then lmax is used
    mmin : int, optional
        Minimum spherical harmonic order to retain, by default 0
    mmax : int, optional
        Maximum spherical harmonic order to retain, by default None
        If None, lmax is used

    Returns
    -------
    glm_filt : array_like
        Array of filtered Gauss coefficients of cos(m*phi)
    hlm_filt : array_like
        Array of filtered Gauss coefficients of sin(m*phi)
    """

    glm_filt = np.copy(glm)
    hlm_filt = np.copy(hlm)

    if lCutMax is None:
        lCutMax = lmax
    if mmax is None:
        mmax = model_mmax

    if larr is not None:
        if max(larr) > lmax:
            raise ValueError("Values in filter array must be <= lmax = %d" %lmax)
        lkeep = np.isin(np.arange(lmax+1), larr)
    else:
        if lCutMax > lmax or lCutMin > lmax:
            raise ValueError("lCutMin/lCutMax must be <= lmax = %d" %lmax)
        lkeep = (np.arange(lmax+1) >= lCutMin) & (np.arange(lmax+1) <= lCutMax)

    if model_mmax > 0:
        if marr is not None:
            if max(marr) > lmax:
                raise ValueError("Values in filter array must be <= lmax = %d" %lmax)
            mkeep = np.isin(np.arange(lmax+1), marr)
        else:
            if mmin > lmax or mmax > lmax:
                raise ValueError("mmin/mmax must be <= lmax = %d" %lmax)
            mkeep = (np.arange(lmax+1) >= mmin) & (np.arange(lmax+1) <= mmax)
    else:
        mkeep = np.ones(lmax+1, dtype=bool)

    for l in range(lmax+1):
        for m in range(l+1):
            if not (lkeep[l] and mkeep[m]):
                glm_filt[idx[l,m]] = 0.
                hlm_filt[idx[l,m]] = 0.

    return glm_filt,hlm_filt


def sphInt(f,g,phi,th2D,theta):
    """
    Utility function for integrating a product of 2D arrays defined
    on a spherical surface.

    Parameters
    ----------
    f : ndarray(float, ndim=2)
        Defined on a (longitude,co-latitude) grid
    g : ndarray(float, ndim=2)
        Defined on a (longitude,co-latitude) grid
    phi : array_like
        Longitude
    th2D : ndarray(float, ndim=2)
        Co-latitude defined on (longitude,co-latitude) grid
    theta : array_like
        Co-latitude

    Returns
    -------
    phiInt : float
        f*g integrated over the spherical surface defined by (phi, theta)
    """
    from scipy.integrate import simpson

    thetaInt = simpson(f * g * np.sin(th2D),theta,axis=1)
    phiInt = simpson(thetaInt,phi)

    return phiInt


def getGauss(lmax,Br,r,phi,theta,th2D,p2D,planetname="earth"):
    """
    Get Gauss coefficients from a surface field defined on the grid produced
    by get_grid. This is the inverse of getB.

    Parameters
    ----------
    lmax : int
        Maximum degree of spherical harmonic
    Br : ndarray(float, ndim=2)
        Radial magnetic field defined on (phi,theta) grid
    r : float
        Radial level (planetary radius = 1)
    phi : array_like
        Longitude
    theta : array_like
        Co-latitude
    th2D : ndarray(float, ndim=2)
        Co-latitude defined on (longitude,co-latitude) grid
    p2D : ndarray(float, ndim=2)
        Longitude defined on (longitude,co-latitude) grid
    planetname : str, optional
        Name of the planet, by default "earth"

    Returns
    -------
    glm : array_like
        Gauss coefficients of cos(m*phi)
    hlm : array_like
        Gauss coefficients of sin(m*phi)
    """
    _, wtheta = roots_legendre(len(theta))
    wphi = np.full(len(phi), phi[1] - phi[0])
    wphi[[0,-1]] *= 0.5
    weights = wphi[:,None] * wtheta[None,:]

    glm = []
    hlm = []

    for l in range(lmax+1):
        for m in range(l+1):
            if l == 0:
                glm.append(0.)
                hlm.append(0.)
                continue

            ylm = sph_harm_y(l, m, th2D, p2D)
            fac = phase_factor(planetname, m) * (l+1) * r**(-l-2) * np.sqrt((4.*np.pi)/(2*l+1))
            if m != 0:
                fac *= np.sqrt(2) * 0.5

            glm.append(np.sum(weights * Br * np.real(ylm)) / fac)
            hlm.append(np.sum(weights * Br * np.imag(ylm)) / fac)

    return np.array(glm), np.array(hlm)


def get_dipole_tilt(glm,hlm,idx,mmax):
    """
    Computes dipole tilt co-latitude and longitude in degrees.

    Returns
    -------
    dipTheta : float
        Dipole tilt co-latitude in degrees
    dipPhi : float
        Dipole longitude in degrees
    """
    if mmax == 0:
        return 0, 0

    g10, g11, h11 = glm[idx[1,0]], glm[idx[1,1]], hlm[idx[1,1]]
    dipTheta = np.arctan(np.sqrt(g11**2 + h11**2)/g10) * 180./np.pi
    dipPhi = np.arctan(h11/g11) * 180./np.pi
    return dipTheta, dipPhi

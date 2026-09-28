#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from .libgauss import get_grid
from .models import phase_factor
from .utils import has_module, require


def _use_shtns():
    if has_module('shtns'):
        return True
    print("SHTns library not found, falling back to scipy!")
    print("Consider installing SHTns from https://bitbucket.org/nschaeff/shtns")
    return False


def _check_shapes(r, theta, phi):
    if not (np.shape(r) == np.shape(theta) == np.shape(phi)):
        raise ValueError("Please make sure all three arrays are of the same shape")
    return (np.asarray(r, dtype=np.float64),
            np.asarray(theta, dtype=np.float64),
            np.asarray(phi, dtype=np.float64))


def _potential_field_scipy(glm, hlm, idx, lmax, mmax, r, theta, phi, planetname=None):
    """Evaluates the potential field at broadcastable arrays r, theta, phi."""
    from scipy.special import sph_harm_y

    br = 0j
    bt = 0j
    bp = 0j

    for l in range(1,lmax+1):
        fac_l = (1/r)**(l+2)
        for m in range(min(mmax+1,l+1)):
            ylm, dylm = sph_harm_y(l,m,theta,phi,diff_n=1)
            if m == 0:
                Nlm = np.sqrt(4 * np.pi / (2*l+1))
            else:
                Nlm = np.sqrt(8 * np.pi / (2*l+1)) * phase_factor(planetname, m)

            coeff = fac_l * (glm[idx[l,m]] - 1j*hlm[idx[l,m]]) * Nlm

            br = br + coeff * (l+1) * ylm
            bt = bt - coeff * dylm[...,0]
            bp = bp - coeff * dylm[...,1] / np.sin(theta)

    return np.real(br), np.real(bt), np.real(bp)


def get_pol_from_Gauss(planetname,glm,hlm,lmax,mmax,idx):

    bpol = np.zeros(len(glm),dtype=np.complex128)

    for l in range(1,lmax+1):
        for m in range(min(mmax+1,l+1)):

            fac_m = phase_factor(planetname, m)

            if m == 0:
                norm = np.sqrt((4*np.pi/(2*l+1)))/l
            else:
                norm = np.sqrt((2*np.pi/(2*l+1)))/l

            bpol[idx[l,m]] = norm*fac_m*(glm[idx[l,m]] - 1j*hlm[idx[l,m]])

    return bpol


def extrapot_scipy(glm, hlm, idx, lmax, mmax, rplanet, rout, nphi=None, planetname=None):

    """
    This function extrapolates a potential field to an array of desired radial
    levels. It uses the scipy library for computing spherical harmonics.

    Parameters
    ----------
    glm : ndarray(float, ndim=1)
        Array of Gauss coefficients g_lm
    hlm : ndarray(float, ndim=1)
        Array of Gauss coefficients h_lm
    idx : ndarray(int, ndim=1)
        Array of indices to map [l,m] to an index
    lmax : int
        Maximum spherical harmonic degree of field model
    mmax : int
        Maximum spherical harmonic order of field model
    rplanet : float
        Radius of the planet at which the magnetic field is measured or defined. This is used to scale the field correctly during extrapolation.
    rout : array_like
        Array of radial levels to which the field should be extrapolated
    nphi : int, optional
        Number of grid points in longitude, can be automatically
        selected, by default None
    planetname : str, optional
        Name of the planet, determines the phase convention, by default None

    Returns
    -------
    brout : ndarray(float, ndim=3)
        3D array of extrapolated radial magnetic field, shape : (nphi,ntheta,nr)
    btout : ndarray(float, ndim=3)
        3D array of extrapolated co-latitudinal magnetic field, shape : (nphi,ntheta,nr)
    bpout : ndarray(float, ndim=3)
        3D array of extrapolated azimuthal magnetic field, shape : (nphi,ntheta,nr)
    """

    if nphi is None:
        nphi   = int(max(256,lmax*3))

    phi2d, theta2d, _, _ = get_grid(nphi,nphi//2)
    rout = np.atleast_1d(rout)/rplanet

    brout, btout, bpout = _potential_field_scipy(glm, hlm, idx, lmax, mmax,
                                                 rout[None,None,:],
                                                 theta2d[...,None], phi2d[...,None],
                                                 planetname=planetname)
    shape = phi2d.shape + (len(rout),)
    return (np.broadcast_to(brout, shape).copy(),
            np.broadcast_to(btout, shape).copy(),
            np.broadcast_to(bpout, shape).copy())

def extrapot_shtns(bpol,idx,lmax,mmax,rplanet,rout,nphi=None):
    """
    This function extrapolates a potential field to an array of desired radial
    levels. It uses the SHTns library (https://bitbucket.org/nschaeff/shtns)
    for spherical harmonic transforms.

    Parameters
    ----------
    bpol : ndarray(complex128, ndim=1)
        Array of poloidal coefficients computed from Gauss coefficients
    idx : ndarray(int, ndim=1)
        Array of indices to map [l,m] to an index
    lmax : int
        Maximum spherical harmonic degree of field model
    mmax : int
        Maximum spherical harmonic order of field model
    rplanet : float
        Radius at which the magnetic field is measured or defined
    rout : array_like
        Array of radial levels to which the field should be extrapolated
    nphi : int, optional
        Number of grid points in longitude, can be automatically
        selected, by default None

    Returns
    -------
    brout : ndarray(float, ndim=3)
        3D array of extrapolated radial magnetic field, shape : (nphi,ntheta,nr)
    btout : ndarray(float, ndim=3)
        3D array of extrapolated co-latitudinal magnetic field, shape : (nphi,ntheta,nr)
    bpout : ndarray(float, ndim=3)
        3D array of extrapolated azimuthal magnetic field, shape : (nphi,ntheta,nr)
    """

    shtns = require('shtns')

    nrout = len(rout)
    polar_opt = 1e-15

    if nphi is None:
        nphi   = int(max(256,lmax*3))
    ntheta = nphi//2

    norm=shtns.sht_orthonormal

    sh = shtns.sht(lmax,mmax=mmax,norm=norm)
    ntheta, nphi = sh.set_grid(ntheta, nphi, polar_opt=polar_opt)

    L = sh.l * (sh.l + 1)

    # Take care of shtns index convention

    bpolcmb = sh.spec_array()

    for l in range(1,lmax+1):
        for m in range(min(mmax+1,l+1)):
            bpolcmb[sh.idx(l,m)] = bpol[idx[l,m]]

    btor = np.zeros_like(bpolcmb)

    brout = np.zeros([ntheta,nphi,nrout])
    btout = np.zeros([ntheta,nphi,nrout])
    bpout = np.zeros([ntheta,nphi,nrout])

    for k,radius in enumerate(rout):
        print(("%d/%d" %(k,nrout)))

        radratio = rplanet/radius
        bpol = bpolcmb * radratio**(sh.l)
        brlm = bpol * L/radius**2
        brout[...,k] = sh.synth(brlm)

        slm = -sh.l/radius**2 * bpol

        btout[...,k], bpout[...,k] = sh.synth(slm,btor)

    brout = np.transpose(brout,(1,0,2))
    btout = np.transpose(btout,(1,0,2))
    bpout = np.transpose(bpout,(1,0,2))

    return brout, btout, bpout

def extrapot(planetname, glm, hlm, idx, lmax, mmax, rplanet, rout, nphi=None):
    """
    This function extrapolates a potential field to an array of desired radial
    levels. It uses the SHTns library (https://bitbucket.org/nschaeff/shtns)
    for spherical harmonic transforms, and falls back to scipy if SHTns is not
    available.

    Parameters
    ----------
    glm : ndarray(float, ndim=1)
        Array of Gauss coefficients g_lm
    hlm : ndarray(float, ndim=1)
        Array of Gauss coefficients h_lm
    idx : ndarray(int, ndim=1)
        Array of indices to map [l,m] to an index
    lmax : int
        Maximum spherical harmonic degree of field model
    mmax : int
        Maximum spherical harmonic order of field model
    rplanet : float
        Radius at which the magnetic field is measured or defined
    rout : array_like
        Array of radial levels to which the field should be extrapolated
    nphi : int, optional
        Number of grid points in longitude, can be automatically
        selected, by default None

    Returns
    -------
    brout : ndarray(float, ndim=3)
        3D array of extrapolated radial magnetic field, shape : (nphi,ntheta,nr)
    btout : ndarray(float, ndim=3)
        3D array of extrapolated co-latitudinal magnetic field, shape : (nphi,ntheta,nr)
    bpout : ndarray(float, ndim=3)
        3D array of extrapolated azimuthal magnetic field, shape : (nphi,ntheta,nr)
    """

    if _use_shtns():
        bpol = get_pol_from_Gauss(planetname, glm, hlm, lmax, mmax, idx)
        return extrapot_shtns(bpol, idx, lmax, mmax, rplanet, rout, nphi=nphi)

    return extrapot_scipy(glm, hlm, idx, lmax, mmax, rplanet, rout,
                          nphi=nphi, planetname=planetname)

def get_field_along_path_scipy(glm, hlm, idx, lmax, r, theta, phi, mmax=None, planetname=None):

    """Gets field along a specific trajectory defined by 1-D
       arrays r, theta, phi. Uses scipy for computing spherical harmonics.

    Parameters
    ----------
    glm : ndarray(float, ndim=1)
        Array of Gauss coefficients g_lm
    hlm : ndarray(float, ndim=1)
        Array of Gauss coefficients h_lm
    idx : ndarray(int, ndim=1)
        Array of indices to map [l,m] to an index
    lmax : int
        Maximum spherical harmonic degree of field model
    r : array_like
        Array of radial distances
    theta : array_like
        Array of co-latitudes in radians
    phi : array_like
        Array of longitudes in radians
    mmax : int, optional
        Maximum spherical harmonic order of field model, by default lmax
    planetname : str, optional
        Name of the planet, determines the phase convention, by default None

    Returns
    -------
    brout : array_like
        Array of extrapolated radial magnetic field values
    btout : array_like
        Array of extrapolated co-latitudinal magnetic field values
    bpout : array_like
        Array of extrapolated azimuthal magnetic field values

    Raises
    ------
    ValueError
        If the shapes of the three arrays r, theta, phi
        are not the same, raises an error.
    """
    r, theta, phi = _check_shapes(r, theta, phi)
    if mmax is None:
        mmax = lmax

    return _potential_field_scipy(glm, hlm, idx, lmax, mmax, r, theta, phi,
                                  planetname=planetname)


def get_field_along_path_shtns(bpol,idx,lmax,mmax,
                         rplanet,r,theta,phi):
    """Gets field along a specific trajectory defined by 1-D
       arrays r, theta, phi

    Parameters
    ----------
    bpol : ndarray(complex128, ndim=1)
        Array of poloidal coefficients computed from Gauss coefficients
    idx : ndarray(int, ndim=1)
        Array of indices to map [l,m] to an index
    lmax : int
        Maximum spherical harmonic degree of field model
    mmax : int
        Maximum spherical harmonic order of field model
    rplanet : float
        Radius at which the magnetic field is measured or defined
    r : array_like
        Array of radial distances
    theta : array_like
        Array of co-latitudes in radians
    phi : array_like
        Array of longitudes in radians

    Returns
    -------
    brout : array_like
        Array of extrapolated radial magnetic field values
    btout : array_like
        Array of extrapolated co-latitudinal magnetic field values
    bpout : array_like
        Array of extrapolated azimuthal magnetic field values

    Raises
    ------
    ValueError
        If the shapes of the three arrays r, theta, phi
        are not the same, raises an error.
    """

    r, theta, phi = _check_shapes(r, theta, phi)
    shtns = require('shtns')

    mmax = lmax
    norm=shtns.sht_orthonormal
    sh = shtns.sht(lmax,mmax=mmax,norm=norm)

    L = sh.l * (sh.l + 1)

    # Take care of shtns index convention

    bpolcmb = sh.spec_array()

    if mmax > 0:
        for l in range(1,lmax+1):
            for m in range(l+1):
                bpolcmb[sh.idx(l,m)] = bpol[idx[l,m]]
    else:
        for l in range(1,lmax+1):
                bpolcmb[sh.idx(l,0)] = bpol[idx[l,0]]

    brout = np.zeros_like(r)
    btout = np.zeros_like(r)
    bpout = np.zeros_like(r)

    # Assuming array of dimension 1 of r, theta, phi
    for k, radius in enumerate(r):
        radratio = rplanet/radius
        bpol = bpolcmb * radratio**(sh.l)
        qlm = bpol * L/radius**2
        slm = -sh.l/radius**2 * bpol
        tlm = np.zeros_like(qlm)
        brout[k], btout[k], bpout[k] = sh.SHqst_to_point(qlm,slm,tlm,
                                                         np.cos(theta[k]),
                                                         phi[k])

    return brout, btout, bpout


def get_field_along_path(planetname, glm, hlm, idx, lmax, mmax, r, theta, phi):
    """Gets field along a specific trajectory defined by 1-D
       arrays r, theta, phi

    Parameters
    ----------
    glm : ndarray(float, ndim=1)
        Array of Gauss coefficients g_lm
    hlm : ndarray(float, ndim=1)
        Array of Gauss coefficients h_lm
    idx : ndarray(int, ndim=1)
        Array of indices to map [l,m] to an index
    lmax : int
        Maximum spherical harmonic degree of field model
    r : array_like
        Array of radial distances
    theta : array_like
        Array of co-latitudes in radians
    phi : array_like
        Array of longitudes in radians

    Returns
    -------
    brout : array_like
        Array of extrapolated radial magnetic field values
    btout : array_like
        Array of extrapolated co-latitudinal magnetic field values
    bpout : array_like
        Array of extrapolated azimuthal magnetic field values

    Raises
    ------
    ValueError
        If the shapes of the three arrays r, theta, phi
        are not the same, raises an error.
    """

    r, theta, phi = _check_shapes(r, theta, phi)

    if _use_shtns():
        bpol = get_pol_from_Gauss(planetname, glm, hlm, lmax, mmax, idx)
        return get_field_along_path_shtns(bpol, idx, lmax, mmax, 1.0, r, theta, phi)

    return get_field_along_path_scipy(glm, hlm, idx, lmax, r, theta, phi,
                                      mmax=mmax, planetname=planetname)


def export_xshells(planet, filename, r=1.0, info=True):
    """
    Writes the potential magnetic field of a planet to a file readable by the
    xshells simulation code. See https://nschaeff.bitbucket.io/xshells

    Parameters
    ----------
    planet : Planet class instance
        Class containing Gauss coefficients,
    filename : string
        Name of file to export to
    r : float, optional
        Radial level for radial field computation, by default 1.0
    info : bool, optional
        Whether to print information about the planet, by default True

    Returns
    -------
    None
    """

    lmax = planet.lmax
    mmax = planet.mmax
    nlm = (mmax+1)*(lmax+1) - (mmax*(mmax+1))//2;
    bpol = np.zeros(nlm, dtype=complex)

    glm, hlm = planet.glm, planet.hlm
    idx = planet.idx

    i=1
    for l in range(1,lmax+1):   # m=0
        f = r**(-l-2)
        bpol[i] = glm[idx[l,0]] * f / l
        i+=1

    for m in range(1,mmax+1):
        for l in range(m,lmax+1):
            f = np.sqrt(0.5) * r**(-l-2)
            ix = idx[l,m]
            bpol[i] = (glm[ix] + 1.j*hlm[ix]) * f / l
            i+=1

    with open(filename,"w") as f:
        f.write("%%XS Pol lmax=%d mmax=%d\n" % (lmax,mmax))
        f.write("%%XS %s surface magnetic field from model %s, exported by planetMagFields, see https://github.com/AnkitBarik/planetMagFields\n" % (planet.name, planet.model))
        for q in bpol:
            f.write("%10.7g %10.7g\n" % (np.real(q),np.imag(q)))

    if info:
        print(("Planet: %s" %planet.name.capitalize()))
        print("Model: %s" %planet.model)
        if planet.name == 'earth':
            print("Year = %d" %planet.year)
        print("To use as an imposed field in Xshells, modify your xshells.par file to set:")
        print("  b = potential(%s)  # imposed from inner boundary" % filename)
        print("or")
        print("  b = potential(%s,out)  # imposed from outer boundary" % filename)

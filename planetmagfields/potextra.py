#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

def extrapot(rcmb,brcmb,rout,lmax=None):
    """
    This function extrapolates a potential field to an array of desired radial
    levels. It uses the SHTns library (https://bitbucket.org/nschaeff/shtns)
    for spherical harmonic transforms.

    Parameters
    ----------
    rcmb : float
        Radius at which the magnetic field is measured or defined
    brcmb : ndarray(float, ndim=2)
        2D array of radial magnetic field defined on a grid of longitude and
        co-latitude
    rout : array_like
        Array of radial levels to which the field should be extrapolated
    lmax : int, optional
        Maximum spherical harmonic degree, if None, automatically chosen
        from the grid, by default None

    Returns
    -------
    brout : ndarray(float, ndim=3)
        3D array of extrapolated radial magnetic field, shape : (nphi,ntheta,nr)
    btout : ndarray(float, ndim=3)
        3D array of extrapolated co-latitudinal magnetic field, shape : (nphi,ntheta,nr)
    bpout : ndarray(float, ndim=3)
        3D array of extrapolated azimuthal magnetic field, shape : (nphi,ntheta,nr)
    """

    nphi, ntheta = brcmb.shape
    nrout = len(rout)

    polar_opt = 1e-10

    if lmax is None:
        lmax = int(nphi/3)
    mmax = lmax

    try:
        import shtns
    except ImportError:
        print("Potential extrapolation requires the SHTns library")
        print("It can be obtained here: https://bitbucket.org/nschaeff/shtns")

    norm=shtns.sht_orthonormal | shtns.SHT_NO_CS_PHASE

    sh = shtns.sht(lmax,mmax=mmax,norm=norm)
    ntheta, nphi = sh.set_grid(ntheta, nphi, polar_opt=polar_opt)


    L = sh.l * (sh.l + 1)

    brlm = sh.analys(brcmb.T)
    bpolcmb = np.zeros_like(brlm)
    bpolcmb[1:] = rcmb**2 * brlm[1:]/L[1:]
    btor = np.zeros_like(brlm)

    brout = np.zeros([ntheta,nphi,nrout])
    btout = np.zeros([ntheta,nphi,nrout])
    bpout = np.zeros([ntheta,nphi,nrout])

    for k,radius in enumerate(rout):
        print(("%d/%d" %(k,nrout)))

        radratio = rcmb/radius
        bpol = bpolcmb * radratio**(sh.l)
        brlm = bpol * L/radius**2
        brout[...,k] = sh.synth(brlm)

        dbpoldr = -sh.l/radius * bpol
        slm = dbpoldr

        btout[...,k], bpout[...,k] = sh.synth(slm,btor)

    brout = np.transpose(brout,(1,0,2))
    btout = np.transpose(btout,(1,0,2))
    bpout = np.transpose(bpout,(1,0,2))

    return brout, btout, bpout

def get_grid(r,theta,phi):
    """
    Produces a 3D grid for storing values of extrapolated field. Used in producing
    vtk file for 3D visualization.

    Parameters
    ----------
    r : array_like
        Array of radial levels
    theta : array_like
        Array of colatitudes, ranging from 0 to pi
    phi : array_like
        Array of longitudes, ranging from 0 to 2*pi

    Returns
    -------
    r3D : ndarray(float, ndim=3)
        3D array of radius values at each grid point, shape (nphi,ntheta,nr)
    th3D : ndarray(float, ndim=3)
        3D array of colatitude values at each grid point, shape (nphi,ntheta,nr)
    p3D : ndarray(float, ndim=3)
        3D array of longitude values at each grid point, shape (nphi,ntheta,nr)
    x : ndarray(float, ndim=3)
        3D array of Cartesian x values at each grid point, shape (nphi,ntheta,nr)
    y : ndarray(float, ndim=3)
        3D array of Cartesian y values at each grid point, shape (nphi,ntheta,nr)
    z : ndarray(float, ndim=3)
        3D array of Cartesian z values at each grid point, shape (nphi,ntheta,nr)
    s : ndarray(float, ndim=3)
        3D array of cylindrical radius values at each grid point, shape (nphi,ntheta,nr)
    """

    nr,ntheta,nphi = len(r),len(theta),len(phi)

    r3D  = np.zeros([nphi,ntheta,nr])
    th3D = np.zeros([nphi,ntheta,nr])
    p3D  = np.zeros([nphi,ntheta,nr])

    for i in range(nr):
        r3D[...,i] = r[i]
    for j in range(ntheta):
        th3D[:,j,:] = theta[j]
    for k in range(nphi):
        p3D[k,...] = phi[k]

    s = r3D * np.sin(th3D)
    x = s *   np.cos(p3D)
    y = s *   np.sin(p3D)
    z = r3D * np.cos(th3D)

    return r3D,th3D,p3D, x,y,z, s

def get_cart(vr,vt,vp,th3D,p3D):
    """
    Converts a set of vector components in spherical coordinate system
    to Cartesian. All components, both at input and output have shape
    (nphi,ntheta,nr), number of points in longitude, colatitude and radial
    directions, respectively.

    Parameters
    ----------
    vr : ndarray(float, ndim=3)
        3D array of radial component of vector
    vt : ndarray(float, ndim=3)
        3D array of colatitude component of vector
    vp : ndarray(float, ndim=3)
        3D array of azimuthal component of vector
    th3D : ndarray(float, ndim=3)
        3D array of colatitude values on the grid
    p3D : ndarray(float, ndim=3)
        3D array of longitude values on the grid

    Returns
    -------
    vx : ndarray(float, ndim=3)
        3D array of vector component in Cartesian x direction
    vy : ndarray(float, ndim=3)
        3D array of vector component in Cartesian x direction
    vz : ndarray(float, ndim=3)
        3D array of vector component in Cartesian x direction
    """

    vs = vr * np.sin(th3D) + vt * np.cos(th3D)
    vz = vr * np.cos(th3D) - vt * np.sin(th3D)

    vx = vs * np.cos(p3D) - vp * np.sin(p3D)
    vy = vs * np.sin(p3D) + vp * np.cos(p3D)

    return vx,vy,vz

def writeVts(name,br,btheta,bphi,r,theta,phi):
    """
    Writes an unstructured vtk file for 3D visualization.

    Parameters
    ----------
    name : str
        Name of the file
    br : ndarray(float, ndim=3)
        3D array of radial component of vector
    btheta : ndarray(float, ndim=3)
        3D array of colatitudinal component of vector
    bphi : ndarray(float, ndim=3)
        3D array of azimuthal component of vector
    r : array_like
        Array of radial levels
    theta : array_like
        Array of colatitudes
    phi : array_like
        Array of longitudes

    Returns
    -------
    None
    """

    r3D,th3D,p3D, x,y,z, s = get_grid(r,theta,phi)

    print("grid shape=",th3D.shape)

    bx,by,bz = get_cart(br, btheta, bphi,th3D,p3D)

    try:
        try: # Version 2 changed naming convention of functions
            from evtk.hl import structuredToVTK
            gridToVTK = structuredToVTK
        except:
            import evtk
            gridToVTK = evtk.hl.gridToVTK
    except:
        try:
            from pyevtk.hl import gridToVTK
        except:
            print("This requires the use of evtk library!")
            print("You can get it from https://github.com/paulo-herrera/PyEVTK")


    br = np.asfortranarray(br)
    bx = np.asfortranarray(bx)
    by = np.asfortranarray(by)
    bz = np.asfortranarray(bz)

    gridToVTK("%s"%name,x,y,z,pointData= {"radius":r3D,
                                          "Radial mag field":br,
                                          "Mag Field":(bx, by,bz)})

    print("Output written to %s!" %name)
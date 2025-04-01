#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

def get_pol_from_Gauss(planetname,glm,hlm,lmax,mmax,idx):

    bpol = np.zeros(len(glm),dtype=np.complex128)

    if mmax > 0:
        for l in range(1,lmax+1):
            for m in range(l+1):

                if planetname in ["earth"]:
                    fac_m = 1.
                else:
                    fac_m = (-1)**m

                if m == 0:
                    norm = np.sqrt((4*np.pi/(2*l+1)))/l
                else:
                    norm = np.sqrt((2*np.pi/(2*l+1)))/l

                bpol[idx[l,m]] = norm*fac_m*(glm[idx[l,m]] - 1j*hlm[idx[l,m]])
    else:
        for l in range(1,lmax+1):

                fac_m = 1
                norm = np.sqrt((4*np.pi/(2*l+1)))/l

                bpol[idx[l]] = norm*fac_m*(glm[idx[l]] - 1j*hlm[idx[l]])

    return bpol


def extrapot(bpol,idx,lmax,model_mmax,rcmb,rout,nphi=None):
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
    model_mmax : int
        Maximum spherical harmonic order of field model
    rcmb : float
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

    # nphi, ntheta = brcmb.shape
    nrout = len(rout)

    polar_opt = 1e-15

    if nphi is None:
        nphi   = int(max(256,lmax*3))
    ntheta = nphi//2

    mmax = lmax

    try:
        import shtns
    except ImportError:
        print("Potential extrapolation requires the SHTns library")
        print("It can be obtained here: https://bitbucket.org/nschaeff/shtns")

    norm=shtns.sht_orthonormal

    sh = shtns.sht(lmax,mmax=mmax,norm=norm)
    ntheta, nphi = sh.set_grid(ntheta, nphi, polar_opt=polar_opt)

    L = sh.l * (sh.l + 1)

    # Take care of shtns index convention

    bpolcmb = sh.spec_array()

    if model_mmax > 0:
        for l in range(1,lmax+1):
            for m in range(l+1):
                bpolcmb[sh.idx(l,m)] = bpol[idx[l,m]]
    else:
        for l in range(1,lmax+1):
                bpolcmb[sh.idx(l,0)] = bpol[idx[l]]

    # brlm = sh.analys(brcmb.T)
    # bpolcmb = np.zeros_like(brlm)
    # bpolcmb[1:] = rcmb**2 * brlm[1:]/L[1:]
    btor = np.zeros_like(bpolcmb)

    brout = np.zeros([ntheta,nphi,nrout])
    btout = np.zeros([ntheta,nphi,nrout])
    bpout = np.zeros([ntheta,nphi,nrout])

    for k,radius in enumerate(rout):
        print(("%d/%d" %(k,nrout)))

        radratio = rcmb/radius
        bpol = bpolcmb * radratio**(sh.l)
        brlm = bpol * L/radius**2
        brout[...,k] = sh.synth(brlm)

        slm = -sh.l/radius**2 * bpol

        btout[...,k], bpout[...,k] = sh.synth(slm,btor)

    brout = np.transpose(brout,(1,0,2))
    btout = np.transpose(btout,(1,0,2))
    bpout = np.transpose(bpout,(1,0,2))

    return brout, btout, bpout

def get_field_along_path(bpol,idx,lmax,model_mmax,
                         rcmb,r,theta,phi):
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
    model_mmax : int
        Maximum spherical harmonic order of field model
    rcmb : float
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

    # Check dimensions
    if ( np.shape(r) != np.shape(theta)  or
         np.shape(theta) != np.shape(phi)  or
         np.shape(r) != np.shape(phi) ):
        raise ValueError("Please make sure all three arrays are of the same shape")

    # Ensure float
    r     = np.float64(r)
    theta = np.float64(theta)
    phi   = np.float64(phi)

    try:
        import shtns
    except ImportError:
        print("Orbit track extrapolation requires the SHTns library")
        print("It can be obtained here: https://bitbucket.org/nschaeff/shtns")

    mmax = lmax
    norm=shtns.sht_orthonormal
    sh = shtns.sht(lmax,mmax=mmax,norm=norm)

    L = sh.l * (sh.l + 1)

    # Take care of shtns index convention

    bpolcmb = sh.spec_array()

    if model_mmax > 0:
        for l in range(1,lmax+1):
            for m in range(l+1):
                bpolcmb[sh.idx(l,m)] = bpol[idx[l,m]]
    else:
        for l in range(1,lmax+1):
                bpolcmb[sh.idx(l,0)] = bpol[idx[l]]

    brout = np.zeros_like(r)
    btout = np.zeros_like(r)
    bpout = np.zeros_like(r)

    # Assuming array of dimension 1 of r, theta, phi
    for k, radius in enumerate(r):
        radratio = rcmb/radius
        bpol = bpolcmb * radratio**(sh.l)
        qlm = bpol * L/radius**2
        slm = -sh.l/radius**2 * bpol
        tlm = np.zeros_like(qlm)
        brout[k], btout[k], bpout[k] = sh.SHqst_to_point(qlm,slm,tlm,
                                                         np.cos(theta[k]),
                                                         phi[k])

    return brout, btout, bpout


###########################
# For writing vts files
###########################

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

def writeVts(name,br,btheta,bphi,r,theta,phi,r_planet=1):
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

    r3D,th3D,p3D, x,y,z, s = get_grid(r*r_planet,theta,phi)

    print("grid shape=",th3D.shape)

    bx,by,bz = get_cart(br, btheta, bphi,th3D,p3D)

    try:
        from pyevtk.hl import gridToVTK
    except:
        print("This requires the use of pyevtk library!")
        print("You can get it from https://github.com/paulo-herrera/PyEVTK")

    br = np.asfortranarray(br)
    bx = np.asfortranarray(bx)
    by = np.asfortranarray(by)
    bz = np.asfortranarray(bz)

    gridToVTK("%s"%name,x,y,z,pointData= {"radius":r3D,
                                          "Radial mag field":br,
                                          "Mag Field":(bx, by,bz)})

    print("Output written to %s!" %name)


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
            f = sqrt(0.5) * r**(-l-2)
            ix = idx[l,m]
            bpol[i] = (glm[ix] + 1.j*hlm[ix]) * f / l
            i+=1

    f = open(filename,"w")
    f.write("%%XS Pol lmax=%d mmax=%d\n" % (lmax,mmax))
    f.write("%%XS %s surface magnetic field from model %s, exported by planetMagFields, see https://github.com/AnkitBarik/planetMagFields\n" % (planet.name, planet.model))
    for q in bpol:
        f.write("%10.7g %10.7g\n" % (real(q),imag(q)))
    f.close()

    if info:
        print(("Planet: %s" %planet.name.capitalize()))
        print("Model: %s" %planet.model)
        if planet.name == 'earth':
            print("Year = %d" %planet.year)
        print("To use as an imposed field in Xshells, modify your xshells.par file to set:")
        print("  b = potential(%s)  # imposed from inner boundary" % filename)
        print("or")
        print("  b = potential(%s,out)  # imposed from outer boundary" % filename)

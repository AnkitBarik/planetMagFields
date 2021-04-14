#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

def extrapot(lmax,rcmb,brcmb,rout):

    nphi, ntheta = brcmb.shape
    nrout = len(rout)

    polar_opt = 1e-10

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

def get_cart(vr,vt,vp,r3D,th3D,p3D):

    vs = vr * np.sin(th3D) + vt * np.cos(th3D)
    vz = vr * np.cos(th3D) - vt * np.sin(th3D)

    vx = vs * np.cos(p3D) - vp * np.sin(p3D)
    vy = vs * np.sin(p3D) + vp * np.cos(p3D)

    return vx,vy,vz

def writeVts(name,br,btheta,bphi,r,theta,phi):

    r3D,th3D,p3D, x,y,z, s = get_grid(r,theta,phi)

    print("grid shape=",th3D.shape)

    bx,by,bz = get_cart(br, btheta, bphi,r3D,th3D,p3D)

    try:
        try: # Version 2 changed naming convention of functions
            #import evtk
            from evtk.hl import structuredToVTK
            #gridToVTK = evtk.hl.structuredToVTK
            gridToVTK = structuredToVTK
        except:
            import evtk
            gridToVTK = evtk.hl.gridToVTK
    except:
        print("movie2vtk requires the use of evtk library!")
        print("You can get it from https://github.com/paulo-herrera/PyEVTK")


    br = np.asfortranarray(br)
    bx = np.asfortranarray(bx)
    by = np.asfortranarray(by)
    bz = np.asfortranarray(bz)

    gridToVTK("%s"%name,x,y,z,pointData= {"radius":r3D,
                                          "Radial mag field":br,
                                          "Mag Field":(bx, by,bz)})


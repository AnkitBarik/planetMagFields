#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import scipy.special as sp

def get_data(datDir,planet="earth"):
    
    datfile = datDir + planet + '.dat'

    if planet == "earth":
        datfile = datDir + 'IGRF13.dat'
        dat = np.loadtxt(datfile,usecols=[1,2,-2])
        gh  = np.genfromtxt(datfile,usecols=[0],dtype='str')

        mask = gh == 'g'
        gDat = dat[mask,:]
        g   = gDat[:,-1]
        gl  = gDat[:,0]

        mask = gh == 'h'
        hDat = dat[mask,:]
        h   = hDat[:,-1]
        hm  = hDat[:,1]
        hIdx = np.where(hm == 1.)[0]
        h   = np.insert(h,hIdx,0.)

        lmax = np.int32(gl.max())

    elif planet in ["mercury","saturn"]:
        dat = np.loadtxt(datfile,usecols=[3])
        g   = dat.flatten()
        lmax = len(g) - 1
        h = np.zeros_like(g)
    
    elif planet == "jupiter":
        dat = np.loadtxt(datfile)
        l_dat = dat[:,-2]
        ghlm = dat[:,1]

        lmax = np.int32(l_dat.max())

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

    elif planet in ['uranus','neptune','ganymede']:
        dat = np.loadtxt(datfile,usecols=[1,2,3])
        gh  = np.genfromtxt(datfile,usecols=[0],dtype='str')

        mask = gh == 'g'
        gDat = dat[mask,:]
        g   = gDat[:,-1]
        gl  = gDat[:,0]

        mask = gh == 'h'
        hDat = dat[mask,:]
        h   = hDat[:,-1]

        gl = np.int32(gl)

        lmax = np.int32(gl.max())

    return g,h,lmax

def get_grid(nphi=256,ntheta=128):
    
    phi    = np.linspace(0.,2*np.pi,nphi)
    x,w    = sp.roots_legendre(ntheta)
    theta  = np.sort(np.arccos(x))

    p2D  = np.zeros([nphi, ntheta])
    th2D = np.zeros([nphi, ntheta])

    for i in range(nphi):
        p2D[i,:] = phi[i]

    for j in range(ntheta):
        th2D[:,j] = theta[j]

    return p2D, th2D

def gen_arr(lmax, l1,m1,mode='g'):

    '''
    Generate Gauss coefficient array for testing
    purposes
    '''

    idx = np.zeros([lmax+1,lmax+1])
    lArr = []
    mArr = []

    count = 0

    glm = []
    hlm = []

    for l in range(lmax+1):
        for m in range(l+1):

            if l in l1 and m in m1:
                if mode == 'g' or mode == 'gh':
                    glm.append(1.)
                else:
                    glm.append(0.)
                if mode == 'h' or mode == 'gh':
                    hlm.append(1.)
                else:
                    hlm.append(0.)
            else:
                glm.append(0.)
                hlm.append(0.)

            idx[l,m] = count
            lArr.append(l)
            mArr.append(m)
            
            count += 1

    glm  = np.array(glm)
    hlm  = np.array(hlm)
    lArr = np.array(lArr)
    mArr = np.array(mArr)
    idx =  np.int32(idx)

    return glm, hlm, lArr, mArr, idx

def gen_idx(lmax):

    '''
    Generate index to convert from (l,m) to [idx]
    '''

    idx = np.zeros([lmax+1,lmax+1])
    count = 0

    for l in range(1,lmax+1):
        for m in range(l+1):

            idx[l,m] = count
            count += 1

    return np.int32(idx)

def getB(lmax,glm,hlm,idx,r,p2D,th2D,planet="earth"):
    
    Br = np.zeros_like(p2D)

    for l in range(1,lmax+1):
        for m in range(l+1):
            ylm = sp.sph_harm(m, l, p2D, th2D)
            if planet in ["earth"]:
                fac_m = 1.
            else:
                fac_m = (-1)**m
            fac = fac_m * (l+1) * r**(-l-2) * np.sqrt((4.*np.pi)/(2*l+1))
            
            G = glm[idx[l,m]] * np.real(ylm)
            H = hlm[idx[l,m]] * np.imag(ylm) 

            Br +=   np.real(fac * (G + H))

    return Br

def getBm0(lmax,g,p2D,th2D):

    Br = np.zeros_like(p2D)

    for l in range(lmax):
        l1 = l+1
        ylm = sp.sph_harm(0, l1, p2D, th2D)
        fac = (l1 + 1) * np.sqrt((4.*np.pi)/(2*l+1))

        Br += fac * g[l] * np.real(ylm)

    return Br

def sphInt(f,g,phi,th2D,theta):
    
    thetaInt = np.trapz(f * g * np.sin(th2D),theta,axis=1)
    phiInt = np.trapz(thetaInt,phi)

    return phiInt


def getGauss(lmax,Br,r,phi,theta,th2D,p2D):

    '''
    Get Gauss coefficients from a surface field
    '''
    glm = []
    hlm = []
    
    comp = complex(0,1)

    for l in range(0,lmax+1):
        for m in range(0,l+1):
           
            ylm = (-1)**m * sp.sph_harm(m, l, p2D, th2D)

            ylm_conj = np.conjugate(ylm)

            fac = r**(l+2)/(l+1) 

            if m==0:
                fac *= 0.5
            
            I1 = sphInt(Br,ylm,phi,th2D,theta)
            I2 = sphInt(Br,ylm_conj,phi,th2D,theta)
          
            g = fac * (I2 + I1) 
            h = comp * fac * (I2 - I1) 
            
            glm.append(g)
            hlm.append(h)
            
    glm = np.array(glm)
    hlm = np.array(hlm)

    return glm, hlm

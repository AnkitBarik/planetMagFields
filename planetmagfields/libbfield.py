#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from .libgauss import get_data,gen_idx,get_grid,getB,getBm0
from .plotlib import *
from .utils import planetlist, stdDatDir

import sys



def getBr(datDir=stdDatDir,planet="earth",r=1,info=True):

    if planet not in planetlist:
        print("Planet must be one of the following!")
        print(planetlist)
        sys.exit()

    if planet == '--help':
        print("Usage: ./magField.py <planet name>")
        print("Example: ./magField.py earth")
        sys.exit()

    g,h,lmax,idx = get_data(datDir,planet)

    p2D,th2D = get_grid()

    if planet in ["mercury","saturn"]:
        Br = getBm0(lmax,g,r,p2D,th2D) * 1e-3
        dipTheta = 0.
        dipPhi = 0.
    else:
        Br = getB(lmax,g,h,idx,r,p2D,th2D,planet=planet) * 1e-3

        dipTheta = np.arctan(np.sqrt(g[idx[1,1]]**2 + h[idx[1,1]]**2)/g[idx[1,0]]) * 180./np.pi
        dipPhi = np.arctan(h[idx[1,1]]/g[idx[1,1]]) * 180./np.pi

    if info:
        print(("Planet: %s" %planet.capitalize()))
        #print(("Depth (fraction of surface radius) = %.2f" %r))
        print(("l_max = %d" %lmax))
        print(("Dipole tilt (degrees) = %f" %dipTheta))

    return p2D, th2D, Br, dipTheta, dipPhi

def plotAllFields(datDir=stdDatDir,r=1.0,levels=30,cmap='RdBu_r',proj='Mollweide'):

    print("")
    print('|=========|======|=======|')
    print(('|%-8s | %-2s| %-5s |' %('Planet','Theta','Phi')))
    print('|=========|======|=======|')

    plt.figure(figsize=(12,12))
    for k, planet in enumerate(planetlist):
        p2D,th2D,Br,dipTheta,dipPhi = getBr(datDir=datDir,planet=planet,r=r,info=False)

        if planet == "ganymede":
            nplot = 8
        else:
            nplot = k+1

        if proj.lower() == 'hammer':
            ax = plt.subplot(3,3,nplot)
        else:
            import cartopy.crs as ccrs
            projection = eval('ccrs.'+proj+'()')
            ax = plt.subplot(3,3,nplot,projection=projection)

        plotB_subplot(p2D,th2D,Br,ax,planet=planet,levels=levels,cmap=cmap,proj=proj)

        if planet in ["mercury","saturn"]:
            print(('|%-8s | %-4.1f | %-5.1f |' %(planet.capitalize(),dipTheta, dipPhi)))
        else:
            print(('|%-8s | %-3.1f | %-5.1f |' %(planet.capitalize(),dipTheta, dipPhi)))

    print('|---------|------|-------|')

    if r == 1:
        plt.suptitle(r'Radial magnetic field ($\mu$T) at surface', fontsize=20)
    else:
        plt.suptitle(r'Radial magnetic field ($\mu$T) at $r/r_{\rm surface} = %.2f$' %r, fontsize=20)




def plotMagField(datDir=stdDatDir,planet="earth",r=1,levels=30,proj='moll',cmap='RdBu_r'):

    p2D, th2D, Br, dum1,dum2 = getBr(datDir=datDir,planet=planet,r=r)
    plt.figure(figsize=(12,6.75))
    plotB(p2D,th2D,Br,planet=planet,levels=levels,proj=proj,cmap=cmap,r=r)

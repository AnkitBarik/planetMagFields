#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from .libgauss import get_data,gen_idx,get_grid,getB,getBm0
from .plotlib import *
import cartopy.crs as ccrs
import sys

planetlist = ["mercury","earth","jupiter","saturn","uranus","neptune","ganymede"]

def getBr(datDir="data/",planet="earth",r=1,info=True):

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
        Br = getBm0(lmax,g,p2D,th2D) * 1e-3
        dipTheta = 0.
        dipPhi = 0.
    else:
        Br = getB(lmax,g,h,idx,r,p2D,th2D,planet=planet) * 1e-3

        dipTheta = np.arctan(np.sqrt(g[idx[1,1]]**2 + h[idx[1,1]]**2)/g[idx[1,0]]) * 180./np.pi
        dipPhi = np.arctan(h[idx[1,1]]/g[idx[1,1]]) * 180./np.pi

    if info:
        print(("Planet: %s" %planet.capitalize()))
        print(("Depth (fraction of surface radius) = %.2f" %r))
        print(("l_max = %d" %lmax))
        print(("Dipole tilt (degrees) = %f" %dipTheta))

    return p2D, th2D, Br, dipTheta, dipPhi

def plotAllFields(datDir="data/",r=1.0):

    print("")
    print('|=========|======|=======|')
    print(('|%-8s | %-2s| %-5s |' %('Planet','Theta','Phi')))
    print('|=========|======|=======|')

    plt.figure(figsize=(12,12))
    for k, planet in enumerate(planetlist):
        p2D,th2D,Br,dipTheta,dipPhi = getBr(datDir=datDir,planet=planet,r=r,info=False)

        if planet == "ganymede":
            ax = plt.subplot(3,3,8,projection=ccrs.Mollweide())
        else:
            ax = plt.subplot(3,3,k+1,projection=ccrs.Mollweide())

        plotB_subplot(p2D,th2D,Br,ax,planet=planet)

        if planet in ["mercury","saturn"]:
            print(('|%-8s | %-4.1f | %-5.1f |' %(planet.capitalize(),dipTheta, dipPhi)))
        else:
            print(('|%-8s | %-3.1f | %-5.1f |' %(planet.capitalize(),dipTheta, dipPhi)))

    print('|---------|------|-------|')



def plotMagField(datDir="data/",planet="earth",r=1):

    p2D, th2D, Br, dum1,dum2 = getBr(datDir=datDir,planet=planet,r=r)
    plt.figure(figsize=(16,9))
    plotB(p2D,th2D,Br,planet=planet)


if __name__=="__main__":

    if len(sys.argv) == 3:
        planet = str(sys.argv[1]).lower()
        r      = np.float32(sys.argv[2])
    elif len(sys.argv) == 2:
        print("Radius not specified, using surface\n")
        planet = str(sys.argv[1]).lower()
        r = 1.
    elif len(sys.argv) > 3:
        print("Too many arguments, exiting ...\n")
        sys.exit()
    else:
        print("Planet or radius not specified, plotting for Earth's surface\n")
        planet="earth"
        r=1.

    if planet == 'all':
        plotAllFields()
    else:
        plotMagField(planet=planet,r=r)

    plt.tight_layout()
    plt.show()

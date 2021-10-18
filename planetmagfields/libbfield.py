#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from .libgauss import get_grid,getB,getBm0
from .plotlib import *
from .utils import planetlist, stdDatDir


def getBr(planet, r=1, info=True):
    p2D,th2D = get_grid()

    if planet.name in ["mercury", "saturn"]:
        Br = getBm0(planet.lmax,
                    planet.glm,
                    r,
                    p2D,
                    th2D) * 1e-3
        dipTheta = 0.
        dipPhi = 0.
    else:
        Br = getB(planet.lmax,
                  planet.glm,
                  planet.hlm,
                  planet.idx,
                  r,
                  p2D,
                  th2D,
                  planet=planet.name) * 1e-3

        dipTheta = np.arctan(np.sqrt(planet.glm[planet.idx[1,1]]**2 + planet.hlm[planet.idx[1,1]]**2)
                                    /planet.glm[planet.idx[1,0]]) * 180./np.pi
        dipPhi = np.arctan(planet.hlm[planet.idx[1,1]]/planet.glm[planet.idx[1,1]]) * 180./np.pi

    if info:
        print(("Planet: %s" %planet.name.capitalize()))
        #print(("Depth (fraction of surface radius) = %.2f" %r))
        print(("l_max = %d" %planet.lmax))
        print(("Dipole tilt (degrees) = %f" %dipTheta))

    return p2D, th2D, Br, dipTheta, dipPhi

def plotAllFields(datDir=stdDatDir,r=1.0,levels=30,cmap='RdBu_r',proj='Mollweide'):

    from .planet import planet as Planet

    print("")
    print('|=========|======|=======|')
    print(('|%-8s | %-2s| %-5s |' %('Planet','Theta','Phi')))
    print('|=========|======|=======|')

    plt.figure(figsize=(12,12))

    for k, name in enumerate(planetlist):
        planet = Planet(name=name,r=r,info=False)

        if name == "ganymede":
            nplot = 8
        else:
            nplot = k+1

        if proj.lower() == 'hammer':
            ax = plt.subplot(3,3,nplot)
        else:
            import cartopy.crs as ccrs
            projection = eval('ccrs.'+proj+'()')
            ax = plt.subplot(3,3,nplot,projection=projection)

        plotB_subplot(planet.p2D,
                      planet.th2D,
                      planet.Br,
                      ax,
                      planet=name,
                      levels=levels,
                      cmap=cmap,
                      proj=proj)

        if name in ["mercury","saturn"]:
            print(('|%-8s | %-4.1f | %-5.1f |' %(name.capitalize(),planet.dipTheta, planet.dipPhi)))
        else:
            print(('|%-8s | %-3.1f | %-5.1f |' %(name.capitalize(),planet.dipTheta, planet.dipPhi)))

    print('|---------|------|-------|')

    if r == 1:
        plt.suptitle(r'Radial magnetic field ($\mu$T) at surface', fontsize=20)
    else:
        plt.suptitle(r'Radial magnetic field ($\mu$T) at $r/r_{\rm surface} = %.2f$' %r, fontsize=20)


def plotMagField(name,r=1,levels=30,proj='moll',cmap='RdBu_r'):
    from .planet import planet as Planet

    planet = Planet(name)
    plt.figure(figsize=(12,6.75))
    plotB(planet.p2D,
          planet.th2D,
          planet.Br,
          planet=planet.name,
          levels=levels,proj=proj,cmap=cmap,r=planet.r)

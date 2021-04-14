#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

def hammer2cart(ttheta, pphi, colat=False):
    """
    This function is used to define the Hammer projection for
    default plotting. Copied from MagIC python plotting script:

    https://github.com/magic-sph/magic/blob/master/python/magic/plotlib.py
    """

    if not colat: # for lat and phi \in [-pi, pi]
        xx = 2.*np.sqrt(2.) * np.cos(ttheta)*np.sin(pphi/2.)\
             /np.sqrt(1.+np.cos(ttheta)*np.cos(pphi/2.))
        yy = np.sqrt(2.) * np.sin(ttheta)\
             /np.sqrt(1.+np.cos(ttheta)*np.cos(pphi/2.))
    else:  # for colat and phi \in [0, 2pi]
        xx = -2.*np.sqrt(2.) * np.sin(ttheta)*np.cos(pphi/2.)\
             /np.sqrt(1.+np.sin(ttheta)*np.sin(pphi/2.))
        yy = np.sqrt(2.) * np.cos(ttheta)\
             /np.sqrt(1.+np.sin(ttheta)*np.sin(pphi/2.))
    return xx, yy

def plotSurf(p2D,th2D,B,levels=60,cmap='RdBu_r',proj='Mollweide'):

    bmax = np.abs(B).max()
    digits = int(np.log10(bmax)) + 1

    if digits > 1:
        bmax = np.round(bmax)
    else:
        bmax = np.round(bmax,decimals=1)

    divnorm = colors.TwoSlopeNorm(vmin=-bmax, vcenter=0, vmax=bmax)
    cs = np.linspace(-bmax,bmax,levels)


    lon2D = p2D - np.pi
    lat2D = np.pi/2 - th2D

    try:
        import cartopy.crs as ccrs
    except:
        print("cartopy library not available, using Hammer projection")
        proj = 'hammer'

    if proj.lower() == 'hammer':
        ax = plt.axes()
        xx,yy = hammer2cart(lat2D,lon2D)
        cont = ax.contourf(xx,yy,B,cs,cmap=cmap,norm=divnorm,extend='both')
    else:
        projection = eval('ccrs.'+proj+'()')

        ax = plt.axes(projection=projection)

        cont = ax.contourf(lon2D*180/np.pi,lat2D*180/np.pi,B,cs,  \
            transform=ccrs.PlateCarree(),cmap=cmap,norm=divnorm,extend='both')

    cbar = plt.colorbar(cont,orientation='horizontal',fraction=0.06, pad=0.04,ticks=[-bmax,0,bmax])

    ax.axis('equal')
    ax.axis('off')

    return ax, cbar

def plotB(p2D,th2D,B,r=1,planet="earth",levels=60,cmap='RdBu_r',proj='Mollweide'):

    planet = planet.lower()

    ax,cbar = plotSurf(p2D,th2D,B,levels=levels,cmap=cmap,proj=proj)
    cbar.ax.set_xlabel(r'Radial magnetic field ($\mu$T)',fontsize=25)
    cbar.ax.tick_params(labelsize=20)

    if proj.lower() != 'hammer' and planet == 'earth':
        ax.coastlines()

    if r==1:
        radLabel = '  Surface'
    else:
        radLabel = r'  $r/r_{\rm surface}=%.2f$' %r

    ax.set_title(planet.capitalize() + radLabel,fontsize=25,pad=20)

def plotB_subplot(p2D,th2D,B,ax,planet="earth",levels=60,cmap='RdBu_r',proj='Mollweide'):
    planet = planet.lower()

    bmax = np.abs(B).max()
    digits = int(np.log10(bmax)) + 1

    if digits > 1:
        bmax = np.round(bmax)
    else:
        bmax = np.round(bmax,decimals=1)

    p2D -= np.pi
    th2D -= np.pi/2
    th2D = -th2D

    cs = np.linspace(-bmax,bmax,levels)
    divnorm = colors.TwoSlopeNorm(vmin=-bmax, vcenter=0, vmax=bmax)

    try:
        import cartopy.crs as ccrs
    except:
        print("cartopy library not available, using Hammer projection")
        proj = 'hammer'

    if proj.lower() == 'hammer':
        xx,yy = hammer2cart(th2D,p2D)
        cont = ax.contourf(xx,yy,B,cs,cmap=cmap,norm=divnorm,extend='both')
    else:
        if planet == "earth":
            ax.coastlines()

        cont = ax.contourf(p2D*180/np.pi,th2D*180/np.pi,B,cs,  \
            transform=ccrs.PlateCarree(),cmap=cmap,norm=divnorm,extend='both')

    cbar = plt.colorbar(cont,orientation='horizontal',fraction=0.06, pad=0.04,ticks=[-bmax,0,bmax])
    cbar.ax.tick_params(labelsize=15)

    ax.set_title(planet.capitalize(),fontsize=20)
    ax.axis('equal')
    ax.axis('off')

def plot_spec(l,E,r,planet):

    plt.semilogy(l[1:],E[1:],'-o',lw=1.2,color='#449c99',mfc='#347b79',mec='#347b79',ms=8)
    plt.xlabel(r'$l$',fontsize=30)
    plt.ylabel(r'$R_l$ (nT$^2$)',fontsize=30)
    plt.xticks(l)

    plt.tick_params(which='major',labelsize=15,length=10)
    plt.tick_params(which='minor',length=5)

    plt.grid(True,alpha=0.5)

    if r==1:
        radLabel = '  Surface'
    else:
        radLabel = r'  $r/r_{\rm surface}=%.2f$' %r

    plt.title('Lowes spectrum, '+ planet.capitalize() + radLabel,fontsize=20,pad=10)

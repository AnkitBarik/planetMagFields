#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cartopy.crs as ccrs

def plotB(p2D,th2D,B,r=1,planet="earth"):

    planet = planet.lower()

    bmax = B.max()
    bmin = B.min()

    projection = ccrs.Mollweide()
    ax = plt.axes(projection=projection)

    if planet == "earth":
        ax.coastlines()

    lon2D = p2D - np.pi
    lat2D = np.pi/2 - th2D

    divnorm = colors.TwoSlopeNorm(vmin=bmin, vcenter=0, vmax=bmax)

    cont = ax.contourf(lon2D*180/np.pi,lat2D*180/np.pi,B,100,  \
           transform=ccrs.PlateCarree(),cmap='RdBu_r',norm=divnorm)

    cbar = plt.colorbar(cont,orientation='horizontal',fraction=0.06, pad=0.04,ticks=[bmin,0,bmax])
    cbar.ax.set_xlabel(r'Radial magnetic field ($\mu$T)',fontsize=30)
    cbar.ax.tick_params(labelsize=20)

    if r==1:
        radLabel = '  Surface'
    else:
        radLabel = r'  $r/r_{\rm surface}=%.2f$' %r

    ax.set_title(planet.capitalize() + radLabel,fontsize=30,pad=20)

def plotB_subplot(p2D,th2D,B,ax,planet="earth"):
    planet = planet.lower()

    bmax = B.max()
    bmin = B.min()

    if planet == "earth":
        ax.coastlines()

    p2D -= np.pi
    th2D -= np.pi/2
    th2D = -th2D

    divnorm = colors.TwoSlopeNorm(vmin=bmin, vcenter=0, vmax=bmax)

    cont = ax.contourf(p2D*180/np.pi,th2D*180/np.pi,B,100,  \
           transform=ccrs.PlateCarree(),cmap='RdBu_r',norm=divnorm)

    cbar = plt.colorbar(cont,orientation='horizontal',fraction=0.06, pad=0.04,ticks=[bmin,0,bmax])
    cbar.ax.set_xlabel(r'Radial magnetic field ($\mu$T)',fontsize=15)
    cbar.ax.tick_params(labelsize=15)

    ax.set_title(planet.capitalize(),fontsize=20,pad=20)

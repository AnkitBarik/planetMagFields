#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import os
import sys
sys.path.append(os.path.abspath('../../'))
from planetmagfields import Planet
from planetmagfields.libgauss import get_grid
import cartopy.crs as ccrs

def plot_surface(ax,p2D,th2D,B,levels=21,cmap='RdBu_r'):

    bmax = np.abs(B).max()
    digits = int(np.log10(bmax)) + 1

    if digits > 1:
        bmax = np.round(bmax)
    else:
        bmax = np.round(bmax,decimals=1)

    divnorm = colors.TwoSlopeNorm(vmin=-bmax, vcenter=0, vmax=bmax)
    cs = np.linspace(-bmax,bmax,levels)
    im = cm.ScalarMappable(norm=divnorm)

    lon2D = p2D - np.pi
    lat2D = np.pi/2 - th2D

    cont = ax.contourf(lon2D*180/np.pi,lat2D*180/np.pi,B,cs,  \
        transform=ccrs.PlateCarree(),cmap=cmap,norm=divnorm,extend='both',
        zorder=-9)

    for c in cont.collections:
        c.set_edgecolor('face')

    ax.axis('equal')
    ax.axis('off')
    ax.set_rasterization_zorder(-1)

    return bmax, cont

p2D,th2D = get_grid()
p = Planet(name='jupiter',r=0.85,nphi=256,info=False,model='jrm33')

proj = ccrs.Mollweide()
fig,ax = plt.subplots(figsize=(16,6),nrows=1,ncols=2,subplot_kw={'projection': proj})

br_ref = np.loadtxt('./Br_reference085.dat')

bmax, im = plot_surface(ax[0],p2D,th2D,br_ref)
bmax, im = plot_surface(ax[1],p2D,th2D,p.Br)

ax[0].set_title(r'jupiterMag $B_r$',fontsize=30)
ax[1].set_title(r'planetMagFields $B_r$',fontsize=30)

fig.suptitle('Jupiter, $r/r_{surface}=0.85$',fontsize=40)

cbar = plt.colorbar(im,ax=ax.ravel().tolist(),
                    orientation='horizontal',fraction=0.06,
                    pad=0.04,ticks=[-bmax,0,bmax])
cbar.ax.set_xlabel(r'Radial magnetic field ($\mu$T)',fontsize=20)
cbar.ax.tick_params(labelsize=20)

plt.tight_layout()
#plt.savefig('../../paper/figures/jup_bench.pdf',
#            dpi=200,bbox_inches='tight')
plt.show()

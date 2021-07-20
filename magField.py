#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from planetmagfields.libbfield import *
import sys

levels=30
cmap='RdBu_r'
proj = 'Mollweide'
r = 1

if len(sys.argv) > 4:
    print("Too many arguments, exiting ...\n")
    sys.exit()
elif len(sys.argv) == 4:
    planet = str(sys.argv[1]).lower()
    r      = np.float32(sys.argv[2])
    proj   = str(sys.argv[3])
elif len(sys.argv) == 3:
    planet = str(sys.argv[1]).lower()
    try:
        r      = np.float32(sys.argv[2])
    except:
        proj   = str(sys.argv[2])
elif len(sys.argv) == 2:
    print("Radius not specified, using surface\n")
    planet = str(sys.argv[1]).lower()
else:
    print("Planet or radius not specified, plotting for Earth's surface\n")
    planet="earth"
    r=1.

if planet == 'all':
    plotAllFields(datDir='./planetmagfields/data/',r=r,levels=levels,cmap=cmap,proj=proj)
    plt.tight_layout()
    plt.subplots_adjust(top=0.895,
                        bottom=0.035,
                        left=0.023,
                        right=0.976,
                        hspace=0.38,
                        wspace=0.109)
else:
    plotMagField(planet=planet,r=r,datDir='./planetmagfields/data/',levels=levels,cmap=cmap,proj=proj)
    plt.tight_layout()

plt.show()

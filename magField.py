#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import argparse
from planetmagfields.libbfield import plotAllFields
from planetmagfields.planet import Planet


parser=argparse.ArgumentParser(
    description='''Script for easy plotting of planetary magnetic field.''')

parser.add_argument('-p','--planet',type=str,default='earth',
                    help='Planet name (default : earth)',
                    dest='planet')
parser.add_argument('-r','--radius',type=float,default=1,
                    help='Radial level scaled to planetary radius (default : 1)',
                    dest='r')
parser.add_argument('-c','--cmap',type=str,default='RdBu_r',
                    help='Colormap of plot (default : RdBu_r)',
                    dest='cmap')
parser.add_argument('-l','--levels',type=int,default=20,
                    help='Number of contour levels (default : 20)',
                    dest='levels')
parser.add_argument('-m','--mapproj',type=str,default='Mollweide',
                    help='Type of map projection (default : Mollweide)',
                    dest='proj')

parser.add_argument('-o','--model',type=str,default=None,
                    help='Model to be used, uses the latest model by default (default : None)',
                    dest='model')

args=parser.parse_args()

planet = args.planet
r      = args.r
cmap   = args.cmap
levels = args.levels
proj   = args.proj
model  = args.model

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
    pl = Planet(name=planet,r=r,model=model)
    pl.plot(r=r,levels=levels,cmap=cmap,proj=proj)

plt.show()

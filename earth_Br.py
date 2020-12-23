#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from libgauss import gen_idx,get_grid,getB
from plotlib import plotB

datDir = 'data/'

dat = np.loadtxt(datDir + 'IGRF13.dat',usecols=[1,2,-2])
gh  = np.genfromtxt(datDir + 'IGRF13.dat',usecols=[0],dtype='str')

r = 1.

mask = gh == 'g'
gDat = dat[mask,:]

g   = gDat[:,-1]
gl  = gDat[:,0]
gm  = gDat[:,1]

mask = gh == 'h'
hDat = dat[mask,:]

h   = hDat[:,-1]
hm  = hDat[:,1]

hIdx = np.where(hm == 1.)[0]

h   = np.insert(h,hIdx,0.)

lmax = np.int32(gl.max())
idx = gen_idx(lmax)

p2D,th2D = get_grid()

Br = getB(lmax,g,h,idx,r,p2D,th2D) * 1e-3

plotB(p2D,th2D,Br)
plt.show()

######################
# Other useful stuff
######################

# gm  = gDat[:,1]
# hl  = hDat[:,0]
# hm  = np.insert(hm,hIdx,0)
# hl  = np.insert(hl,hIdx,hl[hIdx])

# gl = np.int32(gl)
# hl = np.int32(hl)
# gm = np.int32(gm)
# hm = np.int32(hm)

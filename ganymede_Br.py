#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from libgauss import gen_idx,get_grid,getB
from plotlib import plotB

datDir = 'data/'

dat = np.loadtxt(datDir + 'ganymede.dat',usecols=[1,2,3])
gh  = np.genfromtxt(datDir + 'ganymede.dat',usecols=[0],dtype='str')

r = 1.

mask = gh == 'g'
gDat = dat[mask,:]
g   = gDat[:,-1]
gl  = gDat[:,0]

mask = gh == 'h'
hDat = dat[mask,:]
h   = hDat[:,-1]

gl = np.int32(gl)

lmax = np.int32(gl.max())

idx = gen_idx(lmax)

p2D,th2D = get_grid()

Br = getB(lmax,g,h,idx,r,p2D,th2D) * 1e-3

plotB(p2D,th2D,Br,planet='ganymede')
plt.show()
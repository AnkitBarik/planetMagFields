#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from libgauss import gen_idx,get_grid,getB
from plotlib import plotB

datDir = 'data/'

dat = np.loadtxt(datDir + 'jupiter.dat')

l_dat = dat[:,-2]
m_dat = dat[:,-1]
ghlm = dat[:,1]

r = 1.
lmax = np.int32(l_dat.max())

g = []
h = []
gm = []
hm = []
gl = []
hl = []

########################
# Separate glm and hlm
########################

for i in range(1,lmax+1):
    mask = l_dat == i
    n = len(l_dat[mask])
    half = int(n/2)
    
    g.append(ghlm[mask][:half+1])
    h.append(np.concatenate([[0.],ghlm[mask][half+1:]]))

    gm.append(m_dat[mask][:half+1])
    hm.append(np.concatenate([[0.],m_dat[mask][half+1:]]))

    gl.append(l_dat[mask][:half+1])
    hl.append(np.concatenate([[i],l_dat[mask][half+1:]]))

    
g  = np.concatenate(g)
h  = np.concatenate(h)
gm = np.concatenate(gm)
hm = np.concatenate(hm)
gl = np.concatenate(gl)
hl = np.concatenate(hl)


gl = np.int32(gl)
hl = np.int32(hl)
gm = np.int32(gm)
hm = np.int32(hm)

idx = gen_idx(lmax)

p2D,th2D = get_grid()

Br = getB(lmax,g,h,idx,r,p2D,th2D) * 1e-3

Br = np.roll(Br,70,axis=0)

plotB(p2D,th2D,Br,planet='jupiter')
plt.show()
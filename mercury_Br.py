#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from libgauss import gen_idx,get_grid,getBm0
from plotlib import plotB

import sys

datDir = 'data/'

dat = np.loadtxt(datDir + 'mercury.dat')

r = 1.

gl  = dat[:,0]
g   = dat[:,2]

lmax = np.int32(gl.max())
print(lmax)
nphi = 256
nlat = nphi/2


idx = gen_idx(lmax)

p2D,th2D = get_grid()

Br = getBm0(lmax,g,p2D,th2D) * 1e-3

plotB(p2D,th2D,Br,planet="mercury")
plt.show()
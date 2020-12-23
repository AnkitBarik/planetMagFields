#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from libgauss import gen_idx,get_grid,getBm0
from plotlib import plotB

datDir = 'data/'

dat = np.loadtxt(datDir + 'saturn.dat')

r = 1.

g   = dat.flatten()

lmax = len(g) - 1
print(lmax)
nphi = 256
nlat = nphi/2


idx = gen_idx(lmax)

p2D,th2D = get_grid()

Br = getBm0(lmax,g,p2D,th2D) * 1e-3

plotB(p2D,th2D,Br,planet="saturn")
plt.show()
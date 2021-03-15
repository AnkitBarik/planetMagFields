#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from .libgauss import get_data
from .magField import getBr
from .plotlib import plotB

class planet:

    def __init__(self,name='earth',datDir='./data/'):
    
        self.name   = name.lower()
        self.datDir = datDir
        self.glm, self.hlm, self.lmax, self.idx = \
                get_data(self.datDir,planet=self.name)

        self.p2D, self.th2D, self.Br, self.dipTheta, self.dipPhi = \
                getBr(datDir=self.datDir,planet=self.name,r=1,info=True)

        self.phi = self.p2D[:,0]
        self.theta = self.th2D[0,:]

    def plot(self):
        plt.figure(figsize=(16,9))
        plotB(self.p2D,self.th2D,self.Br,planet=self.name)

    def writeVtsFile(self,potExtra=False,ratio_out=2,nrout=32):
            from .potextra import extrapot, writeVts

            rout = np.linspace(1,ratio_out,nrout)
            if potExtra:
                brout, btout, bpout = extrapot(self.lmax,1.,self.Br,rout)
            else:
                brout = self.Br
                btout = np.zeros_like(self.Br)
                bpout = np.zeros_like(self.Br)

            writeVts(self.name,brout,btout,bpout,rout,self.theta,self.phi)

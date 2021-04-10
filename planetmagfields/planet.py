#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from .libgauss import get_data, filt_Gauss, filt_Gaussm0,getB, getBm0, get_spec
from .libbfield import getBr
from .plotlib import plotB, plot_spec

class planet:

    def __init__(self,name='earth',datDir='./planetmagfields/data/'):

        self.name   = name.lower()
        self.datDir = datDir
        self.glm, self.hlm, self.lmax, self.idx = \
                get_data(self.datDir,planet=self.name)

        self.p2D, self.th2D, self.Br, self.dipTheta, self.dipPhi = \
                getBr(datDir=self.datDir,planet=self.name,r=1,info=True)

        self.phi = self.p2D[:,0]
        self.theta = self.th2D[0,:]
        self.r = 1

    def plot(self,r=1,levels=30,cmap='RdBu_r',proj='Hammer'):
        plt.figure(figsize=(12,6.75))

        if r == 1:
            plotB(self.p2D,self.th2D,self.Br,planet=self.name,levels=levels,cmap=cmap,proj=proj)
        else:
            self.p2D, self.th2D, self.Br, self.dipTheta, self.dipPhi = \
                    getBr(datDir=self.datDir,planet=self.name,r=r,info=False)
            self.r = r
            plotB(self.p2D,self.th2D,self.Br,r=self.r,planet=self.name,levels=levels,cmap=cmap,proj=proj)

        plt.tight_layout()

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

    ## Filtered plots

    def plot_filt(self,r=1,larr=None,marr=None,lCutMin=0,lCutMax=None,mmin=0,mmax=None,levels=30,cmap='RdBu_r',proj='hammer'):

        self.larr_filt = larr
        self.marr_filt = marr
        self.lCutMin = lCutMin
        self.lCutMax = lCutMax
        self.mmin_filt = mmin
        self.mmax_filt = mmax

        if self.lCutMax is None:
            self.lCutMax = self.lmax
        if self.mmax_filt is None:
            self.mmax_filt = self.lmax

        self.r_filt = r

        if self.name in ['mercury','saturn']:
            self.glm_filt,self.hlm_filt =\
                    filt_Gaussm0(self.glm,self.hlm,self.lmax,larr=self.larr_filt,\
                        lCutMin=self.lCutMin,lCutMax=self.lCutMax)
            self.Br_filt = 1e-3*getBm0(self.lmax,self.glm_filt,self.r_filt,self.p2D,self.th2D)
        else:
            self.glm_filt,self.hlm_filt =\
                filt_Gauss(self.glm,self.hlm,self.lmax,self.idx,larr=self.larr_filt,\
                    marr=self.marr_filt,lCutMin=self.lCutMin,lCutMax=self.lCutMax,mmin=self.mmin_filt,mmax=self.mmax_filt)

            self.Br_filt = 1e-3*getB(self.lmax,self.glm_filt,self.hlm_filt,self.idx,self.r_filt,self.p2D,self.th2D,planet=self.name)

        plt.figure(figsize=(12,6.75))

        plotB(self.p2D,self.th2D,self.Br_filt,r=self.r_filt,planet=self.name,levels=levels,cmap=cmap,proj=proj)

        if r==1:
            radLabel = '  Surface'
        else:
            radLabel = r'  $r/r_{\rm surface}=%.2f$' %r

        if self.larr_filt is not None:
            elllabel = r', $l = %s$' %np.str(self.larr_filt)
        else:
            if self.lCutMin > 0:
                if self.lCutMax < self.lmax:
                    elllabel = r', $ %d \leq l \leq %d$' %(self.lCutMin,self.lCutMax)
                else:
                    elllabel = r', $l \geq %d$' %self.lCutMin

            elif self.lCutMax < self.lmax:
                elllabel = r', $l \leq %d$' %self.lCutMax

        if self.marr_filt is not None:
            elllabel += r', $m = %s$' %np.str(self.marr_filt)
        else:
            if self.mmin_filt > 0:
                if self.mmax_filt < self.lmax:
                    elllabel += r', $ %d \leq m \leq %d$' %(self.mmin_filt,self.mmax_filt)
                else:
                    elllabel += r', $m \geq %d$' %self.mmin_filt
            elif self.mmax_filt < self.lmax:
                elllabel += r', $m \leq %d$' %self.mmax_filt

        plt.title(self.name.capitalize() + radLabel + elllabel,fontsize=25,pad=20)
        plt.tight_layout()


    def spec(self,r=1,iplot=True):
        self.emag_spec, emag_10 = get_spec(self.glm,self.hlm,self.idx,self.lmax,planet=self.name,r=r)
        l = np.arange(self.lmax+1)

        self.dip_tot = self.emag_spec[1]/sum(self.emag_spec)
        self.dipolarity = emag_10/sum(self.emag_spec)
        if iplot:
            plt.figure(figsize=(7,7))
            plot_spec(l,self.emag_spec,r,self.name)
            plt.tight_layout()
            plt.show()
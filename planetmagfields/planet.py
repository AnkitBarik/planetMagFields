#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import numpy as np
import matplotlib.pyplot as plt
from .libgauss import get_data, filt_Gauss, filt_Gaussm0,getB, getBm0, get_spec
from .libbfield import getBr
from .plotlib import plotB, plotSurf, plot_spec
from .utils import stdDatDir, planetlist
import sys


class planet:

    def __init__(self,name='earth',datDir=stdDatDir):

        self.name   = name.lower()
        if self.name not in planetlist:
            print("Planet must be one of the following!")
            print(planetlist)

        self.datDir = datDir
        self.glm, self.hlm, self.lmax, self.idx = \
                get_data(self.datDir,planet=self.name)

        self.p2D, self.th2D, self.Br, self.dipTheta, self.dipPhi = \
                getBr(datDir=self.datDir,planet=self.name,r=1,info=True)

        self.phi = self.p2D[:,0]
        self.theta = self.th2D[0,:]
        self.r = 1

    def plot(self,r=1,levels=30,cmap='RdBu_r',proj='Mollweide'):
        plt.figure(figsize=(12,6.75))

        if r == 1:
            ax,cbar = plotSurf(self.p2D,self.th2D,self.Br,levels=levels,cmap=cmap,proj=proj)
        else:
            self.p2D, self.th2D, self.Br, self.dipTheta, self.dipPhi = \
                    getBr(datDir=self.datDir,planet=self.name,r=r,info=False)
            self.r = r
            ax,cbar = plotSurf(self.p2D,self.th2D,self.Br,levels=levels,cmap=cmap,proj=proj)

        cbar.ax.set_xlabel(r'Radial magnetic field ($\mu$T)',fontsize=25)
        cbar.ax.tick_params(labelsize=20)

        if r==1:
            radLabel = '  Surface'
        else:
            radLabel = r'  $r/r_{\rm surface}=%.2f$' %r

        if proj.lower() != 'hammer' and self.name == 'earth':
            ax.coastlines()

        if r==1:
            radLabel = '  Surface'
        else:
            radLabel = r'  $r/r_{\rm surface}=%.2f$' %r

        ax.set_title(self.name.capitalize() + radLabel,fontsize=25,pad=20)
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

    def plot_filt(self,r=1,larr=None,marr=None,lCutMin=0,lCutMax=None,mmin=0,mmax=None,levels=30,cmap='RdBu_r',proj='Mollweide'):

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

        ax,cbar = plotSurf(self.p2D,self.th2D,self.Br_filt,levels=levels,cmap=cmap,proj=proj)

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

        cbar.ax.set_xlabel(r'Radial magnetic field ($\mu$T)',fontsize=25)
        cbar.ax.tick_params(labelsize=20)

        if proj.lower() != 'hammer' and self.name == 'earth':
            ax.coastlines()
        ax.set_title(self.name.capitalize() + radLabel + elllabel,fontsize=25,pad=20)
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

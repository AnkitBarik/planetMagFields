#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import numpy as np
import matplotlib.pyplot as plt
from .libgauss import get_data, filt_Gauss, filt_Gaussm0,getB, getBm0, get_spec
from .libbfield import getBr
from .plotlib import plotSurf, plot_spec
from .utils import stdDatDir, planetlist


class Planet:
    """
    Planet class

    The Planet class contains all information about a planet. It contains
    arrays of Gauss coefficients, glm and hlm, the maximum spherical harmonic
    degree lmax to which data is available, and also computes and stores the
    (optionally filtered) radial magnetic field at a surface and the Lowes
    spectrum.
    """

    def __init__(self,name='earth',r=1.0,nphi=256,datDir=stdDatDir,info=True):
        """
        Initialization of the Planet class.

        Parameters
        ----------
        name : str, optional
            Name of the planet, by default 'earth'
        r : float, optional
            Radial level to compute and plot field on, scaled by the planetary
            radius, by default 1.0
        nphi : int, optional
            Number of points in longitude, number of points in colatitude
            are automatically set to half this number, by default 256
        datDir : str, optional
            Data directory, where the Gauss coefficient data is present,
            named as <planetname>.dat, the standard directory is ./data,
            by default stdDatDir
        info : bool, optional
            If True, prints some information about the planet, by default True
        """

        self.name   = name.lower()
        self.nphi   = nphi
        self.ntheta = nphi//2

        if self.name not in planetlist:
            print("Planet must be one of the following!")
            print(planetlist)

        self.datDir = datDir
        self.glm, self.hlm, self.lmax, self.idx = get_data(self.datDir,
                                                           planetname=self.name)

        self.p2D, self.th2D, self.Br, self.dipTheta, self.dipPhi = getBr(self,r=r,
                                                                        nphi=self.nphi,
                                                                        ntheta=self.ntheta,
                                                                        info=info)

        self.phi = self.p2D[:,0]
        self.theta = self.th2D[0,:]
        self.r = r

    def plot(self,r=1.0,levels=30,cmap='RdBu_r',proj='Mollweide'):
        """
        Plots the radial magnetic field of a planet at a radial surface.

        Parameters
        ----------
        r : float, optional
            Radial surface for plot, by default 1
        levels : int, optional
            Number of contour levels, by default 30
        cmap : str, optional
            Colormap for contours, by default 'RdBu_r'
        proj : str, optional
            Map projection, by default 'Mollweide'

        Returns
        -------
        None
        """

        plt.figure(figsize=(12,6.75))

        if r == 1:
            ax,cbar = plotSurf(self.p2D,self.th2D,self.Br,levels=levels,cmap=cmap,proj=proj)
        else:
            self.p2D, self.th2D, self.Br, self.dipTheta, self.dipPhi = \
                    getBr(planet=self,r=r,info=False)
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
        """
        Writes a .vts file for 3D visualization. Uses the SHTns library
        for potential extrapolation and the pyevtk library for writing
        the vts file.

        Parameters
        ----------
        potExtra : bool, optional
            Whether to use potential extrapolation, by default False
        ratio_out : int, optional
            Radial level to which the magnetic field needs to be upward
            continued, scaled to planetary radius, by default 2
        nrout : int, optional
            Number of radial grid levels, by default 32

        Returns
        -------
        None
        """
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

    def plot_filt(self,r=1.0,larr=None,marr=None,lCutMin=0,lCutMax=None,mmin=0,mmax=None,levels=30,cmap='RdBu_r',proj='Mollweide'):
        """
        Plots a filtered radial magnetic field at a radial level. Filters can be
        set using specific values of degree and order of spherical harmonics given
        through the arrays larr and marr or by providing a range using lCutMin,
        lCutMax and mmin, mmax.

        Parameters
        ----------
        r : float, optional
            Radial level for plot, scaled to planetary radius, by default 1
        larr : array_like, optional
            Array of spherical harmonic degrees, if None, uses lmax, by default None
        marr : array_like, optional
            Array of spherical harmonic orders, if None, uses lmax, by default None
        lCutMin : int, optional
            Minimum spherical harmonic degree to retain, by default 0
        lCutMax : int, optional
            Maximum spherical harmonic degree to retain, if None, uses lmax, by default None
        mmin : int, optional
            Minimum spherical harmonic order to retain, by default 0
        mmax : int, optional
            Maximum spherical harmonic degree to retain, if None, uses lmax, by default None
        levels : int, optional
            Number of contour levels, by default 30
        cmap : str, optional
            Colormap for contours, by default 'RdBu_r'
        proj : str, optional
            Map projection, by default 'Mollweide'

        Returns
        -------
        None
        """

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

            self.Br_filt = 1e-3*getB(self.lmax,self.glm_filt,self.hlm_filt,self.idx,self.r_filt,self.p2D,self.th2D,planetname=self.name)

        plt.figure(figsize=(12,6.75))

        ax,cbar = plotSurf(self.p2D,self.th2D,self.Br_filt,levels=levels,cmap=cmap,proj=proj)

        if r==1:
            radLabel = '  Surface'
        else:
            radLabel = r'  $r/r_{\rm surface}=%.2f$' %r

        if self.larr_filt is not None:
            elllabel = r', $l = %s$' %str(self.larr_filt)
        else:
            if self.lCutMin > 0:
                if self.lCutMax < self.lmax:
                    elllabel = r', $ %d \leq l \leq %d$' %(self.lCutMin,self.lCutMax)
                else:
                    elllabel = r', $l \geq %d$' %self.lCutMin

            elif self.lCutMax < self.lmax:
                elllabel = r', $l \leq %d$' %self.lCutMax

        if self.marr_filt is not None:
            elllabel += r', $m = %s$' %str(self.marr_filt)
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


    def spec(self,r=1.0,iplot=True):
        """
        General plot of Lowes spectrum of a planet at a radial level, scaled
        to planetary radius. Also computes dipolarity (energy of axial dipole)
        to total and total dipolarity (dipTot, energy of total dipole to total
        magnetic energy)

        Parameters
        ----------
        r : float, optional
            Radial level scaled to planetary radius, by default 1.0
        iplot : bool, optional
            If True, generates a plot, by default True

        Returns
        -------
        None
        """
        self.emag_spec, emag_10 = get_spec(self.glm,self.hlm,self.idx,self.lmax,planetname=self.name,r=r)
        l = np.arange(self.lmax+1)

        self.dip_tot = self.emag_spec[1]/sum(self.emag_spec)
        self.dipolarity = emag_10/sum(self.emag_spec)
        if iplot:
            plt.figure(figsize=(7,7))
            plot_spec(l,self.emag_spec,r,self.name)
            plt.tight_layout()
            plt.show()

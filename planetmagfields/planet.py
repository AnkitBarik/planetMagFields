#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from .libdata import get_data
from .libgauss import filt_Gauss,getB, get_spec
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

    def __init__(self,name='earth',model=None,year=None,
                 r=1.0,nphi=256,datDir=stdDatDir,unit='muT',info=True):
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
        unit : str, optional
            Units of magnetic field, can be 'nT', 'muT' or 'Gauss' for nanoTeslas,
            microTeslas and Gauss, respectively. By default, 'muT'
        info : bool, optional
            If True, prints some information about the planet, by default True
        """

        self.name   = name.lower()
        self.nphi   = nphi
        self.ntheta = nphi//2
        self.unit   = unit

        if self.unit.lower() == 'mut':
            self.unitfac = 1e-3
            self.unitlabel = '$\mu$T'
        elif self.unit.lower() == 'nt':
            self.unitfac = 1.
            self.unitlabel = 'nT'
        elif self.unit.lower() == 'gauss':
            self.unitfac = 1e-5
            self.unitlabel = 'Gauss'

        #Automatic selection of latest model
        if model is None:
            if self.name =='earth':
                model = 'igrf14'
            elif self.name =='mercury':
                model = 'wardinski2019'
            elif self.name == 'jupiter':
                model = 'jrm33'
            elif self.name == 'saturn':
                model = 'cassini11+'
            elif self.name == 'uranus':
                model = 'connerny1987'
            elif self.name == 'neptune':
                model = 'connerny1991'
            elif self.name == 'ganymede':
                model = 'kivelson2002'

        self.model = model

        if year is None:
            self.year = 2020
        else:
            self.year = year

        if self.name not in planetlist:
            print("Planet must be one of the following!")
            print(planetlist)

        self.datDir = datDir
        self.glm, self.hlm, self.lmax, self.idx, self.mmax = get_data(self.datDir,
                                                           planetname=self.name,
                                                           model = self.model,
                                                           year=self.year)

        self.p2D, self.th2D, self.Br, self.dipTheta, self.dipPhi = getBr(self,r=r,
                                                                        nphi=self.nphi,
                                                                        ntheta=self.ntheta,
                                                                        info=info)

        self.Br *= self.unitfac

        self.phi = self.p2D[:,0]
        self.theta = self.th2D[0,:]
        self.r = r

    def plot(self,r=None,levels=30,cmap='RdBu_r',
             proj='Mollweide',vmin=None,vmax=None):
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
        vmin : float, optional
            Minimum of colorscale, by default None
        vmax : float, optional
            Maximum of colorscale, by default None

        Returns
        -------
        fig : matplotlib.pyplot.figure instance
            Figure handle
        ax  : matplotlib.axes.Axes instance
            Figure axis
        cbar: matplotlib.axes.Axes instance
            Colorbar axis
        """

        fig = plt.figure(figsize=(12,6.75))

        if r is None:
            r = self.r

        if r == self.r:
            ax,cbar,proj = plotSurf(self.p2D,self.th2D,self.Br,
                                    levels=levels,cmap=cmap,proj=proj,
                                    vmin=vmin,vmax=vmax)
        else:
            self.p2D, self.th2D, self.Br, self.dipTheta, self.dipPhi = \
                    getBr(planet=self,r=r,info=False)
            self.r = r
            self.Br *= self.unitfac
            ax,cbar,proj = plotSurf(self.p2D,self.th2D,self.Br,
                                    levels=levels,cmap=cmap,proj=proj,
                                    vmin=vmin,vmax=vmax)

        cbar.ax.set_xlabel(r'Radial magnetic field (%s)' %self.unitlabel,fontsize=25)
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

        title = self.name.capitalize() + radLabel

        if self.name == 'earth':
            title = title + ', %d' %self.year

        ax.set_title(title,fontsize=25,pad=20)
        plt.tight_layout()

        return fig, ax, cbar


    def extrapolate(self,rout):
        """Potential extrapolation of the magnetic field

        Parameters
        ----------
        rout : array_like
            Array of radial levels

        Returns
        -------
        None
            Assigns three arrays self.br_ex,self.btheta_ex,self.bphi_ex to
            the planet class for radial, colatitudinal and azimuthal components
            of the extrapolated field, respectively.
        """
        from .potextra import get_pol_from_Gauss, extrapot

        # Create poloidal potential from glm and glm
        bpol = get_pol_from_Gauss(self.name,self.glm,self.hlm,
                                  self.lmax,self.mmax,self.idx)

        self.br_ex,self.btheta_ex,self.bphi_ex \
            = extrapot(bpol,self.idx,self.lmax,self.mmax,1,rout,self.nphi)

        self.br_ex     *= self.unitfac
        self.btheta_ex *= self.unitfac
        self.bphi_ex   *= self.unitfac

    def orbit_path(self,r,theta,phi):
        """Extrapolates the magnetic field along an orbit trajectory.
           Assigns objects self.br_orb, self.btheta_orb, self.bphi_orb
           which are extrpolated values of radial, co-latitudinal and
           azimuthal components of the magnetic field, respectively.

        Parameters
        ----------
        r : array_like
            Array of radial distances
        theta : array_like
            Array of co-latitudes in radians
        phi : array_like
            Array of longitudes in radians
        """
        from .potextra import get_pol_from_Gauss, get_field_along_path

        bpol = get_pol_from_Gauss(self.name,self.glm,self.hlm,
                            self.lmax,self.mmax,self.idx)

        self.br_orb,self.btheta_orb,self.bphi_orb=\
                get_field_along_path(bpol,self.idx,self.lmax,self.mmax,
                                   1,r,theta,phi)

        self.br_orb     *= self.unitfac
        self.btheta_orb *= self.unitfac
        self.bphi_orb   *= self.unitfac

    def writeVtsFile(self,potExtra=False,ratio_out=2,nrout=32,r_planet=1):
        """
        Writes an unstructured vtk (.vts) file for 3D visualization. Uses the
        SHTns library for potential extrapolation and the pyevtk library for
        writing the vtk file.

        Parameters
        ----------
        potExtra : bool, optional
            Whether to use potential extrapolation, by default False
        ratio_out : int, optional
            Radial level to which the magnetic field needs to be upward
            continued, scaled to planetary radius, by default 2
        nrout : int, optional
            Number of radial grid levels, by default 32
        r_planet : float, optional
            Radius of planet to get coordinates in dimensional units

        Returns
        -------
        None
        """
        from .potextra import writeVts

        rout = np.linspace(1,ratio_out,nrout)
        if potExtra:
            self.extrapolate(rout)
            brout = self.br_ex
            btout = self.btheta_ex
            bpout = self.bphi_ex
        else:
            brout = self.Br
            btout = np.zeros_like(self.Br)
            bpout = np.zeros_like(self.Br)

        writeVts(self.name,brout,btout,bpout,rout,self.theta,self.phi,r_planet)

    ## Filtered plots

    def plot_filt(self,r=1.0,larr=None,marr=None,lCutMin=0,lCutMax=None,mmin=0,mmax=None,
                  levels=30,cmap='RdBu_r',proj='Mollweide',
                  vmin=None,vmax=None,iplot=True):
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
        vmin : float, optional
            Minimum of colorscale, by default None
        vmax : float, optional
            Maximum of colorscale, by default None
        iplot: logical, optional
            Flag for producing a plot, by default True

        Returns
        -------
        fig : matplotlib.pyplot.figure instance
            Figure handle
        ax  : matplotlib.axes.Axes instance
            Figure axis
        cbar: matplotlib.axes.Axes instance
            Colorbar axis
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

        self.glm_filt,self.hlm_filt =\
                filt_Gauss(self.glm,self.hlm,self.lmax,self.mmax,self.idx,larr=self.larr_filt,
                           marr=self.marr_filt,lCutMin=self.lCutMin,lCutMax=self.lCutMax,
                           mmin=self.mmin_filt,mmax=self.mmax_filt)

        self.Br_filt = getB(self.lmax,self.mmax,self.glm_filt,self.hlm_filt,
                            self.idx,self.r_filt,self.p2D,self.th2D,planetname=self.name)
        self.Br_filt *= self.unitfac

        if iplot:
            fig = plt.figure(figsize=(12,6.75))

            ax,cbar,proj = plotSurf(self.p2D,self.th2D,self.Br_filt,levels=levels,
                                    cmap=cmap,proj=proj,vmin=vmin,vmax=vmax)

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

            cbar.ax.set_xlabel(r'Radial magnetic field (%s)' %self.unitlabel,fontsize=25)
            cbar.ax.tick_params(labelsize=20)

            if proj.lower() != 'hammer' and self.name == 'earth':
                ax.coastlines()
            ax.set_title(self.name.capitalize() + radLabel + elllabel,fontsize=25,pad=20)
            plt.tight_layout()

            return fig, ax, cbar


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
        self.emag_spec, emag_10, self.emag_symm, self.emag_antisymm, self.emag_axi \
            = get_spec(self.glm,self.hlm,
                       self.idx,self.lmax,
                       self.mmax,r=r)
        l = np.arange(self.lmax+1)

        self.dip_tot = self.emag_spec[1]/sum(self.emag_spec)
        self.dipolarity = emag_10/sum(self.emag_spec)
        self.emag_tot = sum(self.emag_spec)
        if iplot:
            plt.figure(figsize=(7,7))
            plot_spec(l,self.emag_spec,r,self.name)
            plt.tight_layout()
            plt.show()

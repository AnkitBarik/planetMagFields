#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from .libdata import get_data
from .libgauss import filt_Gauss, getB, get_grid, get_spec, get_dipole_tilt
from .plotlib import (plotSurf, plotB_subplot, plot_spec, radius_label,
                      filter_label, _import_cartopy)
from .models import check_planet, default_model, planetlist
from .utils import stdDatDir, get_unit


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
                 r=1.0,nphi=256,datDir=stdDatDir,units='muT',info=True):
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
        units : str, optional
            Units of magnetic field, can be 'nT', 'muT' or 'Gauss' for nanoTeslas,
            microTeslas and Gauss, respectively. By default, 'muT'
        info : bool, optional
            If True, prints some information about the planet, by default True
        """

        self.name   = check_planet(name)
        self.nphi   = nphi
        self.ntheta = nphi//2
        self.units  = units
        self.unitfac, self.unitlabel = get_unit(units)

        self.model = model if model is not None else default_model(self.name)
        self.year = 2020 if year is None else year

        self.datDir = datDir
        self.glm, self.hlm, self.lmax, self.idx, self.mmax = get_data(self.datDir,
                                                           planetname=self.name,
                                                           model=self.model,
                                                           year=self.year)

        self.p2D, self.th2D, self.phi, self.theta = get_grid(nphi=self.nphi,
                                                             ntheta=self.ntheta)
        self.dipTheta, self.dipPhi = get_dipole_tilt(self.glm, self.hlm,
                                                     self.idx, self.mmax)
        self.r = r
        self.Br = self.get_Br(r)

        if info:
            self.print_info()

    def get_Br(self, r):
        """Radial magnetic field in self.units at radial level r on the
        (p2D, th2D) grid."""
        return getB(self.lmax, self.mmax, self.glm, self.hlm, self.idx, r,
                    self.p2D, self.th2D, planetname=self.name) * self.unitfac

    def print_info(self):
        print("Planet: %s" %self.name.capitalize())
        print("Model: %s" %self.model)
        print("l_max = %d" %self.lmax)
        print("Dipole tilt (degrees) = %f" %self.dipTheta)
        if self.name == 'earth':
            print("Year = %d" %self.year)

    def _plot_map(self, Br, title, levels, cmap, proj, vmin, vmax):
        fig = plt.figure(figsize=(12,6.75))
        ax,cbar,proj = plotSurf(self.p2D,self.th2D,Br,levels=levels,cmap=cmap,
                                proj=proj,vmin=vmin,vmax=vmax)

        cbar.ax.set_xlabel(r'Radial magnetic field (%s)' %self.unitlabel,fontsize=25)
        cbar.ax.tick_params(labelsize=20)

        if proj.lower() != 'hammer' and self.name == 'earth':
            ax.coastlines()

        ax.set_title(title,fontsize=25,pad=20)
        plt.tight_layout()

        return fig, ax, cbar

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

        if r is None:
            r = self.r

        Br = self.Br if r == self.r else self.get_Br(r)

        title = self.name.capitalize() + radius_label(r)
        if self.name == 'earth':
            title = title + ', %d' %self.year

        return self._plot_map(Br, title, levels, cmap, proj, vmin, vmax)


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
        from .potextra import extrapot

        self.br_ex,self.btheta_ex,self.bphi_ex \
            = extrapot(self.name,self.glm,self.hlm,self.idx,self.lmax,self.mmax,1,rout,self.nphi)

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
       from .potextra import get_field_along_path

       self.br_orb,self.btheta_orb,self.bphi_orb=\
               get_field_along_path(self.name,self.glm,self.hlm,self.idx,self.lmax,self.mmax,r,theta,phi)

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

        from .lib3d import writeVts
        writeVts(self.name,brout,btout,bpout,rout,self.theta,self.phi,r_planet)

    def plot3D(self, fieldlines=False,ratio_out=2,nrout=32,r_planet=1):
        """
        Plots the 3D magnetic field of the planet.
        """
        from .lib3d import plot_surface, render_field_lines

        if fieldlines:
            rout = np.linspace(r_planet,ratio_out,nrout)
            pl = render_field_lines(self.name, self.glm, self.hlm, self.idx, self.lmax, self.mmax, 1,
                       rout, nphi=128, surf=True, clim_fac=1.0,
                       units=self.units, bgcolor='white', cmap='seismic',
                       lightweight=False)
            pl.show()
        else:
            pl, _ = plot_surface(self.theta,self.phi,self.Br,fieldname='Br',cmap='seismic',clim_fac=1, bgcolor='white')
            pl.show()


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

        self.Br_filt = self.unitfac * getB(self.lmax,self.mmax,self.glm_filt,self.hlm_filt,
                            self.idx,self.r_filt,self.p2D,self.th2D,planetname=self.name)

        if iplot:
            title = (self.name.capitalize() + radius_label(r)
                     + filter_label(self.lmax, self.larr_filt, self.marr_filt,
                                    self.lCutMin, self.lCutMax,
                                    self.mmin_filt, self.mmax_filt))
            return self._plot_map(self.Br_filt, title, levels, cmap, proj, vmin, vmax)


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


def plotAllFields(datDir=stdDatDir,r=1.0,levels=30,cmap='RdBu_r',
                  proj='Mollweide',units='muT',vmin=None,vmax=None):
    """
    Plots fields of all the planets for which data is available. It's provided in
    models.planetlist.

    Parameters
    ----------
    datDir : str, optional
        Data directory, where the Gauss coefficient data is present,
        named as <planetname>_<modelname>.dat, by default stdDatDir
    r : float, optional
        Radial level to compute and plot field on, scaled by the planetary
        radius, by default 1.0
    levels : int, optional
        Number of contour levels, by default 30
    cmap : str, optional
        Colormap for contours, by default 'RdBu_r'
    proj : str, optional
        Map projection, by default 'Mollweide'
    units : str, optional
        Units of magnetic field, can be 'nT', 'muT' or 'Gauss', by default 'muT'
    vmin : float, optional
        Minimum of colorscale, by default None
    vmax : float, optional
        Maximum of colorscale, by default None
    """

    print("")
    print('|=========|======|=======|')
    print(('|%-8s | %-2s| %-5s |' %('Planet','Theta','Phi')))
    print('|=========|======|=======|')

    plt.figure(figsize=(12,12))

    ccrs = None if proj.lower() == 'hammer' else _import_cartopy()
    if ccrs is None:
        proj = 'hammer'

    for k, name in enumerate(planetlist):
        planet = Planet(name=name,datDir=datDir,r=r,info=False,units=units)

        nplot = 8 if name == "ganymede" else k+1

        if ccrs is None:
            ax = plt.subplot(3,3,nplot)
        else:
            ax = plt.subplot(3,3,nplot,projection=getattr(ccrs, proj)())

        plotB_subplot(ax,planet.p2D,
                      planet.th2D,
                      planet.Br,
                      planetname=name,
                      levels=levels,
                      cmap=cmap,
                      proj=proj,
                      vmin=vmin,
                      vmax=vmax)

        if name in ["mercury","saturn"]:
            print(('|%-8s | %-4.1f | %-5.1f |' %(name.capitalize(),planet.dipTheta, planet.dipPhi)))
        else:
            print(('|%-8s | %-3.1f | %-5.1f |' %(name.capitalize(),planet.dipTheta, planet.dipPhi)))

    print('|---------|------|-------|')

    unitlabel = get_unit(units)[1]
    if r == 1:
        plt.suptitle(r'Radial magnetic field (%s) at surface' %unitlabel, fontsize=20)
    else:
        plt.suptitle(r'Radial magnetic field (%s) at $r/r_{\rm surface} = %.2f$' %(unitlabel,r), fontsize=20)

#!/usr/bin/env python
# -*- coding: utf-8 -*-

from netCDF4 import Dataset
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline, PchipInterpolator, Akima1DInterpolator
from matplotlib.ticker import FixedLocator, FixedFormatter
import matplotlib.colors as colors
#
class Lion_data:
    def __init__(self,data_params):
        default_params = {
            "filename": "dataset.nc",  "time_var_name": "time",
            "temp_var_name": "temp", "salt_var_name": "salt",
            "zr_var_name": "z_r", "time_start": 0, "time_end": 24*30,
            "time_scaling":1./(24.*3600.), "time_shift":0., "z_scaling":1.,"data_title":"ED"
        }
        params = default_params.copy()
        params.update(data_params)
        self.config = params

    def get_data(self,OBS):
        nc = Dataset(self.config["filename"] , mode='r')
        ts = self.config["time_start"]
        tf = self.config["time_end"]
        temp = nc.variables[self.config["temp_var_name"]][ts:tf,:]
        salt = nc.variables[self.config["salt_var_name"]][ts:tf,:]
        zz        = self.config["z_scaling"]*nc.variables[self.config["zr_var_name"]][:]
        tt        = self.config["time_scaling"]*nc.variables[self.config["time_var_name"]][ts:tf] - self.config["time_shift"]
        nc.close()
        self.nt   = len(tt)
        self.dpth = np.tile( zz, (self.nt,1) )
        if OBS:
            # remove anomalous values
            self.temp   = temp
            self.salt   = np.delete(salt     , [0,1, 3, 4, 5, 6, 8, 9, 10, 11, 13, 14, 15, 18, 25], axis=1)
            self.dpth_s = np.delete(self.dpth, [0,1, 3, 4, 5, 6, 8, 9, 10, 11, 13, 14, 15, 18, 25], axis=1)
            self.nz_s   = len(self.salt[0,:])
            print(self.nz_s)
            self.time_s = np.transpose(np.tile( tt,  (self.nz_s,1) ))
        else:
            self.temp = temp
            self.salt = salt
        self.nz   = len(zz)
        self.nt   = len(tt)
        self.time = np.transpose(np.tile( tt,  (self.nz,1) ))
        self.z1d  = zz

    def plot_temp_data(self,ax,tmin,tmax,colorMap,Tlevels,cmap):
        im = ax.pcolormesh(self.time, self.dpth, self.temp, cmap=colorMap, vmin=tmin, vmax=tmax, shading='auto')
        cs10 = ax.contour(self.time, self.dpth, self.temp, colors='k', levels=Tlevels, linewidths=1.5)
        ax.clabel(cs10, inline=True, fontsize=10, fmt = '%1.2f')
        ax.set_title(self.config["data_title"])
        ax.set_ylabel(r'depth $[{\rm m}]$')
        if cmap:
            f.subplots_adjust(right=0.8)
            cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
            f.colorbar(im, cax=cbar_ax)

    def plot_salt_data(self,ax,tmin,tmax,colorMap,Tlevels,cmap):
        if cmap:
            im = ax.pcolormesh(self.time_s, self.dpth_s, self.salt, cmap=colorMap, vmin=tmin, vmax=tmax, shading='auto')
            cs10 = ax.contour(self.time_s, self.dpth_s, self.salt, colors='k', levels=Tlevels, linewidths=1.5)
        else:
            im = ax.pcolormesh(self.time, self.dpth, self.salt, cmap=colorMap, vmin=tmin, vmax=tmax, shading='auto')
            cs10 = ax.contour(self.time, self.dpth, self.salt, colors='k', levels=Tlevels, linewidths=1.5)
        ax.clabel(cs10, inline=True, fontsize=10, fmt = '%1.2f')
        ax.set_title(self.config["data_title"])
        ax.set_ylabel(r'depth $[{\rm m}]$')
        if cmap:
            f.subplots_adjust(right=0.8)
            cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
            f.colorbar(im, cax=cbar_ax)
#=========================================================
# Creating different instance and calling class methods
#=========================================================
evd_params = {'filename': '../data/run_asics1.nc', 'time_var_name': 'ocean_time', 'data_title':'TKE+EVD'}
scm_evd = Lion_data(evd_params)
scm_evd.get_data(OBS=False)

edmf_params = {'filename': '../data/run_asics3.nc', 'time_var_name': 'ocean_time', 'data_title':'EDMF-Energy'}
scm_edmf = Lion_data(edmf_params)
scm_edmf.get_data(OBS=False)

data_params = {'filename': '../data/LION_mooring_data.nc', 'time_var_name': 'ocean_time', 'time_end': 598, 'data_title':'Lion mooring data'}
scm_data = Lion_data(data_params)
scm_data.get_data(OBS=True)
#=========================================================
# set min/max values and colormaps
tmin=12.875; tmax=13.205; colorMapT='twilight'; titleT='Temperature [Celsius]'
smin=38.45 ; smax=38.55 ; colorMapS='viridis' ; titleS='Salinity [Celsius]'
hmax=2300. ; Tlevels=[12.90, 13, 13.10]; Slevels=[38.49, 38.50, 38.52]
#=========================================================
# CREATE THE TEMPERATURE PLOT
#=========================================================
f, axes = plt.subplots(3, 1, sharex=True, sharey=True)
f.set_figheight(9); f.set_figwidth(10)
f.suptitle(titleT, fontsize=16)
ax = axes[0]; ax.set_ylim(-hmax,0.); ax.set_xlim(0.,26.)
scm_evd.plot_temp_data(ax,tmin,tmax,colorMapT,Tlevels,False)
ax = axes[1]
scm_edmf.plot_temp_data(ax,tmin,tmax,colorMapT,Tlevels,False)
ax = axes[2]
scm_data.plot_temp_data(ax,tmin,tmax,colorMapT,Tlevels,True)
ax.set_xticks([1,3,5,7,9,11,13,15,17,19,21,23,25])
ax.set_xticklabels([r'$16\;{\rm Jan}$',r'$18\;{\rm Jan}$',r'$20\;{\rm Jan}$',r'$22\;{\rm Jan}$',r'$24\;{\rm Jan}$',r'$26\;{\rm Jan}$',r'$28\;{\rm Jan}$',r'$30\;{\rm Jan}$',r'$01\;{\rm Feb}$',r'$03\;{\rm Feb}$',r'$05\;{\rm Feb}$',r'$07\;{\rm Feb}$',r'$09\;{\rm Feb}$'], rotation=45)
ax.set_xlabel("Time")
# create figure
plt.savefig('../figures/LION_MOORING_temp.png', dpi=500)
#=========================================================
# CREATE THE SALINITY PLOT
#=========================================================
f, axes = plt.subplots(3, 1, sharex=True, sharey=True)
f.set_figheight(9); f.set_figwidth(10)
f.suptitle(titleS, fontsize=16)
ax = axes[0]; ax.set_ylim(-hmax,0.); ax.set_xlim(0.,26.)
scm_evd.plot_salt_data(ax,smin,smax,colorMapS,Slevels,False)
ax = axes[1]
scm_edmf.plot_salt_data(ax,smin,smax,colorMapS,Slevels,False)
ax = axes[2]
scm_data.plot_salt_data(ax,smin,smax,colorMapS,Slevels,True)
ax.set_xticks([1,3,5,7,9,11,13,15,17,19,21,23,25])
ax.set_xticklabels([r'$16\;{\rm Jan}$',r'$18\;{\rm Jan}$',r'$20\;{\rm Jan}$',r'$22\;{\rm Jan}$',r'$24\;{\rm Jan}$',r'$26\;{\rm Jan}$',r'$28\;{\rm Jan}$',r'$30\;{\rm Jan}$',r'$01\;{\rm Feb}$',r'$03\;{\rm Feb}$',r'$05\;{\rm Feb}$',r'$07\;{\rm Feb}$',r'$09\;{\rm Feb}$'], rotation=45)
ax.set_xlabel("Time")
# create figure
plt.savefig('../figures/LION_MOORING_salt.png', dpi=500)
#=========================================================

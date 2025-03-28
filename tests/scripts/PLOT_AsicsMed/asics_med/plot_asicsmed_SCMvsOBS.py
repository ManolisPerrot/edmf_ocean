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
plt.rcParams['text.usetex'] = True

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

    def get_data(self):
        nc = Dataset(self.config["filename"] , mode='r')
        ts = self.config["time_start"]
        tf = self.config["time_end"]
        self.temp = nc.variables[self.config["temp_var_name"]][ts:tf,:]
        self.salt = nc.variables[self.config["salt_var_name"]][ts:tf,:]
        zz        = self.config["z_scaling"]*nc.variables[self.config["zr_var_name"]][:]
        tt        = self.config["time_scaling"]*nc.variables[self.config["time_var_name"]][ts:tf] - self.config["time_shift"]
        self.nz   = len(zz)
        self.nt   = len(tt)
        nc.close()
        self.dpth = np.tile( zz, (self.nt,1) )
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
evd_params = {'filename': '../data/run_asics_TKE+EVD.nc', 'time_var_name': 'ocean_time', 'data_title':'TKE+EVD'}
scm_evd = Lion_data(evd_params)
scm_evd.get_data()

edmf_energy_params = {'filename': '../data/run_asics_EDMF-Engy.nc', 'time_var_name': 'ocean_time', 'data_title':'\textbf{EDMF-Energy}'}
scm_edmf_energy = Lion_data(edmf_energy_params)
scm_edmf_energy.get_data()

edmf_params = {'filename': '../data/run_asics_EDMF.nc', 'time_var_name': 'ocean_time', 'data_title':'\textbf{EDMF-Inconsistent}'}
scm_edmf = Lion_data(edmf_params)
scm_edmf.get_data()

Keps_params = {'filename': '../data/run_asics_Keps.nc', 'time_var_name': 'ocean_time', 'data_title':'\textbf{k-epsilon}'}
scm_Keps = Lion_data(Keps_params)
scm_Keps.get_data()

data_params = {'filename': '../data/AsicsMed_data.nc', 'time_var_name': 'time_counter',
               'temp_var_name': 'votemper', 'salt_var_name': 'vosaline', 'zr_var_name': 'deptht',
               'time_start': 11, 'time_end': 11+30*2,'time_scaling':1.,'time_shift':15,'z_scaling':-1., 'data_title':'Lion observational data'}
scm_data = Lion_data(data_params)
scm_data.get_data()
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
plt.savefig('../figures/ASICS_MED_temp.png', dpi=500)
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
plt.savefig('../figures/ASICS_MED_salt.png', dpi=500)
#=========================================================
# CREATE THE PLOT for the BIAS
#=========================================================
time_b      = np.zeros((scm_data.nt,scm_evd.nz))
dpth_b      = np.zeros((scm_data.nt,scm_evd.nz))
#==============================================================
# Step1: subsample the scm data in time to get a 12h frequency
#==============================================================
indices   = np.arange(0, scm_edmf.nt, 12, dtype=int)
time_b    = scm_evd.time[indices,:];     dpth_b    = scm_evd.dpth[indices,:]
temp_evd  = scm_evd.temp[indices,:];     salt_evd  = scm_evd.salt[indices,:]
temp_edmf = scm_edmf.temp[indices,:];    salt_edmf = scm_edmf.salt[indices,:]
temp_edmf_energy = scm_edmf_energy.temp[indices,:];    salt_edmf_energy = scm_edmf_energy.salt[indices,:]
temp_Keps = scm_Keps.temp[indices,:];    salt_Keps = scm_Keps.salt[indices,:]
#=============================================================
# Step2: subsample the scm data in time to get a 12h frequency
#==============================================================
temp_data = temp_edmf.copy()
for kt in range(scm_data.nt):
    t1d     = scm_data.temp[kt,:]
    t1d_new = PchipInterpolator(-1.*scm_data.z1d, t1d)(-1.*scm_evd.z1d)
#    t1d_new = Akima1DInterpolator(-1.*scm_data.z1d, t1d)(-1.*scm_evd.z1d)
    temp_data[kt,:] = t1d_new[:]
#
salt_data = salt_edmf.copy()
for kt in range(scm_data.nt):
    t1d     = scm_data.salt[kt,:]
    t1d_new = PchipInterpolator(-1.*scm_data.z1d, t1d)(-1.*scm_evd.z1d)
#    t1d_new = Akima1DInterpolator(-1.*scm_data.z1d, t1d)(-1.*scm_evd.z1d)
    salt_data[kt,:] = t1d_new[:]
#=========================================================
# Create the plot for the temperature/salinity bias
#=========================================================
colorMapE='RdBu_r'; hmax=2300.
boundsT = np.array([-0.1, -0.075, -0.05, -0.025, 0.025, 0.05, 0.075, 0.1])
normT = colors.BoundaryNorm(boundaries=boundsT, ncolors=256)
boundsS = np.array([-0.03, -0.025, -0.02, -0.015, -0.01, -0.005, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03])
normS = colors.BoundaryNorm(boundaries=boundsS, ncolors=256)
#-------------------- Temperature
f, axes = plt.subplots(4, 1, sharex=True, sharey=True)
f.set_figheight(10); f.set_figwidth(10)
axindex = -1
#
axindex +=1;ax = axes[axindex]; ax.set_ylim(-hmax,0.); ax.set_xlim(0.,26.)
ax.set_title(r"Temperature difference [ (\textbf{TKE+EVD}) - OBS ]")
im = ax.pcolormesh(time_b, dpth_b, temp_evd-temp_data, norm=normT, cmap=colorMapE, shading='auto')
ax.set_ylabel(r'depth $[{\rm m}]$')
#
axindex +=1;ax = axes[axindex]
ax.set_title(r"Temperature difference [ \textbf{k-epsilon} - OBS ]")
im = ax.pcolormesh(time_b, dpth_b, temp_Keps-temp_data, norm=normT, cmap=colorMapE, shading='auto')
ax.set_ylabel(r'depth $[{\rm m}]$')
#
axindex +=1;ax = axes[axindex]
ax.set_title(r"Temperature difference [ (\textbf{EDMF-Inconsistent}) - OBS ]")
im = ax.pcolormesh(time_b, dpth_b, temp_edmf-temp_data, norm=normT, cmap=colorMapE, shading='auto')
ax.set_ylabel(r'depth $[{\rm m}]$')
#
axindex +=1;ax = axes[axindex]
ax.set_title(r"Temperature difference [ (\textbf{EDMF-Energy}) - OBS ]")
im = ax.pcolormesh(time_b, dpth_b, temp_edmf_energy-temp_data, norm=normT, cmap=colorMapE, shading='auto')
ax.set_ylabel(r'depth $[{\rm m}]$')
#
ax.set_xticks([1,3,5,7,9,11,13,15,17,19,21,23,25])
ax.set_xticklabels([r'$16\;{\rm Jan}$',r'$18\;{\rm Jan}$',r'$20\;{\rm Jan}$',r'$22\;{\rm Jan}$',r'$24\;{\rm Jan}$',r'$26\;{\rm Jan}$',r'$28\;{\rm Jan}$',r'$30\;{\rm Jan}$',r'$01\;{\rm Feb}$',r'$03\;{\rm Feb}$',r'$05\;{\rm Feb}$',r'$07\;{\rm Feb}$',r'$09\;{\rm Feb}$'], rotation=45)
ax.set_xlabel("Time")
f.subplots_adjust(right=0.8)
cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
cbar = f.colorbar(im, cax=cbar_ax)
cbar.set_label("K")
plt.savefig('../figures/ASICS_MED_temp_diff.png', dpi=500)
#-------------------- Salinity
f, axes = plt.subplots(4, 1, sharex=True, sharey=True)
f.set_figheight(10); f.set_figwidth(10)
axindex = -1
#
axindex +=1;ax = axes[axindex]; ax.set_ylim(-hmax,0.); ax.set_xlim(0.,26.)
ax.set_title(r"Salinity difference [ (\textbf{TKE+EVD}) - OBS ]")
im = ax.pcolormesh(time_b, dpth_b, salt_evd-salt_data, norm=normS, cmap=colorMapE, shading='auto')
ax.set_ylabel(r'depth $[{\rm m}]$')
#
axindex +=1;ax = axes[axindex]
ax.set_title(r"Salinity difference [ \textbf{k-epsilon} - OBS ]")
im = ax.pcolormesh(time_b, dpth_b, salt_Keps-salt_data, norm=normS, cmap=colorMapE, shading='auto')
ax.set_ylabel(r'depth $[{\rm m}]$')
#
axindex +=1;ax = axes[axindex]
ax.set_title(r"Salinity difference [ (\textbf{EDMF-Inconsistent}) - OBS ]")
im = ax.pcolormesh(time_b, dpth_b, salt_edmf-salt_data, norm=normS, cmap=colorMapE, shading='auto')
ax.set_ylabel(r'depth $[{\rm m}]$')
#
axindex +=1;ax = axes[axindex]
ax.set_title(r"Salinity difference [ (\textbf{EDMF-Energy}) - OBS ]")
im = ax.pcolormesh(time_b, dpth_b, salt_edmf_energy-salt_data, norm=normS, cmap=colorMapE, shading='auto')
ax.set_ylabel(r'depth $[{\rm m}]$')
#
ax.set_xticks([1,3,5,7,9,11,13,15,17,19,21,23,25])
ax.set_xticklabels([r'$16\;{\rm Jan}$',r'$18\;{\rm Jan}$',r'$20\;{\rm Jan}$',r'$22\;{\rm Jan}$',r'$24\;{\rm Jan}$',r'$26\;{\rm Jan}$',r'$28\;{\rm Jan}$',r'$30\;{\rm Jan}$',r'$01\;{\rm Feb}$',r'$03\;{\rm Feb}$',r'$05\;{\rm Feb}$',r'$07\;{\rm Feb}$',r'$09\;{\rm Feb}$'], rotation=45)
ax.set_xlabel("Time")
f.subplots_adjust(right=0.8)
cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
cbar = f.colorbar(im, cax=cbar_ax)
cbar.set_label("psu")
plt.savefig('../figures/ASICS_MED_salt_diff.png', dpi=500)
#=========================================================

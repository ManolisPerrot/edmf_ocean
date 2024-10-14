#!/usr/bin/env python
# coding: utf-8
__docformat__ = 'reStructuredText'
###########################################
# Imports
###########################################
from sys import exit
#import add
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl
from netCDF4 import Dataset
import xarray as xr
###########################################
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
###########################################
#===========================================================================
# read data from gotm outputs
#===========================================================================
nc_file = 'wang1_kpp_LR.nc'
fh      = Dataset(nc_file, mode='r')
temp1   = fh.variables['temp'][-1,:,0,0]
akt1    = fh.variables['nuh'][-1,:,0,0]
gam1    = fh.variables['gamh'][-1,:,0,0]
zr1     = fh.variables['z'][0,:,0,0]
zw1     = fh.variables['zi'][0,:,0,0]
fh.close()
#===========================================================================
nz = len(zr1); wt1 = np.zeros(nz-1); dz = zw1[1]-zw1[0]
wt1 = -akt1[1:nz-1]*(temp1[1:nz-1]-temp1[0:nz-2])/dz + gam1[1:nz-1]
#===========================================================================
nc_file  = 'wang1_keps_LR.nc'
fh       = Dataset(nc_file, mode='r')
temp2    = fh.variables['temp'][-1,:,0,0]
tke2     = fh.variables['tke'][-1,:,0,0]
akt2     = fh.variables['nuh'][-1,:,0,0]
zr2      = fh.variables['z'][0,:,0,0]
zw2      = fh.variables['zi'][0,:,0,0]
fh.close()
#===========================================================================
nz = len(zr2); wt2 = np.zeros(nz-1); dz = zw2[1]-zw2[0]
wt2 = -akt2[1:nz-1]*(temp2[1:nz-1]-temp2[0:nz-2])/dz
#===========================================================================
cases = ['WANG1_FR_lat60']
les = {}
LG_MEAN={}
LG_RES={}
LG_SBG={}
BU_KE={}
TH_les={}
U_les={}
V_les={}
WTH={}
WU={}
WV={}
TKE={}
WTKE={}
for case in cases:
    file = 'GN_01.1.OC_01.000.nc'
    path = '../../tests/data/'+case+'/'
    les[case] = xr.open_dataset(path+file)
    LG_MEAN[case]= xr.open_dataset(path+file,group ='/LES_budgets/Mean/Cartesian/Not_time_averaged/Not_normalized/cart')
    LG_RES[case] = xr.open_dataset(path+file,group ='/LES_budgets/Resolved/Cartesian/Not_time_averaged/Not_normalized/cart')    
    LG_SBG[case] = xr.open_dataset(path+file,group ='/LES_budgets/Subgrid/Cartesian/Not_time_averaged/Not_normalized/cart')
    BU_KE[case] = xr.open_dataset(path+file,group ='/LES_budgets/BU_KE/Cartesian/Not_time_averaged/Not_normalized/cart')
    TH_les[case] = (LG_MEAN[case].MEAN_TH - 273.15).data
    U_les [case] = (LG_MEAN[case].MEAN_U).data
    V_les [case] = (LG_MEAN[case].MEAN_V).data
    WTH   [case] = (LG_RES[case].RES_WTH + LG_SBG[case].SBG_WTHL).data
    WU    [case] = (LG_RES[case].RES_WU  + LG_SBG[case].SBG_WU).data
    WV    [case] = (LG_RES[case].RES_WV  + LG_SBG[case].SBG_WV).data
    TKE   [case] = (LG_RES[case].RES_KE  + LG_SBG[case].SBG_TKE).data 
    WTKE  [case] = (LG_RES[case].RES_WKE + LG_SBG[case].SBG_WTKE).data
z_r_les = (les[case].level_les - (les[case].level_les[0] + les[case].level_les[-1])).data
#===========================================================================
edmf = xr.open_dataset('scm_WANG1_FR_lat60.nc')
evd  = xr.open_dataset('scm_WANG1_EVD.nc')
#===========================================================================
fig, axes = plt.subplots(1, 3, sharex=False, sharey=True)
fig.set_figheight(4)
fig.set_figwidth(10)
#===========================================================================
ax = axes[0]
ax.set_title(r"$\overline{\theta}$",fontsize=14)
xmin = 2.98; xmax = 2.99
ax.set_xlim(xmin,xmax)
ax.plot(  temp1, zr1,linewidth=2, color='k',label='KPP')
ax.plot(  temp2, zr2,linewidth=2, color='b',label=r'$k$-$\varepsilon$')
ax.plot(  TH_les[case][-1], z_r_les	    ,label='LES')
ax.plot(  edmf['temp'][-1], edmf['z_r']	  ,label='edmf')
ax.plot(  evd['temp'][-1], edmf['z_r']	  ,label='evd')

ax.legend(loc=4,fontsize=10)
#===========================================================================
ax = axes[1]
ax.set_title(r"$\overline{w'\theta'}$",fontsize=14)
ax.plot(  wt1, zw1[1:nz-1],linewidth=2, color='k',label='KPP')
ax.plot(  wt2, zw2[1:nz-1],linewidth=2, color='b',label=r'$k$-$\varepsilon$')
ax.plot(  WTH[case][-1], z_r_les     	          ,label='LES')
ax.plot(  -edmf['WT'][-1], edmf['z_w']	  ,label='edmf')
ax.plot(  -evd['WT'][-1], edmf['z_w']	  ,label='evd')

#===========================================================================
ax = axes[2]
ax.set_title(r"$k$",fontsize=14)
ax.plot( tke2, zw2,linewidth=2, color='b')
ax.plot(  TKE[case][-1], z_r_les	    ,label='LES')
ax.plot(  edmf['tke'][-1], edmf['z_w']	  ,label='edmf')
ax.plot(  evd['tke'][-1], edmf['z_w']	  ,label='evd')

#========================================================================================================================
plt.tight_layout()
#===============================
plt.show()
##############################################

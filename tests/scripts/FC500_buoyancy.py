#!/usr/bin/env python
# coding: utf-8
__docformat__ = 'reStructuredText'
###########################################
# Imports
###########################################


from sys import exit
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import scipy.signal
from scipy.interpolate import interp1d
import xarray as xr
import time as TIME
import xrft as xrft
from matplotlib.ticker import MaxNLocator
from cases_exhaustive_settings import case_settings

plt.rcParams['text.usetex'] = True
#plt.rcParams.update({'font.size': 18})
plt.rcParams.update({'figure.facecolor':'white'})
plt.rcParams.update({'savefig.facecolor':'white'})
###########################################



cases=['FC500']

case='FC500'
alpha = case_settings[case]['alpha']
cp = case_settings[case]['cpoce']
rho0 = case_settings[case]['rho0']

g = 9.81

B0 = -g*alpha/(rho0*cp)*case_settings[case]['stflx']


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
dz_WU_les ={}
dz_WV_les ={}
time_les={}
hours = {}
z_r_les = {}
mld = {}
B_les = {}
W2_les = {}
B2_les = {}
WB = {}
instant=-1
for case in cases:
    file = 'GN_01.1.OC_01.000_copy.nc'
    path = '../../tests/data/'+case+'/'
    les[case] = xr.open_dataset(path+file)
    print('Opening LES case: ',case)
    LG_MEAN[case]= xr.open_dataset(path+file,group ='/LES_budgets/Mean/Cartesian/Not_time_averaged/Not_normalized/cart')
    LG_RES[case] = xr.open_dataset(path+file,group ='/LES_budgets/Resolved/Cartesian/Not_time_averaged/Not_normalized/cart')    
    LG_SBG[case] = xr.open_dataset(path+file,group ='/LES_budgets/Subgrid/Cartesian/Not_time_averaged/Not_normalized/cart')
    BU_KE[case] = xr.open_dataset(path+file,group ='/LES_budgets/BU_KE/Cartesian/Not_time_averaged/Not_normalized/cart')
    TH_les[case] = (LG_MEAN[case].MEAN_TH - 273.15)
    B_les [case] = g*alpha*TH_les[case]
    U_les [case] = (LG_MEAN[case].MEAN_U)
    V_les [case] = (LG_MEAN[case].MEAN_V)
    # V_les [case] = (LG_MEAN[case].MEAN_U)
    # U_les [case] = (LG_MEAN[case].MEAN_V)
    WTH   [case] = (LG_RES[case].RES_WTH + LG_SBG[case].SBG_WTHL)
    WB   [case] = g*alpha*WTH[case]
    B2_les[case] = (g*alpha)**2*(LG_RES[case]['RES_TH2']+LG_SBG[case]['SBG_THL2'])
    W2_les[case] = (LG_RES[case].RES_W2+LG_SBG[case]['SBG_W2'])
    WU    [case] = (LG_RES[case].RES_WU  + LG_SBG[case].SBG_WU)
    WV    [case] = (LG_RES[case].RES_WV  + LG_SBG[case].SBG_WV)
    # WV    [case] = (LG_RES[case].RES_WU  + LG_SBG[case].SBG_WU)
    # WU    [case] = (LG_RES[case].RES_WV  + LG_SBG[case].SBG_WV)
    TKE   [case] = (LG_RES[case].RES_KE  + LG_SBG[case].SBG_TKE)
    WTKE  [case] = (LG_RES[case].RES_WKE + LG_SBG[case].SBG_WTKE)

    # remap level_les on negative depth values
    z_r_les[case] = -(les[case].level_les - (les[case].level_les[0] + les[case].level_les[-1])).data
    

    mld[case] = (z_r_les[case][(-WTH[case]).argmax(axis=1)])


wstar = (B0*mld[case][instant])**(1/3)
bstar = B0/wstar

#========

################################# Styles #################################
#SCM
styles = {case: '-' for case in cases}
# cocolors=['tab:gray' ,'tab:orange' , 'tab:green' , 'tab:blue']
cocolors=['tab:blue','tab:orange' , 'tab:green' , 'tab:gray']
colors = {case: cocolors[i] for i,case in enumerate(cases)}
alpha  = {case: 1 for case in cases}
# linewidth = [2]*(len(run_label))
#LES
styles_les = {case: 'k+' for case in cases}
colors_les = colors
# colors_les = {case: 'k' }
alpha_les  = {case:1 for case in cases}
marker_les=''
width_les = {case:2 for case in cases}
linestyle_les=':'
################################ Ploting ##################################
fig, axes = plt.subplots(nrows=2, ncols=3, sharex=False,
                        sharey=True, constrained_layout=True)

ax_index=-1


# ===============================================================
ax_index+=1
ax = axes.flat[ax_index]
ax.set_ylabel(r'$z / h $')
ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

# ===============================================================
ax_index+=1
ax = axes.flat[ax_index]
ax.set_xlabel(r'$ (\overline{b}-\overline{b}_{\mathrm{min}}) / b_*$')
ax.set_ylabel(r'$z / h $')

for case in cases:  
    ax.plot(-(B_les[case][instant,:] - B_les[case][instant,:].max())/bstar ,z_r_les[case]/mld[case][instant])
    ax.plot(-(B_les[case][0,:] - B_les[case][instant,:].max())/bstar ,z_r_les[case]/mld[case][instant],'k--', alpha=0.5)   


ax.set_xlim((-5,50))
# ===============================================================
ax_index+=1
ax = axes.flat[ax_index]

ax.set_xlabel(r'$\overline{w^\prime b^\prime}/B_0$')

for case in cases:
    ax.plot(WB[case][instant,:]/B0, z_r_les[case]/mld[case][instant])
    
# ===============================================================
# ===============================================================
ax_index+=1
ax = axes.flat[ax_index]
ax.set_xlabel(r'$\overline{b^{\prime 2}}/b_*^2$')

for case in cases:
    ax.plot(B2_les[case][instant,:]/bstar**2, z_r_les[case]/mld[case][instant])

# ===============================================================
ax_index+=1
ax = axes.flat[ax_index]
ax.set_xlabel(r'$\overline{w^{\prime 2}}/w_*^2$')

for case in cases:
    ax.plot(W2_les[case][instant,:]/wstar**2, z_r_les[case]/mld[case][instant])

# ===============================================================
# ===============================================================
ax_index+=1
ax = axes.flat[ax_index]
ax.remove()
# ===============================================================
# ===============================================================

# adding subplot labels
subplot_label = ['',r'\rm{(a)}', r'\rm{(b)}', r'\rm{(c)}',
                r'\rm{(d)}', r'\rm{(e)}', r'\rm{(f)}',r'\rm{(g)}',r'\rm{(h)}',r'\rm{(i)}',r'\rm{(j)}',r'\rm{(k)}',r'\rm{(l)}']

for i,ax in enumerate(axes.flat):
    # ax.set_ylim((-2, 0))
    ax.set_box_aspect(1)
    ax.text(0.96, 0.12, subplot_label[i], transform=ax.transAxes, bbox=dict(facecolor='1.', edgecolor='none'), fontweight='bold', va='top', ha='right')
    # ax.grid()
    ax.axvline(0,   color='tab:grey', lw=0.5)
    ax.axhline(0.82, color='tab:grey', lw=0.5)
    # ax.axhline(1,   color='tab:grey', lw=0.5)
    ax.axhline(1.15,color='tab:grey', lw=0.5)
    ax.axhline(0.1, color='tab:grey', lw=0.5)
    ax.set_ylim((0,1.3))
plt.savefig('../figures/FC500_buoyancy.png', dpi=600, bbox_inches='tight')

# plt.show()

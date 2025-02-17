#!/usr/bin/env python
# coding: utf-8
__docformat__ = 'reStructuredText'
###########################################
# Imports
###########################################
from sys import exit
import sys
sys.path.insert(0,"../../library/F2PY/")
#import add
import numpy as np
import matplotlib
#matplotlib.use('TkAgg')
#matplotlib.use('pdf')
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scm_class_oce import SCM
#from F2PY import scm_class
import scipy.signal
###########################################
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
###########################################
case        = 'W005_FC500'
saving_path = '../figures/'
saving_name = case+'72h_profile_W005_FC500.png'
###########################################
end_time   = 72
#===========================================================================
compar_les = 1
if compar_les==1:
  les_dir    = "../data/W005_C500_NO_COR/"
  les_file   = les_dir+"GN_01.1.OC_01.000_copy.nc"
  mean_var   =  '/LES_budgets/Mean/Cartesian/Not_time_averaged/Not_normalized/cart'
  subg_var   =  '/LES_budgets/Subgrid/Cartesian/Not_time_averaged/Not_normalized/cart'
  res_var    =  '/LES_budgets/Resolved/Cartesian/Not_time_averaged/Not_normalized/cart'
  ## read temperature reference data from LES results
  fh1       = Dataset(les_file, mode='r')
  indx_les  = end_time-1
  temp_les  = fh1[mean_var].variables['MEAN_TH'][indx_les,:]-273.15
  uu_les    = fh1[mean_var].variables['MEAN_U'][indx_les,:]
  wt_les    = fh1[subg_var].variables['SBG_WTHL'][indx_les,:]
  wu_les    = fh1[subg_var].variables['SBG_WU'][indx_les,:]
  tke1_les  = fh1[res_var].variables['RES_KE'][indx_les,:]
  tke2_les  = fh1[subg_var].variables['SBG_TKE'][indx_les,:]
  we1_les   = fh1[subg_var].variables['SBG_WTKE'][indx_les,:]
  we_les    = fh1[res_var].variables['RES_WKE'][indx_les,:]
  wt_les    = wt_les + fh1[res_var].variables['RES_WTH'][indx_les,:]
  wu_les    = wu_les + fh1[res_var].variables['RES_WU'][indx_les,:]
  level_les = fh1.variables['level_les'][:]
  fh1.close()
  z_r_les = (level_les[:] - (level_les[0] + level_les[-1]))
  nz      = len(z_r_les)
  hinv_les = - z_r_les[np.argmin(wt_les)]
#===========================================================================
# Define parameters specific to each run (overwrite common parameters):
#run_label = ['TKE+EVD', 'K-eps', 'EDMF', 'EDMF-Engy']

config_params = { 'nz': 100, 'lat0':0.0, 'sustr': 0.05/1027.0, 'stflx': -500., 'T0': 2., 'N0': 1.9620001275490499e-6}

run_label = ['TKE+EVD', 'EDMF', 'EDMF-Engy']
runs = [
    {
        'eddy_diff': True,
        'evd': True,
        'mass_flux_tra': False,
        'mass_flux_dyn': False,
        'mass_flux_tke': False,
        'mass_flux_tke_trplCorr': False,
        'mass_flux_small_ap': False,
        'eddy_diff_scheme' : 'TKE',
        'output_filename': 'run1.nc'
    },
    {
        'eddy_diff': True,
        'evd': False,
        'mass_flux_tra': True,
        'mass_flux_dyn': True,
        'mass_flux_tke': False,
        'mass_flux_tke_trplCorr': False,
        'mass_flux_small_ap': False,
        'entr_scheme': 'R10',
        'eddy_diff_scheme' : 'TKE',
        'output_filename': 'run2.nc'
    },
    {
        'eddy_diff': True,
        'evd': False,
        'mass_flux_tra': True,
        'mass_flux_dyn': True,
        'mass_flux_tke': True,
        'mass_flux_tke_trplCorr': True,
        'mass_flux_small_ap': False,
        'entr_scheme': 'R10',
        'eddy_diff_scheme' : 'TKE',
        'output_filename': 'run3.nc'
    }
        ]
#
scm = [0]*len(runs)
# Run the SCM
for i, run_params in enumerate(runs):
    
    scm[i] = SCM(run_params, config_params)
    scm[i].run_direct()
    if scm[i].MF_tra or scm[i].MF_dyn: print("zinv =", scm[i].zinv)
#===========================================================================
vscale = 1./325.
if compar_les==1: vscale = 1./hinv_les
zmax   = -1.2
styles = ['solid', 'dashed', 'dashed', 'solid']
colors = ['0.4', '0.4', 'k', 'k']
linewidth = [2.0, 2.0, 3, 3]
#===========================================================================
fig, axes = plt.subplots(nrows=2, ncols=3, sharex=False,sharey=True, figsize=(10, 8) )
#===========================================================================
ax = axes[0,0]  #== THETA
ax.set_title(r"$\overline{\theta}$",fontsize=18)
for i, label in enumerate(run_label):
  ax.plot(scm[i].t_np1[:, 0], vscale*scm[i].z_r, linestyle=styles[i], color = colors[i],linewidth=linewidth[i], label=label)
if compar_les==1:
  ax.plot(  temp_les, z_r_les*vscale, 'r', linewidth=2, label = 'LES' )
#==
xmin = 1.55; xmax = 1.8
ax.set_xlim(xmin,xmax)
ax.set_ylim(zmax, 0)
ax.set_xlabel(r"$^{\circ}{\rm C}$",fontsize=18)
ax.set_ylabel(r"$z/h$",fontsize=18)
ax.legend(loc=2,fontsize=12)
#============================================================================
ax = axes[1,0]
ax.set_title(r"$\overline{w^{'}\theta^{'}}$",fontsize=18)
for i, label in enumerate(run_label):
  ax.plot(-scm[i].wtmf[:]-scm[i].wted[:],vscale*scm[i].z_w, linestyle=styles[i], color = colors[i],linewidth=linewidth[i])
if compar_les==1:
  ax.plot(  wt_les[:], z_r_les*vscale, 'r', linewidth=2)
#==
ax.set_xlabel(r"$^{\circ}{\rm C}\;{\rm m}\;{\rm s}^{-1}$",fontsize=14)
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#============================================================================
ax = axes[0,2]
ax.set_title(r"$k$",fontsize=18)
for i, label in enumerate(run_label):
  ax.plot(scm[i].tke_np1[:],vscale*scm[i].z_w, linestyle=styles[i], color = colors[i],linewidth=linewidth[i])
if compar_les==1:
  ax.plot(  tke1_les[:]+tke2_les-0.000075,z_r_les*vscale, linewidth=2, color='r'  )
#==
ax.set_xlabel(r"${\rm m}^2\;{\rm s}^{-2}$",fontsize=18)
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#============================================================================
ax = axes[1,2]
ax.set_title(r"$\overline{w' \mathbf{u}'\cdot \mathbf{u}'/2}$",fontsize=18)
for i, label in enumerate(run_label):
  ax.plot(scm[i].wtke[:],vscale*scm[i].z_r, linestyle=styles[i], color = colors[i],linewidth=linewidth[i])
if compar_les==1:
  ax.plot(  we1_les[:]+we_les[:],z_r_les*vscale, linewidth=2, color='r'  )
#==
ax.set_xlabel(r"${\rm m}^3\;{\rm s}^{-3}$",fontsize=18)
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#============================================================================
ax = axes[0,1]
ax.set_title(r"$\overline{u}$",fontsize=18)
xmin=-2e-5; xmax = 0.1
ax.set_xlim(xmin,xmax)
for i, label in enumerate(run_label):
  ax.plot(scm[i].u_np1[:],vscale*scm[i].z_r, linestyle=styles[i], color = colors[i],linewidth=linewidth[i])
if compar_les==1:
  ax.plot(  uu_les[:],z_r_les*vscale, linewidth=2, color='r'  )
#==
ax.set_xlabel(r"${\rm m}\;{\rm s}^{-1}$",fontsize=18)
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#===========================================================================
ax = axes[1,1]
ax.set_title(r"$\overline{w^{'}u^{'}}$",fontsize=18)
for i, label in enumerate(run_label):
  ax.plot((scm[i].wued[:]+scm[i].wumf[:]),vscale*scm[i].z_w, linestyle=styles[i], color = colors[i],linewidth=linewidth[i])
if compar_les==1:
  ax.plot( -wu_les , z_r_les*vscale, 'r', linewidth=2)
#==
ax.set_xlabel(r"${\rm m}^2\;{\rm s}^{-2}$",fontsize=18)
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#ax = axes[4]
#ax.set_title(r"$K_\theta$",fontsize=14)
#for i, label in enumerate(run_label):
#  ax.plot(scm[i].akt[:],vscale*scm[i].z_w, linestyle=styles[i], color = colors[i],linewidth=linewidth[i])
#==
#ax.set_xlabel(r"${\rm m}^2\;{\rm s}^{-1}$",fontsize=14)
#ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#===============================
plt.tight_layout()
#plt.savefig('W005_FC500_case.pdf', bbox_inches='tight')
#===============================
#plt.subplot_tool()
plt.show()

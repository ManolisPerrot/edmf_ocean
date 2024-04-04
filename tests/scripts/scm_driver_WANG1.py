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
from scm_class import SCM
#from F2PY import scm_class
import scipy.signal
###########################################
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
###########################################
case        = 'WANG1_FR'
saving_path = '../figures/'
saving_name = case+'72h_profile_WANG1.png'
###########################################
end_time   = 276
#===========================================================================
compar_les = 1
if compar_les==1:
  les_dir    = "/mnt/c/Users/flori/OneDrive/Bureau/Codes/kets_and_edmf/3d_runs/WANG1_FR/"
  les_file   = les_dir+"GN_01.1.OC_01.000.nc"
  mean_var   =  '/LES_budgets/Mean/Cartesian/Not_time_averaged/Not_normalized/cart'
  subg_var   =  '/LES_budgets/Subgrid/Cartesian/Not_time_averaged/Not_normalized/cart'
  res_var    =  '/LES_budgets/Resolved/Cartesian/Not_time_averaged/Not_normalized/cart'
  ## read temperature reference data from LES results
  fh1       = Dataset(les_file, mode='r')
  indx_les  = end_time-1
  temp_les  = fh1[mean_var].variables['MEAN_TH'][indx_les,:]-273.15
  uu_les    = fh1[mean_var].variables['MEAN_U'][indx_les,:]
  vv_les    = fh1[mean_var].variables['MEAN_V'][indx_les,:]
  wt_les    = fh1[subg_var].variables['SBG_WTHL'][indx_les,:]
  wu_les    = fh1[subg_var].variables['SBG_WU'][indx_les,:]
  wv_les    = fh1[subg_var].variables['SBG_WV'][indx_les,:]
  tke1_les  = fh1[res_var].variables['RES_KE'][indx_les,:]
  tke2_les  = fh1[subg_var].variables['SBG_TKE'][indx_les,:]
  we1_les   = fh1[subg_var].variables['SBG_WTKE'][indx_les,:]
  we_les    = fh1[res_var].variables['RES_WKE'][indx_les,:]
  wt_les    = wt_les + fh1[res_var].variables['RES_WTH'][indx_les,:]
  wu_les    = wu_les + fh1[res_var].variables['RES_WU'][indx_les,:]
  wv_les    = wv_les + fh1[res_var].variables['RES_WV'][indx_les,:]
  level_les = fh1.variables['level_les'][:]
  fh1.close()
  z_r_les = (level_les[:] - (level_les[0] + level_les[-1]))
  nz      = len(z_r_les)
  hinv_les = - z_r_les[np.argmin(wt_les)]
#  hinv_les = - zw_les[np.argmin(wt_les)]
#===========================================================================
common_params = {
    'nz': 450,
    'dt': 360.,
    'h0': 4500.,
    'thetas': 6.5,
    'hc': 1000000000000.,
    'nbhours': end_time,
    'outfreq': 1,
    'T0': 3.,
    'N0': 1.5865490613891073e-08,
    'cpoce': 3985.,
    'mld_ini': -1000.,
    'mld_iniS': 0.,
    'Tcoef': 0.2054,
    'Scoef': 0.8216,
    'Tref': 3.,
    'Sref': 32.6,
    'SaltCst': 32.6,
    'lat0': 60.,
    'svstr': 0.,
    'sustr': 0.,
    'stflx': -111.12982,
    'srflx': 0.,
    'ssflx': 0.,
    'diurnal_cycle': False,
    'btflx': 'no flux',
    'eddy_diff': True,
    'evd': False,
    'mass_flux_tra': False,
    'mass_flux_dyn': False,
    'mass_flux_tke': False,
    'mass_flux_tke_trplCorr': False,
    'mass_flux_small_ap': True,
    'lin_eos': True,
    'extrap_ak_surf': True,
    'tke_sfc_dirichlet': False,
    'eddy_diff_tke_const': 'NEMO',
    'entr_scheme': 'R10',
    'Cent': 0.9,
    'Cdet': 1.95,       # 'Cdet': 2.5,
    'wp_a': 1.,
    'wp_b': 1.,      # 'wp_b': 1.
    'wp_bp': 0.005*250,   #*300,     #      0.002,
    'up_c': 0.25,
    'vp_c': 0.25,
    'bc_ap': 0.2,    #0.3,
    'delta_bkg': 0.5*0.005*250,   # 0.006,
    'wp0': -1.0e-08,
    'akvmin': 1.e-4,
    'aktmin': 1.e-5,
    'mxlmin': 1.0,
    'output_filename': 'run',
    'write_netcdf': True,
    'avg_last_hour': False,
}
#===========================================================================
# Define parameters specific to each run (overwrite common parameters):
run_label = ['ED', 'EDMF-Engy', 'EDMF-Engy-Cor']
runs = [
    {
        'eddy_diff': True,
        'evd': False,
        'mass_flux_tra': False,
        'mass_flux_dyn': False,
        'mass_flux_tke': False,
        'mass_flux_tke_trplCorr': False,
        'output_filename': 'run1.nc'
    },
    {
        'eddy_diff': True,
        'evd': False,
        'mass_flux_tra': True,
        'mass_flux_dyn': True,
        'mass_flux_tke': True,
        'mass_flux_tke_trplCorr': True,
        'entr_scheme': 'R10',
        'output_filename': 'run2.nc'
    },
        {
        'eddy_diff': True,
        'evd': False,
        'mass_flux_tra': True,
        'mass_flux_dyn': True,
        'mass_flux_tke': True,
        'mass_flux_tke_trplCorr': True,
        'entr_scheme': 'R10corNT',
        'output_filename': 'run3.nc'
    }
        ]
#
scm = [0]*len(runs)
# Run the SCM
for i, run_params in enumerate(runs):
    params = common_params.copy()  # Create a copy of common_params
    params.update(run_params)  # Update with run_params
    scm[i] = SCM(**params)
    scm[i].run_direct()
    print("zinv =", scm[i].zinv)
#===========================================================================
vscale = 1./3750.
if compar_les==1: vscale = 1./hinv_les
zmax   = -1.2
styles = ['solid', 'dashed', 'solid']
colors = ['0.5', 'k', 'k']
linewidth = [3]*(len(run_label))
#===========================================================================
fig, axes = plt.subplots(nrows=2, ncols=3, sharex=False,sharey=True, figsize=(10, 8) )
#===========================================================================
ax = axes[0,0]
ax.set_title(r"$\overline{\theta}$",fontsize=14)
for i, label in enumerate(run_label):
  ax.plot(scm[i].t_np1[:, 0], vscale*scm[i].z_r, linestyle=styles[i], color = colors[i],linewidth=linewidth[i], label=label)
if compar_les==1:
  ax.plot(  temp_les, z_r_les*vscale, 'r', linewidth=2, label = 'LES' )
#==
xmin = 2.95; xmax = 2.99
ax.set_xlim(xmin,xmax)
ax.set_ylim(zmax, 0)
ax.set_xlabel(r"$^{\circ}{\rm C}$",fontsize=14)
ax.set_ylabel(r"$z/h$",fontsize=14)
ax.legend(loc=2,fontsize=10)
#============================================================================
ax = axes[0,1]
ax.set_title(r"$\overline{w^{'}\theta^{'}}$",fontsize=14)
for i, label in enumerate(run_label):
  ax.plot(-scm[i].wtmf[:]-scm[i].wted[:],vscale*scm[i].z_w, linestyle=styles[i], color = colors[i],linewidth=linewidth[i])
if compar_les==1:
  ax.plot(  wt_les[:], z_r_les*vscale, 'r', linewidth=2)
#==
ax.set_xlabel(r"$^{\circ}{\rm C}\;{\rm m}\;{\rm s}^{-1}$",fontsize=14)
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#============================================================================
ax = axes[0,2]
ax.set_title(r"$k$",fontsize=14)
for i, label in enumerate(run_label):
  ax.plot(scm[i].tke_np1[:],vscale*scm[i].z_w, linestyle=styles[i], color = colors[i],linewidth=linewidth[i])
if compar_les==1:
  ax.plot(  tke1_les[:]+tke2_les,z_r_les*vscale, linewidth=2, color='r'  )
#==
ax.set_xlabel(r"${\rm m}^2\;{\rm s}^{-2}$",fontsize=14)
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#============================================================================
ax = axes[1,0]
ax.set_title(r"$\overline{u}$",fontsize=14)
for i, label in enumerate(run_label):
  ax.plot(scm[i].u_np1[:],vscale*scm[i].z_r, linestyle=styles[i], color = colors[i],linewidth=linewidth[i])
if compar_les==1:
  ax.plot(  uu_les, z_r_les*vscale, 'r', linewidth=1)
#==
ax.set_ylabel(r"$z/h$",fontsize=14)
ax.set_xlabel(r"${\rm m}^2\;{\rm s}^{-2}$",fontsize=14)
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#============================================================================
ax = axes[1,1]
ax.set_title(r"$\overline{v}$",fontsize=14)
for i, label in enumerate(run_label):
  ax.plot(scm[i].v_np1[:],vscale*scm[i].z_r, linestyle=styles[i], color = colors[i],linewidth=linewidth[i])
if compar_les==1:
  ax.plot(  vv_les, z_r_les*vscale, 'r', linewidth=1)
#==
ax.set_xlabel(r"${\rm m}\;{\rm s}^{-1}$",fontsize=14)
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#============================================================================
ax = axes[1,2]
ax.set_title(r"$\overline{w' \mathbf{u}'\cdot \mathbf{u}'/2}$",fontsize=14)
for i, label in enumerate(run_label):
  ax.plot(scm[i].wtke[:],vscale*scm[i].z_r, linestyle=styles[i], color = colors[i],linewidth=linewidth[i])
if compar_les==1:
  ax.plot(  we1_les[:]+we_les[:],z_r_les*vscale, linewidth=2, color='r'  )
#==
ax.set_xlabel(r"${\rm m}^3\;{\rm s}^{-3}$",fontsize=14)
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#===============================
plt.tight_layout()
#plt.savefig('convect_VR18_avg.pdf', bbox_inches='tight')
#===============================
plt.show()

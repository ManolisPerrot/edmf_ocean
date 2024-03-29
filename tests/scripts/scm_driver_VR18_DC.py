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
case        = 'Van_Roekel_2018_DC'
saving_path = '../figures/'
saving_name = case+'72h_profile_VR18_DC.png'
###########################################
end_time   = 72
#===========================================================================
compar_les = 0
if compar_les==1:
  les_dir    = "/mnt/c/Users/flori/Data/lesNetcdf/"
  les_file   = les_dir+"noMLdiurnal_av.nc"
  ## read temperature reference data from LES results
  fh1       = Dataset(les_file, mode='r')
  indx_les  = end_time-1
  temp_les  = fh1.variables['temp'][indx_les,:]  # zm
  salt_les  = fh1.variables['salt'][indx_les,:]  # zm
  wt_les    = fh1.variables['wt'][indx_les,:]    # zt
  tke_les   = 0.5*(fh1.variables['uu'][indx_les,:] +
                 fh1.variables['vv'][indx_les,:] +
                 fh1.variables['ww'][indx_les,:]) # zt
  zr_les = fh1.variables['zm'][:]
  zw_les = fh1.variables['zt'][:]
  fh1.close()
#  hinv_les = - zw_les[np.argmin(wt_les)]
#===========================================================================
common_params = {
    'nz': 50,
    'dt': 60.,
    'h0': 150.,
    'thetas': 6.5,
    'hc':     75.,
    'nbhours': 72,
    'outfreq': 1,
    'output_filename': "scm_output.nc",
    'T0': 20.,
    'N0': 1.967748046875e-5,
    'cpoce': 4200.,
    'mld_ini': 0.,
    'mld_iniS': 0.,
    'Tcoef': 0.2054,
    'Scoef': 0.8216,
    'Tref': 25.,
    'Sref': 35.,
    'SaltCst': 35.,
    'lat0': 43.2885,
    'svstr': 0.,
    'sustr': 0.,
    'stflx': -75.,
    'srflx': 300.,
    'ssflx': 1.37*35. / (1000.*86400.),
    'diurnal_cycle': True,
    'btflx': 'linear_continuation',
    'eddy_diff': True,
    'evd': False,
    'mass_flux_tra': False,
    'mass_flux_dyn': False,
    'mass_flux_tke': False,
    'mass_flux_tke_trplCorr': False,
    'mass_flux_small_ap': False,
    'lin_eos': True,
    'extrap_ak_surf': True,
    'tke_sfc_dirichlet': False,
    'eddy_diff_tke_const': 'NEMO',
    'entr_scheme': 'R10',
    'Cent': 0.99,
    'Cdet': 1.99,       # 'Cdet': 2.5,
    'wp_a': 1.,
    'wp_b': 1.,      # 'wp_b': 1.
    'wp_bp': 0.005*300,     #      0.002,
    'up_c': 0.25,
    'vp_c': 0.25,
    'bc_ap': 0.2,    #0.3,
    'delta_bkg': 0.5*0.005*300,   # 0.006,
    'wp0': -1.0e-08,
    'akvmin': 1.e-5,
    'aktmin': 1.e-6,
    'mxlmin': 5.0,
    'output_filename': 'run',
    'write_netcdf': True,
    'avg_last_hour': True,
}
#===========================================================================
# Define parameters specific to each run (overwrite common parameters):
run_label = ['ED', 'EDMF', 'EDMF-Energy']
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
        'mass_flux_tke': False,
        'mass_flux_tke_trplCorr': False,
        'output_filename': 'run2.nc'
    },
        {
        'eddy_diff': True,
        'evd': False,
        'mass_flux_tra': True,
        'mass_flux_dyn': True,
        'mass_flux_tke': True,
        'mass_flux_tke_trplCorr': True,
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
vscale = 1./30.
zmax   = -1.2
styles = ['solid', 'dashed', 'solid']
colors = ['0.5', 'k', 'k']
linewidth = [3]*(len(run_label))
#===========================================================================
fig, axes = plt.subplots(nrows=1, ncols=4, sharex=False,sharey=True, figsize=(10, 4) )
#===========================================================================
ax = axes[0]
ax.set_title(r"$\overline{\theta}$",fontsize=14)
for i, label in enumerate(run_label):
  ax.plot(scm[i].t_avg[:, 0], vscale*scm[i].z_r, linestyle=styles[i], color = colors[i],linewidth=linewidth[i], label=label)
if compar_les==1:
  ax.plot(  temp_les, zr_les*vscale, 'r', linewidth=2, label = 'LES' )
#==
xmin = 19.6; xmax = 19.95
ax.set_xlim(xmin,xmax)
ax.set_ylim(zmax, 0)
ax.set_xlabel(r"$^{\circ}{\rm C}$",fontsize=14)
ax.set_ylabel(r"$z/h$",fontsize=14)
ax.legend(loc=2,fontsize=10)
#===========================================================================
ax = axes[1]
ax.set_title(r"$\overline{S}$",fontsize=14)
for i, label in enumerate(run_label):
  ax.plot(scm[i].t_avg[:, 1], vscale*scm[i].z_r, linestyle=styles[i], color = colors[i],linewidth=linewidth[i])
if compar_les==1:
  ax.plot(  salt_les, zr_les*vscale, 'r', linewidth=2)
#==
xmin = 34.995; xmax = 35.02
ax.set_xlim(xmin,xmax)
ax.set_ylim(zmax, 0)
ax.set_xlabel(r"${\rm psu}$",fontsize=14)
#============================================================================
ax = axes[2]
ax.set_title(r"$\overline{w^{'}\theta^{'}}$",fontsize=14)
for i, label in enumerate(run_label):
  ax.plot(scm[i].wt_avg[:],vscale*scm[i].z_w, linestyle=styles[i], color = colors[i],linewidth=linewidth[i])
if compar_les==1:
  ax.plot(  wt_les[:], zw_les*vscale, 'r', linewidth=2)
#==
ax.set_xlabel(r"$^{\circ}{\rm C}\;{\rm m}\;{\rm s}^{-1}$",fontsize=14)
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#============================================================================
ax = axes[3]
ax.set_title(r"$k$",fontsize=14)
for i, label in enumerate(run_label):
  ax.plot(scm[i].tke_avg[:],vscale*scm[i].z_w, linestyle=styles[i], color = colors[i],linewidth=linewidth[i])
if compar_les==1:
  ax.plot(  tke_les[:],zw_les*vscale, linewidth=2, color='r'  )
#==
ax.set_xlabel(r"${\rm m}^2\;{\rm s}^{-2}$",fontsize=14)
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#===============================
plt.tight_layout()
#plt.savefig('convect_VR18_avg.pdf', bbox_inches='tight')
#===============================
plt.show()

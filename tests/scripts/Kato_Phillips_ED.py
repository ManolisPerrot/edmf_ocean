

#!/usr/bin/env python
# coding: utf-8

import sys  # to put the SCM into the PYTHONPATH

sys.path.append('../../library/F2PY')


###########################################
# Imports
###########################################
from sys import exit
import time as TIME
import xarray as xr
from scipy.interpolate import interp1d
import scipy.signal
from scm_class import SCM
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np


###################################################
#plt.rcParams['text.usetex'] = True
#plt.rcParams.update({'font.size': 18})
#plt.rcParams.update({'figure.facecolor': 'white'})
#plt.rcParams.update({'savefig.facecolor': 'white'})
###########################################
###########################################
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
###########################################
# colors
blue, orange, magenta, grey, green = '#0db4c3', '#eea021', '#ff0364', '#606172', '#3fb532'

# ===========================================================================
case = 'Kato_Phillips_NOCOR'


saving_path = '../figures/'
saving_name = case+'30h_profiles.png'
# ===========================================================================

# Define the common parameters:
common_params = {
    'nz': 50,
    'dt': 30.,
    'h0': 50.,
    'thetas': 6.5,
    'hc': 1000000000000,
    'nbhours': 30,
    'outfreq': 1,
    'output_filename': "scm_output_KP_NOCOR.nc",
    'T0': 16.,
    'N0': 1.e-4,
    'Tcoef': 0.2054,
    'SaltCst': 35.,
    'lat0': 0.,
    'sustr': 1.e-4,
    'svstr': 0.,
    'stflx': 0.,
    'srflx': 0.,
    'ssflx': 0.,
    'btflx': 'no flux',
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
    'wp_bp': 0.003*250,     #      0.002,
    'up_c': 0.5,
    'vp_c': 0.5,
    'bc_ap': 0.2,    #0.3,
    'delta_bkg': 0.005*250,   # 0.006,
    'wp0' : -0.5e-03,
    'akvmin': 1.e-5,
    'aktmin': 1.e-6,
    'mxlmin': 1.,
    'output_filename': 'run',
    'write_netcdf': True
}

# Define parameters specific to each run (overwrite common parameters):
#run_label = ['ED']

runs = [
    {
        'btflx': 'linear_continuation',
        'eddy_diff': True,
        'evd': False,
        'mass_flux_tra': False,
        'mass_flux_dyn': False,
        'mass_flux_tke': False,
        'mass_flux_tke_trplCorr': False,
        'output_filename': 'run1.nc',
        'mxlmin': 1.
    },
    {
        'btflx': 'linear_continuation',
        'eddy_diff': True,
        'evd': False,
        'mass_flux_tra': False,
        'mass_flux_dyn': False,
        'mass_flux_tke': False,
        'mass_flux_tke_trplCorr': False,
        'output_filename': 'run2.nc',
        'mxlmin': 0.05
    }
        ]

scm = [0]*len(runs)

# Run the SCM
for i, run_params in enumerate(runs):
    params = common_params.copy()  # Create a copy of common_params
    params.update(run_params)  # Update with run_params
    scm[i] = SCM(**params)
    scm[i].run_direct()
################################# PLOTTING
fig, axes = plt.subplots(nrows=1, ncols=3, sharex=False,sharey=True, figsize=(7, 4))
#=================
ax = axes[0]
ax.set_title(r"$\overline{\theta}$",fontsize=14)
xmin = 13; xmax = 16.5
ax.set_xlim(xmin,xmax)
ax.plot(  scm[0].t_n[:,scm[0].itemp], scm[0].z_r,linewidth=2, color='k')
ax.plot(  scm[1].t_n[:,scm[1].itemp], scm[0].z_r,linewidth=2, color='b')
# theoretical mixed layer depth
hmxl2  = - 1.05*np.sqrt(scm[0].ustr_sfc)*np.sqrt(30.*3600./0.01)
ax.axhline(y = hmxl2, color = 'r', linestyle = '-')
#=================
ax = axes[1]
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax.set_title(r"$k$",fontsize=14)
ax.plot(  scm[0].tke_np1[:], scm[0].z_w,linewidth=2, color='k')
ax.plot(  scm[1].tke_np1[:], scm[1].z_w,linewidth=2, color='b')
#=================
ax = axes[2]
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax.set_title(r"$K_m$",fontsize=14)
ax.plot(  scm[0].akv[:], scm[0].z_w,linewidth=2, color='k')
ax.plot(  scm[1].akv[:], scm[1].z_w,linewidth=2, color='b')
#=================
plt.tight_layout()
plt.show()
#

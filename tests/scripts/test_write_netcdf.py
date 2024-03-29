

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
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 18})
plt.rcParams.update({'figure.facecolor': 'white'})
plt.rcParams.update({'savefig.facecolor': 'white'})
###########################################


# colors
blue, orange, magenta, grey, green = '#0db4c3', '#eea021', '#ff0364', '#606172', '#3fb532'

# ===========================================================================
case = 'W005_C500_NO_COR'


saving_path = '../figures/'
saving_name = case+'72h_profile_LES_vs_EDMF.png'


# loading LES output once before using the function

file = 'GN_01.1.OC_01.000.nc'


path = '../data/'+case+'/'
les = xr.open_dataset(path+file)
LG_MEAN = xr.open_dataset(
    path+file, group='/LES_budgets/Mean/Cartesian/Not_time_averaged/Not_normalized/cart')
TH_les = (LG_MEAN.MEAN_TH - 273.15).data
U_les = (LG_MEAN.MEAN_U).data
V_les = (LG_MEAN.MEAN_V).data

LG_RES = xr.open_dataset(
    path+file, group='/LES_budgets/Resolved/Cartesian/Not_time_averaged/Not_normalized/cart')
LG_SBG = xr.open_dataset(
    path+file, group='/LES_budgets/Subgrid/Cartesian/Not_time_averaged/Not_normalized/cart')
WTH = (LG_RES.RES_WTH + LG_SBG.SBG_WTHL).data
WU = (LG_RES.RES_WU + LG_SBG.SBG_WU).data
WV = (LG_RES.RES_WV + LG_SBG.SBG_WV).data
corr_tke_les = 0. #-2e-04+1.0e-08
TKE = (LG_RES.RES_KE + LG_SBG.SBG_TKE).data + corr_tke_les
WTKE = (LG_RES.RES_WKE + LG_SBG.SBG_WTKE).data

time_les = les.time_les

# numpy array of integer hours, starting at inital time + 1h
time = ((time_les - time_les[0]) / np.timedelta64(1, 'h')).data.astype(int) + 1

# remap level_les on negative depth values
z_r_les = (les.level_les - (les.level_les[0] + les.level_les[-1])).data
# ===========================================================================

# Define the common parameters:
common_params = {
    'nz': 100,
    'dt': 50.,
    'h0': 2000.,
    'thetas': 6.5,
    'hc': 400,
    'nbhours': 72,
    'outfreq': 1,
    'output_filename': "scm_output.nc",
    'T0': 2.,
    'N0': 1.9620001275490499e-6,
    'Tcoef': 0.2048,
    'SaltCst': 35.,
    'lat0': 0.,
    'sustr': 0.,
    'svstr': 0.,
    'stflx': -500.,
    'srflx': 0.,
    'ssflx': 0.,
    'eddy_diff': True,
    'evd': False,
    'mass_flux_tra': True,
    'mass_flux_dyn': True,
    'mass_flux_tke': True,
    'mass_flux_tke_trplCorr': True,
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
    'output_filename': 'run',
    'write_netcdf': False
}





if case == 'W05_C500':
    common_params['sustr'] = 0.5/1027

if case == 'W005_C500_NO_COR':
    common_params['sustr'] = 0.05/1027

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




scm = [0]*len(runs)

# Run the SCM
for i, run_params in enumerate(runs):
    params = common_params.copy()  # Create a copy of common_params
    params.update(run_params)  # Update with run_params
    scm[i] = SCM(**params)
    scm[i].run_direct()
    print("zinv =", scm[i].zinv)
    print('write_netcdf',params['write_netcdf'])
    instant=5
    plt.plot(scm[i].t_history[:,instant], scm[i].z_r)
    plt.plot(TH_les[instant], z_r_les)
    plt.plot(scm[i].u_history[:,instant], scm[i].z_r)
    plt.plot(U_les[instant], z_r_les)

    plt.plot(scm[i].v_history[:,instant], scm[i].z_r)
    plt.plot(V_les[instant], z_r_les)
    plt.show()

    level= 90

    plt.plot(scm[i].t_history[level,:])
    plt.plot(TH_les[:,level])
    plt.plot(scm[i].u_history[level,:])
    plt.plot(U_les[:,level])

    plt.plot(scm[i].v_history[level,:])
    plt.plot(V_les[:,level])
    plt.show()
exit()




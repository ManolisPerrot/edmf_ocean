

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
# from scm_class import SCM
from scm_class_oce import SCM
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from case_configs import case_params
import subprocess
from unit_tests import is_in_range
from WANG1_plot_instant_panel import plot_instant_panel
from WANG1_plot_velocities_mean_only import plot_mean_velocities
from WANG1_plot_condsamp import plot_condsamp_panel
# ===========================================================================
# cases = ['WANG1_FR_lat30','WANG1_FR_lat60','WANG1_FR_lat90']
# cases = ['WANG1_NR_new','WANG1_FR_lat60','WANG1_FR_lat30']
cases = ['WANG1_FR_lat30','WANG1_FR_lat60',]
# cases = ['WANG1_FR_lat30',]
# ====================================Define configurations=======================
# Define the common parameters (attention some of them will be overwritten by case_configurations.py):
edmf_params = {
    'mass_flux_tra': True,
    'mass_flux_dyn': True,
    'mass_flux_tke': True,
    'mass_flux_tke_trplCorr': True,
    'mass_flux_small_ap': False,
    'lin_eos': True,
    'extrap_ak_surf': True,
    'tke_sfc_dirichlet': False,
    'eddy_diff_tke_const': 'NEMO',
    # 'entr_scheme': 'R10corNT',
    # 'entr_scheme': 'R10',
    'trad_coriolis_mod': True,
    'Cent': 0.99,    
    'Cdet': 1.99,       # 'Cdet': 2.5,
    'wp_a': 1.,
    'wp_b': 1.,      # 'wp_b': 1.
    'wp_bp': 0.003*250,     #      0.002,
    'up_c': 0.5,
    'vp_c': 0.5,
    'bc_ap': 0.18,#     0.35,    #0.3,
    'delta_bkg': 0.0045*250,   # 0.02,
    'output_filename': 'run',
    'wp0':-1e-02,
    'write_netcdf': True
}

run_label = cases

# runs = [
#     {
#         'entr_scheme': 'R10',
#         'output_filename': 'scm_'+cases[0]+'.nc',
#     },
#         {
#         'output_filename': 'scm_'+cases[1]+'.nc',
#     },
#             {
#         'output_filename': 'scm_'+cases[2]+'.nc',
#     },
#             {
#         'output_filename': 'scm_'+cases[3]+'.nc',
#     },
#     # {
#     #     'output_filename': 'scm_'+cases[4]+'.nc',
#     # },
#         ]
# runs = [
#     {
#         # 'entr_scheme': 'R10',
#         'output_filename': 'scm_'+cases[0]+'.nc',
#     },]
# for case in cases[1:]:
#     runs.append({'output_filename': 'scm_'+case+'.nc'})

runs = [ {'output_filename': 'scm_'+case+'.nc'} for case in cases]

scm = [0]*len(runs)

# Run the SCM
for i, run_params in enumerate(runs):
    print('Running case', cases[i])
    params = edmf_params.copy()  # Create a copy of common_params
    params.update(run_params)  # Update with run_params
    params.update(case_params[cases[i]])
    # !!!!!!!!!!!!!!!!!!! Modifying dt to speed up preliminary tests, to remove !!!!!!!!!!!!!!!!!!
    params.update({'dt':1000})
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # scm[i] = SCM(**params)
    scm[i] = SCM(params,cfg_params={})
    scm[i].run_direct()
    # test zinv
    # if scm[i].MF_tra or scm[i].MF_dyn: 
    #     print(run_label[i])
    #     reference=mld
    #     is_in_range(value=scm[i].zinv, value_name='zinv', reference=reference,tolerance=10 )

# cases = cases+['WANG1_NR_new']
plot_mean_velocities(cases=cases,savefig=True,)
plot_instant_panel(cases=cases,savefig=True,avg_start=-10,avg_stop=-8)
plot_condsamp_panel(cases=cases,plot_vorticity=True)
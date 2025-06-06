import sys  # to put the SCM into the PYTHONPATH

sys.path.append('../../library_old/F2PY')

#!/usr/bin/env python
# coding: utf-8

###########################################
# Imports
###########################################
from sys import exit
import time as TIME
import xarray as xr
from scipy.interpolate import interp1d
import scipy.signal
#old code
from scm_class import SCM
#new code
# from scm_class_oce import SCM
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from unit_tests import is_in_range
###################################################
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 18})
plt.rcParams.update({'figure.facecolor': 'white'})
plt.rcParams.update({'savefig.facecolor': 'white'})
###########################################
blue, orange, magenta, grey, green = '#0db4c3', '#eea021', '#ff0364', '#606172', '#3fb532'

# ===========================================================================
#case = 'FC500'
case = 'W005_C500_NO_COR'


saving_path = '../figures/'
saving_name = case+'72h_profile_LES_vs_EDMF.png'
# ===========================================================================
# Define the common parameters:
common_params = {
    'nz': 2000,
    'dt': 50.,
    'h0': 2000.,
    'thetas': 6.5,
    'hc': 1.e+16,
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
    'mass_flux_small_ap': True,
    'lin_eos': True,
    'extrap_ak_surf': True,
    'tke_sfc_dirichlet': False,
    'eddy_diff_tke_const': 'NEMO',
    'entr_scheme': 'R10',
    'Cent': 0.99,
    'Cdet': 1.99,       # 'Cdet': 2.5,
    'wp_a': 1.,   #1
    'wp_b': 1.,     #1.
    'wp_bp': 0.003*250,    
    'up_c': 0.5, 
    'vp_c': 0.5,
    'bc_ap': 0.2,    #0.3,
    'delta_bkg': 0.0045*250,   # 0.005,
    'wp0' : -0.5e-08,
    'output_filename': 'run',
    'beta_bc_P09': 0.3,
    'write_netcdf': True
}

if case == 'W05_C500':
    common_params['sustr'] = 0.5/1027

if case == 'W005_C500_NO_COR':
    common_params['sustr'] = 0.05/1027

# Define parameters specific to each run (overwrite common parameters):

run_label = [r'EDMF-TKE-inconsistent',  r'EDMF-boundary-inconsistent', r'EDMF-Energy',r'EDMF-boundary-consistent']
runs = [
    {
        'eddy_diff': True,
        'evd': False,
        'mass_flux_tra': True,
        'mass_flux_dyn': True,
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
        'bc_P09': 'inconsistent',
        'output_filename': 'run2.nc'
    },
    {
        'eddy_diff': True,
        'evd': False,
        'mass_flux_tra': True,
        'mass_flux_dyn': True,
        'mass_flux_tke': True,
        'mass_flux_tke_trplCorr': True,
        'bc_P09': 'false',
        'output_filename': 'run3.nc'
    },
    {
        'eddy_diff': True,
        'evd': False,
        'mass_flux_tra': True,
        'mass_flux_dyn': True,
        'mass_flux_tke': True,
        'mass_flux_tke_trplCorr': True,
        'bc_P09': 'consistent',
        'output_filename': 'run4.nc'
    },
        ]


scm = [0]*len(runs)

# Run the SCM
for i, run_params in enumerate(runs):
    params = common_params.copy()  # Create a copy of common_params
    params.update(run_params)  # Update with run_params
    #old code
    scm[i] = SCM(params)
    #new code
    # scm[i] = SCM(params,cfg_params={}) #empty cfg_params because its is already contained in params...
    scm[i].run_direct()

    # do unit tests
    if scm[i].MF_tra or scm[i].MF_dyn: 
        zinv_ref = -320.
        tolerance=10.
        print(run_label[i])
        is_in_range(value=scm[i].zinv, value_name='zinv', reference=zinv_ref,tolerance=tolerance )

# LOAD outputs

out = [0]*len(runs)

for i, run_params in enumerate(runs):
    print('opening '+run_params['output_filename'])
    out[i] = xr.open_dataset(run_params['output_filename'])

instant = 71

#mld = (-z_r_les[(-WTH[instant]).argmax()]).data



# PLOTTING
styles = ['r-', 'b-', 'k-']
alpha = [0.75]*(len(run_label))
linewidth = [3]*(len(run_label))

style_les = 'ko'
alpha_les = 0.5
linewidth_les = 4


#================================================================
#================================================================
fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True,
                         sharey=True , figsize=(9,6))
# ===============================================================
ax = axes
ax.set_title(r'\rm{Vertically integrated total energy budget}')
ax.set_xlabel(r'\rm{time} (hours)')
ax.set_ylabel(r'${\rm m}^{3}\;{\rm s}^{-3}$')
ax.plot( (out[0]['Etot'] )[1:]  ,color='tab:red'  , linewidth=3 , alpha=1, linestyle = '-', label=run_label[0])
ax.plot(0.5*2.456e-7*out[0]['zinv'], color='tab:gray', linestyle='--', label=r'$B_0 h/2$')
ax.plot( (out[2]['Etot'] )[1:]  ,color='tab:green', linewidth=3 , alpha=0.75, linestyle = '-', label=run_label[2])
ax.plot( (out[3]['Etot'] )[1:]  ,color='tab:green', linewidth=3 , alpha=1, linestyle = '--', label=run_label[3])
ax.plot( (out[1]['Etot'] )[1:]  ,color='tab:green', linewidth=3 , alpha=1, linestyle = ':', label=run_label[1])

ax.set_yscale('symlog',linthresh=1e-13)
ax.axhline(y=0, color='k', linestyle='--', linewidth=1)
ax.legend(fontsize=12)
# plt.plot()
fig.tight_layout()
plt.savefig(saving_path+'W005_C500_engy_budgets', bbox_inches = 'tight', dpi=300)
print('fig saved at '+saving_path+'W005_C500_engy_budgets')
# plt.show()
#================================================================
#================================================================

# fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True,
#                          sharey=True, figsize=(6, 4) )
# # ===============================================================
# ax = axes

# ax.set_title(r'\rm{Vertically integrated total energy budget}')
# ax.set_xlabel(r'\rm{time} (hours)')
# ax.set_ylabel(r'${\rm m}^{3}\;{\rm s}^{-3}$')
# ax.plot( (out[2]['Etke'] )[1:]  , linewidth=3 , alpha=1, linestyle = '-', label='Etke')
# ax.plot( (out[2]['Ekin'] )[1:]  , linewidth=3 , alpha=1, linestyle = '-', label='Ekin')
# ax.plot( (out[2]['Epot'] )[1:]  , linewidth=3 , alpha=1, linestyle = '-', label='Epot')
# ax.plot( (out[2]['Etot'] )[1:]  , linewidth=3 , alpha=1, linestyle = '-', label='Etot')

# # ===============================================================

# ax.legend()
# saving_name=case+'_energy_diag'
# plt.savefig(saving_path+saving_name, bbox_inches = 'tight', dpi=300)
# print('fig saved at'+saving_path+saving_name)


# plt.show()
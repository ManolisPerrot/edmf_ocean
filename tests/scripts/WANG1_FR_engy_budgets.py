import sys  # to put the SCM into the PYTHONPATH

sys.path.append('../../library/F2PY')

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
from scm_class import SCM
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from case_configs import case_params
###################################################
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 18})
plt.rcParams.update({'figure.facecolor': 'white'})
plt.rcParams.update({'savefig.facecolor': 'white'})
###########################################

# ===========================================================================
case = 'WANG1_FR'

# ====================================Define configurations=======================
# Define the common parameters (attention some of them will be overwritten by case_configurations.py):
common_params = {
    'nz': 450,
    'dt': 50.,
    'h0': 4500.,
    'thetas': 6.5,
    'hc': 4000000000000,
    'nbhours': 276,
    'outfreq': 1,
    'output_filename': "scm_output.nc",
    'T0': 3.,
    'N0': 1.5865490613891073e-08, #equals to f(60°)
    'mld_ini': -1000.,
    'Tcoef': 0.2048,
    'SaltCst': 32.6,
    'lat0': 60.,
    'sustr': 0.,
    'svstr': 0.,
    'stflx': -111.12982, #Wm-2
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
    'eddy_diff_tke_const': 'MNH',
    'entr_scheme': 'R10',
    'Cent': 0.9,
    'Cdet': 1.95,       # 'Cdet': 2.5,
    'wp_a': 1.,
    'wp_b': 1.,      # 'wp_b': 1.
    'wp_bp': 0.005*250,     #      0.002,
    'up_c': 0.25,
    'vp_c': 0.25,
    'bc_ap': 0.2,    #0.3,
    'delta_bkg': 0.0025*250,   # 0.006,
    'output_filename': 'run',
    'write_netcdf': True
}

common_params.update(case_params[case])  # Update with the specific case configuration in case_params[case]


# Define parameters specific to each run (overwrite common parameters):

run_label = ['ED',  'EDMF-Energy','EDMF-Energy-cor' ]
runs = [
    {
        'eddy_diff': True,
        'evd': False,
        'mass_flux_tra': False,
        'mass_flux_dyn': False,
        'mass_flux_tke': False,
        'mass_flux_tke_trplCorr': False,
        'output_filename': 'run1.nc',
    'write_netcdf': False
    },
            {
        'eddy_diff': True,
        'evd': False,
        'mass_flux_tra': True,
        'mass_flux_dyn': True,
        'mass_flux_tke': True,
        'mass_flux_tke_trplCorr': True,
        'entr_scheme': 'R10',
        'output_filename': 'run2.nc',
    'write_netcdf': False
    },
        {
        'eddy_diff': True,
        'evd': False,
        'mass_flux_tra': True,
        'mass_flux_dyn': True,
        'mass_flux_tke': True,
        'mass_flux_tke_trplCorr': True,
        'entr_scheme': 'R10corNT',
        'output_filename': 'scm_WANG1_FR.nc',
    'write_netcdf': True
    }
        ]




scm = [0]*len(runs)

# Run the SCM
for i, run_params in enumerate(runs):
    params = common_params.copy()  # Create a copy of common_params
    params.update(run_params)  # Update with run_params
    scm[i] = SCM(**params)
    scm[i].run_direct()
    if params['mass_flux_tra']==True:
        print("zinv =", scm[i].zinv)

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

fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True,
                         sharey=True, figsize=(6, 4) )
# ===============================================================
ax = axes

ax.set_title(r'\rm{Vertically integrated total energy budget}')
ax.set_xlabel(r'\rm{time} (hours)')
ax.set_ylabel(r'${\rm m}^{3}\;{\rm s}^{-3}$')
ax.plot( (out[2]['Etot'] )[1:] , color='tab:orange' , linewidth=3 , alpha=1, linestyle = '-', label='EDMF-Energy')
ax.plot( (out[1]['Etot'] )[1:] , color='tab:blue'   , linewidth=3 , alpha=1, linestyle = '-', label='EDMF')
# ===============================================================

ax.legend()

saving_name = '../figures/'+case+'_engy_budgets'
plt.savefig(saving_name, bbox_inches = 'tight', dpi=300)
print('fig saved at'+saving_name)


#================================================================
#================================================================

fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True,
                         sharey=True, figsize=(6, 4) )
# ===============================================================
ax = axes

ax.set_title(r'\rm{Vertically integrated total energy budget}')
ax.set_xlabel(r'\rm{time} (hours)')
ax.set_ylabel(r'${\rm m}^{3}\;{\rm s}^{-3}$')
ax.plot( (out[2]['Etke'] )[1:]  , linewidth=3 , alpha=1, linestyle = '-', label='Etke')
ax.plot( (out[2]['Ekin'] )[1:]  , linewidth=3 , alpha=1, linestyle = '-', label='Ekin')
ax.plot( (out[2]['Epot'] )[1:]  , linewidth=3 , alpha=1, linestyle = '-', label='Epot')
ax.plot( (out[2]['Etot'] )[1:]  , linewidth=3 , alpha=1, linestyle = '-', label='Etot')

# ===============================================================

ax.legend()
saving_name='../figures/'+case+'_energy_diag'
plt.savefig(saving_name, bbox_inches = 'tight', dpi=300)
print('fig saved at '+saving_name)


plt.show()
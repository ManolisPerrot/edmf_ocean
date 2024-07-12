

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
from unit_tests import is_in_range

###################################################
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 18})
plt.rcParams.update({'figure.facecolor': 'white'})
plt.rcParams.update({'savefig.facecolor': 'white'})
###########################################


# colors
blue, orange, magenta, grey, green = '#0db4c3', '#eea021', '#ff0364', '#606172', '#3fb532'

# ===========================================================================


# saving_path = '../figures/'
# saving_name = case+'72h_profile_LES_vs_EDMF.png'



# ===========================================================================

Tcoef= 0.20495


# Define the common parameters:
common_params = {
    'nz': 512,
    'dt': 360.,
    'h0': 100.,
    'thetas': 6.5,         # /!\ ATTENTION pas la même def de la grille stretchée!!
    'hc': 100, # /!\ ATTENTION pas la même def de la grille stretchée!!
    'gridType': 'ORCA75',
    'nbhours': 7*24+1,
    'outfreq': 1,
    'output_filename': "scm_output.nc",
    'T0': 20.,
    'N0': 1.9620001275490499e-5,
    'Tcoef': 0.20495, #= alpha*rho0
    'Tref':20.,
    'SaltCst': 35.,
    'Sref': 35.,
    'Scoef':0.8198,   # beta*rho0
    'lat0': 0.,
    'cpoce': 3985,
    'rho0': 1024.763,
    'sustr': 0.,
    'svstr': 0.,
    'stflx': -75.,
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
    'wp_a': 1.3,
    'wp_b': 1.3,      # 'wp_b': 1.
    'wp_bp': 0.003*250,     #      0.002,
    'up_c': 0.5,
    'vp_c': 0.5,
    'bc_ap': 0.2,    #0.3,
    'delta_bkg': 0.006*250,   # 0.006,
    'wp0' : -0.5e-08,
    'output_filename': 'run',
    'write_netcdf': True
}



# Define parameters specific to each run (overwrite common parameters):

run_label = ['ED+EVD', 'EDMF', 'EDMF-Energy']
runs = [
    {
        'eddy_diff': True,
        'evd': True,
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
        'output_filename': 'FC_alexandre_edmf_ocean.nc'
    }
        ]




scm = [0]*len(runs)

# Run the SCM
for i, run_params in enumerate(runs):
    params = common_params.copy()  # Create a copy of common_params
    params.update(run_params)  # Update with run_params
    scm[i] = SCM(**params)
    scm[i].run_direct()

    # test zinv
    if scm[i].MF_tra or scm[i].MF_dyn: 
        print(run_label[i])
        reference=-320
        #reference=-mld #LES value
        is_in_range(value=scm[i].zinv, value_name='zinv', reference=reference,tolerance=10 )

# LOAD outputs

# TH_scm = [0]*len(runs)
# U_scm = [0]*len(runs)
# V_scm = [0]*len(runs)

# #interpolate scm output on LES #TODO do the converse to reduce computation cost?


# for i, run_params in enumerate(runs):
#     TH_scm = scm[i].t_history
#     U_scm  = scm[i].u_history
#     V_scm  = scm[i].v_history

# instant = 71


# mld = (-z_r_les[(-WTH[instant]).argmax()]).data


# ################################# PLOTTING
# styles = ['-', '-', '-']
# #colors = ['k',blue,orange]
# colors = ['k','tab:blue','tab:orange']
# alpha = [0.5,1,1]
# linewidth = [4]*(len(run_label))

# style_les = 'ko'
# alpha_les = 1
# linewidth_les = 4





# ###########################################################################################
# #============================================ FC ==========================================
# ###########################################################################################

# if case == 'FC500':
#     fig, axes = plt.subplots(nrows=1, ncols=4, sharex=False,
#                          sharey=True, figsize=(15, 5))

#     # ===============================================================
#     ax = axes.flat[0]
#     ax.set_xlabel(r'$ C$')
#     ax.set_ylabel(r'$z / h $')
#     ax.set_title(r'$\overline{\theta}$')


#     ax.plot(TH_les[instant], z_r_les/mld, style_les,
#             alpha=alpha_les, linewidth=linewidth_les,  label='LES')

#     for i, label in enumerate(run_label):
#         ax.plot(scm[i].t_np1[:, 0], scm[i].z_r/mld, linestyle=styles[i], color = colors[i],
#                 alpha=alpha[i], linewidth=linewidth[i], label=label)

#     ax.set_xlim((1.65, 1.78))
#     ax.set_ylim((-1.3, 0))

#     # ===============================================================


#     ax = axes.flat[1]
#     ax.set_title(r'$\overline{w^\prime \theta^\prime}$')

#     ax.plot(WTH[instant], z_r_les/mld, style_les,
#             alpha=alpha_les, linewidth=linewidth_les, label='LES')

#     for i, label in enumerate(run_label):
#         if run_label == 'ED':
#             ax.plot(-(scm[i].wted), scm[i].z_w/mld, linestyle=styles[i], color = colors[i],
#                     alpha=alpha[i], linewidth=linewidth[i], label=label)
#         else:
#             ax.plot(-(scm[i].wted + scm[i].wtmf), scm[i].z_w/mld, linestyle=styles[i], color = colors[i],
#                     alpha=alpha[i], linewidth=linewidth[i], label=label)
            
#     ax.plot( -scm[1].wted, scm[1].z_w/mld, color ='tab:blue'  , linestyle ='--', alpha=1.0 , linewidth=3 )
#     ax.plot( -scm[2].wted, scm[2].z_w/mld, color ='tab:orange', linestyle ='--', alpha=1.0 , linewidth=3 )

#     ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

#     #ax.plot( -out[1]['WT_ED'][-1,:], out[1].z_w/mld, color ='tab:blue'  , linestyle ='--', alpha=1.0 , linewidth=3 )
#     #ax.plot( -out[2]['WT_ED'][-1,:], out[2].z_w/mld, color ='tab:orange', linestyle ='--', alpha=1.0 , linewidth=3 )

#     ax.set_ylim((-1.3, 0))

#     ax.set_xlabel(r'$K.m.s^{-1}$')


#     # ===============================================================
#     # ===============================================================
#     ax = axes.flat[2]
#     ax.set_title(r'$k$')

#     ax.plot(TKE[instant], z_r_les/mld, style_les,
#             alpha=alpha_les, linewidth=linewidth_les, label='LES')

#     for i, label in enumerate(run_label):
#         ax.plot(scm[i].tke_np1, scm[i].z_w/mld, linestyle=styles[i], color = colors[i],
#                 alpha=alpha[i], linewidth=linewidth[i], label=label)
#     ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

#     ax.set_xlim((-0.0001, 0.0009))
#     ax.set_ylim((-1.3, 0))

#     ax.set_xlabel(r'$m^2.s^{-2}$')



#     # ===============================================================
#     ax = axes.flat[3]
#     ax.set_title(r'$\overline{w^\prime \frac{\mathbf{u}^{\prime 2}}{2}  } + \frac{1}{\rho_0} \overline{w^\prime p^{\dagger \prime} }$')

#     #add velocity-pressure correlation
#     cond_samp = xr.open_dataset(path+    'FC500_object_diags_Cw_m1_72h.nc')  

#     ax.plot(WTKE[instant]-cond_samp['TOT_intra_WPHI_over_RHO_0'][-1], z_r_les/mld, style_les,
#             alpha=alpha_les, linewidth=linewidth_les, label='LES')


#     for i, label in enumerate(run_label):
#         ax.plot((scm[i].wtke), scm[i].z_r/mld, linestyle=styles[i], color = colors[i],
#                 alpha=alpha[i], linewidth=linewidth[i], label=label)


#     ax.set_xlim((- 1.3e-5, 0))
#     ax.set_ylim((-1.3, 0))

#     ax.set_xlabel(r'$m^3.s^{-3}$')

#     handles, labels = ax.get_legend_handles_labels()
#     fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(
#         0.5, -0.05), fancybox=False, shadow=False, ncol=4)

#     for ax in axes.flat:
#         ax.set_box_aspect(1)



#     fig.tight_layout()

#     # ===============================================================






# plt.savefig(saving_path+saving_name, bbox_inches='tight', dpi=300)


# print('figure saved at'+saving_path+saving_name)

# #plt.show()
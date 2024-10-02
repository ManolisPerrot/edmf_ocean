

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
from case_configs import case_params
import subprocess
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
case = 'WANG1_FR'


saving_path = '../figures/'
saving_name = case+'_profile_LES_vs_EDMF.png'


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
BU_KE = xr.open_dataset(path+file,group ='/LES_budgets/BU_KE/Cartesian/Not_time_averaged/Not_normalized/cart')

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

# choose intant to plot    
instant = -1 
### ATTENTION: scm contains time 0 + 276 hours = 277 instants
###            les contains hour 1 to 276 = 276 instants
### so SAME END but different start in the lsit

# compute LES mixed layer depth
# to non-dimensionalize depths
mld = (-z_r_les[(-WTH[instant]).argmax()])

# ===========================================================================

# ====================================Define configurations=======================
# Define the common parameters (attention some of them will be overwritten by case_configurations.py):
common_params = {
    'nz': 4500,
    'dt': 200.,
    'h0': 4500.,
    'thetas': 6.5,
    'hc': 4000000000000,
    'nbhours': 276,
    'outfreq': 1,
    'output_filename': "scm_output.nc",
    'T0': 3.,
    'N0': 1.5865490613891073e-08, #equals to f(60°)
    'mld_ini_temp': -1000.,
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
    'mass_flux_small_ap': False,
    'lin_eos': True,
    'extrap_ak_surf': True,
    'tke_sfc_dirichlet': False,
    'eddy_diff_tke_const': 'MNH',
    'entr_scheme': 'R10',
    'Cent': 0.99,    
    'Cdet': 1.99,       # 'Cdet': 2.5,
    'wp_a': 1.,
    'wp_b': 1.,      # 'wp_b': 1.
    'wp_bp': 0.003*250,     #      0.002,
    'up_c': 0.25,
    'vp_c': 0.25,
    'bc_ap': 0.35,    #0.3,
    'delta_bkg': 0.006*250,   # 0.02,
    'output_filename': 'run',
    'wp0':-1.e-02,
    'write_netcdf': True
}

# Good for theta and velocities: 
    # 'Cent': 0.99,
    # 'Cdet': 1.99,       # 'Cdet': 2.5,
    # 'wp_a': 1.,
    # 'wp_b': 1.,      # 'wp_b': 1.
    # 'wp_bp': 0.003*250,     #      0.002,
    # 'up_c': 0.2,
    # 'vp_c': 0.2,
    # 'bc_ap': 0.2,    #0.3,
    # 'delta_bkg': 0.02*250,   # 0.02,
    # 'output_filename': 'run',
    # 'wp0':-1.e-08,

    # 'Cent': 0.3,
    # 'Cdet': 1.99,       # 'Cdet': 2.5,
    # 'wp_a': 1.,
    # 'wp_b': 1.,      # 'wp_b': 1.
    # 'wp_bp': 0.003*250,     #      0.002,
    # 'up_c': 0.5,
    # 'vp_c': 0.5,
    # 'bc_ap': 0.35,    #0.3,
    # 'delta_bkg': 0.02*250,   # 0.02,
    # 'output_filename': 'run',
    # 'wp0':-1.e-02,


common_params.update(case_params[case])  # Update with the specific case configuration in case_params[case]


# Define parameters specific to each run (overwrite common parameters):

# run_label = ['ED+EVD',  'EDMF-Energy','EDMF-Energy-cor' ]
# runs = [
#     {
#         'eddy_diff': True,
#         'evd': True,
#         'mass_flux_tra': False,
#         'mass_flux_dyn': False,
#         'mass_flux_tke': False,
#         'mass_flux_tke_trplCorr': False,
#         'output_filename': 'run1.nc',
#     'write_netcdf': False
#     },
#             {
#         'eddy_diff': True,
#         'evd': False,
#         'mass_flux_tra': True,
#         'mass_flux_dyn': True,
#         'mass_flux_tke': True,
#         'mass_flux_tke_trplCorr': True,
#         'entr_scheme': 'R10',
#         'output_filename': 'run2.nc',
#     'write_netcdf': False
#     },
#         {
#         'eddy_diff': True,
#         'evd': False,
#         'mass_flux_tra': True,
#         'mass_flux_dyn': True,
#         'mass_flux_tke': True,
#         'mass_flux_tke_trplCorr': True,
#         'entr_scheme': 'R10corNT',
#         'output_filename': 'scm_WANG1_FR.nc',
#     'write_netcdf': True
#     }
#         ]

run_label = ['EDMF-Energy-cor' ]
runs = [
        {
        'eddy_diff': True,
        'evd': False,
        'mass_flux_tra': True,
        'mass_flux_dyn': True,
        'mass_flux_tke': True,
        'mass_flux_tke_trplCorr': True,
        'entr_scheme': 'R10corNT',
        # 'entr_scheme': 'R10',
        'output_filename': 'scm_WANG1_FR.nc',
        'write_netcdf': True,
        'trad_coriolis_mod': True
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
        reference=mld
        is_in_range(value=scm[i].zinv, value_name='zinv', reference=reference,tolerance=10 )

# LOAD outputs
out = [0]*len(runs)

for i, run_params in enumerate(runs):
    print('opening '+run_params['output_filename'])
    out[i] = xr.open_dataset(run_params['output_filename'])

# choose intant to plot    
instant = -1 
### ATTENTION: scm contains time 0 + 276 hours = 277 instants
###            les contains hour 1 to 276 = 276 instants
### so SAME END but different start in the lsit

# compute LES mixed layer depth
# to non-dimensionalize depths
mld = (-z_r_les[(-WTH[instant]).argmax()]).data
mld=1

################################# PLOTTING #################################
styles = ['-', '-', '-']
#colors = ['k',blue,orange]
colors = ['k','tab:blue','tab:orange']
alpha = [0.5,1,1]
linewidth = [4]*(len(run_label))

style_les = 'ko'
alpha_les = 1
linewidth_les = 4


def plot_intant_panel(instant=instant):
    #============================================ WC ===============================================
    fig, axes = plt.subplots(nrows=4, ncols=3, sharex=False,
                            sharey=True, figsize=(15, 12))

    ax_index=-1

    # ===============================================================
    ax_index+=1
    ax = axes.flat[ax_index]
    ax.set_xlabel(r'$^{\circ}{\rm C}$')
    ax.set_ylabel(r'$z / h $')
    ax.set_title(r'$\overline{\theta}$')


    ax.plot(TH_les[instant], z_r_les/mld, style_les,
            alpha=alpha_les, linewidth=linewidth_les,  label='LES')

    for i, label in enumerate(run_label):
        # ax.plot(scm[i].t_np1[:, 0], scm[i].z_r/mld, linestyle=styles[i], color = colors[i],
        ax.plot(out[i].temp[instant], out[i].z_r/mld, linestyle=styles[i], color = colors[i],

                alpha=alpha[i], linewidth=linewidth[i], label=label)

    #ax.set_xlim((1.55, 1.8))
    #ax.set_ylim((-1.3, 0))


    # ===============================================================

    ax_index+=1
    ax = axes.flat[ax_index]
    ax.set_title(r'$\overline{u}$')

    ax.set_xlabel(r'${\rm m}\;{\rm s}^{-1}$')

    ax.plot(U_les[instant], z_r_les/mld, style_les,
            alpha=alpha_les, linewidth=linewidth_les, label='LES')

    for i, label in enumerate(run_label):
        # ax.plot((scm[i].u_np1), scm[i].z_r/mld, linestyle=styles[i], color = colors[i],
        #         alpha=alpha[i], linewidth=linewidth[i], label=label)
        # ax.plot((scm[i].up), scm[i].z_w/mld, linestyle='--', color = colors[i],
        #         alpha=alpha[i], linewidth=linewidth[i], label=label)
        # if i != 0: 
        #     ax.plot((out[i].w_p[instant]), out[i].z_w/mld, linestyle=':', color = colors[i],
        #         alpha=alpha[i], linewidth=linewidth[i], label=label)
        ax.plot((out[i]['u'][instant]), out[i].z_r/mld, linestyle=styles[i], color = colors[i],
                alpha=alpha[i], linewidth=linewidth[i], label=label)
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    #ax.set_ylim((-1.3, 0))

    # ===============================================================

    ax_index+=1
    ax = axes.flat[ax_index]
    ax.set_title(r'$\overline{v}$')

    ax.set_xlabel(r'${\rm m}\;{\rm s}^{-1}$')

    ax.plot(V_les[instant], z_r_les/mld, style_les,
            alpha=alpha_les, linewidth=linewidth_les, label='LES')

    for i, label in enumerate(run_label):
        # ax.plot((scm[i].u_np1), scm[i].z_r/mld, linestyle=styles[i], color = colors[i],
        #         alpha=alpha[i], linewidth=linewidth[i], label=label)
        # ax.plot((scm[i].up), scm[i].z_w/mld, linestyle='--', color = colors[i],
        #         alpha=alpha[i], linewidth=linewidth[i], label=label)
        # if i != 0: 
        #     ax.plot((out[i].w_p[instant]), out[i].z_w/mld, linestyle=':', color = colors[i],
        #         alpha=alpha[i], linewidth=linewidth[i], label=label)
        ax.plot((out[i]['v'][instant]), out[i].z_r/mld, linestyle=styles[i], color = colors[i],
                alpha=alpha[i], linewidth=linewidth[i], label=label)
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    #ax.set_ylim((-1.3, 0))
    # ===============================================================
    ax_index+=1
    ax = axes.flat[ax_index]

    ax.set_title(r'$\overline{w^\prime \theta^\prime}$')

    ax.plot(WTH[instant], z_r_les/mld, style_les,
            alpha=alpha_les, linewidth=linewidth_les, label='LES')
    
    for i, label in enumerate(run_label):
        ax.plot( -out[i]['WT'][instant,:], out[i].z_w/mld, color ='tab:blue'  , linestyle ='-', alpha=1.0 , linewidth=3 )
        # ax.plot( -out[i]['WT_ED'][instant,:], out[i].z_w/mld, color ='tab:blue'  , linestyle =':', alpha=1.0 , linewidth=3 )
        # ax.plot( -out[i]['WT_MF'][instant,:], out[i].z_w/mld, color ='tab:blue'  , linestyle ='--', alpha=1.0 , linewidth=3 )

    ax.set_xlabel(r'${\rm K}\;{\rm m}\;{\rm s}^{-1}$')

    # ===============================================================
    # ===============================================================
    ax_index+=1
    ax = axes.flat[ax_index]
    ax.set_title(r'$k$')

    ax.plot(TKE[instant], z_r_les/mld, style_les,
            alpha=alpha_les, linewidth=linewidth_les, label='LES')

    for i, label in enumerate(run_label):
        #ax.plot(scm[i].tke_np1, scm[i].z_w/mld, linestyle=styles[i], color = colors[i],
        #        alpha=alpha[i], linewidth=linewidth[i], label=label)
        
        ax.plot(out[i]['tke'][instant], out[i].z_w/mld, linestyle=styles[i], color = colors[i],
                alpha=alpha[i], linewidth=linewidth[i], label=label)
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

    #ax.set_xlim((-0.0001, 0.001))
    #ax.set_ylim((-1.3, 0))

    ax.set_xlabel(r'${\rm m}^2\;{\rm s}^{-2}$')
    # ===============================================================


    ax_index+=1
    ax = axes.flat[ax_index]
    ax.set_title(r'$\overline{w^\prime \frac{u^{\prime 2}}{2}  }$')

    ax.plot(WTKE[instant], z_r_les/mld, style_les,
            alpha=alpha_les, linewidth=linewidth_les, label='LES')

    for i, label in enumerate(run_label):
        # ax.plot((scm[i].wtke), scm[i].z_r/mld, linestyle=styles[i], color = colors[i],
        #         alpha=alpha[i], linewidth=linewidth[i], label=label)
        ax.plot(out[i]['WTKE'][instant], out[i].z_r/mld, linestyle=styles[i], color = colors[i],
                alpha=alpha[i], linewidth=linewidth[i], label=label)
        
    # if i==1:
    #     ax.plot((out[i]['we'][instant]), out[i].z_w/mld, linestyle=':', color = colors[i],
    #             alpha=alpha[i], linewidth=linewidth[i], label=label)


    #ax.set_xlim((- 1e-5, 0))
    #ax.set_ylim((-1.3, 0))

    ax.set_xlabel(r'${\rm m}^3\;{\rm s}^{-3}$')

    # ===============================================================
    # ax_index+=1
    # ax = axes.flat[ax_index]

    # ax.set_title(r'$\overline{w^\prime b^\prime}$')

    # coef = 9.81*out[0].attrs['alpha']
    # ax.plot((BU_KE['RES_TP'] + BU_KE['SBG_TP'])[instant], z_r_les/mld, style_les,
    #         alpha=alpha_les, linewidth=linewidth_les, label='LES')
    
    # for i, label in enumerate(run_label):
    #     ax.plot( -coef*out[i]['WT'][-1,:], out[i].z_w/mld, color ='tab:blue'  , linestyle ='--', alpha=1.0 , linewidth=3 )
    # ax.set_xlabel(r'${\rm K}\;{\rm m}\;{\rm s}^{-1}$')

    # # ===============================================================

    # ===============================================================
    ax_index+=1
    ax = axes.flat[ax_index]

    ax.set_title(r'$\overline{w^\prime u^\prime}$')

    ax.plot(-(WU)[instant], z_r_les/mld, style_les,
    # ax.plot((WV)[instant], z_r_les/mld, style_les,
            alpha=alpha_les, linewidth=linewidth_les, label='LES')
    
    for i, label in enumerate(run_label):
        ax.plot( out[i]['WU'][instant,:], out[i].z_w/mld, color ='tab:blue'  , linestyle ='-', alpha=1.0 , linewidth=3 )
        # ax.plot( out[i]['WU_MF'][instant,:], out[i].z_w/mld, color ='tab:blue'  , linestyle ='--', alpha=1.0 , linewidth=3 )
        # ax.plot( out[i]['WU_ED'][instant,:], out[i].z_w/mld, color ='tab:blue'  , linestyle =':', alpha=1.0 , linewidth=3 )

    ax.set_xlabel(r'${\rm m^2}\;{\rm s}^{-2}$')

    # ===============================================================
        # ===============================================================
    ax_index+=1
    ax = axes.flat[ax_index]

    ax.set_title(r'$\overline{w^\prime v^\prime}$')

    ax.plot(-(WV)[instant], z_r_les/mld, style_les,
    # ax.plot((WU)[instant], z_r_les/mld, style_les,
            alpha=alpha_les, linewidth=linewidth_les, label='LES')
    
    for i, label in enumerate(run_label):
        ax.plot( out[i]['WV'][instant,:], out[i].z_w/mld, color ='tab:blue'  , linestyle ='-', alpha=1.0 , linewidth=3 )
        # ax.plot( out[i]['WV_MF'][instant,:], out[i].z_w/mld, color ='tab:blue'  , linestyle ='--', alpha=1.0 , linewidth=3 )
        # ax.plot( out[i]['WV_ED'][instant,:], out[i].z_w/mld, color ='tab:blue'  , linestyle =':', alpha=1.0 , linewidth=3 )

    ax.set_xlabel(r'${\rm m^2}\;{\rm s}^{-2}$')

    # ===============================================================

    # ax_index+=1
    # ax = axes.flat[ax_index]

    # ax.set_xlabel(r'$s^{-1}$')
    # ax.set_title(r'$\omega_p$')


    # for i, label in enumerate(run_label):
    #     if label=='EDMF-Energy-cor':
    #         # ax.plot((scm[i].vortp), scm[i].z_w/mld, linestyle=styles[i], color = colors[i],
    #         #     alpha=alpha[i], linewidth=linewidth[i], label=label)
    #         ax.plot((out[i].vort_p)[instant], out[i].z_w/mld, linestyle=styles[i], color = colors[i],
    #             alpha=alpha[i], linewidth=linewidth[i], label=label)

    # #ax.set_ylim((-1.3, 0))

    # ===============================================================

    ax_index+=1
    ax = axes.flat[ax_index]

    ax.set_xlabel(r'$m s^{-1}$')
    ax.set_title(r'$w_p$')


    for i, label in enumerate(run_label):
        if (label=='EDMF-Energy-cor') or (label=='EDMF-Energy'):
            True
            # ax.plot((scm[i].wp), scm[i].z_w/mld, linestyle=styles[i], color = colors[i],
            #     alpha=alpha[i], linewidth=linewidth[i], label=label)
            
            # ax.plot((out[i].w_p)[instant], out[i].z_w/mld, linestyle=styles[i], color = colors[i],
                # alpha=alpha[i], linewidth=linewidth[i], label=label)

    #ax.set_ylim((-1.3, 0))

    # ===============================================================
    # ===============================================================

    ax_index+=1
    ax = axes.flat[ax_index]

    ax.set_xlabel(r'$m s^{-1}$')
    ax.set_title(r'$w_p$')

    for i, label in enumerate(run_label):
        if (label=='EDMF-Energy-cor') or (label=='EDMF-Energy'):
            True
            # ax.plot((scm[i].ent), scm[i].z_r/mld, linestyle='-', color = colors[i],
            #     alpha=alpha[i], linewidth=linewidth[i], label=label)
            # ax.plot((scm[i].det), scm[i].z_r/mld, linestyle='--', color = colors[i],
            #     alpha=alpha[i], linewidth=linewidth[i], label=label)
            # ax.plot((scm[i].ent-scm[i].det), scm[i].z_r/mld, linestyle=':', color = colors[i],
            # alpha=alpha[i], linewidth=linewidth[i], label=label)
            # ax.plot(out[i]['Ent'][instant], out[i].z_r/mld, linestyle=styles[i], color = colors[i],
            #     alpha=alpha[i], linewidth=linewidth[i], label=label)
            # ax.plot(out[i]['Det'][instant], out[i].z_r/mld, linestyle=styles[i], color = colors[i],
            #     alpha=alpha[i], linewidth=linewidth[i], label=label)
            # ax.plot((out[i]['Ent'][instant]-out[i]['Det'][instant]), out[i].z_r/mld, linestyle=styles[i], color = colors[i],
            #     alpha=alpha[i], linewidth=linewidth[i], label=label)
    # ax.set_xlim(-0.01,0.01)
    #ax.set_ylim((-1.3, 0))

    # ===============================================================

    # adding subplot labels
    subplot_label = [r'\rm{(a)}', r'\rm{(b)}', r'\rm{(c)}',
                    r'\rm{(d)}', r'\rm{(e)}', r'\rm{(f)}',r'\rm{(g)}',r'\rm{(h)}',r'\rm{(i)}',r'\rm{(j)}',r'\rm{(k)}',r'\rm{(l)}']

    for i,ax in enumerate(axes.flat):
        # ax.set_ylim((-2, 0))
        ax.set_box_aspect(1)
        ax.text(0.15, 0.98, subplot_label[i], transform=ax.transAxes,
                fontsize=16, bbox=dict(facecolor='1.', edgecolor='none', pad=3.0), fontweight='bold', va='top', ha='right')




    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(
        0.5, -0.05), fancybox=False, shadow=False, ncol=4)

    fig.tight_layout()

    plt.savefig(saving_path+saving_name, bbox_inches='tight', dpi=300)
    print('figure saved at'+saving_path+saving_name)

plot_intant_panel()

### Plot velocity panel
subprocess.run(["python", "WANG1_LES_vs_SCM_velocities.py"])

#plt.show()


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
# plt.rcParams.update({'font.size': 18})
plt.rcParams.update({'figure.facecolor': 'white'})
plt.rcParams.update({'savefig.facecolor': 'white'})
###########################################


# colors
blue, orange, magenta, grey, green = '#0db4c3', '#eea021', '#ff0364', '#606172', '#3fb532'

# ===========================================================================
cases = ['WANG1_FR_Q2000_dT1e-3_lat60_120h']
les = {}
LG_MEAN={}
LG_RES={}
LG_SBG={}
BU_KE={}
TH_les={}
U_les={}
V_les={}
WTH={}
WU={}
WV={}
TKE={}
WTKE={}
for case in cases:
    file = 'GN_01.1.OC_01.000.nc'
    path = '../../tests/data/'+case+'/'
    les[case] = xr.open_dataset(path+file)
    LG_MEAN[case]= xr.open_dataset(path+file,group ='/LES_budgets/Mean/Cartesian/Not_time_averaged/Not_normalized/cart')
    LG_RES[case] = xr.open_dataset(path+file,group ='/LES_budgets/Resolved/Cartesian/Not_time_averaged/Not_normalized/cart')    
    LG_SBG[case] = xr.open_dataset(path+file,group ='/LES_budgets/Subgrid/Cartesian/Not_time_averaged/Not_normalized/cart')
    BU_KE[case] = xr.open_dataset(path+file,group ='/LES_budgets/BU_KE/Cartesian/Not_time_averaged/Not_normalized/cart')
    TH_les[case] = (LG_MEAN[case].MEAN_TH - 273.15).data
    U_les [case] = (LG_MEAN[case].MEAN_U).data
    V_les [case] = (LG_MEAN[case].MEAN_V).data
    WTH   [case] = (LG_RES[case].RES_WTH + LG_SBG[case].SBG_WTHL).data
    WU    [case] = (LG_RES[case].RES_WU  + LG_SBG[case].SBG_WU).data
    WV    [case] = (LG_RES[case].RES_WV  + LG_SBG[case].SBG_WV).data
    TKE   [case] = (LG_RES[case].RES_KE  + LG_SBG[case].SBG_TKE).data 
    WTKE  [case] = (LG_RES[case].RES_WKE + LG_SBG[case].SBG_WTKE).data

time_les = les[case].time_les

# numpy array of integer hours, starting at inital time + 1h
time = ((time_les - time_les[0]) / np.timedelta64(1, 'h')).data.astype(int) + 1

# remap level_les on negative depth values
z_r_les = (les[case].level_les - (les[case].level_les[0] + les[case].level_les[-1])).data

# choose intant to plot    
instant = -1 
### ATTENTION: scm contains time 0 + 276 hours = 277 instants
###            les contains hour 1 to 276 = 276 instants
### so SAME END but different start in the lsit

# compute LES mixed layer depth
# to non-dimensionalize depths
mld = (-z_r_les[(-WTH[cases[0]][instant]).argmax()])

# ===========================================================================

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
    'entr_scheme': 'R10corNT',
    'trad_coriolis_mod': True,
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
    'wp0':-1e-2,
    'write_netcdf': True
}

run_label = cases

runs = [ {'output_filename': 'scm_'+case+'.nc'} for case in cases ]

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
#     {
#         'output_filename': 'scm_'+cases[4]+'.nc',
#     },
#         ]



scm = [0]*len(runs)

# Run the SCM
for i, run_params in enumerate(runs):
    params = edmf_params.copy()  # Create a copy of common_params
    params.update(run_params)  # Update with run_params
    params.update(case_params[cases[i]])
    scm[i] = SCM(**params)
    scm[i].run_direct()
    # test zinv
    if scm[i].MF_tra or scm[i].MF_dyn: 
        print(run_label[i])
        reference=mld
        is_in_range(value=scm[i].zinv, value_name='zinv', reference=reference,tolerance=10 )

# LOAD outputs
out = [0]*len(runs)
out_dic={}
for i, run_params in enumerate(runs):
    print('opening '+run_params['output_filename'])
    out[i] = xr.open_dataset(run_params['output_filename'])
    out_dic[run_label[i]] = out[i]  
# choose intant to plot    
instant = -1 
### ATTENTION: scm contains time 0 + 276 hours = 277 instants
###            les contains hour 1 to 276 = 276 instants
### so SAME END but different start in the lsit

# compute LES mixed layer depth
# to non-dimensionalize depths
mld = (-z_r_les[(-WTH[cases[0]][instant]).argmax()])
# mld=1

################################# PLOTTING #################################
#SCM
styles = {case: '-' for case in cases}
colors = {cases[0]: 'k'}
# colors = {cases[0]:'tab:gray' ,cases[1]:'tab:orange' , cases[2]:'tab:green' , cases[3]:'tab:blue'}
alpha  = {case: 1 for case in cases}
# linewidth = [2]*(len(run_label))

#LES

styles_les = {case: '--' for case in cases}
colors_les = colors
alpha_les  = {case:0.3 for case in cases}
marker_les='o'
width_les = {case:3 for case in cases}
 





def plot_intant_panel(instant=instant):
    #============================================ WC ===============================================
    fig, axes = plt.subplots(nrows=1, ncols=3, sharex=False,
                            sharey=True, constrained_layout=True)

    ax_index=-1

    # ===============================================================
    ax_index+=1
    ax = axes.flat[ax_index]
    ax.set_xlabel(r'$^{\circ}{\rm C}$')
    ax.set_ylabel(r'$z / h $')
    ax.set_title(r'$\overline{\theta}$')


    for case in cases:
        
        ax.plot(out_dic[case].temp[instant], out_dic[case].z_r/mld, linestyle=styles[case], color = colors[case], alpha=alpha[case], label=case)
        
        ax.plot(TH_les[case][instant], z_r_les/mld, colors_les[case], 
            alpha=alpha_les[case],linewidth=width_les[case])

    ax.set_xlim((9.65,9.85))
    ax.set_ylim((-1.3, 0))
    # ===============================================================
    ax_index+=1
    ax = axes.flat[ax_index]

    ax.set_title(r'$\overline{w^\prime \theta^\prime}$')

    for case in cases:
    
        ax.plot(-out_dic[case]['WT'][instant], out_dic[case].z_w/mld, linestyle=styles[case], color = colors[case], alpha=alpha[case], label=case)
    
        ax.plot(WTH[case][instant], z_r_les/mld, colors_les[case], 
        alpha=alpha_les[case],linewidth=width_les[case])

    ax.set_xlabel(r'${\rm K}\;{\rm m}\;{\rm s}^{-1}$')

    # ===============================================================
    # ===============================================================
    ax_index+=1
    ax = axes.flat[ax_index]
    ax.set_title(r'$k$')

    for case in cases:
    
        ax.plot(out_dic[case]['tke'][instant], out_dic[case].z_w/mld, linestyle=styles[case], color = colors[case], alpha=alpha[case], label=case)
    
        ax.plot(TKE[case][instant], z_r_les/mld, colors_les[case], 
        alpha=alpha_les[case],linewidth=width_les[case])

    ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

    #ax.set_xlim((-0.0001, 0.001))
    #ax.set_ylim((-1.3, 0))

    ax.set_xlabel(r'${\rm m}^2\;{\rm s}^{-2}$')
    # ===============================================================
    # ===============================================================

    # adding subplot labels
    subplot_label = [r'\rm{(a)}', r'\rm{(b)}', r'\rm{(c)}',
                    r'\rm{(d)}', r'\rm{(e)}', r'\rm{(f)}',r'\rm{(g)}',r'\rm{(h)}',r'\rm{(i)}',r'\rm{(j)}',r'\rm{(k)}',r'\rm{(l)}']

    for i,ax in enumerate(axes.flat):
        # ax.set_ylim((-2, 0))
        ax.set_box_aspect(1)
        ax.text(0.15, 0.95, subplot_label[i], transform=ax.transAxes, bbox=dict(facecolor='1.', edgecolor='none'), fontweight='bold', va='top', ha='right')


        # ax.text(0.15, 0.95, subplot_label[i], transform=ax.transAxes, bbox=dict(facecolor='1.', edgecolor='none', pad=3.0), fontweight='bold', va='top', ha='right')

    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(
        0.5, 0.2), fancybox=False, shadow=False, ncol=4)

    # fig.tight_layout()

    saving_path = '../figures/'
    saving_name = 'NT_new_test.png'    
    plt.savefig(saving_path+saving_name, bbox_inches='tight', dpi=300)
    print('figure saved at'+saving_path+saving_name)

# mld=1
plot_intant_panel(instant=-1)

def several_instants_temperature(nb_plots=4):
    fig,axs = plt.subplots(nrows=1,ncols=1,constrained_layout=True)

    ax=axs
    for k in range(nb_plots):
        instant=len(time)//(nb_plots)*(k+1) -1    
        # p=ax.plot(TH_les['WANG1_NR'][instant],z_r_les,alpha=0.5)
        # color = p[0].get_color()
        # scm='No rotation'
        # ax.plot(out_dic[scm]['temp'][instant+1], out_dic[scm]['z_r'],color=color,linestyle='-')
        
        p=ax.plot(TH_les['WANG1_FR'][instant],z_r_les,alpha=0.5)
        color = p[0].get_color()
        scm='Traditional Coriolis'
        ax.plot(out_dic[scm]['temp'][instant+1], out_dic[scm]['z_r'],color=color,linestyle='-',label='EDMF')

        scm='EVD'
        ax.plot(out_dic[scm]['temp'][instant+1], out_dic[scm]['z_r'],color='k',linestyle='-',label='EVD')
        # ax.plot(out['EVD']['temp'][instant], out['EVD']['z_r'],color=color,linestyle='-',alpha=0.5)
    
    instant=-1
    p=ax.plot(TH_les['WANG1_NR'][instant],z_r_les,alpha=0.5)
    color = p[0].get_color()
    scm='No rotation'
    ax.plot(out_dic[scm]['temp'][instant], out_dic[scm]['z_r'],color=color,linestyle='--')
    plt.show()
# several_instants_temperature()

# run_label = ['No rotation','Traditional Coriolis','EVD' ]
### Plot velocity panel
subprocess.run(["python", "WANG1_LES_vs_SCM_velocities_mean_only.py"])

plt.show()
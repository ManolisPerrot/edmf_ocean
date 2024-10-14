

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
from cases_exhaustive_settings import case_settings
from matplotlib.axes import Axes

###################################################
plt.rcParams['text.usetex'] = True
# plt.rcParams.update({'font.size': 18})
plt.rcParams.update({'figure.facecolor': 'white'})
plt.rcParams.update({'savefig.facecolor': 'white'})
###########################################


# colors
blue, orange, magenta, grey, green = '#0db4c3', '#eea021', '#ff0364', '#606172', '#3fb532'

# ===========================================================================
case = 'W005_C500_NO_COR'


saving_path = '../figures/'
saving_name = case+'scheme_comparison.png'


# loading LES output once before using the function

file = 'GN_01.1.OC_01.000_copy.nc'


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
instant=-1
mld = (-z_r_les[(-WTH[instant]).argmax()])

# ===========================================================================

#load physical parameters of the case
physical_params = case_settings[case].copy()
# Define SCM and mass-flux default parameters (overwrites physical_params):
scm_params = {
    'nz': 100,
    'dt': 50.,
    'h0': 2000.,
    'thetas': 6.5,
    'hc': 400,
    'nbhours': 72,
    'outfreq': 1,
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
    'wp_a': 1.3,   #1
    'wp_b': 1.3,     #1.
    'wp_bp': 0.003*250,    
    'up_c': 0.5, 
    'vp_c': 0.5,
    'bc_ap': 0.2,    #0.3,
    'delta_bkg': 0.009*250,   # 0.005,
    'wp0' : -0.5e-08,
    'output_filename': 'run',
    'beta_bc_P09': 0.3,
    'write_netcdf': True
}


# Define parameters specific to each run (overwrite previous parameters):
run_label = [r'TKE+EVD',r'EDMF',r'$k-\epsilon$ (GOTM)']
runs = [
        {
        'eddy_diff': True,
        'evd': True,
        'mass_flux_tra': False,
        'mass_flux_dyn': False,
        'mass_flux_tke': False,
        'mass_flux_tke_trplCorr': False,
        'output_filename': case+'_EVD.nc'

    },
    {
        'eddy_diff': True,
        'evd': False,
        'mass_flux_tra': True,
        'mass_flux_dyn': True,
        'mass_flux_tke': True,
        'mass_flux_tke_trplCorr': True,
        'output_filename': case+'_EDMF.nc'

    },
        {
        'eddy_diff': True,
        'evd': False,
        'eddy_diff_scheme': 'Keps',
        'mass_flux_tra': False,
        'mass_flux_dyn': False,
        'mass_flux_tke': False,
        'mass_flux_tke_trplCorr': False,
        'output_filename': case+'_Keps.nc'
    },
        ]

scm = [0]*len(runs)

# Run the SCM
for i, run_params in enumerate(runs):
    params = physical_params.copy()  # Create a copy of common_params
    params.update(scm_params)  # Update with run_params
    params.update(run_params)  # Update with run_params
    scm[i] = SCM(**params)
    scm[i].run_direct()

    # test zinv
    if scm[i].MF_tra or scm[i].MF_dyn: 
        print(run_label[i])
        reference=-320
        #reference=-mld #LES value
        is_in_range(value=scm[i].zinv, value_name='zinv', reference=reference,tolerance=10 )

#========= Opening KPP outputs
nc_file='../data/W005_C500_NO_COR/W005_FC500_kpp_cvmix.nc'
fh1      = Dataset(nc_file, mode='r')
temp_kpp  = fh1.variables['temp'][:,:,0,0]
akt_kpp   = fh1.variables['nuh'] [:,:,0,0]
gam_kpp   = fh1.variables['gamh'][:,:,0,0]
u_kpp     = fh1.variables['u']   [:,:,0,0]
akv_kpp   = fh1.variables['num'] [:,:,0,0]
gamu_kpp   = fh1.variables['gamu'][:,:,0,0]
zr_kpp    = fh1.variables['z'][0,:,0,0]
zw_kpp    = fh1.variables['zi'][0,:,0,0]
nz = len(zr_kpp); dz = zw_kpp[1]-zw_kpp[0]
wt_kpp = np.zeros((len(temp_kpp[:,0]),nz-1))
wt_kpp = -akt_kpp[:,1:nz-1]*(temp_kpp[:,1:nz-1]-temp_kpp[:,0:nz-2])/dz + gam_kpp[:,1:nz-1]
wu_kpp = np.zeros((len(temp_kpp[:,0]),nz-1))
wu_kpp = -akv_kpp[:,1:nz-1]*(u_kpp[:,1:nz-1]-u_kpp[:,0:nz-2])/dz + gamu_kpp[:,1:nz-1]
fh1.close()
#====================================
# LOAD outputs

# TH_scm = [0]*len(runs)
# U_scm = [0]*len(runs)
# V_scm = [0]*len(runs)

# #interpolate scm output on LES #TODO do the converse to reduce computation cost?


# for i, run_params in enumerate(runs):
#     TH_scm = scm[i].t_history
#     U_scm  = scm[i].u_history
#     V_scm  = scm[i].v_history

instant = 71
mld = (-z_r_les[(-WTH[instant]).argmax()]).data


################################# PLOTTING
styles = ['-', '-', '-','--',':']
#colors = ['k',blue,orange]
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
alpha = [1]*len(run_label)
linewidth = [2]*(len(run_label))

style_les = 'k'
alpha_les = 0.5
linewidth_les = 0
marker_les='+'


# run_label = [r'EDMF']
plot_kpp=False
#============================================ WC ===============================================
if case == 'W005_C500_NO_COR':
    fig, axes = plt.subplots(nrows=2, ncols=3, sharex=False,
                         sharey=True, constrained_layout=True)
# ===============================================================
    ax = axes.flat[0]
    ax.set_xlabel(r'$^{\circ}{\rm C}$')
    ax.set_ylabel(r'$z / h $')
    ax.set_title(r'$\overline{\theta}$')


    ax.plot(TH_les[instant], z_r_les/mld, style_les,
            alpha=alpha_les, linewidth=linewidth_les, marker=marker_les,  label='LES')

    for i, label in enumerate(run_label):
        if True:         ax.plot(scm[i].t_np1[:, 0], scm[i].z_r/mld, linestyle=styles[i], color = colors[i],
                alpha=alpha[i], linewidth=linewidth[i], label=label)
    if plot_kpp:
        ax.plot(temp_kpp[-1,:], zr_kpp/mld, linestyle=styles[i], color = 'tab:purple',
                alpha=alpha[i], linewidth=linewidth[i], label=r'KPP (CVMix)')
    ax.set_xlim((1.71, 1.755))
    # ax.set_ylim((-1.3, 0))


    # ===============================================================

    ax = axes.flat[1]
    ax.set_title(r'$\overline{u}$')

    ax.set_xlabel(r'${\rm m}\;{\rm s}^{-1}$')

    ax.plot(U_les[instant], z_r_les/mld, style_les,
            alpha=alpha_les, linewidth=linewidth_les, marker=marker_les, label='LES')

    for i, label in enumerate(run_label):
        if True:         ax.plot((scm[i].u_np1), scm[i].z_r/mld, linestyle=styles[i], color = colors[i],
                alpha=alpha[i], linewidth=linewidth[i], label=label)
    
    if plot_kpp: 
        ax.plot(u_kpp[-1,:], zr_kpp/mld, linestyle=styles[i], color = 'tab:purple',
                alpha=alpha[i], linewidth=linewidth[i], label=r'KPP (CVMix)')
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    ax.set_ylim((-1.3, 0))


    # ===============================================================


    ax = axes.flat[3]
    ax.set_title(r'$\overline{w^\prime \theta^\prime}$')

    ax.plot(WTH[instant], z_r_les/mld, style_les,
            alpha=alpha_les, linewidth=linewidth_les, marker=marker_les, label='LES')

    for i, label in enumerate(run_label):
        if True:        
            if run_label == 'EVD':
                ax.plot(-(scm[i].wted), scm[i].z_w/mld, styles[i], color = colors[i],
                   alpha=alpha[i], linewidth=linewidth[i], label=label)
            else:
                ax.plot(-(scm[i].wted + scm[i].wtmf), scm[i].z_w/mld, styles[i], color = colors[i],
                   alpha=alpha[i], linewidth=linewidth[i], label=label)
    if plot_kpp:
        ax.plot(wt_kpp[-1,:], zw_kpp[1:nz-1]/mld,
            linestyle=styles[i], color = 'tab:purple',
                alpha=alpha[i], linewidth=linewidth[i], label=r'KPP (CVMix)')
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

    # for i, label in enumerate(run_label[1:]):
    #     ax.plot(-(scm[i].wted ), scm[i].z_w/mld, linestyle='--', color = colors[i],
    #                alpha=alpha[i], linewidth=linewidth[i], label=label)
    # ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

    # ax.plot( -scm[1].wted, scm[1].z_w/mld, color ='tab:red'  , linestyle ='-.', alpha=1.0 , linewidth=2 )
    # ax.plot( -scm[2].wted, scm[2].z_w/mld, color ='tab:green', linestyle ='-.', alpha=1.0 , linewidth=2 )
    # ax.plot( -out[1]['WT_ED'][-1,:], out[1].z_w/mld, color ='tab:blue'  , linestyle ='--', alpha=1.0 , linewidth=3 )
    # ax.plot( -out[2]['WT_ED'][-1,:], out[2].z_w/mld, color ='tab:orange', linestyle ='--', alpha=1.0 , linewidth=3 )


    ax.set_ylim((-1.3, 0))

    ax.set_xlabel(r'${\rm K}\;{\rm m}\;{\rm s}^{-1}$')


    # ===============================================================
    # ===============================================================
    ax = axes.flat[2]
    ax.set_title(r'$k$')

    ax.plot(TKE[instant], z_r_les/mld, style_les,
            alpha=alpha_les, linewidth=linewidth_les, marker=marker_les, label='LES')

    for i, label in enumerate(run_label):
        if True:         ax.plot(scm[i].tke_np1, scm[i].z_w/mld, linestyle=styles[i], color = colors[i],
                alpha=alpha[i], linewidth=linewidth[i], label=label)
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

    ax.set_xlim((-0.0001, 0.001))
    ax.set_ylim((-1.3, 0))

    ax.set_xlabel(r'${\rm m}^2\;{\rm s}^{-2}$')
    # ===============================================================


    ax = axes.flat[5]
    ax.set_title(r'$\overline{w^\prime \frac{\mathbf{u}^{\prime 2}}{2}  } + \frac{1}{\rho_0} \overline{w^\prime p^{\dagger \prime} }$')
#add velocity-pressure correlation
    cond_samp = xr.open_dataset(path+    'W005_C500_NO_COR_object_diags_Cw_m1_72h.nc')  

    ax.plot(WTKE[instant]-cond_samp['TOT_intra_WPHI_over_RHO_0'][-1], z_r_les/mld, style_les,
            alpha=alpha_les, linewidth=linewidth_les, marker=marker_les, label='LES')
    # ax.plot(WTKE[instant], z_r_les/mld, style_les,
    #         alpha=alpha_les, linewidth=linewidth_les, marker=marker_les, label='LES')

    for i, label in enumerate(run_label):
        if True:         ax.plot((scm[i].wtke), scm[i].z_r/mld, linestyle=styles[i], color = colors[i],
                alpha=alpha[i], linewidth=linewidth[i], label=label)


    ax.set_xlim((- 1.5e-5, 2e-6))
    ax.set_ylim((-1.3, 0))

    ax.set_xlabel(r'${\rm m}^3\;{\rm s}^{-3}$')

    # ===============================================================
    ax = axes.flat[4]

    ax.set_xlabel(r'${\rm m}^2\;{\rm s}^{-2}$')
    ax.set_title(r'$\overline{w^\prime u^\prime}$')

    ax.plot(-WU[instant], z_r_les/mld, style_les,alpha=alpha_les, linewidth=linewidth_les, marker=marker_les, label='LES')


    for i, label in enumerate(run_label):
        if True:             ax.plot((scm[i].wued + scm[i].wumf), scm[i].z_w/mld, linestyle=styles[i], color = colors[i],
                    alpha=alpha[i], linewidth=linewidth[i], label=label)
    if plot_kpp:
        ax.plot(-wu_kpp[-1,:], zw_kpp[1:nz-1]/mld,
             linestyle=styles[i], color = 'tab:purple',
                alpha=alpha[i], linewidth=linewidth[i], label=r'KPP (CVMix)')
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

    # ax.plot( scm[1].wued, scm[1].z_w/mld,  color ='tab:red'  , linestyle ='-.', alpha=1.0 , linewidth=2 ) 
    # ax.plot( scm[2].wued, scm[2].z_w/mld,  color ='tab:green', linestyle ='-.', alpha=1.0 , linewidth=2 )

    ax.set_ylim((-1.3, 0))



    # ===============================================================




    # adding subplot labels
    subplot_label = [r'\rm{(a)}', r'\rm{(b)}', r'\rm{(c)}',
                    r'\rm{(d)}', r'\rm{(e)}', r'\rm{(f)}']

    for i,ax in enumerate(axes.flat):
        # ax.set_box_aspect(1)
        ax.text(0.15, 0.98, subplot_label[i], transform=ax.transAxes,
                 bbox=dict(facecolor='1.', edgecolor='none', pad=3.0), fontweight='bold', va='top', ha='right')




    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(
        0.5, -0.008), fancybox=False, shadow=False, ncol=5)

    # fig.tight_layout()

plt.savefig(saving_path+saving_name, bbox_inches='tight', dpi=300)
print('figure saved at '+saving_path+saving_name)
# plt.show()
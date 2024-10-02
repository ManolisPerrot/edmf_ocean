#!/usr/bin/env python
# coding: utf-8
__docformat__ = 'reStructuredText'
###########################################
# Imports
###########################################


from sys import exit
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import scipy.signal
from scipy.interpolate import interp1d
import xarray as xr
import time as TIME
import xrft as xrft
from matplotlib.ticker import MaxNLocator


plt.rcParams['text.usetex'] = True
#plt.rcParams.update({'font.size': 18})
plt.rcParams.update({'figure.facecolor':'white'})
plt.rcParams.update({'savefig.facecolor':'white'})
###########################################

import numpy as np

## define Coriolis frequency
lat0=60
omega          = 7.292116e-05; rad            = np.pi/180.
fcor      = 2.*omega*np.sin(rad*lat0)
ecor      = 2.*omega*np.cos(rad*lat0)
g=9.809
alpha=2.6e-4
rho0=1027
cp=4178
stflx =111
B0 = g*alpha/(rho0*cp)*stflx
# B0=5e-8*10
wstar = (B0*1000)**(1/3)

# # Ro = np.maximum((np.abs(B0)/fcor)**(0.5)/(fcor*np.abs(zinv)) , (np.abs(B0)/fcor)**(0.5)/(fcor*1000.))
# Ro = (np.abs(B0)/fcor)**(0.5)/(fcor*500.)

# cff= np.tanh(Ro**0.37)
# # print(cff, Ro)
# # mf_params[-1] = delta_bkg/cff
# dTdz = 1/(g*alpha)*(16*fcor)**2 





# ===========================================================================
# cases = ['WANG1_FR_lat60']
# cases = ['WANG1_FR_Q2000_lat60_120h']
cases = ['WANG1_FR_Q2000_dT1e-3_lat60_120h']
# cases = ['WANG1_FR_lat30','WANG1_FR_lat60']


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
dz_WU_les ={}
dz_WV_les ={}
for case in cases:
    file = 'GN_01.1.OC_01.000.nc'
    path = '../../tests/data/'+case+'/'
    les[case] = xr.open_dataset(path+file)
    LG_MEAN[case]= xr.open_dataset(path+file,group ='/LES_budgets/Mean/Cartesian/Not_time_averaged/Not_normalized/cart')
    LG_RES[case] = xr.open_dataset(path+file,group ='/LES_budgets/Resolved/Cartesian/Not_time_averaged/Not_normalized/cart')    
    LG_SBG[case] = xr.open_dataset(path+file,group ='/LES_budgets/Subgrid/Cartesian/Not_time_averaged/Not_normalized/cart')
    BU_KE[case] = xr.open_dataset(path+file,group ='/LES_budgets/BU_KE/Cartesian/Not_time_averaged/Not_normalized/cart')
    TH_les[case] = (LG_MEAN[case].MEAN_TH - 273.15)
    V_les [case] = (LG_MEAN[case].MEAN_U)
    U_les [case] = (LG_MEAN[case].MEAN_V)
    WTH   [case] = (LG_RES[case].RES_WTH + LG_SBG[case].SBG_WTHL).data
    WV    [case] = (LG_RES[case].RES_WU  + LG_SBG[case].SBG_WU)
    WU    [case] = (LG_RES[case].RES_WV  + LG_SBG[case].SBG_WV)
    TKE   [case] = (LG_RES[case].RES_KE  + LG_SBG[case].SBG_TKE).data 
    WTKE  [case] = (LG_RES[case].RES_WKE + LG_SBG[case].SBG_WTKE).data
    z_r_les = (les[case].level_les - (les[case].level_les[0] + les[case].level_les[-1])).data
    dz_WU_les[case] = (WU[case][:].mean(dim='time_les').data[1:]-WU[case][:].mean(dim='time_les').data[:-1])/(z_r_les[1:]-z_r_les[:-1])
    dz_WV_les[case] = (WV[case][:].mean(dim='time_les').data[1:]-WV[case][:].mean(dim='time_les').data[:-1])/(z_r_les[1:]-z_r_les[:-1])

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




run_label = cases

runs = [ {'output_filename': 'scm_'+case+'.nc'} for case in cases ]


# runs = [
#     {
#         # 'entr_scheme': 'R10',
#         'output_filename': 'scm_'+cases[0]+'.nc',
#     },
#     #     {
#     #     'output_filename': 'scm_'+cases[1]+'.nc',
#     # },
#     #         {
#     #     'output_filename': 'scm_'+cases[2]+'.nc',
#     # },
#     #         {
#     #     'output_filename': 'scm_'+cases[3]+'.nc',
#     # },
#         ]


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

dz_WU_scm={}
dz_WV_scm={}
for case in cases:
    dz_WU_scm[case] = -(out_dic[case]['WU'][:].mean(dim='time').data[1:]-out_dic[case]['WU'][:].mean(dim='time').data[:-1])/(out_dic[case]['z_w'].data[1:]-out_dic[case]['z_w'].data[:-1])
    dz_WV_scm[case] = -(out_dic[case]['WV'][:].mean(dim='time').data[1:]-out_dic[case]['WV'][:].mean(dim='time').data[:-1])/(out_dic[case]['z_w'].data[1:]-out_dic[case]['z_w'].data[:-1])



################################ PLOTTING #################################
#SCM
styles = {case: '-' for case in cases}
colors = {cases[0]:'k'}

# colors = {cases[0]:'tab:gray' ,cases[1]:'tab:orange'}
# colors = {cases[0]:'tab:gray' ,cases[1]:'tab:orange' , cases[2]:'tab:green' , cases[3]:'tab:blue'}

# colors = { cases[0]:'tab:blue'}
alpha = {case: 1 for case in cases}
# linewidth = [2]*(len(run_label))

#LES

styles_les = {case: '-' for case in cases}
colors_les = colors
alpha_les  = {case:0.3 for case in cases}
marker_les='o'
width_les = {case:3 for case in cases}




def plot_mean_velocities():
    #============================================ WC ===============================================
    fig, axes = plt.subplots(nrows=1, ncols=2, sharex=False,
                            sharey=True, constrained_layout=True)

    ax_index=-1

    # ===============================================================
    ax_index+=1
    ax = axes.flat[ax_index]
    # ax.set_xlabel(r'$\mathrm{m s^{-1}}$')
    ax.set_ylabel(r'$z / h $')
    ax.set_title(r'Zonal velocity')

    for case in cases:
        
        ax.plot(out_dic[case].u[:].mean(dim='time')/wstar, out_dic[case].z_r/mld, linestyle=styles[case], color = colors[case], alpha=alpha[case], label=r'$\frac{\langle \overline{u} \rangle_{t}}{w_*}$, SCM')
        ax.plot(-1/fcor*dz_WV_scm[case]/wstar, out_dic[case].z_w[1:]/mld, linestyle=':', color = colors[case], alpha=alpha[case], label=r'$-\frac{1}{f} \partial_z \frac{\langle \overline{w^\prime v^\prime} \rangle_t}{w_*}$, SCM')
        
        ax.plot(U_les[case][:].mean(dim='time_les')/wstar, z_r_les/mld, colors_les[case], 
            alpha=alpha_les[case],linewidth=width_les[case], label=r'$\frac{ \langle \overline{u} \rangle_{t} }{w_*}$, LES')
        ax.plot(1/fcor*dz_WV_les[case]/wstar, z_r_les[1:]/mld, colors_les[case], 
            alpha=alpha_les[case],linewidth=width_les[case],linestyle=':',label=r'$-\frac{1}{f} \partial_z \frac{\langle \overline{w^\prime v^\prime} \rangle_t}{w_*}$, LES')

        # var = -( 1/(fcor*wstar)*out_dic[case].a_p*(out_dic[case].w_p)**2 ).diff(dim='z_w')[:].mean(dim='time')
        # ax.plot(var, out_dic[case].z_r/mld)
    # ax.set_xlim((-0.01, 0.06))
    ax.set_ylim((-1.3, 0))
    ax.legend()

    # ===============================================================
    ax_index+=1
    ax = axes.flat[ax_index]
    # ax.set_xlabel(r'$\mathrm{m s^{-1}}$')
    # ax.set_ylabel(r'$z / h $')
    ax.set_title(r'Meridonal velocity')

    for case in cases:
        
        ax.plot(out_dic[case].v[:].mean(dim='time')/wstar, out_dic[case].z_r/mld, linestyle=styles[case], color = colors[case], alpha=alpha[case], label=r'$\frac{\langle \overline{v} \rangle_{t}}{w_*}$, SCM')
        ax.plot(1/fcor*dz_WU_scm[case]/wstar, out_dic[case].z_w[1:]/mld, linestyle=':', color = colors[case], alpha=alpha[case], label=r'$\frac{1}{f} \partial_z \frac{\langle \overline{w^\prime u^\prime} \rangle_t}{w_*} $, SCM')

        
        ax.plot(V_les[case][:].mean(dim='time_les')/wstar, z_r_les/mld, colors_les[case], 
            alpha=alpha_les[case],linewidth=width_les[case],label=r'$\frac{\langle \overline{v} \rangle_{t}}{w_*}$, LES')
        ax.plot(-1/fcor*dz_WU_les[case]/wstar, z_r_les[1:]/mld, colors_les[case], 
            alpha=alpha_les[case],linewidth=width_les[case],linestyle=':',label=r'$\frac{1}{f} \partial_z \frac{ \langle \overline{w^\prime u^\prime} \rangle_t}{w_*} $, LES')
    # ax.set_xlim((2.93, 2.98))
    ax.set_ylim((-1.3, 0))
    ax.legend()



    # ===============================================================

    # adding subplot labels
    subplot_label = [r'\rm{(a)}', r'\rm{(b)}', r'\rm{(c)}',
                    r'\rm{(d)}', r'\rm{(e)}', r'\rm{(f)}',r'\rm{(g)}',r'\rm{(h)}',r'\rm{(i)}',r'\rm{(j)}',r'\rm{(k)}',r'\rm{(l)}']

    for i,ax in enumerate(axes.flat):
        # ax.set_ylim((-2, 0))
        # ax.set_box_aspect(1)
        ax.text(0.1, 0.97, subplot_label[i], transform=ax.transAxes, bbox=dict(facecolor='1.', edgecolor='none'), fontweight='bold', va='top', ha='right')


        # ax.text(0.15, 0.95, subplot_label[i], transform=ax.transAxes, bbox=dict(facecolor='1.', edgecolor='none', pad=3.0), fontweight='bold', va='top', ha='right')

    # handles, labels = axes.flat[0].get_legend_handles_labels()
    # fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(
    #     0.5, 0.2), fancybox=False, shadow=False, ncol=4)

    # fig.tight_layout()

    saving_path = '../figures/'
    saving_name = cases[0]+'_mean_velocities.png'    
    plt.savefig(saving_path+saving_name, bbox_inches='tight', dpi=300)
    print('figure saved at'+saving_path+saving_name)
    plt.show()

plot_mean_velocities()


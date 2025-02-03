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
from load_MNH import load_MNH
import matplotlib as mpl


plt.rcParams['text.usetex'] = True
#plt.rcParams.update({'font.size': 18})
plt.rcParams.update({'figure.facecolor':'white'})
plt.rcParams.update({'savefig.facecolor':'white'})
###########################################

import numpy as np

def plot_mean_velocities(cases=['WANG1_FR_lat60'], savefig=False):
#==========================================
    ## define Coriolis frequency
    # lat0={'WANG1_FR_lat30':30.,'WANG1_FR_lat60':60.,'WANG1_FR_Q2000_dT1e-3_lat60_120h':60.}
    # omega          = 7.292116e-05; rad            = np.pi/180.
    # fcor,ecor={},{}
    # for case in cases:
    #     fcor[case]      = 2.*omega*np.sin(rad*lat0[case])
    #     ecor[case]      = 2.*omega*np.cos(rad*lat0[case])
    g=9.809
    alpha=2.6e-4
    rho0=1027
    cp=4178
    stflx =111
    B0 = g*alpha/(rho0*cp)*stflx
    # B0=5e-8*10
    wstar = (B0*1000)**(1/3)

    #======= Open LES
    les,LG_MEAN,LG_RES,LG_SBG,BU_KE,TH_les,U_les,V_les,WTH,WU,WV,TKE,WTKE,dz_WU_les ,dz_WV_les ,time_les,hours,z_r_les,MLD = load_MNH(cases)
    #======= Open SCM + removing time=0
    out_dic={case: xr.open_dataset('scm_'+case+'.nc').isel(time=slice(1, None)) for case in cases}
    instant = -1 
    mld = MLD[cases[0]][instant]
    #========
    dz_WU_scm={}
    dz_WV_scm={}
    for case in cases:
        dz_WU_scm[case] = -(out_dic[case]['WU'].mean(dim='time').data[1:]-out_dic[case]['WU'].mean(dim='time').data[:-1])/(out_dic[case]['z_w'].data[1:]-out_dic[case]['z_w'].data[:-1])
        dz_WV_scm[case] = -(out_dic[case]['WV'].mean(dim='time').data[1:]-out_dic[case]['WV'].mean(dim='time').data[:-1])/(out_dic[case]['z_w'].data[1:]-out_dic[case]['z_w'].data[:-1])

    ################################ PLOTTING #################################
    #SCM
    styles = {case: '-' for case in cases}
    cocolors = plt.rcParams['axes.prop_cycle']
    colors = {case: cocolors.by_key()['color'][i] for i,case in enumerate(cases)}
    # colors = {cases[0]:'tab:gray' ,cases[1]:'tab:orange'}
    # colors = {cases[0]:'tab:gray' ,cases[1]:'tab:orange' , cases[2]:'tab:green' , cases[3]:'tab:blue'}

    # colors = { cases[0]:'tab:blue'}
    alpha = {case: 1 for case in cases}
    # linewidth = [2]*(len(run_label))

    #LES

    styles_les = {case: '-' for case in cases}
    colors_les = colors
    alpha_les  = {case:1 for case in cases}
    marker_les='o'
    linestyle_les=':'
    # width_les = {case:3 for case in cases}
    #============================================ WC ===============================================
    fig, axes = plt.subplots(nrows=1, ncols=2, sharex=False,
                            sharey=True, constrained_layout=True)

    ax_index=-1

    # ===============================================================
    ax_index+=1
    ax = axes.flat[ax_index]
    # ax.set_xlabel(r'$\mathrm{m s^{-1}}$')
    ax.set_ylabel(r'$z / h $')
    ax.set_title(r'Normalized U')

    for case in cases:
        
        ax.plot(out_dic[case].u[:].mean(dim='time')/wstar, out_dic[case].z_r/mld, linestyle=styles[case], color = colors[case], alpha=alpha[case], label=r'$\frac{\langle \overline{u} \rangle_{t}}{w_*}$, SCM '+case)
        # ax.plot(-1/fcor[case]*dz_WV_scm[case]/wstar, out_dic[case].z_w[1:]/mld, linestyle=':', color = colors[case], alpha=alpha[case], label=r'$-\frac{1}{f} \partial_z \frac{\langle \overline{w^\prime v^\prime} \rangle_t}{w_*}$, SCM '+case)
        
        ax.plot(U_les[case][:].mean(dim='time_les')/wstar, z_r_les[case]/mld, colors_les[case], 
            alpha=alpha_les[case],linestyle=linestyle_les, label=r'$\frac{ \langle \overline{u} \rangle_{t} }{w_*}$, LES '+case)
        # ax.plot(1/fcor[case]*dz_WV_les[case]/wstar, z_r_les[case][1:]/mld, colors_les[case], 
        #     alpha=alpha_les[case],linestyle=linestyle_les,linestyle=':',label=r'$-\frac{1}{f} \partial_z \frac{\langle \overline{w^\prime v^\prime} \rangle_t}{w_*}$, LES '+case)

        # var = -( 1/(fcor[case]*wstar)*out_dic[case].a_p*(out_dic[case].w_p)**2 ).diff(dim='z_w')[:].mean(dim='time')
        # ax.plot(var, out_dic[case].z_r/mld)
    # ax.set_xlim((-0.01, 0.06))
    ax.set_ylim((-1.3, 0))
    # handles, labels = ax.get_legend_handles_labels()
    # labels=[r'EDMF, lat=30',r'LES, lat=30',r'EDMF, lat=60',r'LES, lat=60']
    # ax.legend(handles, labels, fontsize=14)

    # ===============================================================
    ax_index+=1
    ax = axes.flat[ax_index]
    # ax.set_xlabel(r'$\mathrm{m s^{-1}}$')
    # ax.set_ylabel(r'$z / h $')
    ax.set_title(r'Normalized V')

    for case in cases:
        
        ax.plot(out_dic[case].v[:].mean(dim='time')/wstar, out_dic[case].z_r/mld, linestyle=styles[case], color = colors[case], alpha=alpha[case], label=r'$\frac{\langle \overline{v} \rangle_{t}}{w_*}$, SCM')

        # ax.plot(1/fcor[case]*dz_WU_scm[case]/wstar, out_dic[case].z_w[1:]/mld, linestyle=':', color = colors[case], alpha=alpha[case], label=r'$\frac{1}{f} \partial_z \frac{\langle \overline{w^\prime u^\prime} \rangle_t}{w_*} $, SCM')

        
        ax.plot(V_les[case][:].mean(dim='time_les')/wstar, z_r_les[case]/mld, colors_les[case], 
            alpha=alpha_les[case],linestyle=linestyle_les,label=r'$\frac{\langle \overline{v} \rangle_{t}}{w_*}$, LES')

        # ax.plot(-1/fcor[case]*dz_WU_les[case]/wstar, z_r_les[case][1:]/mld, colors_les[case], 
        #     alpha=alpha_les[case],linestyle=linestyle_les,linestyle=':',label=r'$\frac{1}{f} \partial_z \frac{ \langle \overline{w^\prime u^\prime} \rangle_t}{w_*} $, LES')
    # ax.set_xlim((2.93, 2.98))
    ax.set_ylim((-1.3, 0))
    handles, labels = ax.get_legend_handles_labels()
    labels=[r'lat=30, EDMF',r'lat=30, LES',r'lat=60, EDMF',r'lat=60, LES']
    ax.legend(handles, labels, fontsize=14)

    # ===============================================================

    # adding subplot labels
    subplot_label = [r'\rm{(a)}', r'\rm{(b)}', r'\rm{(c)}',
                    r'\rm{(d)}', r'\rm{(e)}', r'\rm{(f)}',r'\rm{(g)}',r'\rm{(h)}',r'\rm{(i)}',r'\rm{(j)}',r'\rm{(k)}',r'\rm{(l)}']

    for i,ax in enumerate(axes.flat):
        ax.text(0.1, 0.97, subplot_label[i], transform=ax.transAxes, bbox=dict(facecolor='1.', edgecolor='none'), fontweight='bold', va='top', ha='right')
        ax.set_box_aspect(1)
    if savefig:
        saving_path = '../figures/'
        saving_name = cases[0]+'_mean_velocities.png'    
        plt.savefig(saving_path+saving_name, bbox_inches='tight', dpi=300)
        print('figure saved at '+saving_path+saving_name)
    plt.show()


if __name__ == '__main__':
    cases=['WANG1_FR_lat30','WANG1_FR_lat60']
    plot_mean_velocities(cases=cases,savefig=True,)


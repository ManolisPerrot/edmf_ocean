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

plt.rcParams['text.usetex'] = True
#plt.rcParams.update({'font.size': 18})
plt.rcParams.update({'figure.facecolor':'white'})
plt.rcParams.update({'savefig.facecolor':'white'})
###########################################

def plot_instant_panel(avg_start=-2,avg_stop=-1,cases=['WANG1_FR_lat60'],savefig=False,saving_name = 'WANG1_NR_FR306090_profiles.png' ):

    #======= Open LES
    les,LG_MEAN,LG_RES,LG_SBG,BU_KE,TH_les,U_les,V_les,WTH,WU,WV,TKE,WTKE,dz_WU_les ,dz_WV_les ,time_les,hours,z_r_les,MLD = load_MNH(cases)
    #======= Open SCM + removing time=0
    out_dic={case: xr.open_dataset('scm_'+case+'.nc').isel(time=slice(1, None)) for case in cases}
    mld = MLD[cases[0]][avg_start]
    #========

    ################################# Styles #################################
    #SCM
    styles = {case: '-' for case in cases}
    # cocolors=['tab:gray' ,'tab:orange' , 'tab:green' , 'tab:blue']
    cocolors=['tab:blue','tab:orange' , 'tab:green' , 'tab:gray']
    colors = {case: cocolors[i] for i,case in enumerate(cases)}
    alpha  = {case: 1 for case in cases}
    # linewidth = [2]*(len(run_label))
    #LES
    styles_les = {case: 'k+' for case in cases}
    colors_les = colors
    # colors_les = {case: 'k' }
    alpha_les  = {case:1 for case in cases}
    marker_les=''
    width_les = {case:2 for case in cases}
    linestyle_les=':'
    ################################ Ploting ##################################
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
        # ax.plot(TH_les[case][avg_start], z_r_les[case]/mld,marker=marker_les, color=colors_les[case], 
        #     alpha=alpha_les[case],linewidth=width_les[case],linestyle=linestyle_les,label=r'LES')   
        # ax.plot(out_dic[case].temp[avg_start], out_dic[case].z_r/mld, linestyle=styles[case], color = colors[case], alpha=alpha[case], label=r'EDMF')

        ax.plot(TH_les[case][avg_start:avg_stop].mean(dim='time_les'), z_r_les[case]/mld,marker=marker_les, color=colors_les[case], 
            alpha=alpha_les[case],linewidth=width_les[case],linestyle=linestyle_les,label=r'LES')   
        ax.plot(out_dic[case].temp[avg_start:avg_stop].mean(dim='time'), out_dic[case].z_r/mld, linestyle=styles[case], color = colors[case], alpha=alpha[case], label=r'EDMF')


    handles, labels = ax.get_legend_handles_labels()
    labels=[r'LES',r'EDMF']
    ax.legend(handles[-2:], labels,loc='center left', fancybox=False, shadow=False,)    
    # ax.set_xlim((2.9820, 2.9885))
    ax.set_ylim((-1.4, 0))
    # ===============================================================
    ax_index+=1
    ax = axes.flat[ax_index]

    ax.set_title(r'$\overline{w^\prime \theta^\prime}$')

    for case in cases:
        # ax.plot(WTH[case][avg_start], z_r_les[case]/mld, marker=marker_les,color=colors_les[case], 
        # alpha=alpha_les[case],linewidth=width_les[case],linestyle=linestyle_les)
        # ax.plot(-out_dic[case]['WT'][avg_start], out_dic[case].z_w/mld, linestyle=styles[case], color = colors[case], alpha=alpha[case], label=case)

        ax.plot(WTH[case][avg_start:avg_stop].mean(dim='time_les'), z_r_les[case]/mld, marker=marker_les,color=colors_les[case], 
        alpha=alpha_les[case],linewidth=width_les[case],linestyle=linestyle_les)
        ax.plot(-out_dic[case]['WT'][avg_start:avg_stop].mean(dim='time'), out_dic[case].z_w/mld, linestyle=styles[case], color = colors[case], alpha=alpha[case], label=case)

    ax.set_xlabel(r'${\rm K}\;{\rm m}\;{\rm s}^{-1}$')

    # ===============================================================
    # ===============================================================
    ax_index+=1
    ax = axes.flat[ax_index]
    ax.set_title(r'$k$')

    for case in cases:
        ax.plot(TKE[case][avg_start:avg_stop].mean(dim='time_les'), z_r_les[case]/mld, marker=marker_les,color=colors_les[case], 
        alpha=alpha_les[case],linewidth=width_les[case],linestyle=linestyle_les)    
        ax.plot(out_dic[case]['tke'][avg_start:avg_stop].mean(dim='time'), out_dic[case].z_w/mld, linestyle=styles[case], color = colors[case], alpha=alpha[case], label=case)
    


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


    #     # ax.text(0.15, 0.95, subplot_label[i], transform=ax.transAxes, bbox=dict(facecolor='1.', edgecolor='none', pad=3.0), fontweight='bold', va='top', ha='right')

    handles, labels = axes.flat[0].get_legend_handles_labels()
    cases_name = {'WANG1_FR_lat30': r'lat=30°N','WANG1_FR_lat60':r'lat=60°N','WANG1_FR_lat90':r'lat=90°N','WANG1_NR_new':r'No rotation'}
    labels = [cases_name[case] for case in cases]
    fig.legend(handles[1::2], labels, loc='upper center', bbox_to_anchor=(
        0.5, 0.2), fancybox=False, shadow=False, ncol=4)

    # fig.tight_layout()

    if savefig:
        saving_path = '../figures/'   
        plt.savefig(saving_path+saving_name, bbox_inches='tight', dpi=300)
        print('figure saved at '+saving_path+saving_name)
    plt.show()

if __name__ == '__main__':
    cases = ['WANG1_FR_lat30','WANG1_FR_lat60','WANG1_FR_lat90','WANG1_NR_new']
    plot_instant_panel(cases=cases,savefig=True,avg_start=-10,avg_stop=-8)



















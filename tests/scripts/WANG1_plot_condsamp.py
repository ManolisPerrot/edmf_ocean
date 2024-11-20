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

def plot_condsamp_panel(instant=-1,cases=['WANG1_FR_lat60'],savefig=False,saving_name = 'WANG1_FR_condsamp.png',plot_vorticity=False ):

    #== Opening LES and SCM
    les,scm = {},{}
    for case in cases:
        file = case+'_object_diags_Cw_m1_276h.nc'
        path = '../data/'+case+'/'
        les[case] = xr.open_dataset(path+file)
        scm[case] = xr.open_dataset('scm_'+case+'.nc')
#----------------------------------------

    #======== Compute dz_(ap wp)
    MF_div = {}
    for case in cases:
        MF_div[case] = (les[case]['DW_FRAC']*les[case]['DW_WT']*les[case]['DW_THT']).diff('level')

    def E_D_LES(les):
    # ------------------------------------------------------------
    # Computing E,D on LES from passive tracer SVT001
    # ------------------------------------------------------------
    # $$ E\_ m\_ D = \partial_z (a_p w_p)$$
    # $$\tilde{E} = \frac{a_p w_p}{\phi_e - \phi_p} \partial_z \phi_p
    # $$
        dz = les['level'].data[1] - les['level'].data[0]
        regul = -0
        # interpolate on the level_w grid
        UP_SVT_interp = les['UP_SVT001'].interp(level=les.level_w).data
        DW_SVT_interp = les['DW_SVT001'].interp(level=les.level_w).data
        DW_FRAC_interp = les['DW_FRAC'].interp(level=les.level_w).data
        DW_WT_interp = les['DW_WT'].interp(level=les.level_w).data

        E_minus_D = 1/dz * (les['DW_FRAC'][:, 1:].data * les['DW_WT'][:, 1:].data -
                            les['DW_FRAC'][:, :-1].data * les['DW_WT'][:, :-1].data)

        Etilde = (DW_FRAC_interp[:, 1:-1] * DW_WT_interp[:, 1:-1]) / (UP_SVT_interp[:, 1:-1] -
                                                                    DW_SVT_interp[:, 1:-1]) * 1/dz * (les['DW_SVT001'][:, 1:].data - les['DW_SVT001'][:, :-1].data)

        # Etilde est le proxy pour E calculer avec traceur ou température
        # Astuce pour éviter d'avoir D < 0:
        # D = max (Etilde - (E-D), 0)  
        # Ce qui imlique que E = max(Etilde, E-D) (car quand E-D > Etilde, on impose D=0 donc ça implique que E = E-D)

        D = np.maximum(Etilde - E_minus_D, np.zeros_like(Etilde))
        E = np.maximum(Etilde, E_minus_D)

        return E,D

    E_les = { case:E_D_LES(les[case])[0] for case in cases }
    D_les = { case:E_D_LES(les[case])[1] for case in cases }

    zadim = les[cases[0]].level/(-les[cases[0]].MLD_FLUX[instant])
    mld = (-les[cases[0]].MLD_FLUX[instant])

    #==== Plotting Styles ====
    masks = ['TOT', 'DW']
    # color =     {'TOT': 'k','DW': 'k','UP': 'tab:orange'}
    cocolors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    cocolors[3] = 'tab:gray'
    colors = {case: cocolors[i] for i,case in enumerate(cases)}
    linewidth = {'TOT': 2,  'DW': 2, 'UP': 2}
    linestyle = {case: ':' for case in cases}
    linestyle_scm = {case: '-' for case in cases}

    alpha = {case: 1 for case in cases}
    alpha_les = {case: 1 for case in cases}
    #=========================

    subplot_label = [r'\rm{(a)}', r'\rm{(b)}', r'\rm{(c)}',
                    r'\rm{(d)}', r'\rm{(e)}', r'\rm{(f)}', r'\rm{(g)}',r'\rm{(h)}',r'\rm{(i)}']


    fig, axs = plt.subplots(nrows=2, ncols=3, sharey=True,constrained_layout=True)
    i_ax = -1
    #-----------------------------------------------------------------------
    i_ax+=1
    ax=axs.flat[i_ax]
    ax.set_title(r'$\theta_p$')
    ax.set_xlabel(r'$\mathrm{^oC}$')
    ax.set_ylabel(r'$z/h$')
    ax.set_xlim(2.9825,2.9880)
    for case in cases:
        mask='DW'
        ax.plot(les[case][mask+'_THT'][instant]-273.15 , zadim, color=colors[case], linewidth=linewidth[mask], linestyle=linestyle[case],alpha=alpha_les[case],label=case+'LES')
        #
        ax.plot(scm[case]['temp_p'][instant] , scm[case]['z_w']/mld, color=colors[case], linewidth=linewidth[mask], linestyle=linestyle_scm[case],alpha=alpha[case],label=case+'SCM')
    
    # handles, labels = ax.get_legend_handles_labels()
    # labels=[r'LES',r'EDMF']
    # ax.legend(handles[-2:], labels,loc='center left', fancybox=False, shadow=False,)  
    # handles, labels = ax.get_legend_handles_labels()
    # ax.legend(loc='lower left')

    # labels=[r'$\overline{\theta}$ (NR)',r'$\theta_p$ (NR)', r'$\overline{\theta}$ (FR)',r'$\theta_p$ (NR)']
    # ax.legend(handles,labels,fancybox=False)
    # #-----------------------------------------------------------------------
    i_ax+=1
    ax=axs.flat[i_ax]
    ax.set_title(r'$w_p$')
    ax.set_xlabel(r'$\mathrm{m s^{-1}}$')

    for case in cases:
        mask='DW'
        ax.plot(les[case][mask+'_WT'][instant] , zadim, color=colors[case], linewidth=linewidth[mask], linestyle=linestyle[case],alpha=alpha_les[case],label=case)
        ax.plot(scm[case]['w_p'][instant] , scm[case]['z_w']/mld, color=colors[case], linewidth=linewidth[mask], linestyle=linestyle_scm[case],alpha=alpha[case],label=case+mask)
    #-----------------------------------------------------------------------
    i_ax+=1
    ax=axs.flat[i_ax]
    ax.set_title(r'$a_p$')

    for case in cases:
        mask='DW'
        ax.plot(les[case][mask+'_FRAC'][instant] , zadim, color=colors[case], linewidth=linewidth[mask], linestyle=linestyle[case],alpha=alpha_les[case])
        ax.plot(scm[case]['a_p'][instant] , scm[case]['z_w']/mld, color=colors[case], linewidth=linewidth[mask], linestyle=linestyle_scm[case],alpha=alpha[case],label=case+mask)
    # # #-----------------------------------------------------------------------
    # i_ax+=1
    # ax=axs.flat[i_ax]
    # ax.set_title(r'$\partial_z(a_p w_p \theta_p)$')
    # ax.set_xlabel(r'$\mathrm{K s^{-1}}$')

    # for case in cases:
    #     mask='DW'
    #     ax.plot(MF_div[case][instant] , zadim[1:], color=colors[case], linewidth=linewidth[mask], linestyle=linestyle[case],alpha=alpha[case])

    # ax.set_xlim(-0.06,0.2)
    # # #-----------------------------------------------------------------------

    # #-----------------------------------------------------------------------
    if plot_vorticity:
        i_ax+=1
        ax=axs.flat[i_ax]
        ax.set_title(r'$\omega_p$')
        ax.set_xlabel(r'$\mathrm{s^{-1}}$')
        ax.set_ylabel(r'$z/h$')

        for case in cases:
            if case=='WANG1_NR_new':
                continue
            else:
                mask='DW'
                ax.plot(les[case][mask+'_VORT_z'][instant] , zadim, color=colors[case], linewidth=linewidth[mask], linestyle=linestyle[case],alpha=alpha_les[case])
                ax.plot(scm[case]['vort_p'][instant] , scm[case]['z_w']/mld, color=colors[case], linewidth=linewidth[mask], linestyle=linestyle_scm[case],alpha=alpha[case],label=case+mask)
    
    #--------------------------
    i_ax+=1
    while i_ax < len(axs.flat):
        axs.flat[i_ax].remove(); i_ax+=1
    # #-----------------------------------------------------------------------
    # # #-----------------------------------------------------------------------
    # i_ax+=1
    # ax=axs.flat[i_ax]
    # ax.set_title(r'$E-D$')
    # ax.set_xlabel(r'$\mathrm{s^{-1}}$')

    # for case in cases:
    #     mask='DW'
    #     ax.plot(E_les[case][instant]-D_les[case][instant] , zadim[1:], color=colors[case], linewidth=linewidth[mask], linestyle=linestyle[case],alpha=alpha_les[case])
    #     ax.plot(scm[case]['Ent'][instant]-scm[case]['Det'][instant] , scm[case]['z_r']/mld, color=colors[case], linewidth=linewidth[mask], linestyle=linestyle_scm[case],alpha=alpha[case],label=case+mask)
    #     # ax.plot(E_les[case][instant] , zadim[1:], color=colors[case], linewidth=linewidth[mask], linestyle=linestyle[case],alpha=alpha_les[case])
    #     # ax.plot(scm[case]['Ent'][instant] , scm[case]['z_r']/mld, color=colors[case], linewidth=linewidth[mask], linestyle=linestyle_scm[case],alpha=alpha[case],label=case+mask)
    # ax.set_xlim(-1e-5,1e-5)
    # # ax.set_xscale('log')
    # # #-----------------------------------------------------------------------
    # # #-----------------------------------------------------------------------
    # i_ax+=1
    # ax=axs.flat[i_ax]
    # ax.set_title(r'$E$')
    # ax.set_xlabel(r'$\mathrm{s^{-1}}$')

    # for case in cases:
    #     mask='DW'
    #     # ax.plot(E_les[case][instant]-D_les[case][instant] , zadim[1:], color=colors[case], linewidth=linewidth[mask], linestyle='-',alpha=alpha[case])
    #     # ax.plot(scm[case]['Ent'][instant]-scm[case]['Det'][instant] , scm[case]['z_r']/mld, color=colors[case], linewidth=linewidth[mask], linestyle=linestyle_scm[case],alpha=alpha[case],label=case+mask)
    #     ax.plot(E_les[case][instant] , zadim[1:], color=colors[case], linewidth=linewidth[mask], linestyle=linestyle[case],alpha=alpha_les[case])
    #     ax.plot(scm[case]['Ent'][instant] , scm[case]['z_r']/mld, color=colors[case], linewidth=linewidth[mask], linestyle=linestyle_scm[case],alpha=alpha[case],label=case+mask)
    # # ax.set_xlim(-1e-5,1e-4)
    # ax.set_xscale('log')
    # #-----------------------------------------------------------------------
    # # #-----------------------------------------------------------------------
    # i_ax+=1
    # ax=axs.flat[i_ax]
    # ax.set_title(r'$D$')
    # ax.set_xlabel(r'$\mathrm{s^{-1}}$')

    # for case in cases:
    #     mask='DW'
    #     ax.plot(D_les[case][instant] , zadim[1:], color=colors[case], linewidth=linewidth[mask], linestyle=linestyle[case],alpha=alpha_les[case])
    #     ax.plot(scm[case]['Det'][instant] , scm[case]['z_r']/mld, color=colors[case], linewidth=linewidth[mask], linestyle=linestyle_scm[case],alpha=alpha[case],label=case+mask)
    # ax.set_xscale('log')

    # # ax.set_xlim(-1e-5,1e-4)
    # ax.set_xscale('log')
    # i_ax+=1
    # ax=axs.flat[i_ax]
    # ax.remove()
    # # # #-----------------------------------------------------------------------
    # i_ax+=1
    # ax=axs.flat[i_ax]
    # ax.set_title(r'$D$')
    # ax.set_xlabel(r'$\mathrm{s^{-1}}$')

    # for case in cases:
    #     mask='DW'
    #     ax.plot(D_les[case][instant] , zadim[1:], color=colors[case], linewidth=linewidth[mask], linestyle='--',alpha=alpha[case])
    # ax.set_xlim(-3e-6,5e-5)

    subplot_label = [r'\rm{(a)}', r'\rm{(b)}', r'\rm{(c)}',
                        r'\rm{(d)}', r'\rm{(e)}', r'\rm{(f)}',r'\rm{(g)}',r'\rm{(h)}',r'\rm{(i)}',r'\rm{(j)}',r'\rm{(k)}',r'\rm{(l)}']

    for i,ax in enumerate(axs.flat):
        # ax.set_box_aspect(1)
        ax.grid(alpha=0.5)
        ax.set_ylim(-1.2,0)
        # adding subplot labels
        ax.text(0.15, 0.97, subplot_label[i], transform=ax.transAxes, bbox=dict(facecolor='1.', edgecolor='none',alpha=0.5), fontweight='bold', va='top', ha='right')


    # handles, labels = axs.flat[0].get_legend_handles_labels()
    # # labels = [r'NR, LES',r'NR, SCM',r'FR90, LES',r'FR90, SCM']

    # # fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(1.01, 0.8),
    # #            fancybox=False, shadow=False, ncols=2)
    # axs.flat[0].legend(handles, labels,fancybox=False)


    handles, labels = axs.flat[0].get_legend_handles_labels()
    cases_name = {'WANG1_FR_lat30': r'lat=30°N','WANG1_FR_lat60':r'lat=60°N','WANG1_FR_lat90':r'lat=90°N','WANG1_NR_new':r'No rotation'}
    labels = [cases_name[case] for case in cases]

    fig.legend(handles[1::2], labels, 
               bbox_to_anchor=(0.95, 0.3), 
               fancybox=False, shadow=False, ncol=2)
    fig.suptitle(r'Plume variables')

    if savefig:
        saving_path = '../figures/'   
        plt.savefig(saving_path+saving_name, bbox_inches='tight', dpi=300)
        print('figure saved at '+saving_path+saving_name)
    plt.show()

if __name__ == '__main__':
    cases = ['WANG1_FR_lat30','WANG1_FR_lat60','WANG1_FR_lat90','WANG1_NR_new']
    plot_condsamp_panel(cases=cases,plot_vorticity=True,savefig=True)



















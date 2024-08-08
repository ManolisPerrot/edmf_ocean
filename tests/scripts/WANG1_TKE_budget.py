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


#cases = ['WANG1_NR','WANG1_FR','WANG1_lat_90','W005_C500_NO_COR','WANG1_FR_domain_div_2','WANG1_FR_10f2']
cases = ['WANG1_FR']

les = {}
LG_MEAN={}
LG_RES={}
LG_SBG={}
BU_KE={}

for case in cases:
    file = 'GN_01.1.OC_01.000.nc'
    path = '../data/'+case+'/'
    print('opening '+path+file)
    les[case] = xr.open_dataset(path+file)
    LG_MEAN[case]= xr.open_dataset(path+file,group ='/LES_budgets/Mean/Cartesian/Not_time_averaged/Not_normalized/cart')
    LG_RES[case] = xr.open_dataset(path+file,group ='/LES_budgets/Resolved/Cartesian/Not_time_averaged/Not_normalized/cart')    
    LG_SBG[case] = xr.open_dataset(path+file,group ='/LES_budgets/Subgrid/Cartesian/Not_time_averaged/Not_normalized/cart')
    BU_KE[case] = xr.open_dataset(path+file,group ='/LES_budgets/BU_KE/Cartesian/Not_time_averaged/Not_normalized/cart')

instant=-1
Z = les[case]['level_les'] - les[case]['level_w'][-1]
Z_w_les = 0.5*(Z.data[1:]+Z.data[:-1])

mld = Z[(-LG_MEAN[case]['MEAN_TH'][instant]).argmax()]

WU_les = LG_RES[case]['RES_WU'] + LG_SBG[case]['SBG_WU']
WV_les = LG_RES[case]['RES_WV'] + LG_SBG[case]['SBG_WV']    
dz_WU_les = (WU_les.mean(dim='time_les').data[1:]-WU_les.mean(dim='time_les').data[:-1])/(Z.data[1:]-Z.data[:-1])
dz_WV_les = (WV_les.mean(dim='time_les').data[1:]-WV_les.mean(dim='time_les').data[:-1])/(Z.data[1:]-Z.data[:-1])

zlim = mld - 50
zadim = Z/(-mld)

#open scm
run_params ={'output_filename': 'scm_WANG1_FR.nc'}
print('opening '+run_params['output_filename'])
scm = xr.open_dataset(run_params['output_filename'])

dz_WU_scm = (scm['WU'].mean(dim='time').data[1:]-scm['WU'].mean(dim='time').data[:-1])/(scm['z_w'].data[1:]-scm['z_w'].data[:-1])
dz_WV_scm = (scm['WV'].mean(dim='time').data[1:]-scm['WV'].mean(dim='time').data[:-1])/(scm['z_w'].data[1:]-scm['z_w'].data[:-1])








    ############ Plotting
def plot_waves_panel():
    vmin=-0.004;vmax=0.004
    levels = MaxNLocator(nbins=15).tick_values(vmin,vmax)
    tt = (les[case]['time_les'] - les[case]['time_les'][0])/np.timedelta64(1,'h')+2 #array of hours


    time_les_mesh, z_les_mesh = np.meshgrid(tt, Z)
    time_scm_r_mesh, z_scm_r_mesh= np.meshgrid(scm['ocean_time']/np.timedelta64(1,'h'), scm['z_r'])
    time_scm_w_mesh, z_scm_w_mesh= np.meshgrid(scm['ocean_time']/np.timedelta64(1,'h'), scm['z_w'])

    fig, axs = plt.subplots(nrows=4, ncols=3, sharey=True, layout='constrained', figsize=(10,13))
    i=-1
    #-----------------------------------------------------
    i+=1; ax = axs.flat[i]
    ax.set_title('LES: Total u')
    ax.set_xlabel('time (h)')
    ax.set_ylabel(r'm')
    var=LG_MEAN[case]['MEAN_U']
    ax.contourf( time_les_mesh, z_les_mesh, var.T,vmin=vmin,vmax=vmax,levels=levels,cmap=plt.cm.RdBu_r)

    #-----------------------------------------------------
    i+=1; ax = axs.flat[i]
    ax.set_title('LES: Waves u')
    ax.set_xlabel('time (h)')
    var=LG_MEAN[case]['MEAN_U'] - LG_MEAN[case]['MEAN_U'].mean(dim='time_les')
    ax.contourf( time_les_mesh, z_les_mesh, var.T,vmin=vmin,vmax=vmax,levels=levels,cmap=plt.cm.RdBu_r)

    #-----------------------------------------------------
    i+=1; ax = axs.flat[i]
    ax.set_title('LES: Mean flow u')
    ax.set_xlabel(r'$m s^{-1}$')
    var=LG_MEAN[case]['MEAN_U'].mean(dim='time_les')
    ax.plot( var.T, Z,  label=r'$\overline{u}^{x,y,t}$')
    ax.plot( -1/(scm.attrs['fcorSIN'])*dz_WV_les , Z_w_les,  color='black', linestyle=':', label=r'$-\frac{1}{f} \partial_z \overline{w^\prime v^\prime}^{x,y,t}$')
    ax.set_xlim(-0.004,0.004)
    ax.legend()


    #-----------------------------------------------------
    i+=1; ax = axs.flat[i]
    ax.set_title('SCM: Total u')
    ax.set_xlabel('time (h)')
    ax.set_ylabel(r'm')
    var=scm['u']
    ax.contourf( time_scm_r_mesh, z_scm_r_mesh, var.T,vmin=vmin,vmax=vmax,levels=levels,cmap=plt.cm.RdBu_r)

    #-----------------------------------------------------
    i+=1; ax = axs.flat[i]
    ax.set_title('SCM: Waves u')
    ax.set_xlabel('time (h)')
    var=scm['u'] - scm['u'].mean(dim='time')
    ax.contourf( time_scm_r_mesh, z_scm_r_mesh, var.T,vmin=vmin,vmax=vmax,levels=levels,cmap=plt.cm.RdBu_r)

    #-----------------------------------------------------
    i+=1; ax = axs.flat[i]
    ax.set_title('SCM: Mean flow u')
    ax.set_xlabel(r'$m s^{-1}$')
    ax.plot( scm['u'][1:].mean(dim='time'), scm['z_r'], label=r'$\overline{u}^{x,y,t}$')
    ax.plot( 1/(scm.attrs['fcorSIN'])*dz_WV_scm , scm['z_r'],  color='black', linestyle=':', label=r'$\frac{1}{f} \partial_z \overline{w^\prime v^\prime}^{x,y,t}$')
    ax.set_xlim(-0.004,0.004)
    ax.legend()


    #-----------------------------------------------------
    i+=1; ax = axs.flat[i]
    ax.set_title('LES: Total v')
    ax.set_xlabel('time (h)')
    ax.set_ylabel(r'm')
    var=LG_MEAN[case]['MEAN_V']
    ax.contourf( time_les_mesh, z_les_mesh, var.T,vmin=vmin,vmax=vmax,levels=levels,cmap=plt.cm.RdBu_r)

    #-----------------------------------------------------
    i+=1; ax = axs.flat[i]
    ax.set_title('LES: Waves v')
    ax.set_xlabel('time (h)')
    var=LG_MEAN[case]['MEAN_V'] - LG_MEAN[case]['MEAN_V'].mean(dim='time_les')
    im=ax.contourf( time_les_mesh, z_les_mesh, var.T,vmin=vmin,vmax=vmax,levels=levels,cmap=plt.cm.RdBu_r)

    #-----------------------------------------------------
    i+=1; ax = axs.flat[i]
    ax.set_title('LES: Mean flow v')
    ax.set_xlabel(r'$m s^{-1}$')
    var=LG_MEAN[case]['MEAN_V'].mean(dim='time_les')
    ax.plot( var.T, Z,  label=r'$\overline{v}^{x,y,t}$')
    ax.plot( 1/(scm.attrs['fcorSIN'])*dz_WU_les , Z_w_les,  color='black', linestyle=':', label=r'$\frac{1}{f} \partial_z \overline{w^\prime u^\prime}^{x,y,t}$')
    ax.set_xlim(-0.004,0.004)
    ax.legend()


    #-----------------------------------------------------
    i+=1; ax = axs.flat[i]
    ax.set_title('SCM: Total v')
    ax.set_xlabel('time (h)')
    ax.set_ylabel(r'm')
    var=scm['v']
    ax.contourf( time_scm_r_mesh, z_scm_r_mesh, var.T,vmin=vmin,vmax=vmax,levels=levels,cmap=plt.cm.RdBu_r)

    #-----------------------------------------------------
    i+=1; ax = axs.flat[i]
    ax.set_title('SCM: Waves v')
    ax.set_xlabel('time (h)')
    var=scm['v'] - scm['v'].mean(dim='time')
    ax.contourf( time_scm_r_mesh, z_scm_r_mesh, var.T,vmin=vmin,vmax=vmax,levels=levels,cmap=plt.cm.RdBu_r)

    #-----------------------------------------------------
    i+=1; ax = axs.flat[i]
    ax.set_title('SCM: Mean flow v')
    ax.set_xlabel(r'$m s^{-1}$')
    ax.plot( scm['v'][1:].mean(dim='time'), scm['z_r'], label=r'$\overline{v}^{x,y,t}$')
    ax.plot( -1/(scm.attrs['fcorSIN'])*dz_WU_scm , scm['z_r'],  color='black', linestyle=':', label=r'$\frac{1}{f} \partial_z \overline{w^\prime u^\prime}^{x,y,t}$')
    ax.set_xlim(-0.004,0.004)
    ax.legend()

    #-----------------------------------------------------
    for ax in axs.flat:
        ax.set_ylim(-3300,0)
        ax.set_box_aspect(1)

    #fig.colorbar(im, orientation="vertical" , ticks=np.linspace(vmin,vmax,9))
    plt.suptitle(case)

    saving_path = '../figures/WANG1_LES_vs_SCM_velocities.png'
    plt.savefig(saving_path, bbox_inches='tight', dpi=300)
    print('figure saved at '+saving_path)

    #plt.show()
plot_waves_panel()



def plot_mean_fluxes():
    WU_les_mean=WU_les.mean(dim='time_les')
    WV_les_mean=WV_les.mean(dim='time_les')

    WU_scm_mean=-scm['WU'].mean(dim='time')
    WV_scm_mean=-scm['WV'].mean(dim='time')



    fig, axs = plt.subplots(nrows=1, ncols=2, sharey=True, layout='constrained')
    i=-1
    #-----------------------------------------------------
    i+=1; ax = axs.flat[i]
    ax.set_title(r'$\overline{w^\prime u^\prime}$')
    ax.set_ylabel(r'$m$')
    ax.plot(WU_les_mean, Z,label='LES')
    ax.plot(WU_scm_mean, scm['z_w'],label='SCM')
    ax.axvline(0.,color='k',alpha=1,linestyle='--',linewidth=1)
    ax.legend()
    #-----------------------------------------------------
    i+=1; ax = axs.flat[i]
    ax.set_title(r'$\overline{w^\prime v^\prime}$')
    ax.set_ylabel(r'$m$')
    ax.plot(WV_les_mean, Z,label='LES')
    ax.plot(WV_scm_mean, scm['z_w'],label='SCM')
    ax.axvline(0.,color='k',alpha=1,linestyle='--',linewidth=1)

    ax.legend()

    saving_path = '../figures/WANG1_LES_vs_SCM_momentum_flux_averaged.png'
    plt.savefig(saving_path, bbox_inches='tight', dpi=300)
    print('figure saved at '+saving_path)

#plot_mean_fluxes()
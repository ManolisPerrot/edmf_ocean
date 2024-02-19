#!/usr/bin/env python
# coding: utf-8
__docformat__ = 'reStructuredText'
###########################################
# Imports
###########################################
from sys import exit
import sys
sys.path.insert(0,"../library/F2PY/")
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
#matplotlib.use('pdf')
#matplotlib.use('TkAgg')
from netCDF4 import Dataset
from scm_class_asicsmed import SCM
from scmoce import scm_oce
from scipy.ndimage import gaussian_filter
###########################################
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
###########################################
def single_column(config):
    """
      Run the single column model.
    """
    #=============================== TKE only
    nc_file2  = 'run1'
    nc_file3  = 'run2'
    data_file = '../data/asicsmed/LION_data_ASICSMED.nc'
    #=============================== TKE + EVD
    scm2 = SCM( output_filename=nc_file2,evd=True  ,
                   mass_flux_tra          = False  ,
                   mass_flux_dyn          = False  ,
                   mass_flux_tke          = False  ,
                   mass_flux_tke_trplCorr = False  ,
                   mass_flux_small_ap     = False  ,
                   eddy_diff_tke_const    = 'NEMO' ,
                   tke_sfc_dirichlet      = False,
                   Cent                   = 0.7,
                   Cdet                   = 1.99,
                   wp_bp                  = 0.00005,
                   wp_b                   = 1.0,
                   wp_a                   = 1.0,
                   bc_ap                  = 0.2,
                   delta_bkg              = 0.0001,
                   entr_scheme = 'R10' )
    scm2.run_direct()
    #=============================== EDMF tra+dyn+tke
    scm3 = SCM( output_filename=nc_file3,evd=False,
                   mass_flux_tra          = True   ,
                   mass_flux_dyn          = True  ,
                   mass_flux_tke          = True  ,
                   mass_flux_tke_trplCorr = True  ,
                   mass_flux_small_ap     = False  ,
                   eddy_diff_tke_const    = 'NEMO' ,
                   tke_sfc_dirichlet      = False,
                   Cent                   = 0.7,
                   Cdet                   = 1.99,
                   wp_bp                  = 0.00005,
                   wp_b                   = 1.0,
                   wp_a                   = 1.0,
                   bc_ap                  = 0.2,
                   delta_bkg              = 0.0001,
                   entr_scheme = 'R10' )
    scm3.run_direct()
    #===============================
    nd = 4
    tdays = ["25-jan-2013","29-jan-2013","04-feb-2013","09-feb-2013"]
    tscm  = np.array([10,14,20,25.5]);tscm  = tscm*24
    tdtm  = np.array([230,322,460,620])
    tmin  = np.zeros(nd); tmin[:]=12.875
    tmax  = np.zeros(nd); tmax[:]=13.25
    Dmax  = np.array([1000,1500,2000,2300])
    #===============================
    #### create figure 1
    fig, axes = plt.subplots(nrows = 2, ncols = 4, sharex=False, sharey=False)
    fig.set_figheight(8)
    fig.set_figwidth(12)
    #===============================
    for x in range(nd):
      ntime     = tscm[x]
      ntime_dta = tdtm[x]
      hmax      = Dmax[x]
      #===
      fh        = Dataset(nc_file2, mode='r')
      depth2    = fh.variables['z_r'][:]
      temp2     = fh.variables['temp'][ntime,:]
      salt2     = fh.variables['salt'][ntime,:]
      fh.close()
      #===
      fh        = Dataset(nc_file3, mode='r')
      depth3    = fh.variables['z_r'][:]
      temp3     = fh.variables['temp'][ntime,:]
      salt3     = fh.variables['salt'][ntime,:]
      fh.close()
      #===
      nc_dta    = Dataset(data_file, mode='r')
      print(nc_dta.variables['ocean_time'][ntime_dta]/(3600*24))
      dta_depth = nc_dta.variables['z_r'][:]
      dta_temp  = nc_dta.variables['temp'][ntime_dta,:]
      dta_salt  = nc_dta.variables['salt'][ntime_dta,:]
      nc_dta.close()
      # remove anomalous values
      maskT = (dta_temp[:] > 100.)
      dta_temp[:][maskT] = np.nan
      dta_temp[18] = np.nan
      new_temp   = dta_temp[~np.isnan(dta_temp)]
      new_depthT = dta_depth[~np.isnan(dta_temp)]
      #
      maskS = (dta_salt[:] > 100.)
      dta_salt[:][maskS] = np.nan
      new_salt   = dta_salt[~np.isnan(dta_salt)]
      new_depthS = dta_depth[~np.isnan(dta_salt)]
      #
      ax = axes[0,x]
      ax.set_ylim(-hmax,-150.)
      if x==0:
        ax.plot(temp2, depth2,color='0.5',label='ED+EVD',linewidth=3)
        ax.plot(temp3, depth3,color='k',label='EDMF-Energy',linewidth=3)
        ax.plot(new_temp, new_depthT,marker='o',label='mooring data',linewidth=1.5,linestyle="dashed",color='r')
      else:
        ax.plot(temp2, depth2,color='0.5',linewidth=3)
        ax.plot(temp3, depth3,color='k',linewidth=3)
        ax.plot(new_temp, new_depthT,marker='o',linewidth=1.5,linestyle="dashed",color='r')
      ax.set_title('Temperature on '+tdays[x])
      if x==0: ax.set_ylabel(r'${\rm depth}\;[{\rm m}]$')
      ax.set_xlim(tmin[x],tmax[x])

      ax = axes[1,x]
      ax.set_ylim(-hmax,-150.)
      ax.plot(salt2, depth2,color='0.5',linewidth=3)
      ax.plot(salt3, depth3,color='k',linewidth=3)
      ax.plot(new_salt, new_depthS, marker='o',linewidth=1.5,linestyle="dashed",color='r')
      ax.set_title('Salinity on '+tdays[x])
      if x==0: ax.set_ylabel(r'${\rm depth}\;[{\rm m}]$')
      ax.set_xlim(38.45,38.55)
    #
    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    fig.legend(lines, labels, loc='lower center', bbox_to_anchor=(0.5, -0.04), ncol=4, fontsize="14")
    #===============================
    plt.tight_layout()
    plt.savefig('../figures/asicsmed_results.pdf', bbox_inches='tight')
    #===============================
    #### create figure 2
    #======================================================
    # compute B0 and ustar
    #======================================================
    ndays     = 30; time1     = 0; time2     = 613
    fh        = Dataset(nc_file2, mode='r')
    ustar     = np.sqrt(np.sqrt( np.power((1027./1.22)*fh.variables['taux'][time1:time2],2)
                               + np.power((1027./1.22)*fh.variables['tauy'][time1:time2],2) ))
    time      = fh.variables['ocean_time'][time1:time2] / (3600*24)
    Qsol      = fh.variables['Qs'][time1:time2]/(1027.*3985.)    #K m/s
    Qnsol     = fh.variables['Qns'][time1:time2]/(1027.*3985.)   #K m/s
    Fw        = fh.variables['Fw'][time1:time2]    #psu m / s
    Tt        = fh.variables['temp'][time1:time2,-1] # sst
    Ts        = fh.variables['salt'][time1:time2,-1] # sss
    hmxl1     = fh.variables['hmxl_300m'][time1:time2]
    fh.close()
    #
    fh        = Dataset(nc_file3, mode='r')
    hmxl2     = fh.variables['hmxl_300m'][time1:time2]
    fh.close()
    #=================================================================
    # Estimate surface values of alpha and beta using nonlinear EOS
    #=================================================================
    Q00=+999.842594 ; Q01=+6.793952e-2; Q02=-9.095290e-3
    Q03=+1.001685e-4; Q04=-1.120083e-6; Q05=+6.536332e-9
    U00=+0.824493   ; U01=-4.08990e-3 ; U02=+7.64380e-5
    U03=-8.24670e-7 ; U04=+5.38750e-9 ; V00=-5.72466e-3
    V01=+1.02270e-4 ; V02=-1.65460e-6 ; W00=+4.8314e-4
    sqrtTs=np.sqrt(Ts)
    cff   = 1./1027.
    alpha = -cff * ( Q01+Tt*(2.*Q02+Tt*(3.*Q03+Tt*(4.*Q04+Tt*5.*Q05)))
            +Ts*(U01+Tt*(2.*U02+Tt*(3.*U03+Tt*4.*U04))+sqrtTs*(V01+Tt*2.*V02)))
    beta  = cff*( U00+Tt*(U01+Tt*(U02+Tt*(U03+Tt*U04)))+1.5*(V00+Tt*(V01+Tt*V02))*sqrtTs+2.*W00*Ts)
    Buoy_sfc = 9.81*(alpha*(Qsol+Qnsol)-beta*Fw)
    #==================================================================
    tscm  = np.array([20,28,40,51]);tscm  = tscm*12  # draw vertical blue lines on the plot
    #===============================================
    # Computation of hmxl from observations
    #===============================================
    indx   = np.array([25,24,23,22,21,20,19,17,16,15,14,13,11,10,9,8,7,5,4,3,2])
    nc_dta    = Dataset(data_file, mode='r')
    dta_depth =     nc_dta.variables['z_r'][:]
    time_dta  = nc_dta.variables['ocean_time'][:]/(3600*24)
    ntime     = len(time_dta)
    nc_dta.close()
    nz     = len(indx)
    z_r    = np.zeros(nz)
    z_r[:] = dta_depth[indx[:]]
    z_w    = np.zeros(nz+1)
    z_w[1:nz] = 0.5*(z_r[1:nz] + z_r[0:nz-1]);z_w[0] = z_r[0] - (z_w[1]-z_r[0]);z_w[nz] = z_r[nz-1] + (z_r[nz-1]-z_w[nz-1])
    #
    hmxl_dta  = np.zeros(ntime)
    temp      = np.zeros(nz); salt = np.zeros(nz); salt[:]   = 38.5
    #======================================================
    for tt in range(ntime):
      nc_dta    = Dataset(data_file, mode='r')
      dta_temp  = nc_dta.variables['temp'][tt,:]
      #dta_salt  = nc_dta.variables['salt'][tt,:]
      nc_dta.close()
      temp[:] = dta_temp[indx[:]]    #; salt[:] = dta_salt[indx[:]]
      rho,bvf = scm_oce.rho_eos(temp,salt,z_r,z_w,nz)
      drho    = 0.01; zref = z_r[-5]
      hmxl_dta[tt] = scm_oce.compute_mxl2(bvf,0.01,z_r,nz)
    #===============================
    #### create figure
    f, axes = plt.subplots(nrows = 3, ncols = 1, sharex=True, sharey=False)
    f.set_figheight(5)
    f.set_figwidth(5)
    #
    ax = axes[0]
    ax.set_title(r"$u_\star\;[{\rm m}\;{\rm s}^{-1}]$",fontsize=12)
    ax.axvline(x = time[tscm[0]], color = 'b', linestyle = '-', linewidth=1.0)
    ax.axvline(x = time[tscm[1]], color = 'b', linestyle = '-', linewidth=1.0)
    ax.axvline(x = time[tscm[2]], color = 'b', linestyle = '-', linewidth=1.0)
    ax.axvline(x = time[tscm[3]], color = 'b', linestyle = '-', linewidth=1.0)
    ax.plot(time, ustar,color='k',linewidth=1.5)
    ax.set_xlim(7.5,25.75)
    #
    ax = axes[1]
    ax.set_title(r"$B_0 \; [{\rm m}^2\;{\rm s}^{-3}]$",fontsize=12)
    ax.plot(time, Buoy_sfc,color='k',linewidth=1.5)
    ax.axhline(y = 0., color = 'r', linestyle = '-', linewidth=1.0)
    ax.axvline(x = time[tscm[0]], color = 'b', linestyle = '-', linewidth=1.0)
    ax.axvline(x = time[tscm[1]], color = 'b', linestyle = '-', linewidth=1.0)
    ax.axvline(x = time[tscm[2]], color = 'b', linestyle = '-', linewidth=1.0)
    ax.axvline(x = time[tscm[3]], color = 'b', linestyle = '-', linewidth=1.0)
    #
    ax = axes[2]
    ax.set_title(r"$h_{\rm mxl} \; [{\rm m}]$",fontsize=12)
    ax.plot(time, hmxl1,color='0.5',linewidth=1.5,label='ED+EVD')
    ax.plot(time, hmxl2,color='k',linewidth=1.5,label='EDMF-Energy')
    ax.plot(time_dta, gaussian_filter(hmxl_dta,sigma=3),color='r',linewidth=1.5,label='mooring data')
    ax.axvline(x = time[tscm[0]], color = 'b', linestyle = '-', linewidth=1.0)
    ax.axvline(x = time[tscm[1]], color = 'b', linestyle = '-', linewidth=1.0)
    ax.axvline(x = time[tscm[2]], color = 'b', linestyle = '-', linewidth=1.0)
    ax.axvline(x = time[tscm[3]], color = 'b', linestyle = '-', linewidth=1.0)
    ax.set_xticks([7,9,11,13,15,17,19,21,23,25])
    ax.set_xticklabels([r'$22\;{\rm Jan}$',r'$24\;{\rm Jan}$',r'$26\;{\rm Jan}$',r'$28\;{\rm Jan}$',r'$30\;{\rm Jan}$',r'$01\;{\rm Feb}$',r'$03\;{\rm Feb}$',r'$05\;{\rm Feb}$',r'$07\;{\rm Feb}$',r'$09\;{\rm Feb}$'])
    ax.legend(loc=3,fontsize='6')
    ax.set_ylim(-2300.,0.)
    plt.tight_layout()
    plt.savefig('../figures/asicsmed_forcing_mxl.pdf', bbox_inches='tight')
###########################################
# Main
###########################################
if __name__ == "__main__":
    config = "AsicsMed"
    single_column(config = config)

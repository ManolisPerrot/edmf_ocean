#!/usr/bin/env python
# coding: utf-8
__docformat__ = 'reStructuredText'
###########################################
# Imports
###########################################
from sys import exit
import sys
sys.path.insert(0,"../../library/F2PY/")
#import add
import numpy as np
import matplotlib
matplotlib.rcParams.update(matplotlib.rcParamsDefault)
#matplotlib.use('TkAgg')
#matplotlib.use('pdf')
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scm_class_asicsmed import SCM
from scmoce import scm_oce
#from F2PY import scm_class
import scipy.signal
from scipy.ndimage import gaussian_filter
###########################################
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
###########################################
case        = 'asicsmed'
#===========================================================================
common_params = {
    'nz': 75,
    'dt': 1200.,
    'initial_filename': '../data/asicsmed/init_ASICS_m01d15.nc',
    'sfcforc_filename': '../data/asicsmed/forc_ASICS_y2013.nc',
    'nbhours': 720,
    'outfreq': 1,
    'output_filename': 'run',
    'lat0': 42.04,
    'eddy_diff': True,
    'evd': False,
    'mass_flux_tra': False,
    'mass_flux_dyn': False,
    'mass_flux_tke': False,
    'mass_flux_tke_trplCorr': False,
    'mass_flux_small_ap': False,
    'lin_eos': False,
    'extrap_ak_surf': True,
    'tke_sfc_dirichlet': False,
    'eddy_diff_scheme': 'TKE',
    'eddy_diff_tke_const': 'NEMO',
    'entr_scheme': 'R10',
    'Cent': 0.7,
    'Cdet': 1.99,
    'wp_a': 1.0,
    'wp_b': 1.0,
    'wp_bp': 0.125,
    'up_c': 0.25,
    'vp_c': 0.25,
    'bc_ap': 0.20,
    'delta_bkg': 0.0625,
    'wp0': -1.0e-08,
    'akvmin': 1.e-5,
    'aktmin': 1.e-6,
    'mxlmin': 1.0,
}
#===========================================================================
# Define parameters specific to each run (overwrite common parameters):
run_label = ['ED+EVD', 'EDMF-Energy']
runs = [
    {
        'eddy_diff': True,
        'evd': True,
        'eddy_diff_scheme' : 'TKE',
        'output_filename': 'run_asics1.nc'
    },
        {
        'eddy_diff': True,
        'evd': False,
        'eddy_diff_scheme' : 'TKE',
        'mass_flux_tra': True,
        'mass_flux_dyn': True,
        'mass_flux_tke': True,
        'mass_flux_tke_trplCorr': True,
        'output_filename': 'run_asics3.nc'
    },
        ]
#
colors    = ['tab:gray','tab:green']
linewidth = [3]*(len(run_label))
linestyle = ['solid','solid']
#
scm       = [0]*len(runs)
# Run the SCM
for i, run_params in enumerate(runs):
    params = common_params.copy()  # Create a copy of common_params
    params.update(run_params)  # Update with run_params
    scm[i] = SCM(**params)
    scm[i].run_direct()
    print('results stored in ',scm[i].output)
#===========================================================================
data_file = '../data/asicsmed/LION_data_ASICSMED.nc'
#===========================================================================
nd = 3
tdays = ["25-jan-2013","04-feb-2013","09-feb-2013"]
tscm  = np.array([10,20,25.5]);tscm  = tscm*24
tdtm  = np.array([230,460,620])
tmin  = np.zeros(nd); tmin[:]=12.875
tmax  = np.zeros(nd); tmax[:]=13.25
Dmax  = np.array([1000,2000,2300])
# Create figure 1
fig, axes = plt.subplots(nrows = 3, ncols = 2, sharex=False, sharey=True)
fig.set_figheight(9)
fig.set_figwidth(6)
#===============================
out = [0]*len(runs)

for i in range(len(runs)):
  out[i] = Dataset(scm[i].output, mode='r')
#
for x in range(nd):
  ntime     = tscm[x]
  ntime_dta = tdtm[x]
  hmax      = 2500. #Dmax[x]
  #===
  nc_dta    = Dataset(data_file, mode='r')
  print(nc_dta.variables['ocean_time'][ntime_dta]/(3600*24))
  dta_depth = nc_dta.variables['z_r'][:]
  dta_temp  = nc_dta.variables['temp'][ntime_dta,:]
  dta_salt  = nc_dta.variables['salt'][ntime_dta,:]
  nc_dta.close()
  # remove anomalous values
  maskT = (dta_temp[:] > 100.); dta_temp[:][maskT] = np.nan; dta_temp[18] = np.nan
  new_temp   = dta_temp[~np.isnan(dta_temp)]
  new_depthT = dta_depth[~np.isnan(dta_temp)]
  #
  maskS = (dta_salt[:] > 100.); dta_salt[:][maskS] = np.nan
  new_salt   = dta_salt[~np.isnan(dta_salt)]
  new_depthS = dta_depth[~np.isnan(dta_salt)]
  #=============
  ax = axes[x,0]
  ax.set_ylim(-hmax,-150.)
  #
  if x==0:
    for i in range(len(runs)):
      ax.plot(out[i].variables['temp'][ntime,:], out[i].variables['z_r'][:],
                color = colors[i],linewidth=linewidth[i],linestyle=linestyle[i],label=run_label[i])
    ax.plot(new_temp, new_depthT,marker='o',label='mooring data',linewidth=1.5,linestyle="dashed",color='k')
  else:
    for i in range(len(runs)):
      ax.plot(out[i].variables['temp'][ntime,:], out[i].variables['z_r'][:],
                color = colors[i],linewidth=linewidth[i],linestyle=linestyle[i])
    ax.plot(new_temp, new_depthT,marker='o',linewidth=1.5,linestyle="dashed",color='k')
  #
  for tick in ax.xaxis.get_major_ticks():
                tick.label1.set_fontsize(12)
  for tick in ax.yaxis.get_major_ticks():
                tick.label1.set_fontsize(12)
  #
  ax.set_title('Temperature on '+tdays[x])
  ax.set_xlim(12.875,13.205)
  #=============
  ax = axes[x,1]
  ax.set_ylim(-hmax,-150.)
  #
  for i in range(len(runs)):
    ax.plot(out[i].variables['salt'][ntime,:], out[i].variables['z_r'][:],
                color = colors[i],linewidth=linewidth[i],linestyle=linestyle[i])
  ax.plot(new_salt, new_depthS,marker='o',linewidth=1.5,linestyle="dashed",color='k')
  #
  for tick in ax.xaxis.get_major_ticks():
                tick.label1.set_fontsize(12)
  for tick in ax.yaxis.get_major_ticks():
                tick.label1.set_fontsize(12)
  #
  ax.set_title('Salinity on '+tdays[x])
  ax.set_xlim(38.45,38.55)
#
lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
fig.legend(lines, labels, loc='lower center', bbox_to_anchor=(0.5, -0.04), ncol=4, fontsize="14")
#fig.legend(lines, labels, loc='lower center', ncol=4, fontsize="14")
#==============
for i in range(len(runs)):
  out[i].close()
#==============
plt.tight_layout()
plt.savefig('../figures/asicsmed_tra.pdf', bbox_inches='tight')
# plt.show()
#===============================
#### create figure 2
#======================================================
# compute B0 and ustar
#======================================================
ndays     = 30; time1     = 0; time2     = 613
out[0]    = Dataset(scm[0].output, mode='r')
ustar     = np.sqrt(np.sqrt( np.power((scm[0].rho0/1.22)*out[0].variables['taux'][time1:time2],2)
                           + np.power((scm[0].rho0/1.22)*out[0].variables['tauy'][time1:time2],2) ))
time      = out[0].variables['ocean_time'][time1:time2] / (3600*24)
# data to compute B0
Qsol      = out[0].variables['Qs'][time1:time2] /(scm[0].rho0*scm[0].cp)   #K m/s
Qnsol     = out[0].variables['Qns'][time1:time2]/(scm[0].rho0*scm[0].cp)   #K m/s
Fw        = out[0].variables['Fw'][time1:time2]    #psu m / s
Tt        = out[0].variables['temp'][time1:time2,-1] # sst
Ts        = out[0].variables['salt'][time1:time2,-1] # sss
out[0].close()
#=================================================================
# Estimate surface values of alpha and beta using nonlinear EOS
#=================================================================
Q00=+999.842594 ; Q01=+6.793952e-2; Q02=-9.095290e-3
Q03=+1.001685e-4; Q04=-1.120083e-6; Q05=+6.536332e-9
U00=+0.824493   ; U01=-4.08990e-3 ; U02=+7.64380e-5
U03=-8.24670e-7 ; U04=+5.38750e-9 ; V00=-5.72466e-3
V01=+1.02270e-4 ; V02=-1.65460e-6 ; W00=+4.8314e-4
sqrtTs=np.sqrt(Ts)
cff   = 1./scm[0].rho0
alpha = -cff * ( Q01+Tt*(2.*Q02+Tt*(3.*Q03+Tt*(4.*Q04+Tt*5.*Q05)))
        +Ts*(U01+Tt*(2.*U02+Tt*(3.*U03+Tt*4.*U04))+sqrtTs*(V01+Tt*2.*V02)))
beta  = cff*( U00+Tt*(U01+Tt*(U02+Tt*(U03+Tt*U04)))+1.5*(V00+Tt*(V01+Tt*V02))*sqrtTs+2.*W00*Ts)
Buoy_sfc = scm[0].g*(alpha*(Qsol+Qnsol)-beta*Fw)
#==================================================================
tscm  = np.array([20,40,51]);tscm  = tscm*12  # draw vertical blue lines on the plot
#===============================================
# Computation of hmxl from observations
#===============================================
indx   = np.array([25,24,23,22,21,20,19,17,16,15,14,13,11,10,9,8,7,5,4,3,2])
nc_dta    = Dataset(data_file, mode='r')
dta_depth = nc_dta.variables['z_r'][:]
time_dta  = nc_dta.variables['ocean_time'][:]/(3600*24)
ntime     = len(time_dta)
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
    rho,bvf = scm_oce.rho_eos(temp,salt,z_r,z_w,scm[0].rho0,nz)
    drho    = 0.01; zref = z_r[-5]
    hmxl_dta[tt] = scm_oce.compute_mxl2(bvf,0.01,z_r,-300.,scm[0].rho0,nz)
#========================================================
#### create figure
f, axes = plt.subplots(nrows = 3, ncols = 1, sharex=True, sharey=False)
f.set_figheight(5)
f.set_figwidth(5)
#
ax = axes[0]
ax.set_title(r"$u_\star\;[{\rm m}\;{\rm s}^{-1}]$",fontsize=12)
for i in range(len(tscm)):
    ax.axvline(x = time[tscm[i]], color = 'b', linestyle = '-', linewidth=1.0)
ax.plot(time, ustar,color='k',linewidth=1.5)
ax.set_xlim(7.5,25.75)
#
ax = axes[1]
ax.set_title(r"$B_0 \; [{\rm m}^2\;{\rm s}^{-3}]$",fontsize=12)
ax.plot(time, Buoy_sfc,color='k',linewidth=1.5)
ax.axhline(y = 0., color = 'r', linestyle = '-', linewidth=1.0)
for i in range(len(tscm)):
    ax.axvline(x = time[tscm[i]], color = 'b', linestyle = '-', linewidth=1.0)
#
ax = axes[2]
ax.set_title(r"$h_{\rm mxl} \; [{\rm m}]$",fontsize=12)
ax.plot(time_dta, gaussian_filter(hmxl_dta,sigma=3),color='k',linewidth=1.5,label='mooring data')
#
for i in range(len(runs)):
  out[i] = Dataset(scm[i].output, mode='r')
  ax.plot( out[i].variables['ocean_time'][time1:time2] / (3600*24), out[i].variables['hmxl_300m'][time1:time2],
                                              color = colors[i],linewidth=1.5,linestyle=linestyle[i],label=run_label[i])
  out[i].close()
#
for i in range(len(tscm)):
    ax.axvline(x = time[tscm[i]], color = 'b', linestyle = '-', linewidth=1.0)
ax.set_xticks([7,9,11,13,15,17,19,21,23,25])
ax.set_xticklabels([r'$22\;{\rm Jan}$',r'$24\;{\rm Jan}$',r'$26\;{\rm Jan}$',r'$28\;{\rm Jan}$',r'$30\;{\rm Jan}$',r'$01\;{\rm Feb}$',r'$03\;{\rm Feb}$',r'$05\;{\rm Feb}$',r'$07\;{\rm Feb}$',r'$09\;{\rm Feb}$'])
ax.legend(loc=3,fontsize='6')
ax.set_ylim(-2300.,0.)
plt.tight_layout()
plt.savefig('../figures/asicsmed_hmxl.pdf', bbox_inches='tight')
# plt.show()

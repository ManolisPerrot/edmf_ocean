#!/usr/bin/env python
# coding: utf-8
__docformat__ = 'reStructuredText'
###########################################
# Imports
###########################################
import numpy as np
import matplotlib
matplotlib.rcParams.update(matplotlib.rcParamsDefault)
import matplotlib.pyplot as plt
from netCDF4 import Dataset
###########################################
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
################################################################################
#
class Lion_data:
    def __init__(self,data_params):
        default_params = {
            "filename": "dataset.nc",  "time_var_name": "time",
            "temp_var_name": "temp", "salt_var_name": "salt",
            "zr_var_name": "z_r", "time_start": 0, "time_end": 24*30,
            "time_scaling":1./(24.*3600.), "time_shift":0., "z_scaling":1.,"data_title":"ED"
        }
        params = default_params.copy()
        params.update(data_params)
        self.config = params

    def get_data(self):
        nc = Dataset(self.config["filename"] , mode='r')
        ts = self.config["time_start"]
        tf = self.config["time_end"]
        self.temp = nc.variables[self.config["temp_var_name"]][ts:tf,:]
        self.salt = nc.variables[self.config["salt_var_name"]][ts:tf,:]
        zz        = self.config["z_scaling"]*nc.variables[self.config["zr_var_name"]][:]
        tt        = self.config["time_scaling"]*nc.variables[self.config["time_var_name"]][ts:tf] - self.config["time_shift"]
        self.nz   = len(zz)
        self.nt   = len(tt)
        self.z1d  = zz
        # for the asics_med data flip arrays vertically
        if self.config["z_scaling"] < 0.:
            self.z1d  = np.flip(zz)
            self.temp = np.flip(self.temp, axis=1)
            self.salt = np.flip(self.salt, axis=1)
        nc.close()
        self.time = tt
        self.rho0 = 1027.
#
    def get_sfc_forcing(self):
        nc = Dataset(self.config["filename"] , mode='r')
        ts = self.config["time_start"]
        tf = self.config["time_end"]
        rhoA = 1.22;  cp = 3985.
        self.ustar = np.sqrt(np.sqrt( np.power((self.rho0/rhoA)*nc.variables['taux'][ts:tf],2)
                                    + np.power((self.rho0/rhoA)*nc.variables['tauy'][ts:tf],2) ))
        self.qsol  = nc.variables['Qs'][ts:tf] /(self.rho0*cp)   #K m/s
        self.qnsol = nc.variables['Qns'][ts:tf]/(self.rho0*cp)   #K m/s
        self.Fw    = nc.variables['Fw'][ts:tf]    #psu m / s
        self.sst   = nc.variables['temp'][ts:tf,-1] # sst
        self.sss   = nc.variables['salt'][ts:tf,-1] # sss
        self.g     = 9.81
        nc.close()
#
    def get_surface_buoyancy(self):
        Q00=+999.842594 ; Q01=+6.793952e-2; Q02=-9.095290e-3
        Q03=+1.001685e-4; Q04=-1.120083e-6; Q05=+6.536332e-9
        U00=+0.824493   ; U01=-4.08990e-3 ; U02=+7.64380e-5
        U03=-8.24670e-7 ; U04=+5.38750e-9 ; V00=-5.72466e-3
        V01=+1.02270e-4 ; V02=-1.65460e-6 ; W00=+4.8314e-4
        sqrtTs=np.sqrt(self.sss)
        Tt    = self.sst
        Ts    = self.sss
        cff   = 1./self.rho0
        alpha = -cff * ( Q01+Tt*(2.*Q02+Tt*(3.*Q03+Tt*(4.*Q04+Tt*5.*Q05)))
              +Ts*(U01+Tt*(2.*U02+Tt*(3.*U03+Tt*4.*U04))+sqrtTs*(V01+Tt*2.*V02)))
        beta  = cff*( U00+Tt*(U01+Tt*(U02+Tt*(U03+Tt*U04)))+1.5*(V00+Tt*(V01+Tt*V02))*sqrtTs+2.*W00*Ts)
        self.buoy_sfc = self.g*(alpha*(self.qsol+self.qnsol)-beta*self.Fw)
#
    def get_hmxl(self):
        self.hmxl = np.zeros(self.nt)
        for kt in range(self.nt):
            self.hmxl[kt] = compute_mxl_temp(self.temp[kt,:],0.05,self.z1d,-300.)
#
def compute_mxl_temp(temp, dtemp, zr, zref):
    nz     = len(temp)
    # Find the starting index kstart where depth is just below zref
    kstart = nz-1
    while zr[kstart] > zref and kstart > 0:
        kstart -= 1
    # Compute temperature gradient dTdz
    dTdz = np.zeros(nz)
    for k in range(nz - 1):
        dTdz[k] = (temp[k+1]-temp[k]) / (zr[k+1] - zr[k])
    # Initialize values
    t_c = dtemp
    hmxl = zr[0]  # Initialize at the near bottom value
    cff_old = 0.0
    # Compute mixed layer depth
    for k in range(kstart, 1, -1):
        cff_new = cff_old + max(dTdz[k-1], 0.) * (zr[k] - zr[k-1])
        if cff_new >= t_c:
            hmxl = ((cff_new - t_c) * zr[k] + (t_c - cff_old) * zr[k - 1]) / (cff_new - cff_old)
            break
        cff_old = cff_new

    return hmxl
#=========================================================
# Creating different instance and calling class methods
#=========================================================
evd_params = {'filename': '../data/run_asics1.nc', 'time_var_name': 'ocean_time', 'data_title':'TKE+EVD'}
scm_evd = Lion_data(evd_params)
scm_evd.get_data()
scm_evd.get_sfc_forcing()
scm_evd.get_surface_buoyancy()
scm_evd.get_hmxl()

edmf_params = {'filename': '../data/run_asics3.nc', 'time_var_name': 'ocean_time', 'data_title':'EDMF-Energy'}
scm_edmf = Lion_data(edmf_params)
scm_edmf.get_data()
scm_edmf.get_hmxl()

data_params = {'filename': '../data/AsicsMed_data.nc', 'time_var_name': 'time_counter',
               'temp_var_name': 'votemper', 'salt_var_name': 'vosaline', 'zr_var_name': 'deptht',
               'time_start': 11, 'time_end': 11+30*2,'time_scaling':1.,'time_shift':15,'z_scaling':-1., 'data_title':'Lion observational data'}
scm_data = Lion_data(data_params)
scm_data.get_data()
scm_data.get_hmxl()
#=================================================================
# create the figure
#=================================================================
f, axes = plt.subplots(nrows = 3, ncols = 1, sharex=True, sharey=False)
f.set_figheight(5); f.set_figwidth(5)
### Ustar
ax = axes[0]
ax.set_title(r"$u_\star\;[{\rm m}\;{\rm s}^{-1}]$",fontsize=12)
ax.plot(scm_evd.time, scm_evd.ustar,color='k',linewidth=1.5)
ax.set_xlim(0.0,25.)
### Surface buoyancy
ax = axes[1]
ax.set_title(r"$B_0 \; [{\rm m}^2\;{\rm s}^{-3}]$",fontsize=12)
ax.plot(scm_evd.time, scm_evd.buoy_sfc,color='k',linewidth=1.5)
ax.axhline(y = 0., color = 'r', linestyle = '-', linewidth=1.0)
### hmxl
ax = axes[2]
ax.set_title(r"$h_{\rm mxl} \; [{\rm m}]$",fontsize=12)
ax.set_ylim(-2300.,0.)
ax.plot(scm_evd.time, scm_evd.hmxl,color='tab:gray',linewidth=1.5,label=scm_evd.config["data_title"])
ax.plot(scm_edmf.time, scm_edmf.hmxl,color='tab:green',linewidth=1.5,label=scm_edmf.config["data_title"])
ax.plot(scm_data.time, scm_data.hmxl,color='k',linewidth=1.5,label=scm_data.config["data_title"])
ax.set_xticks([1,3,5,7,9,11,13,15,17,19,21,23,25])
ax.set_xticklabels([r'$16\;{\rm Jan}$',r'$18\;{\rm Jan}$',r'$20\;{\rm Jan}$',r'$22\;{\rm Jan}$',r'$24\;{\rm Jan}$',r'$26\;{\rm Jan}$',r'$28\;{\rm Jan}$',r'$30\;{\rm Jan}$',r'$01\;{\rm Feb}$',r'$03\;{\rm Feb}$',r'$05\;{\rm Feb}$',r'$07\;{\rm Feb}$',r'$09\;{\rm Feb}$'], rotation=45)
ax.set_xlabel("Time")
ax.legend(loc=3,fontsize='6')
plt.tight_layout()
plt.savefig('../figures/ASICS_MED_frc_hmxl.pdf', bbox_inches='tight')



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
from scm_class_oce import SCM
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from unit_tests import is_in_range

###################################################
plt.rcParams['text.usetex'] = True
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['mathtext.it'] = 'STIXGeneral:italic'
plt.rcParams['mathtext.bf'] = 'STIXGeneral:italic:bold'
plt.rcParams.update({'font.size': 18})
plt.rcParams.update({'figure.facecolor': 'white'})
plt.rcParams.update({'savefig.facecolor': 'white'})
###########################################


# colors
blue, orange, magenta, grey, green = '#0db4c3', '#eea021', '#ff0364', '#606172', '#3fb532'

# ===========================================================================
case = 'W005_C500_NO_COR'


saving_path = '../figures/'
saving_name = case+'_dt_dz_MLD_sensitivity.png'


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
# mld = (-z_r_les[(-WTH[instant]).argmax()])

# ===========================================================================

# Define the common parameters:
common_params = {
    'nz': 100,
    'dt': 50.,
    'h0': 2000.,
    'thetas': 6.5,
    'hc': 1e16,
    'nbhours': 72,
    'outfreq': 1,
    'output_filename': "scm_output.nc",
    'T0': 2., 
    'N0': 1.9620001275490499e-6,
    'Tcoef': 0.2048,
    'SaltCst': 35.,
    'lat0': 0.,
    'sustr': 0.,
    'svstr': 0.,
    'stflx': -500.,
    'srflx': 0.,
    'ssflx': 0.,
    'eddy_diff': True,
    'evd': False,
    'mass_flux_tra': True,
    'mass_flux_dyn': True,
    'mass_flux_tke': True,
    'mass_flux_tke_trplCorr': True,
    'mass_flux_small_ap': True,
    'lin_eos': True,
    'extrap_ak_surf': True,
    'tke_sfc_dirichlet': False,
    'eddy_diff_tke_const': 'NEMO',
    'entr_scheme': 'R10',
    'Cent': 0.99,
    'Cdet': 1.99,       # 'Cdet': 2.5,
    'wp_a': 1.,   #1
    'wp_b': 1.,     #1.
    'wp_bp': 0.003*250,    
    'up_c': 0.5, 
    'vp_c': 0.5,
    'bc_ap': 0.2,    #0.3,
    'delta_bkg': 0.0045*250,   # 0.005,
    'wp0' : -0.5e-08,
    'output_filename': 'run',
    'beta_bc_P09': 0.3,
    'write_netcdf': True
}



if case == 'W05_C500':
    common_params['sustr'] = 0.5/1027

if case == 'W005_C500_NO_COR':
    common_params['sustr'] = 0.05/1027

# Define parameters specific to each run (overwrite common parameters):

# run_label = [r'$\textit{\textbf{u}}_{h,p} = \overline{\textbf{\textit{u}}}$', r'$\textbf{\textit{u}}_{h,p} \neq \overline{\textbf{\textit{u}}}$']
run_label = [r'dz=40m', r'dz=20m',r'dz=2m',r'dt=50s',r'dt=1800s']
runs = [
    {
        'nz': 50,
        'output_filename': 'run1.nc',
    },
    {
        'nz': 100,
        'output_filename': 'run2.nc'
    },
        {
        'nz': 1000,
        'output_filename': 'run3.nc'
    },
    {
        'dt': 50,
        'output_filename': 'run4.nc'
    },
        {
        'dt': 1800,
        'output_filename': 'run5.nc'
    }

        ]



scm = [0]*len(runs)

# Run the SCM
for i, run_params in enumerate(runs):
    params = common_params.copy()  # Create a copy of common_params
    params.update(run_params)  # Update with run_params
    scm[i] = SCM(params,cfg_params={}) #empty cfg_params because its is already contained in params...
    scm[i].run_direct()

    # test zinv
    if scm[i].MF_tra or scm[i].MF_dyn: 
        print(run_label[i])
        reference=-320
        #reference=-mld #LES value
        is_in_range(value=scm[i].zinv, value_name='zinv', reference=reference,tolerance=10 )

# LOAD outputs

out = [0]*len(runs)

for i, run_params in enumerate(runs):
    print('opening '+run_params['output_filename'])
    out[i] = xr.open_dataset(run_params['output_filename'])
# TH_scm = [0]*len(runs)
# U_scm = [0]*len(runs)
# V_scm = [0]*len(runs)

# #interpolate scm output on LES #TODO do the converse to reduce computation cost?


# for i, run_params in enumerate(runs):
#     TH_scm = scm[i].t_history
#     U_scm  = scm[i].u_history
#     V_scm  = scm[i].v_history




###################################################################################
# PLOTTING
###################################################################################
### in the fortran routine, does not corresponds to a grind point but to an analtical reconstruction!!

mld = (-z_r_les[(-WTH).argmax(axis=1)])
time = 3600*out[0]['time']

###Analytical expression

#fixed parameters
alpha = 1.9*10**(-4)
g = 9.81
rho0 = 1027.
cp = 3985.
Q_h = 500. #heat flux (Wm-2)
DT = 1.
DZ = 1000.
N = np.sqrt(alpha*g*DT/DZ)
#N = np.sqrt(7.65*10**(-7))
Q_b = g*alpha*Q_h/(rho0*cp)
#Q_b = 1*10**(-7)
# print('N=',N ,'N^2=',N**2 ,'Q_b=', Q_b)
t = np.linspace(0,3,72) #in days
h3 = -np.sqrt(3*Q_b*t*24*3600/N**2)
h = -np.sqrt(Q_b*t*24*3600/N**2)

marker =    ['o'    ,'o'    ,'o'    ,'D'    ,'o']
fillstyle = ['none' ,'full' ,'none' ,'none' ,'full']
markersize= [10     ,10     ,10     ,10     ,5    ]
fig, ax = plt.subplots(nrows=1, ncols=1, sharex=False,
                         sharey=True,figsize=(12,6))
start=0

# ax.plot(out[0]['time'][start:-1],-mld[start:], 'k',linewidth=3, label='LES')
# ax.plot(out[0]['time'][start:-1],np.sqrt(2.8)*h[start:], color='tab:gray',linewidth=3, label=r'$-\sqrt{2.8 \frac{B_0}{N^2} t}$')


for i, label in enumerate(run_label):
    ax.plot(out[0]['time'][start:-1], (out[i]['z_w'][(out[i]['WT']).argmax(axis=1)] )[start:-1] ,marker=marker[i], fillstyle=fillstyle[i], markersize=markersize[i],linestyle='None',label=label)

ax.set_ylim(-320,-0)
ax.legend()
ax.set_ylabel(r'$m$')
ax.set_xlabel(r'time (hours)')
fig.suptitle(r'Mixed Layer Depth')
plt.savefig(saving_path+saving_name, bbox_inches='tight', dpi=300)
print('fig saved at '+saving_path+saving_name)
# plt.show()

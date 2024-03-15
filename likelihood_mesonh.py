#!/usr/bin/env python
# coding: utf-8

#!/usr/bin/env python
# coding: utf-8

###########################################
# Imports
###########################################

import sys  # to put the SCM into the PYTHONPATH
sys.path.append('edmf/library/F2PY')
from sys import exit
import time as TIME
import xarray as xr
from scipy.interpolate import interp1d
import scipy.signal
from scm_class import SCM
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from case_configs import case_params
from multiprocess import Pool #multiprocessING cannot handle locally defined functions, multiprocess can

# ===================================Choose cases/datasets========================================

cases = ['FC500', 'W005_C500_NO_COR']

# ===================================Choose hyperparameters of the calibration===========
# dimensional error tolerance for L2 norm 
model_error_t,data_error_t=0.001,0.001 #°C
model_error_u,data_error_u=0.001,0.001 #ms-1 
model_error_v,data_error_v=0.001,0.001 #ms-1
# importance of each field in the cost function (non-dimensional)
weight_t=1.
weight_u=1.
weight_v=1.
# Bayesian beta hyperparameter
beta_t = weight_t * (model_error_t**2 + data_error_t**2) 
beta_u = weight_u * (model_error_u**2 + data_error_u**2) 
beta_v = weight_v * (model_error_v**2 + data_error_v**2) 

# use H1 Sobolev norm
sobolev=True
#If sobolev=True, the following hyperparameters are
# dimensional error tolerance for H1 norm
model_error_dz_t,data_error_dz_t=0.001,0.001 #°Cm-1
model_error_dz_u,data_error_dz_u=0.001,0.001 #s-1
model_error_dz_v,data_error_dz_v=0.001,0.001 #s-1
# importance of each field in the cost function (non-dimensional)
weight_dzt=1.
weight_dzu=1.
weight_dzv=1.
#TODO find the form of beta for H1
beta_t_h1 = 1.
beta_u_h1 = 1.
beta_v_h1 = 1.
# ===================================Load LES========================================
TH_les = {}
U_les  = {}
V_les  = {}
time = {}
z_r_les = {}
z_w_les = {}
dz_TH_les= {}
dz_U_les = {}
dz_V_les = {}
for case in cases:
    file = 'GN_01.1.OC_01.000_copy.nc'
    path= 'edmf/data/'+case+'/'
    les = xr.open_dataset(path+file)
    LG_MEAN = xr.open_dataset(path+file,group ='/LES_budgets/Mean/Cartesian/Not_time_averaged/Not_normalized/cart')
    TH_les[case] = (LG_MEAN.MEAN_TH - 273.15).data.transpose() #transpose to have coordinates as level then time, 
                                                         #as in the SCM  
    U_les[case]=LG_MEAN.MEAN_U.data.transpose()
    V_les[case]=LG_MEAN.MEAN_V.data.transpose()
    time_les = les.time_les 
    time[case] = ((time_les - time_les[0]) / np.timedelta64(1, 'h')).data.astype(int) + 1 #numpy array of integer hours, starting at inital time + 1h
    z_r_les[case] = (les.level_les - (les.level_les[0] + les.level_les[-1])).data #remap level_les on negative depth values
    z_w_les[case] = (les['level_w'][1:] - (les.level_les[0] + les.level_les[-1])).data
    if sobolev==True:
        dz_TH_les[case] = (np.divide( ((TH_les[case][1:,:]-TH_les[case][0:-1,:])).T , (z_r_les[case][1:]-z_r_les[case][0:-1]))).T
        dz_U_les [case] = (np.divide( ((U_les[case][1:,:]-U_les[case][0:-1,:])  ).T , (z_r_les[case][1:]-z_r_les[case][0:-1]))).T
        dz_V_les [case] = (np.divide( ((V_les[case][1:,:]-V_les[case][0:-1,:])  ).T , (z_r_les[case][1:]-z_r_les[case][0:-1]))).T

# ====================================Define configurations=======================
# Define the common parameters (attention some of them will be overwritten by case_configurations.py):
common_params = {
    'nz': 100,
    'dt': 50.,
    'h0': 2000.,
    'thetas': 6.5,
    'hc': 400,
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
    'mass_flux_small_ap': False,
    'lin_eos': True,
    'extrap_ak_surf': True,
    'tke_sfc_dirichlet': False,
    'eddy_diff_tke_const': 'NEMO',
    'entr_scheme': 'R10',
    'Cent': 0.99,
    'Cdet': 1.99,       
    'wp_a': 1.,
    'wp_b': 1.25,      
    'wp_bp': 0.003,    
    'up_c': 0.5,
    'vp_c': 0.5,
    'bc_ap': 0.2,    
    'delta_bkg': 0.005,
    'wpmin'    : 1.e-08,
    'output_filename': 'run'
}


def likelihood_mesonh(    
    Cent      = 0.99,
    Cdet      = 1.99,
    wp_a      = 1.,
    wp_b      = 1.25,
    wp_bp     = 0.003,
    up_c      = 0.5,
    bc_ap     = 0.2,
    delta_bkg = 0.005,
    wpmin     = 1.e-08):

    # Load the case specific parameters
    # ATTENTION, any parameter entered in case params will 
    # OVERWRITE common params. Double-check scm_configs before running

    params_to_estimate = {
                        'Cent': Cent,
                        'Cdet': Cdet,
                        'wp_a': wp_a,
                        'wp_b': wp_b,
                        'wp_bp': wp_bp,
                        'up_c': up_c,
                        'vp_c': up_c,
                        'bc_ap': bc_ap,
                        'delta_bkg': delta_bkg,
                        'wpmin' : wpmin}

    scm = {}

    def likelihood_of_one_case(case_index):
        case=cases[case_index]
        # ====================================Run the SCM cases=======================================
        params = common_params.copy()  # Create a copy of common_params
        params.update(case_params[case])  # Update with the specific case hyperparameters in case_params[case]
        params.update(params_to_estimate) # Update with the parameters to estimate 
        scm[case] = SCM(**params)
        scm[case].run_direct()            # Run the SCM
        print(case+": zinv =", scm[case].zinv)
        
        #interpolate scm output on LES #TODO do the converse to reduce computation cost?
        TH_scm = scm[case].t_history
        U_scm = scm[case].u_history
        V_scm = scm[case].v_history
        z_r_scm = scm[case].z_r
        TH_scm_int = np.zeros( TH_les[case].shape )
        TH_scm_int = interp1d(z_r_scm, TH_scm, kind = 'linear', axis=0)(z_r_les[case])
        U_scm_int = np.zeros( U_les[case].shape )
        U_scm_int = interp1d(z_r_scm, U_scm, kind = 'linear', axis=0)(z_r_les[case])
        V_scm_int = np.zeros( V_les[case].shape )
        V_scm_int = interp1d(z_r_scm, V_scm, kind = 'linear', axis=0)(z_r_les[case])

        if sobolev==True:
            dz_TH_scm_int = np.divide((TH_scm_int[1:,:]-TH_scm_int[0:-1,:]).T, z_r_scm[1:]-z_r_scm[0:-1] ).T
            dz_U_scm_int =  np.divide((U_scm_int[1:,:]-U_scm_int[0:-1,:]  ).T, z_r_scm[1:]-z_r_scm[0:-1] ).T
            dz_V_scm_int =  np.divide((V_scm_int[1:,:]-V_scm_int[0:-1,:]  ).T, z_r_scm[1:]-z_r_scm[0:-1] ).T

        # compute the space-time L2 average,
        # divided by total depth and duration
        #trapz is a trapezoidal integral 

        metric_t = np.trapz( np.trapz( (TH_scm_int - TH_les[case])**2, z_r_les[case], axis=0) , time[case]) * 1/(z_r_les[case][-1] - z_r_les[case][0]) * 1 / (time[case][-1] - time[case][0]) 

        metric_u = np.trapz( np.trapz( (U_scm_int - U_les[case])**2, z_r_les[case], axis=0) , time[case]) * 1/(z_r_les[case][-1] - z_r_les[case][0]) * 1 / (time[case][-1] - time[case][0]) 

        metric_v = np.trapz( np.trapz( (V_scm_int - V_les[case])**2, z_r_les[case], axis=0) , time[case]) * 1/(z_r_les[case][-1] - z_r_les[case][0]) * 1 / (time[case][-1] - time[case][0]) 

        likelihood = np.exp(-beta_t*metric_t - beta_u*metric_u - beta_v*metric_v)

        if sobolev==True:
            metric_t_h1 = np.trapz( np.trapz( (dz_TH_scm_int - dz_TH_les[case])**2, z_w_les[case][1:-1], axis=0) , time[case]) * 1/(z_w_les[case][1:-1][-1] - z_w_les[case][1:-1][0]) * 1 / (time[case][-1] - time[case][0])  

            metric_u_h1 = np.trapz( np.trapz( (dz_U_scm_int - dz_U_les[case])**2, z_w_les[case][1:-1], axis=0) , time[case]) * 1/(z_w_les[case][1:-1][-1] - z_w_les[case][1:-1][0]) * 1 / (time[case][-1] - time[case][0])  
            
            metric_v_h1 = np.trapz( np.trapz( (dz_V_scm_int - dz_V_les[case])**2, z_w_les[case][1:-1], axis=0) , time[case]) * 1/(z_w_les[case][1:-1][-1] - z_w_les[case][1:-1][0]) * 1 / (time[case][-1] - time[case][0])  

            likelihood = likelihood * np.exp(-beta_t_h1*metric_t_h1 - beta_u_h1*metric_u_h1 - beta_v_h1*metric_v_h1)
        
        return likelihood

    likelihoods=np.zeros(len(cases))
    # parrallelized for-loop
    if __name__ == '__main__':
        with Pool() as p:
            likelihoods = p.map(likelihood_of_one_case, range(likelihoods.size))

    # total likelihood is the product of likelihood of each case
    return np.prod(likelihoods)

likelihood_mesonh()

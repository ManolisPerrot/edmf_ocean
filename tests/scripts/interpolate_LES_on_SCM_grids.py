#!/usr/bin/env python
# coding: utf-8

###########################################
# Imports
###########################################

import sys  # to put the SCM into the PYTHONPATH
from sys import exit
import time as TIME
import xarray as xr
from scipy.interpolate import interp1d
import scipy.signal
from scipy.interpolate import griddata
from scm_class import SCM
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from case_configs import case_params, default_params
from multiprocess import Pool #multiprocessING cannot handle locally defined functions, multiprocess can
import subprocess


# ===================================Functions========================================

def get_SCM_grid(case):
    ''' Run one hour of the SCM to get vertical grids z_r_scm and z_w_scm '''
    params = default_params  # Create empty custom parameters
    params.update(case_params[case])  # Update with the specific case hyperparameters in case_params[case]
                                      # containing grid specification
    #print(params)
    params.update(  {'nbhours': 1}) #specify only one hour of simulation 
    #print(params)
    scm = SCM(**params)
    scm.run_direct()            # Run the SCM
    #print(case+": zinv =", scm.zinv)
    return scm.z_r, scm.z_w, {key: params[key] for key in ['nz','h0','thetas']}


def interpolate_les_on_scm(X_les, z_les, z_scm):
    # Interpolate les on scm grid. 
    # If les points are out out scm grid, NaNs are put 
    X_les_int = griddata(z_les, X_les, z_scm, method='linear')

    # Find indices of non-NaN values in X_les_int
    valid_indices = ~np.isnan(X_les_int[:,0])

    # Filter X_les_int and z_scm based on valid indices
    X_les_int_filtered = X_les_int[valid_indices]
    z_scm_filtered = z_scm[valid_indices] 

    return X_les_int_filtered, z_scm_filtered, valid_indices

#===========================================================================================
default_cases=[key for key in case_params]

def regrid_and_save(cases=default_cases):
  
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
          file = 'GN_01.1.OC_01.000.nc'
          path= './data/'+case+'/'
          les = xr.open_dataset(path+file)
          print('opening', path+file)
          LG_MEAN = xr.open_dataset(path+file,group ='/LES_budgets/Mean/Cartesian/Not_time_averaged/Not_normalized/cart')
          TH_les[case] = (LG_MEAN.MEAN_TH - 273.15).data.transpose() #transpose to have coordinates as level then time, 
                                                                  #as in the SCM  
          U_les[case]=LG_MEAN.MEAN_U.data.transpose()
          V_les[case]=LG_MEAN.MEAN_V.data.transpose()
          time_les = les.time_les 
          time[case] = (((time_les - time_les[0]) / np.timedelta64(1, 'h')).data.astype(int) + 1)*np.timedelta64(1, 'h') #numpy array of integer hours, starting at inital time + 1h
          z_r_les[case] = (les.level_les - (les.level_les[0] + les.level_les[-1])).data #remap level_les on negative depth values
          z_w_les[case] = (les['level_w'][1:] - (les.level_les[0] + les.level_les[-1])).data
          dz_TH_les[case] = (np.divide( ((TH_les[case][1:,:]-TH_les[case][0:-1,:])).T , (z_r_les[case][1:]-z_r_les[case][0:-1]))).T
          dz_U_les [case] = (np.divide( ((U_les[case][1:,:]-U_les[case][0:-1,:])  ).T , (z_r_les[case][1:]-z_r_les[case][0:-1]))).T
          dz_V_les [case] = (np.divide( ((V_les[case][1:,:]-V_les[case][0:-1,:])  ).T , (z_r_les[case][1:]-z_r_les[case][0:-1]))).T

          z_r_scm, z_w_scm, grid_params = get_SCM_grid(case)
          print('Interpolate LES on SCM with grid spec', grid_params)

          TH_les_int,z_r_scm_filtered, z_r_boolean_filter = interpolate_les_on_scm(TH_les[case], z_r_les[case], z_r_scm )
          U_les_int ,z_r_scm_filtered, z_r_boolean_filter = interpolate_les_on_scm(U_les[case], z_r_les[case], z_r_scm )
          V_les_int ,z_r_scm_filtered, z_r_boolean_filter = interpolate_les_on_scm(V_les[case], z_r_les[case], z_r_scm )
          
          # z_w_boolean filter is not z_w_scm 
          # but the z_w that we would obtain computing dz_X_scm
          dz_TH_les_int,z_w_scm_filtered, z_w_boolean_filter = interpolate_les_on_scm(dz_TH_les[case], .5*(z_r_les[case][1:]+z_r_les[case][:-1]), .5*(z_r_scm[1:]+z_r_scm[:-1]) )
          dz_U_les_int ,z_w_scm_filtered, z_w_boolean_filter = interpolate_les_on_scm(dz_U_les[case] , .5*(z_r_les[case][1:]+z_r_les[case][:-1]), .5*(z_r_scm[1:]+z_r_scm[:-1]) )
          dz_V_les_int ,z_w_scm_filtered, z_w_boolean_filter = interpolate_les_on_scm(dz_V_les[case] , .5*(z_r_les[case][1:]+z_r_les[case][:-1]), .5*(z_r_scm[1:]+z_r_scm[:-1]) )

          dsout = xr.Dataset( 
                  data_vars={     'TH_les':    (['time','z_r'], TH_les_int.T ),
                                  'U_les':     (['time','z_r'], U_les_int.T ),
                                  'V_les':     (['time','z_r'], V_les_int.T ), 
                                  'dz_TH_les': (['time','z_w'], dz_TH_les_int.T ),
                                  'dz_U_les':  (['time','z_w'], dz_U_les_int.T ),
                                  'dz_V_les':  (['time','z_w'], dz_V_les_int.T )
                        }, 
                  coords={        'z_r': z_r_scm_filtered,
                                  'z_w': z_w_scm_filtered,
                                  'time': time[case],
                                  'z_r_boolean_filter': z_r_boolean_filter,
                                  'z_w_boolean_filter': z_w_boolean_filter
                          },
                  attrs={         'Description': 'LES interpolated on the SCM grid defined in case_configs',
                                  'Case': case,
                                  'SCM grid nz':  grid_params['nz'],
                                  'SCM grid h0':  grid_params['h0'],
                                  'SCM grid thetas':  grid_params['thetas'],
                          })
          
          dsout.to_netcdf(path+case+'_interpolated_on_SCM.nc', mode='w', format="NETCDF4")
          print('Interpolated LES saved at', path+case+'_interpolated_on_SCM.nc' )

if __name__ == '__main__':
     regrid_and_save()
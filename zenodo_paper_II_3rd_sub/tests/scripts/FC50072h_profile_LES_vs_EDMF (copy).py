#!/usr/bin/env python
# coding: utf-8

###########################################
# Imports
###########################################
from sys import exit
import time as TIME
import xarray as xr
from scipy.interpolate import interp1d
import scipy.signal
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
# ===========================================================================
case = 'FC500'
saving_path = '../figures/'
saving_name = case+'72h_profile_LES_vs_EDMF.png'
# ===========================================================================
# loading LES output 
file = 'GN_01.1.OC_01.000_copy.nc'
path = '../data/'+case+'/'
les = xr.open_dataset(path+file)
LG_MEAN = xr.open_dataset(
    path+file, group='/LES_budgets/Mean/Cartesian/Not_time_averaged/Not_normalized/cart')
TH_les = (LG_MEAN.MEAN_TH - 273.15).data
U_les = (LG_MEAN.MEAN_U).data
V_les = (LG_MEAN.MEAN_V).data
#
LG_RES = xr.open_dataset(
    path+file, group='/LES_budgets/Resolved/Cartesian/Not_time_averaged/Not_normalized/cart')
LG_SBG = xr.open_dataset(
    path+file, group='/LES_budgets/Subgrid/Cartesian/Not_time_averaged/Not_normalized/cart')
# Tubulent vertical fluxes
WTH = (LG_RES.RES_WTH + LG_SBG.SBG_WTHL).data
WU = (LG_RES.RES_WU + LG_SBG.SBG_WU).data
WV = (LG_RES.RES_WV + LG_SBG.SBG_WV).data
corr_tke_les = 0. #-2e-04+1.0e-08
TKE = (LG_RES.RES_KE + LG_SBG.SBG_TKE).data + corr_tke_les
WTKE = (LG_RES.RES_WKE + LG_SBG.SBG_WTKE).data
#
time_les = les.time_les
# numpy array of integer hours, starting at inital time + 1h
time = ((time_les - time_les[0]) / np.timedelta64(1, 'h')).data.astype(int) + 1
# remap level_les on negative depth values
z_r_les = (les.level_les - (les.level_les[0] + les.level_les[-1])).data
instant=-1
mld = (-z_r_les[(-WTH[instant]).argmax()])
# ===========================================================================



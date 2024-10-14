#!/usr/bin/env python
# coding: utf-8
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

def load_MNH(cases):
#========================================
    les = {}
    LG_MEAN={}
    LG_RES={}
    LG_SBG={}
    BU_KE={}
    TH_les={}
    U_les={}
    V_les={}
    WTH={}
    WU={}
    WV={}
    TKE={}
    WTKE={}
    dz_WU_les ={}
    dz_WV_les ={}
    time_les={}
    hours = {}
    z_r_les = {}
    mld = {}
    for case in cases:
        file = 'GN_01.1.OC_01.000.nc'
        path = '../../tests/data/'+case+'/'
        les[case] = xr.open_dataset(path+file)
        print('Opening LES case: ',case)
        LG_MEAN[case]= xr.open_dataset(path+file,group ='/LES_budgets/Mean/Cartesian/Not_time_averaged/Not_normalized/cart')
        LG_RES[case] = xr.open_dataset(path+file,group ='/LES_budgets/Resolved/Cartesian/Not_time_averaged/Not_normalized/cart')    
        LG_SBG[case] = xr.open_dataset(path+file,group ='/LES_budgets/Subgrid/Cartesian/Not_time_averaged/Not_normalized/cart')
        BU_KE[case] = xr.open_dataset(path+file,group ='/LES_budgets/BU_KE/Cartesian/Not_time_averaged/Not_normalized/cart')
        TH_les[case] = (LG_MEAN[case].MEAN_TH - 273.15)
        U_les [case] = (LG_MEAN[case].MEAN_U)
        V_les [case] = (LG_MEAN[case].MEAN_V)
        # V_les [case] = (LG_MEAN[case].MEAN_U)
        # U_les [case] = (LG_MEAN[case].MEAN_V)
        WTH   [case] = (LG_RES[case].RES_WTH + LG_SBG[case].SBG_WTHL)
        WU    [case] = (LG_RES[case].RES_WU  + LG_SBG[case].SBG_WU)
        WV    [case] = (LG_RES[case].RES_WV  + LG_SBG[case].SBG_WV)
        # WV    [case] = (LG_RES[case].RES_WU  + LG_SBG[case].SBG_WU)
        # WU    [case] = (LG_RES[case].RES_WV  + LG_SBG[case].SBG_WV)
        TKE   [case] = (LG_RES[case].RES_KE  + LG_SBG[case].SBG_TKE).data 
        WTKE  [case] = (LG_RES[case].RES_WKE + LG_SBG[case].SBG_WTKE).data
        # remap level_les on negative depth values
        z_r_les[case] = (les[case].level_les - (les[case].level_les[0] + les[case].level_les[-1])).data
        dz_WU_les[case] = (WU[case][:].mean(dim='time_les').data[1:]-WU[case][:].mean(dim='time_les').data[:-1])/(z_r_les[case][1:]-z_r_les[case][:-1])
        dz_WV_les[case] = (WV[case][:].mean(dim='time_les').data[1:]-WV[case][:].mean(dim='time_les').data[:-1])/(z_r_les[case][1:]-z_r_les[case][:-1])
        time_les[case] = les[case].time_les
        hours[case] = ((time_les[case] - time_les[case][0]) / np.timedelta64(1, 'h')).data.astype(int) + 1

        mld[case] = (-z_r_les[case][(-WTH[case]).argmax(axis=1)])
    return les,LG_MEAN,LG_RES,LG_SBG,BU_KE,TH_les,U_les,V_les,WTH,WU,WV,TKE,WTKE,dz_WU_les ,dz_WV_les ,time_les,hours,z_r_les,mld


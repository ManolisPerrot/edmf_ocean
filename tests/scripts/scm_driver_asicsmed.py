#!/usr/bin/env python
# coding: utf-8
__docformat__ = 'reStructuredText'
############################################################################
# Imports
############################################################################
from sys import exit
import sys
sys.path.insert(0,"../../library/F2PY/")
#import add
import numpy as np
import matplotlib
#matplotlib.use('TkAgg')
#matplotlib.use('pdf')
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scm_class_oce_asicsmed import SCM
#from F2PY import scm_class
import scipy.signal
############################################################################
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
#===========================================================================
config_params = { 'nz': 75, 'dt':1200, 'lat0':42.04}

run_label = ['TKE+EVD', 'EDMF', 'EDMF-Engy', 'Keps','EDMF-boundary',]
out_files = ['run_asics_'+run_label[i]+'.nc' for i in range(len(run_label))]

runs = [
    {
        'eddy_diff': True,
        'evd': True,
        'eddy_diff_scheme' : 'TKE',
        'output_filename': out_files[0]
    },
    {
        'eddy_diff': True,
        'evd': False,
        'mass_flux_tra': True,
        'mass_flux_dyn': True,
        'mass_flux_tke': False,
        'mass_flux_tke_trplCorr': False,
        'mass_flux_small_ap': True,
        'entr_scheme': 'R10',
        'eddy_diff_scheme' : 'TKE',
        'output_filename': out_files[1]
    },
    {
        'eddy_diff': True,
        'evd': False,
        'mass_flux_tra': True,
        'mass_flux_dyn': True,
        'mass_flux_tke': True,
        'mass_flux_tke_trplCorr': True,
        'mass_flux_small_ap': True,
        'entr_scheme': 'R10',
        'eddy_diff_scheme' : 'TKE',
        'output_filename': out_files[2]
    },
    {
        'eddy_diff': True,
        'evd': False,
        'eddy_diff_scheme' : 'Keps',
        'output_filename': out_files[3]
    },
        {
        'eddy_diff': True,
        'evd': False,
        'mass_flux_tra': True,
        'mass_flux_dyn': True,
        'mass_flux_tke': True,
        'mass_flux_tke_trplCorr': True,
        'mass_flux_small_ap': True,
        'entr_scheme': 'R10',
        'eddy_diff_scheme' : 'TKE',
        'output_filename': out_files[4],
        'bc_P09': 'inconsistent',
    },
        ]
#
scm = [0]*len(runs)
# Run the SCM
for i, run_params in enumerate(runs):
    scm[i] = SCM(run_params, config_params)
    scm[i].run_direct()
#===========================================================================

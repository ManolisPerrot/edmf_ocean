

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
from scm_class import SCM
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from case_configs import case_params

case ='WANG1_FR'

common_params = {'a': 1
}
print('common_params', common_params)


print('case params', case_params[case])

common_params = common_params.update(case_params[case])  # Update with the specific case configuration in case_params[case]
print('common_params', common_params)



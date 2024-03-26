#!/usr/bin/env python
# coding: utf-8


### File describing each idealized runs

#The comprehensive list of parameters is (see also )
#    {
#     'nz': 100,                                                             
#     'dt': 50.,                           
#     'h0': 2000.,                         
#     'thetas': 6.5,                           
#     'hc': 400,                           
#     'nbhours': 72,                           
#     'outfreq': 1,                        
#     'output_filename': "scm_output.nc",                          
#     'T0': 2.,                        
#     'N0': 1.9620001275490499e-6,                         
#     'Tcoef': 0.2048,                         
#     'SaltCst': 35.,                          
#     'lat0': 0.,                          
#     'sustr': 0.,                         
#     'svstr': 0.,                         
#     'stflx': -500.,                          
#     'srflx': 0.,                         
#     'ssflx': 0.,                         
#     'eddy_diff': True,                           
#     'evd': False,                        
#     'mass_flux_tra': True,                           
#     'mass_flux_dyn': True,                           
#     'mass_flux_tke': True,                           
#     'mass_flux_tke_trplCorr': True,                          
#     'mass_flux_small_ap': False,                         
#     'lin_eos': True,                         
#     'extrap_ak_surf': True,                          
#     'tke_sfc_dirichlet': False,                          
#     'eddy_diff_tke_const': 'NEMO',                           
#     'entr_scheme': 'R10',                        
#     'Cent': 0.99,
#     'Cdet': 2.5,
#     'wp_a': 1.,
#     'wp_b': 1.,
#     'wp_bp': 0.002,
#     'up_c': 0.5,
#     'vp_c': 0.5,
#     'bc_ap': 0.3,
#     'delta_bkg': 0.006,
#     'output_filename': 'run'
# }


"""[summary]
        Args:
            nz: Number of grid points. Defaults to 100.
            h0: Depth of the water column. Defaults to 1000.
            thetas: stretching parameter toward the surface for the vertical grid. Defaults to 6.5.
            hc: Stretching parameter for the vertical grid. Resolution is almost constant between 0 and hc. Defaults to 400.
            T0: Initial surface temperature. Defaults to 2 (Celsius).
            N0: Initial stratification. Defaults to 1.9620001275490499E-6.
            Tcoef: Thermal expansion coefficient to define initial temperature profile. Defaults is 0.2048
            SaltCst: constant salinity value for the initial conditions. Default is 35.
            lat0: Latitude of the water column (used to compute Coriolis frequency). Defaults to 45 (degree).
            sustr, svstr: Horizontal components of the wind stress. Defaults to 0 (m^2.s^-2).
            stflx: non-penetrative surface heat flux. Defaults to -500 (W.m^-2).
            srflx: penetrative surface heat flux. Defaults to 0 (W.m^-2).
            ssflx: freshwater flux. Defaults to 0 (psu.m.s^-1).
            nbhours: length of the simulation in hours. Defaults to 72 (hours).
            outfres: frequency of the outputs in hours. Defaults is 1 (hour).
            output_filename (str): File in which the output will be store. Defaults to "scm_output.nc".
            eddy_diff (bool): Do we apply the eddy-diffusion term ? Defaults to True.
            evd (bool): Do we apply enhanced vertical diffusion ? Defaults to False
            mass_flux_tra (bool): Do we apply the mass flux to the tracer equations ? Defaults to False.
            mass_flux_dyn (bool): Do we apply the mass flux to the horizontal momentum equations ? Defaults to False.
            mass_flux_tke (bool): Do we apply the transfer of KE/PE to TKE due to MF ? Defaults to False.
            mass_flux_tke_trplCorr (bool): Do we apply the modification of the w'e turbulents flux due to MF ? Defaults to False.
            mass_flux_small_ap (bool): Do we apply the small area approximation to the mass-flux equations ? Defaults to True
            lin_eos (bool): use a linear equation of state. Defaults to True  (NOTE: this is the only working options so far)
            extrap_ak_surf (bool): extrapolate the values of eddy-diffusivity/viscosity toward the surface. Defaults to True
            tke_sfc_dirichlet (bool): nature of the surface boundary condition for TKE. Defaults to True
            eddy_diff_tke_const: constants to be used in the TKE scheme ('NEMO', 'MNH' or 'RS81').    Defaults to 'NEMO'
"""

### ATTENTION, any parameter entered here will 
### OVERWRITE the parameter entered in the SCM simulation

case_params = {
    'WANG1_FR':{
        'h0': 4500.,
        'T0': 3.,
        'N0': 1.5865490613891073e-08, #equals to f(60°)
        'SaltCst': 32.6,
        'lat0': 60.,
        'sustr': 0.,
        'svstr': 0.,
        'stflx': -111.12982, #Wm-2
        'srflx': 0.,
        'ssflx': 0.,
                },
    'W005_C500_NO_COR':{
        'h0': 2000.,
        'T0': 2.,
        'N0': 1.9620001275490499e-6,
        'SaltCst': 32.6,
        'lat0': 0.,
        'sustr': 0.05/1027,
        'svstr': 0.,
        'stflx': -500.,
        'srflx': 0.,
        'ssflx': 0.,
    },
    'FC500':{
        'h0': 2000.,
        'T0': 2.,
        'N0': 1.9620001275490499e-6,
        'SaltCst': 32.6,
        'lat0': 0.,
        'sustr': 0.,
        'svstr': 0.,
        'stflx': -500.,
        'srflx': 0.,
        'ssflx': 0.,
    }
                }



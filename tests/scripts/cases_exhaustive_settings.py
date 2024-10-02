#!/usr/bin/env python
# coding: utf-8


### File describing each idealized runs

"""[summary]
        Args:
            nz: Number of grid points. Defaults to 100.
            dt: Dynamical time-step. Defaults to 30.
            h0: Depth of the water column. Defaults to 1000.
            thetas: stretching parameter toward the surface for the vertical grid. Defaults to 6.5.
            hc: Stretching parameter for the vertical grid. Resolution is almost constant between 0 and hc. Defaults to 400.
            T0: Initial surface temperature. Defaults to 2 (Celsius).
            N0: Initial stratification. Defaults to 1.9620001275490499E-6.
            gridType: which type of vertical grid? Defaults to 'croco_new'; variants: 'croco_old', 'ORCA75'
            mld_ini: initial mixed layer depth (m, NEGATIVE). Default to 0. m.
            mld_iniS: initial mixed layer depth for Salinity (m, NEGATIVE). Default to 0. m.
            lin_eos (bool): use a linear equation of state. Defaults to True
            Tcoef: Thermal expansion coefficient to define initial temperature profile. Defaults is 0.2048
            Scoef: Thermal expansion coefficient to define initial temperature profile. Defaults is 0.2048
            SaltCst: constant salinity value for the initial conditions. Default is 35.
            Tref: reference temperature in linear EOS
            Sref: reference salinity in linear EOS
            lat0: Latitude of the water column (used to compute Coriolis frequency). Defaults to 45 (degree).
            sustr, svstr: Horizontal components of the wind stress. Defaults to 0 (m^2.s^-2).
            stflx: non-penetrative surface heat flux. Defaults to -500 (W.m^-2).
            srflx: penetrative surface heat flux. Defaults to 0 (W.m^-2).
            ssflx: freshwater flux. Defaults to 0 (psu.m.s^-1).
            diurnal_cycle (bool):
            btflx: tracer bottom noudary condition ('no flux' or 'linear_continuation')
            nbhours: length of the simulation in hours. Defaults to 72 (hours).
            outfreq: frequency of the outputs in hours. Defaults is 1 (hour).
            output_filename (str): File in which the output will be store. Defaults to "scm_output.nc".
            eddy_diff (bool): Do we apply the eddy-diffusion term ? Defaults to True.
            evd (bool): Do we apply enhanced vertical diffusion ? Defaults to False
            mass_flux_tra (bool): Do we apply the mass flux to the tracer equations ? Defaults to False.
            mass_flux_dyn (bool): Do we apply the mass flux to the horizontal momentum equations ? Defaults to False.
            mass_flux_tke (bool): Do we apply the transfer of KE/PE to TKE due to MF ? Defaults to False.
            mass_flux_tke_trplCorr (bool): Do we apply the modification of the w'e turbulents flux due to MF ? Defaults to False.
            mass_flux_small_ap (bool): Do we apply the small area approximation to the mass-flux equations ? Defaults to True
            extrap_ak_surf (bool): extrapolate the values of eddy-diffusivity/viscosity toward the surface. Defaults to True
            tke_sfc_dirichlet (bool): nature of the surface boundary condition for TKE. Defaults to True
            eddy_diff_tke_const: constants to be used in the TKE scheme ('NEMO', 'MNH' or 'RS81').    Defaults to 'NEMO'
            write_netcdf (bool): Do we write a netcdf file (containing detailed energy diagnostics, but slower to run)? Default to False
            bc_P09 (str): Do you want to use Pergaud (2009) temperature boundary condition? Defaults to 'false'; options are 'inconsistent' (original version) or 'consistent' (fixed version) 
"""


case_settings={
###############
#Alexandre cases        
###############        
'kato_phillips':{   
    'dt': 360.,
    'gridType': 'croco_new',                           
    'nz': 256,
    'h0': 50.,                         
    'thetas': 6.5,                           
    'hc': 50,                           
    'nbhours': 31,                           
    'outfreq': 1,                        
    'T0': 20.,       
    'mld_ini_temp': -0., #initial MLD in temperature                
    'N0': 9.81e-5,
    'SaltCst': 35.,      # constant initial profile of Salt
    'cpoce': 3985,
    'lin_eos': True,
    'rho0': 1024.763,
    'alpha': 2e-4,    
    'Tref':  20.,         #reference temperature 
    'beta': 8e-4 ,           
    'Sref':  35.,        #reference salinity                 
    'lat0': 0.,                          
    'sustr': 9.758353835241997e-5,         #w'u'_ocean_sfc                
    'svstr': 0.,         #w'v'_ocean_sfc       
    'stflx': -0.,      #non-penetrative surface heat flux                  
    'srflx': 0.,         #penetrative surface heat flux      
    'ssflx': 0.,         #freshwater flux
    'diurnal_cycle': False,
    'btflx': 'no flux' }
    ,
    'FC':{   
    'dt': 360.,
    'gridType': 'croco_new',                           
    'nz': 512,
    'h0': 100.,                         
    'thetas': 6.5,                           
    'hc': 100,                           
    'nbhours': 7*24+1,                           
    'outfreq': 1,                        
    'T0': 20.,       
    'mld_ini_temp': -0., #initial MLD in temperature                
    'N0': 1.9620001275490499e-5,
    'SaltCst': 35.,      # constant initial profile of Salt
    'cpoce': 3985,
    'lin_eos': True,
    'rho0': 1024.763,
    'alpha': 2e-4,    
    'Tref':  20.,         #reference temperature 
    'beta': 8e-4 ,           
    'Sref':  35.,        #reference salinity                 
    'lat0': 0.,                          
    'sustr': 0,         #w'u'_ocean_sfc                
    'svstr': 0.,         #w'v'_ocean_sfc       
    'stflx': -75.,      #non-penetrative surface heat flux                  
    'srflx': 0.,         #penetrative surface heat flux      
    'ssflx': 0.,         #freshwater flux
    'diurnal_cycle': False,
    'btflx': 'no flux' }
    ,
    #TODO
    # 'DC':{   
    # 'dt': 360.,
    # 'gridType': 'croco_old',                           
    # 'nz': 256,
    # 'h0': 50.,                         
    # 'thetas': 6.5,                           
    # 'hc': 50,                           
    # 'nbhours': 7*24+1,                           
    # 'outfreq': 1,                        
    # 'T0': 20.,       
    # 'mld_ini_temp': -0., #initial MLD in temperature                
    # 'N0': 1.9620001275490499e-5,
    # 'SaltCst': 35.,      # constant initial profile of Salt
    # 'cpoce': 3985,
    # 'lin_eos': True,
    # 'rho0': 1024.763,
    # 'alpha': 2e-4,    
    # 'Tref':  20.,         #reference temperature 
    # 'beta': 8e-4 ,           
    # 'Sref':  35.,        #reference salinity                 
    # 'lat0': 0.,                          
    # 'sustr': 9.758353835241997e-5,         #w'u'_ocean_sfc                
    # 'svstr': 0.,         #w'v'_ocean_sfc       
    # 'stflx': -75.,      #non-penetrative surface heat flux                  
    # 'srflx': 0.,         #penetrative surface heat flux      
    # 'ssflx': 0.,         #freshwater flux
    # 'diurnal_cycle': False,  # TODO : implémenter le type3 / IB dans edmf, max{cos[2π*(t/86400 – 0.5)], 0} * 235.62, [W/m2]
    # 'btflx': 'no flux' },
    #
    #
###############
#Manolis cases        
###############
    # MesoNH constants
    # cp = 4178 # hard coded in set_rsou.f90 in FINMANOLIS
    # ##### MésoNH ini_cst.f90 constants ###
    # XRH00OCEAN = 1027.
    # XTH00OCEAN = 286.65
    # XSA00OCEAN = 32.6
    # XP00OCEAN = 201.E5
    # XG = 9.80665  #----> PB! si on veut mettre g en argument du SCM, il faut changer toutes les subroutines... 
    # XALPHAOC = 1.9E-4  # Coefficient of thermal expansion of water (K-1)
    # XBETAOC = 7.7475E-4  # Coeff of Haline contraction coeff (S-1)
    'FC500':{   
    'dt': 50.,
    'gridType': 'croco_new',                           
    'nz': 100,
    'h0': 2000.,                         
    'thetas': 6.5,                           
    'hc': 400,                           
    'nbhours': 72,                           
    'outfreq': 1,                        
    'T0': 2.,       
    'mld_ini_temp': -0., #initial MLD in temperature                
    'N0': 1.9620001275490499e-6,
    'SaltCst': 32.6,      # constant initial profile of Salt
    'cpoce': 3985.,
    'lin_eos': True,
    'rho0': 1027.,
    'alpha': 2e-4,      
    'Tref':  13.5,         #reference temperature 
    'beta': 7.7475E-4,      
    'Sref':  32.6,        #reference salinity                 
    'lat0': 0.,                          
    'sustr': 0.,         #w'u'_ocean_sfc                
    'svstr': 0.,         #w'v'_ocean_sfc       
    'stflx': -500.,      #non-penetrative surface heat flux                  
    'srflx': 0.,         #penetrative surface heat flux      
    'ssflx': 0.,         #freshwater flux
    'diurnal_cycle': False,  
    'btflx': 'no flux' }
    ,
    'W05_C500':{   
    'dt': 50.,
    'gridType': 'croco_new',                           
    'nz': 100,
    'h0': 2000.,                         
    'thetas': 6.5,                           
    'hc': 400,                             
    'nbhours': 72,                           
    'outfreq': 1,                        
    'T0': 2.,       
    'mld_ini_temp': -0., #initial MLD in temperature                
    'N0': 1.9620001275490499e-6,
    'SaltCst': 32.6,      # constant initial profile of Salt
    'cpoce': 3985,
    'lin_eos': True,
    'rho0': 1027.,
    'alpha': 1.9e-4,      
    'Tref':  13.5,         #reference temperature 
    'beta': 7.7475E-4,      
    'Sref':  32.6,        #reference salinity                  
    'lat0': 0.,                          
    'sustr': 0.5/1027,         #w'u'_ocean_sfc                
    'svstr': 0.,         #w'v'_ocean_sfc       
    'stflx': -500.,      #non-penetrative surface heat flux                  
    'srflx': 0.,         #penetrative surface heat flux      
    'ssflx': 0.,         #freshwater flux
    'diurnal_cycle': False,  
    'btflx': 'no flux' }
    ,
    'W005_C500_NO_COR':{   
    'dt': 50.,
    'gridType': 'croco_new',                           
    'nz': 100,
    'h0': 2000.,                         
    'thetas': 6.5,                           
    'hc': 400,                             
    'nbhours': 72,                           
    'outfreq': 1,                        
    'T0': 2.,       
    'mld_ini_temp': -0., #initial MLD in temperature                
    'N0': 1.9620001275490499e-6,
    'SaltCst': 32.6,      # constant initial profile of Salt
    'cpoce': 3985,
    'lin_eos': True,
    'rho0': 1027.,
    'alpha': 2e-4,      
    'Tref':  13.5,         #reference temperature 
    'beta': 7.7475E-4,      
    'Sref':  32.6,        #reference salinity                  
    'lat0': 0.,                          
    'sustr': 0.05/1027,         #w'u'_ocean_sfc                
    'svstr': 0.,         #w'v'_ocean_sfc       
    'stflx': -500.,      #non-penetrative surface heat flux                  
    'srflx': 0.,         #penetrative surface heat flux      
    'ssflx': 0.,         #freshwater flux
    'diurnal_cycle': False,  
    'btflx': 'no flux' }
}


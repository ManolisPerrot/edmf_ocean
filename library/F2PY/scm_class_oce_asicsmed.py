#!/usr/bin/env python
# coding: utf-8
__docformat__ = 'reStructuredText'
## Imports
from scmoce_new import scm_oce_new
from initialize import (
    initialize_physical_constants,
    initialize_sfc_boundary_conditions,
    initialize_eddy_diffusion,
    initialize_mass_flux,
    initialize_vertical_grid,
    initialize_model_from_file,
    initialize_ed_arrays,
    initialize_mf_arrays,
    initialize_output,
    initialize_forcing_from_file,
    validate_config
)
from edmf import do_TKE, do_MF
from keps import do_KEPS
from output_asicsmed import output_init, output_state
from datetime import datetime
from netCDF4 import Dataset
import numpy  as np
################################################################################
# Class
################################################################################
class SCM:
    """Single Column model class"""
    ############################################################################
    def __init__(self, run_params, cfg_params):
        # Load and validate configuration
        self.config = self._load_config(run_params, cfg_params)
        validate_config(self)
        # Set grid and tracer parameters
        self.nz        = self.config["nz"]                             ## number of grid points
        self.ntra      = 2                                             ## number of tracers
        self.itemp, self.isalt = 0, 1                                  ## Tracer indices
        self.dt        = self.config["dt"]                             ## time-step of the simulation [seconds]
        dt1            = datetime(1850,1, 1,00,00)
        dt2            = datetime(2013,1,15,00,00)  # initial time for the simulation
        dt3            = datetime(2013,2,15,00,00)  # final time for the simulation
        timediff1      = dt2-dt1
        timediff2      = dt3-dt2
        self.inifrctime= int(timediff1.total_seconds()/3600.)
        self.nbsteps   = int( ( timediff2.total_seconds() )/self.dt )           ## number of time steps
        self.outfreq   = int(  self.config["outfreq"]*3600.)           ## frequency of outputs in seconds
        self.frcname   = self.config["sfcforc_filename"]
        self.kfrcStr   = 336; self.kfrcEnd   = 1681
        self.rhoc10    = 0.01         # density criterion for mixed layer depth (with reference level at 10m)
        self.rhoc300   = 0.01-0.00001 # density criterion for mixed layer depth (with reference level at 300m)
        # Initialize model components
        self._initialize_model()
        self.hmxl10    = np.array(-self.Hz[-1])
        self.hmxl300   = np.array(-self.Hz[-1])
        self.bc_P09    = 'false'
        self.trad_coriolis_mod = False

    ############################################################################
    def _initialize_model(self):
        """Initialize all model components."""
        initialize_physical_constants(self)       ## Physical constants
        initialize_sfc_boundary_conditions(self)
        initialize_eddy_diffusion(self)           ## Initialize eddy diffusion scheme
        initialize_mass_flux(self)                ## Initialize mass flux scheme
        initialize_vertical_grid(self)            ## Initialize vertical grid
        initialize_model_from_file(self)
        initialize_ed_arrays(self)                ## Define arrays for eddy-diffusion
        if self.config["mass_flux_tra"]:          ## Define arrays for mass-flux
            initialize_mf_arrays(self)
        initialize_output(self)
        initialize_forcing_from_file(self)

    ############################################################################
    def _load_config(self, run_params, cfg_params):
        """Return default configuration."""
        default_config = {
            "nz": 75,  "dt": 1200.0, "h0": 2600.0, "gridType": "ORCA75",
            "cpoce": 3985.0, "lin_eos": False, "rho0": 1027.0,"lat0": 42.04,
            "initial_filename":"../data/init_ASICS_m01d15.nc",
            "sfcforc_filename":"../data/forc_ASICS_y2013.nc",
            "nbhours": 720, "outfreq": 1.0, "output_filename": "run_asics.nc",
            "eddy_diff": True, "evd": False,
            "eddy_diff_scheme": "TKE", "mass_flux_tra": False,
            "mass_flux_dyn": False, "mass_flux_tke": False,
            "mass_flux_tke_trplCorr": False, "mass_flux_small_ap": True,
            "extrap_ak_surf": True, "eddy_diff_tke_const": "NEMO",
            "tke_sfc_dirichlet": False, "akvmin": 1.3e-06, "aktmin": 1.4e-07,
            "mxlmin": 1.0, "Cent": 0.7, "Cdet": 1.99, "wp_a": 1.0,
            "wp_b": 1.0, "wp_bp": 0.125, "up_c": 0.25,
            "vp_c": 0.25, "bc_ap": 0.20, "delta_bkg": 0.0625,
            "wp0":-1.0e-08, "entr_scheme": "R10"
        }
        config = default_config.copy()
        config.update(run_params)
        config.update(cfg_params)
        return config

    ############################################################################
    def run_direct(self):
        kout=0
        output_init(self) #to initialize .nc dataset
        #
        sigma    = np.zeros(self.nz)
        rr       = -(self.z_r[:]+self.config["h0"])/self.Hz[:]
        sigma    = 0.5*(1.+np.tanh(2.*(rr-0.64)))
        traLS    = np.zeros(2)
        #
        swr_frac = scm_oce_new.lmd_swfrac(self.Hz,self.nz)   ## Compute fraction of solar shortwave flux penetrating to specified depth
        if self.ED_scheme=='TKE': do_TKE( self )             ## Initialize eddy-diffusion scheme to compute eddy diffusivity/eddy viscosity for the first time-step
        if self.ED_scheme=='Keps': do_KEPS( self )
        MF_sfc_flux = 0.                                     ## Initialize ED sfc flux correction (since Fmass is not rigouroulsy zero at sfc)
        stflx0             = np.zeros(2)
        stflx0[self.isalt] = self.stflx[self.isalt]
        ####################################
        # Time integration loop
        #-----------------------------------
        for kt in range(self.nbsteps):
            # update current time
            traLS[self.itemp] = self.t_n  [17,self.itemp]
            traLS[self.isalt] = self.t_n  [17,self.isalt]
            time = self.dt*float(kt+1)                 ## time in seconds at the end of the current time-step
            self._compute_surface_forcing(time)  ## Compute solar radiation (possibly with diurnal cycle)
            self._advance_tracers(swr_frac)   ## advance tracers to n+1,star (vertical diffusion only)
            self._apply_penalization_tra(sigma,traLS)
            self._advance_dynamics()                  ## advance dynamics to n+1,star (Coriolis + vertical viscosity only)
            self._apply_penalization_dyn(sigma)

            if self.MF_tra or self.MF_dyn:
                self._compute_mass_flux()

            if self.MF_tra:
                self._advance_tracers_MF()
            else:
                self.t_np1[:,self.itemp] = self.t_np1_star[:,self.itemp]
                self.t_np1[:,self.isalt] = self.t_np1_star[:,self.isalt]

            if self.MF_dyn:
                self._advance_dynamics_MF()  ## apply the MF term to the velocity equation + compute the transfer of KE to subgrid scales stored in shearMF
            else:
                self.u_np1[:] = self.u_np1_star[:]
                self.v_np1[:] = self.v_np1_star[:]

            if self.ED_scheme=='TKE' : do_TKE( self )     ## Compute eddy-viscosity / eddy-diffusivity  (TKE scheme)
            if self.ED_scheme=='Keps': do_KEPS( self )    ## Compute eddy-viscosity / eddy-diffusivity  (k-epsilon scheme)
            #==========================================================
            if  time % self.outfreq == 0:
                self._do_turb_fluxes (  )                    ## compute turbulent fluxes to be stored in output file
                # Compute mixed layer depth from bvf
                salt = 0.*self.t_n[:,self.isalt] + 38.5     #salt[:] = self.t_np1[:,self.isalt]
                rho,bvf_loc = scm_oce_new.rho_eos(self.t_np1_star[:,self.itemp],salt,self.z_r,self.z_w,self.rho0,self.nz)
                self.hmxl10,self.hmxl300  = scm_oce_new.compute_mxl(bvf_loc,self.rhoc10,self.rhoc300,self.z_r,self.rho0,self.nz)
                kout+=1; output_state(self,TimeInSecs=time,kout=kout) ## Write outputs in .nc file
            #===================================================================
            self._swap_arrays() ## swap arrays
            #===================================================================
            if kt<self.nbsteps-1: self.wtke[:]  = 0.                            ## diagnostics : reset the array containing w'e

    ############################################################################
    def _compute_surface_forcing(self, time):
        """Advance tracers to the next time step."""
        if time >= self.frc_time[self.kfrc]:   # we need to read a new time
            self.kfrc = self.kfrc + 1
            self.taux1 = self.taux2; self.tauy1 = self.tauy2;
            self.Qs1   = self.Qs2  ; self.Qnet1 = self.Qnet2; self.Emp1 = self.Emp2
            indx=self.kfrcStr+self.kfrc
            frc01      = Dataset(self.frcname, mode='r',format="NETCDF4")
            self.taux2 = np.squeeze(frc01.variables['TAUX2'][indx,0,0])
            self.tauy2 = np.squeeze(frc01.variables['TAUY2'][indx,0,0])
            self.Qs2   = np.squeeze(frc01.variables['FSOL2'][indx,0,0])
            self.Qnet2 = np.squeeze(frc01.variables['FNET2'][indx,0,0])
            self.Emp2  = np.squeeze(frc01.variables['EMP2'][indx,0,0]) * 1.e+6
            frc01.close()
            #======= surface forcings
        if time < self.frc_time[0]:
            wght = 0.
        else:
            wght = (time-self.frc_time[self.kfrc-1])/(self.frc_time[self.kfrc]-self.frc_time[self.kfrc-1])
        irho0cp = 1. / (self.cp*self.rho0)
        self.ustr_sfc  = ((1.-wght)*self.taux1+wght*self.taux2) / self.rho0         # in m2/s2
        self.vstr_sfc  = ((1.-wght)*self.tauy1+wght*self.tauy2) / self.rho0         # in m2/s2
        self.srflx     = ((1.-wght)*self.Qs1 + wght*self.Qs2)*irho0cp
        self.stflx[self.itemp] = ((1.-wght)*(self.Qnet1-self.Qs1)+wght*(self.Qnet2-self.Qs2))*irho0cp
        coeffEMP               = self.t_n[-1,self.isalt] / self.rho0
        self.stflx[self.isalt] = coeffEMP*( (1.-wght)*self.Emp1 + wght*self.Emp2 )

    ############################################################################
    def _compute_solar_radiation(self,time):
        modulation = 1.
        if self.DC: modulation = max(np.cos(2.*np.pi*(time/86400. - 0.5)),0.)
        self.srflx = modulation*self.Qswmax

    ############################################################################
    def _advance_tracers(self, swr_frac):
        """Advance tracers to the next time step."""
        # Implement tracer advancement logic here
        self.t_np1_star, self.tFlx = scm_oce_new.advance_tra_ed(
                                        self.t_n, self.stflx, self.srflx,
                                        swr_frac, self.btflx, self.Hz   , self.akt  ,
                                        self.z_r, 0.*self.eps_n, self.eos_params[1],
                                        self.dt , self.cp, self.nz   , self.ntra  )

    ############################################################################
    def _advance_dynamics(self):
        """Advance dynamics to the next time step."""
        self.u_np1_star, self.v_np1_star, self.uFlx, self.vFlx = scm_oce_new.advance_dyn_cor_ed(
                                        self.u_n, self.v_n, self.ustr_sfc, self.vstr_sfc,
                                        self.ustr_bot, self.vstr_bot,
                                        self.Hz , self.akv, self.z_r, self.fcor, self.dt, self.nz  )

    ############################################################################
    def _compute_mass_flux(self):
        """ Solve plume equations and compute mass flux"""
        do_MF( self )
        if self.MF_small_ap:
            self.Fmass = -(self.ap[:]*self.wp[:])
        else:
            self.Fmass = -(self.ap[:]*self.wp[:])/(1.-self.ap[:])

    ############################################################################
    def _advance_tracers_MF(self):
        """Advance tracers to the next time step adding the mass flux contribution"""
        self.t_np1 = scm_oce_new.advance_tra_mf(
                          self.t_np1_star, self.tp, self.Fmass, self.tFlx, self.Hz,
                          self.dt, self.nz, self.ntra )

    ############################################################################
    def _advance_dynamics_MF(self):
        """Advance dynamics to the next time step adding the mass flux contribution"""
        self.u_np1, self.v_np1 = scm_oce_new.advance_dyn_mf(
                                      self.u_np1_star, self.v_np1_star, self.shearMF,
                                      self.u_n  , self.v_n  , self.up , self.vp ,
                                      self.Fmass, self.uFlx, self.vFlx, self.Hz   , self.dt , self.nz  )
    ############################################################################
    def _apply_penalization_tra(self,sigma,traLS):
        """Apply penalization near the bottom """
        self.t_np1_star[:,self.itemp] = (self.t_np1_star[:,self.itemp]+sigma[:]*self.dt*traLS[self.itemp])/(1.+sigma[:]*self.dt)
        self.t_np1_star[:,self.isalt] = (self.t_np1_star[:,self.isalt]+sigma[:]*self.dt*traLS[self.isalt])/(1.+sigma[:]*self.dt)

    ############################################################################
    def _apply_penalization_dyn(self,sigma):
        """Apply penalization near the bottom """
        self.u_np1_star[:] = np.exp(-self.dt*sigma[:])*self.u_np1_star[:]
        self.v_np1_star[:] = np.exp(-self.dt*sigma[:])*self.v_np1_star[:]

    ############################################################################
    def _swap_arrays(self):
        self.u_n[:]    = self.u_np1[:]
        self.v_n[:]    = self.v_np1[:]
        self.t_n[:,:]  = self.t_np1[:,:]
        self.tke_n[:]  = self.tke_np1[:]
        self.eps_n [:] = self.eps_np1 [:]

    ############################################################################
    def _do_turb_fluxes(self):
        """Turbulent fluxes diagnostics"""
        self.wted[  self.nz] = self.stflx[self.itemp]                            ## (w'theta')_ED
        self.wted[1:self.nz] = self.akt[1:self.nz] * (  self.t_np1_star[1:self.nz  ,self.itemp]
                                                      - self.t_np1_star[0:self.nz-1,self.itemp] ) / self.Hz[0:self.nz-1]
        self.wued[  self.nz] = self.ustr_sfc                                     ## (w'theta')_ED
        self.wved[  self.nz] = self.vstr_sfc

        self.wued[1:self.nz] = self.akv[1:self.nz] * (  self.u_np1_star[1:self.nz] - self.u_np1_star[0:self.nz-1] ) / self.Hz[0:self.nz-1]
        self.wved[1:self.nz] = self.akv[1:self.nz] * (  self.v_np1_star[1:self.nz]
                                                      - self.v_np1_star[0:self.nz-1] ) / self.Hz[0:self.nz-1]
        if self.MF_tra:                                                         ## (w'theta')_MF
          self.wtmf[1:self.nz] =  self.Fmass[1:self.nz]*(
                                              self.tp   [1:self.nz,  self.itemp]
                                            - self.t_np1_star[0:self.nz-1,self.itemp] )

        if self.MF_dyn:                                                         ## (w'u')_MF
          self.wumf[1:self.nz] =  self.Fmass[1:self.nz]*(
                                               self.up   [1:self.nz]
                                             - self.u_np1_star[0:self.nz-1] )
          self.wvmf[1:self.nz] =  self.Fmass[1:self.nz]*(
                                               self.vp   [1:self.nz]
                                             - self.v_np1_star[0:self.nz-1] )

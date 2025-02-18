#!/usr/bin/env python
# coding: utf-8
__docformat__ = 'reStructuredText'
## Imports
from scmoce_new import scm_oce_new
from initialize import (
    initialize_physical_constants,
    initialize_boundary_conditions,
    initialize_eddy_diffusion,
    initialize_mass_flux,
    initialize_vertical_grid,
    initialize_initial_state,
    initialize_ed_arrays,
    initialize_mf_arrays,
    initialize_output_and_averaging,
    initialize_diagnostics,
    validate_config
)
from edmf import do_TKE, do_MF
from keps import do_KEPS
from output import output_init, output_state
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
        self.nbsteps   = int( (self.config["nbhours"]*3600.)/self.dt ) ## number of time steps
        self.outfreq   = int(  self.config["outfreq"]*3600.)           ## frequency of outputs in seconds
        # Initialize model components
        self._initialize_model()

    ############################################################################
    def _initialize_model(self):
        """Initialize all model components."""
        initialize_physical_constants(self)       ## Physical constants
        initialize_boundary_conditions(self)      ## Boundary conditions
        initialize_eddy_diffusion(self)           ## Initialize eddy diffusion scheme
        initialize_mass_flux(self)                ## Initialize mass flux scheme
        initialize_vertical_grid(self)            ## Initialize vertical grid
        initialize_initial_state(self)            ## Define initial state
        initialize_ed_arrays(self)                ## Define arrays for eddy-diffusion
        if self.config["mass_flux_tra"]:          ## Define arrays for mass-flux
            initialize_mf_arrays(self)
        initialize_output_and_averaging(self)     ## Prepare output storage
        initialize_diagnostics(self)              ## Diagnostics

    ############################################################################
    def _load_config(self, run_params, cfg_params):
        """Return default configuration."""
        default_config = {
            "nz": 100,  "dt": 120.0, "h0": 1000.0, "thetas": 6.5,
            "hc": 400.0, "gridType": "croco_new", "T0": 2.0,
            "N0": 1.962e-6, "mld_ini_temp": -0.0, "mld_iniS": -0.0,
            "cpoce": 3985.0, "lin_eos": True, "rho0": 1027.0,
            "alpha": 2e-4, "beta": 8e-4, "SaltCst": 35.0,
            "Tref": 2.0, "Sref": 35.0, "lat0": 0.0,
            "sustr": 0.0, "svstr": 0.0, "stflx": 0.0,
            "srflx": 0.0, "ssflx": 0.0, "diurnal_cycle": False,
            "nbhours": 72, "outfreq": 1.0, "output_filename": "scm_output.nc",
            "eddy_diff": True, "evd": False,
            "eddy_diff_scheme": "TKE", "mass_flux_tra": False,
            "mass_flux_dyn": False, "mass_flux_tke": False,
            "mass_flux_tke_trplCorr": False, "mass_flux_small_ap": True,
            "extrap_ak_surf": True, "eddy_diff_tke_const": "MNH",
            "tke_sfc_dirichlet": False, "akvmin": 1.3e-06, "aktmin": 1.4e-07,
            "mxlmin": 1.0, "Cent": 0.99, "Cdet": 1.99, "wp_a": 1.0,
            "wp_b": 1.0, "wp_bp": 0.75, "up_c": 0.5,
            "vp_c": 0.5, "bc_ap": 0.2, "delta_bkg": 1.125,
            "wp0":-0.5e-08, "entr_scheme": "R10", "bc_P09": False,
            "beta_bc_P09": 0.3, "trad_coriolis_mod": False,
            "write_netcdf": True, "avg_last_hour": False, "btflx": "no flux",
        }
        config = default_config.copy()
        config.update(run_params)
        config.update(cfg_params)
        return config

    ############################################################################
    def run_direct(self):
        kout=0
        if self.write_netcdf:
            output_init(self) #to initialize .nc dataset
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
            time = self.dt*float(kt+1)                 ## time in seconds at the end of the current time-step
            self._compute_solar_radiation(time)        ## Compute solar radiation (possibly with diurnal cycle)
            stflx0[self.itemp] = self.stflx[self.itemp] - MF_sfc_flux  ## Correct ED sfc flux from the MF contribution

            self._advance_tracers(swr_frac, stflx0)   ## advance tracers to n+1,star (vertical diffusion only)
            self._advance_dynamics()                  ## advance dynamics to n+1,star (Coriolis + vertical viscosity only)

            if self.MF_tra or self.MF_dyn:
                self._compute_mass_flux()

            if self.MF_tra:
                self._advance_tracers_MF()
                if self.bc_P09 != 'inconsistent': MF_sfc_flux = self.Fmass[-1]*(self.tp[-1,self.itemp]-self.t_np1_star[-1,self.itemp]) # MF flux at sfc
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
            if time >= self.start_avg and self.do_avg:
              self._do_turb_fluxes (  )
              self._store_for_average (  )
            #==========================================================
            if  time % self.outfreq == 0:
              self._store_in_history(kout);kout+=1
              if self.write_netcdf:
                if self.eddy_diff: self._do_diags_energy(  )  ## compute diagnostics
                self._do_turb_fluxes (  )                    ## compute turbulent fluxes to be stored in output file
                self.hmxl = scm_oce_new.compute_mxl2(self.bvf,0.01,self.z_r,-0.5,self.rho0,self.nz) # Compute mixed layer depth from bvf
                output_state(self,TimeInSecs=time,kout=kout) ## Write outputs in .nc file
            #===================================================================
            self._swap_arrays() ## swap arrays
            #===================================================================
            if kt<self.nbsteps-1: self.wtke[:]  = 0.                            ## diagnostics : reset the array containing w'e
        if self.do_avg: self._compute_average ( )

    ############################################################################
    def _compute_solar_radiation(self,time):
        modulation = 1.
        if self.DC: modulation = max(np.cos(2.*np.pi*(time/86400. - 0.5)),0.)
        self.srflx = modulation*self.Qswmax

    ############################################################################
    def _advance_tracers(self, swr_frac, stflx0):
        """Advance tracers to the next time step."""
        # Implement tracer advancement logic here
        self.t_np1_star, self.tFlx = scm_oce_new.advance_tra_ed(
                                        self.t_n, stflx0, self.srflx,
                                        swr_frac, self.btflx, self.Hz   , self.akt  ,
                                        self.z_r, self.eps_n, self.eos_params[1],
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
    def _store_for_average(self):
        """Store in _avg arrays"""
        self.u_avg[:]   += self.u_np1[:]; self.v_avg[:]   += self.v_np1[:]
        self.t_avg[:]   += self.t_np1[:]; self.tke_avg[:] += self.tke_np1[:]
        self.wt_avg[:]  +=  - (self.wted[:]+self.wtmf[:])
        self.akv_avg[:] += self.akv[:]; self.akt_avg[:] += self.akt[:]
        self.navg       += 1.

    ############################################################################
    def _compute_average(self):
        """Compute the average at the end of the simulation"""
        cff = 1./self.navg
        self.u_avg[:]   = cff*self.u_avg[:]
        self.v_avg[:]   = cff*self.v_avg[:]
        if self.eddy_diff: self.tke_avg[:] = cff*self.tke_avg[:]
        self.wt_avg[:]  = cff*self.wt_avg[:]
        self.t_avg[:]   = cff*self.t_avg[:]
        self.akv_avg[:] = cff*self.akv_avg[:]
        self.akt_avg[:] = cff*self.akt_avg[:]

    ############################################################################
    def _store_in_history(self,kout):
        self.t_history[:,kout] = self.t_np1[:,self.itemp]
        self.u_history[:,kout] = self.u_np1[:]
        self.v_history[:,kout] = self.v_np1[:]
        self.tke_history[:,kout] = self.tke_np1[:]
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

    ############################################################################
    def _do_diags_energy(self):
        """ Energy diagnostics """
        vint_Sh2_TKE   = 0.; vint_Bprod_TKE = 0.; vint_eps_TKE    = 0.
        self.vint_TKE  = 0.; vint_Bprod_Epot= 0.; vint_Sh2_KE     = 0.
        kk = 1
        while kk < self.nz:
          vint_Sh2_TKE    = vint_Sh2_TKE  +(self.z_r[kk]-self.z_r[kk-1])*(self.shear[kk]+self.shearMF[kk])
          dudz            = 0.5*(self.u_n[kk]+self.u_np1[kk])-0.5*(self.u_n[kk-1]+self.u_np1[kk-1])
          vint_Sh2_KE     = vint_Sh2_KE  - dudz*self.uFlx[kk]
          vint_Bprod_TKE  = vint_Bprod_TKE+(self.z_r[kk]-self.z_r[kk-1])*(self.Bprod[kk]+self.buoyMF[kk])
          vint_Bprod_Epot = vint_Bprod_Epot+(self.z_r[kk]-self.z_r[kk-1])*self.alpha*self.g*self.tFlx[kk] #Buoyancy flux
          vint_eps_TKE    = vint_eps_TKE  +(self.z_r[kk]-self.z_r[kk-1])*self.eps_np1[kk]
          self.vint_TKE   = self.vint_TKE + (self.z_r[kk]-self.z_r[kk-1])*(self.tke_np1[kk]-self.tke_n[kk])/self.dt
          kk+=1
        #================================================
        self.vint_Ekin = 0.; self.vint_Epot = 0.; vint_eps_Epot = 0.
        for k in range(self.nz):
          self.vint_Epot = self.vint_Epot + self.Hz[k]*(1.-(self.alpha*self.g/self.cp)*self.z_r[k]) * (
                                                      self.t_np1[k,self.itemp]-self.t_n[k,self.itemp] )/self.dt
          self.vint_Ekin = self.vint_Ekin + 0.5*self.Hz[k] * (
                                     self.u_np1[k]*self.u_np1[k]
                                  -  self.u_n  [k]*self.u_n  [k]
                                  +  self.v_np1[k]*self.v_np1[k]
                                  -  self.v_n  [k]*self.v_n  [k] ) / self.dt   # m3/s3
          vint_eps_Epot = vint_eps_Epot + 0.5*self.Hz[k]*(self.eps_n[k+1]+self.eps_n[k])
        sfc_KE   = 0.5*( (self.u_np1[-1]+self.u_n[-1])*self.uFlx[-1]
                       + (self.v_np1[-1]+self.v_n[-1])*self.vFlx[-1] )
#
        sfc_Epot = (1.+0.5*(self.alpha*self.g/self.cp)*self.Hz[-1])*self.stflx[self.itemp]
        #
        self.vint_Epot = self.cp*( self.vint_Epot - sfc_Epot ) - vint_eps_Epot
        self.vint_Ekin = self.vint_Ekin - sfc_KE
        self.vint_TKE  = self.vint_TKE  - (-vint_eps_TKE-self.wtke[-1]+self.wtke[0]+self.residual)
        self.vint_Etot = self.vint_Epot + self.vint_TKE + self.vint_Ekin
#

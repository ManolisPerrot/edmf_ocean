#!/usr/bin/env python
# coding: utf-8
__docformat__ = 'reStructuredText'
###########################################
## Imports
###########################################
from sys import exit
import numpy as np
import matplotlib.pyplot as plt
from scmoce import scm_oce
from scmtke import scm_tke
from scmmfc import scm_mfc
#from scipy.io import netcdf
from netCDF4 import Dataset
###########################################
# Classes
###########################################
class SCM:
    """Single Column model class"""

    def __init__(self,  nz     = 100 ,   dt =  20. , h0    =  1000.,
                        thetas = 6.5 ,   hc = 400. , T0    =     2.,
                        N0 =1.9620001275490499e-6  , Tcoef = 0.2048,
                        SaltCst = 35., lat0 =  45.  , sustr =     0.,
                        svstr   =  0., stflx=-500. , srflx =     0.,
                        ssflx   =  0., nbhours = 72, outfreq  =  1.,
                        output_filename="scm_output.nc" , eddy_diff = True ,
                        evd = False  , mass_flux_tra          = False  ,
                                       mass_flux_dyn          = False  ,
                                       mass_flux_tke          = False  ,
                                       mass_flux_tke_trplCorr = False  ,
                                       mass_flux_small_ap     = True   ,
                        lin_eos           = True , extrap_ak_surf      = True ,
                        tke_sfc_dirichlet = True , eddy_diff_tke_const = 'NEMO',
                        Cent  = 0.55  , Cdet = -1 , wp_a =  1 , wp_b  = 1   ,
                        wp_bp = 0.0002, up_c = 0.5, vp_c = 0.5, bc_ap = 0.1 ,
                        delta_bkg = 0., wpmin=1.e-08,  entr_scheme = 'P09'  ):
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
################################################################################################################
        ## simulation parameters
        #-----------------------------------
        ## define Coriolis frequency
        omega          = 7.292116e-05; rad            = np.pi/180.
        self.fcor      = 2.*omega*np.sin(rad*lat0)
        self.ecor      = 2.*omega*np.cos(rad*lat0)
        self.nz        = nz                                                     ## number of grid points
        self.ntra      = 2; self.itemp = 0; self.isalt = 1; self.itke = 2       ## number of tracers
        self.ntraMF    = 3
        self.dt        = dt                                                     ## time-step of the simulation [seconds]
        self.nbsteps   = int( (nbhours*3600.)/self.dt )                         ## number of time steps
        self.outfreq   = int( outfreq*3600. )                                   ## frequency of outputs in seconds
        # should check if self.outfreq is a multiple of the time-step
        self.output    = output_filename
        # physical constants
        self.rho0      = 1027; self.cp        = 3985.; self.g         = 9.81
        self.lineos    = True                                                   ## no other option implemented yet
        ####################################
        # surface boundary conditions
        #-----------------------------------
        self.taux      = sustr                                                  ## zonal wind stress       [m2/s2]
        self.tauy      = svstr                                                  ## meridional wind stress  [m2/s2]
        self.stflx     = np.zeros(2)
        cff            = 1./(self.rho0*self.cp)
        self.stflx[self.itemp]  = stflx*cff                                     ## non-penetrative heat flux  [K m / s]
        self.srflx     = srflx*cff                                              ## solar radiation            [K m / s]
        self.stflx[self.isalt]  = ssflx                                         ## freshwater flux         [psu m / s]
        ####################################
        # Eddy-diffusion parameters
        #-----------------------------------
        self.tkemin      = 1.e-8; self.akvmin    = 1.e-4; self.aktmin    = 1.e-5
        self.ED          = eddy_diff     ; self.ED_evd          = evd
        #
        self.ED_tke_sfc_dirichlet = tke_sfc_dirichlet
        self.ED_extrap_sfc        = extrap_ak_surf
        #if self.ED_tke_sfc_dirichlet: self.ED_extrap_sfc    = False
        self.ED_tke_const    = 0  ; self.inv_schmidt     = 1.; cm = 0.1
        if eddy_diff_tke_const=='MNH':
          self.ED_tke_const    = 1; cm = 0.126; self.inv_schmidt = 0.34/cm
        if eddy_diff_tke_const=='RS81':
          self.ED_tke_const    = 2; cm = 0.0667; self.inv_schmidt= 0.4/cm
        ####################################
        # Mass flux options and parameters
        #-----------------------------------
        self.MF_tra      = mass_flux_tra  ; self.MF_dyn          = mass_flux_dyn
        self.MF_tke      = mass_flux_tke  ; self.MF_tke_trplCorr = mass_flux_tke_trplCorr
        self.mass_flux_entr = entr_scheme ; self.MF_small_ap     = mass_flux_small_ap
        self.mf_params  = np.array([Cent,Cdet,wp_a,wp_b,wp_bp,up_c,vp_c,bc_ap,delta_bkg,wpmin])
        #self.mf_params  = np.array([Cent,2.*Cent,wp_a,wp_b,wp_bp,up_c,vp_c,bc_ap,0.5*wp_bp])
        ####################################
        # define vertical grid
        #-----------------------------------
        Sc_r  = np.arange(-nz+0.5, 0.5, 1) / float(nz)
        Sc_w  = np.arange(-nz, 1, 1) / float(nz)
        Cs_r = (1.-np.cosh(thetas*Sc_r))/(np.cosh(thetas)-1.)
        Cs_w = (1.-np.cosh(thetas*Sc_w))/(np.cosh(thetas)-1.)
        self.z_w   = (hc*Sc_w+Cs_w*h0)*h0/(h0+hc)
        self.z_r   = (hc*Sc_r+Cs_r*h0)*h0/(h0+hc)
        self.Hz    = self.z_w[1:]-self.z_w[:-1]
        ####################################
        # define the initial state
        #-----------------------------------
        self.u_n      = np.zeros(self.nz); self.v_n      = np.zeros(self.nz)
        self.u_np1    = np.zeros(self.nz); self.v_np1    = np.zeros(self.nz)
        self.t_n      = np.zeros((self.nz,self.ntra), order='F')
        self.t_np1    = np.zeros((self.nz,self.ntra), order='F')
        self.saltini = SaltCst;
        self.t_n[:,self.isalt] = self.saltini                                   # constant salinity
        self.t_np1[:,self.isalt] = self.t_n[:,self.isalt]                       # constant salinity
        #
        self.T0      = T0;  self.N0      = N0; self.Tcoef = Tcoef
        self.alpha   = self.Tcoef/self.rho0
        self.beta    = 0.
        strat        = N0 / (self.g * self.alpha )
        self.t_n  [:,self.itemp] = T0 + 0.5*strat*( self.z_w[1:]+self.z_w[:-1] )
        self.t_np1[:,self.itemp] = T0 + 0.5*strat*( self.z_w[1:]+self.z_w[:-1] )
        ####################################
        # initialize arrays for EDDY-DIFFUSION (tke, bvf, lup,ldwn,tke,AkT,Akv)
        #-----------------------------------
        self.mxlmin = (self.akvmin / cm) / np.sqrt(self.tkemin)
        self.tke_n   = np.zeros(self.nz+1); self.tke_n[:] = self.tkemin
        self.tke_env = np.zeros(self.nz+1); self.tke_env[:] = self.tkemin
        self.tke_np1 = np.zeros(self.nz+1); self.tke_np1[:] = self.tkemin
        self.lupw  = np.zeros(self.nz+1); self.lupw[:] = self.mxlmin
        self.ldwn  = np.zeros(self.nz+1); self.ldwn[:] = self.mxlmin
        self.bvf   = np.zeros(self.nz+1); self.Prdtl = np.zeros(self.nz+1)
        self.akt   = np.zeros(self.nz+1); self.akt[:] = self.aktmin
        self.akv   = np.zeros(self.nz+1); self.akv[:] = self.akvmin
        self.eps   = np.zeros(self.nz+1)  # TKE dissipation
        self.shear = np.zeros(self.nz+1)
        ####################################
        # initialize arrays for MASS-FLUX
        #-----------------------------------
        self.rhoc   = 0.01    # density criterion for mixed layer depth (consistent with NEMO)
        self.zinv   = np.array(-10.)
        self.hmxl   = np.array(-10.)
        self.ap     = np.zeros(self.nz+1); self.wp     = np.zeros(self.nz+1)
        self.tp     = np.zeros((self.nz+1,self.ntraMF), order='F');
        self.Bp     = np.zeros(self.nz+1); self.Fmass  = np.zeros(self.nz+1)
        self.ent    = np.zeros(self.nz  ); self.det    = np.zeros(self.nz  )
        self.up     = np.zeros(self.nz+1); self.vp     = np.zeros(self.nz+1)
        self.epsPlume = np.zeros(self.nz)  # TKE plume dissipation
        ####################################
        # store MASS-FLUX solution (temperature, velocity)
        # in order to avoid writing a netcdf file
        #-----------------------------------
        dimt = nbhours // outfreq 
        self.t_history = np.zeros((self.nz,dimt), order='F')
        self.u_history = np.zeros((self.nz,dimt), order='F')
        self.v_history = np.zeros((self.nz,dimt), order='F')
        ####################################
        # edmf_diags
        ####################################
        # self.wted   = np.zeros(self.nz+1); self.wtmf   = np.zeros(self.nz+1)
        # self.wued   = np.zeros(self.nz+1); self.wumf   = np.zeros(self.nz+1)
        self.buoyMF = np.zeros(self.nz+1); self.shearMF = np.zeros(self.nz+1)
        self.wtke   = np.zeros(self.nz  )
        self.triple_corr     = np.zeros(self.nz+1);
        ####################################
        # Energy diagnostics
        #-----------------------------------
        # self.vint_Etot = np.array(1e-14); self.vint_Ekin = np.array(1e-14)
        # self.vint_Epot = np.array(1e-14); self.vint_TKE  = np.array(1e-14)
        # self.vint_Eps  = np.array(1e-14)
################################################################################################################



# store MASS-FLUX solution (temperature)
#-----------------------------------
#
    def run_direct(self):
        #=============================
        #self.output_init( ); #to initialize .nc dataset 
        kout = 0
        # save initial temperature and velocity
        # self.t_history[:,kout] = self.t_n[:,self.itemp]
        # self.u_history[:,kout] = self.u_n[:]
        # self.v_history[:,kout] = self.v_n[:]
        swr_frac = scm_oce.lmd_swfrac(self.Hz,self.nz)   ## Compute fraction of solar shortwave flux penetrating to specified depth
        if self.ED: self.do_ED( )                        ## Initialize eddy-diffusion scheme to compute eddy diffusivity/eddy viscosity for the first time-step
        ####################################
        # Time LOOP
        #-----------------------------------
        for kt in range(self.nbsteps):
            # update current time
            time = self.dt*float(kt+1)                   ## time in seconds at the end of the current time-step
            nn = 1 +  kt    % 2                          ## used to alternate the direction of integration for Coriolis term
            self.wtke[:] = 0.                            ## diagnostics : reset the array containing w'e
            #===================================================
            # Advance tracers to n+1 (vertical diffusion only)
            #===================================================
            self.t_np1  = scm_oce.advance_tra_ed(
                                        self.t_n, self.stflx, self.srflx,
                                        swr_frac, self.Hz   , self.akt  ,
                                        self.z_w, self.eps  , self.alpha,
                                        self.dt , self.nz   , self.ntra  )
            #==============================================================
            # advance dynamics to n+1 (Coriolis + vertical viscosity only)
            #==============================================================
            self.u_np1, self.v_np1 = scm_oce.advance_dyn_cor_ed(
                                        self.u_n, self.v_n, self.taux, self.tauy,
                                        self.Hz , self.akv, self.fcor, self.dt  ,
                                        nn, self.nz  )
            #===================
            # Compute mass flux
            #===================
            if self.MF_tra or self.MF_dyn:
                self.do_MF( )
                if self.MF_small_ap:
                  self.Fmass = -(self.ap[:]*self.wp[:])
                else:
                  self.Fmass = -(self.ap[:]*self.wp[:])/(1.-self.ap[:])
            # apply the MF term to the tracer equation
            if self.MF_tra: scm_oce.advance_tra_mf(
                                      self.t_np1, self.tp, self.Fmass, self.Hz,
                                      self.dt, self.nz, self.ntra )
            # apply the MF term to the velocity equation + compute the transfer of KE to subgrid scales stored in shearMF
            if self.MF_dyn: scm_oce.advance_dyn_mf(
                                      self.u_np1, self.v_np1, self.shearMF,
                                      self.u_n  , self.v_n  , self.up , self.vp ,
                                      self.Fmass, self.Hz   , self.dt , self.nz  )
            #==========================================================
            # Compute eddy-viscosity / eddy-diffusivity  (TKE scheme)
            #==========================================================
            if self.ED:
                self.do_ED( )
            #==========================================================
            # Check for outputs & diagnostics
            #==========================================================
            # Turbulent fluxes diagnostics
            # self.wted[  self.nz] = self.stflx[self.itemp]                   ## (w'theta')_ED
            # self.wted[1:self.nz] = self.akt[1:self.nz] * (  self.t_np1[1:self.nz  ,self.itemp]
            #                         - self.t_np1[0:self.nz-1,self.itemp] ) / self.Hz[0:self.nz-1]
            # self.wued[  self.nz] = self.taux                 ## (w'theta')_ED
            # self.wued[1:self.nz] = self.akv[1:self.nz] * (  self.u_np1[1:self.nz]
            #                         - self.u_np1[0:self.nz-1] ) / self.Hz[0:self.nz-1]
            # if self.MF_tra:                                                 ## (w'theta')_MF
            #   self.wtmf[1:self.nz] =  self.Fmass[1:self.nz]*(
            #                               self.tp   [1:self.nz,  self.itemp]
            #                             - self.t_np1[0:self.nz-1,self.itemp] )
            # if self.MF_dyn:                                                 ## (w'u')_MF
            #   self.wumf[1:self.nz] =  self.Fmass[1:self.nz]*(
            #                               self.up   [1:self.nz]
            #                             - self.u_np1[0:self.nz-1] )
            if  time % self.outfreq == 0:
                kout+=1
                #=======================================
                #
                # #Compute b_n and b_np1 for diagnostics purposes
                # if self.lineos:
                #   rho_n,bvf = scm_oce.rho_eos_lin(
                #                  self.t_n[:,self.itemp], self.t_n[:,self.isalt],
                #                  self.z_r , self.alpha , self.beta , self.nz  )
                #   rho_np1,bvf = scm_oce.rho_eos_lin(
                #                  self.t_np1[:,self.itemp], self.t_np1[:,self.isalt],
                #                  self.z_r , self.alpha , self.beta , self.nz  )
                # else:
                #   rho_n,bvf = scm_oce.rho_eos(
                #                  self.t_n[:,self.itemp], self.t_n[:,self.isalt],
                #                  self.z_r, self.z_w, self.nz )
                #   rho_np1,bvf = scm_oce.rho_eos(
                #                  self.t_np1[:,self.itemp], self.t_np1[:,self.isalt],
                #                  self.z_r, self.z_w, self.nz )
                # # Compute mixed layer depth from bvf
                # self.hmxl  = scm_oce.compute_mxl2(bvf,self.rhoc,self.z_r,10.,self.nz)
                # # Compute buoyancy
                # b_n   = -(self.g/self.rho0)*rho_n[:]
                # b_np1 = -(self.g/self.rho0)*rho_np1[:]
                # Energy diagnostics
                # self.vint_Ekin = 0.
                # cor_KE         = 0.
                # for k in range(self.nz):
                #     self.vint_Ekin = self.vint_Ekin + 0.5*self.Hz[k] * (
                #                            self.u_np1[k]*self.u_np1[k]
                #                           -self.u_n  [k]*self.u_n  [k]
                #                           +self.v_np1[k]*self.v_np1[k]
                #                           -self.v_n  [k]*self.v_n  [k] ) / self.dt
                #     cor_KE         = cor_KE + 0.5*self.Hz[k]*self.fcor * (
                #                            self.u_np1[k]*self.v_np1[k]
                #                           -self.u_n  [k]*self.v_n  [k] )
                # sfc_KE = 0.5*( (self.u_np1[-1]+self.u_n[-1])*self.taux
                #              + (self.v_np1[-1]+self.v_n[-1])*self.tauy )
                # #
                # self.vint_Ekin = self.vint_Ekin - cor_KE - sfc_KE
                # self.vint_Eps  = sum( (self.z_r[1:self.nz]
                #                       -self.z_r[0:self.nz-1] ) * self.eps[1:self.nz] )
                # self.vint_Epot = sum( self.Hz[:]*(-self.z_r[:]) *
                #                                       (b_np1[:]-b_n[:])/self.dt
                #                     + self.Hz[:]*self.cp *
                #                                      (self.t_np1[:,self.itemp]-self.t_n[:,self.itemp])/self.dt )
                # sfc_Epot       = -self.cp*self.stflx[self.itemp]
                # self.vint_Epot = self.vint_Epot + sfc_Epot
                # self.vint_Etot = self.vint_Epot + self.vint_TKE + self.vint_Ekin
                #=======================================
                # Write outputs in .nc file
                #self.output_state(TimeInSecs=time,kout=kout)
            # save advanced temperature and velocity
            self.t_history[:,kout-1] = self.t_np1[:,self.itemp]
            self.u_history[:,kout-1] = self.u_np1[:]
            self.v_history[:,kout-1] = self.v_np1[:]
            # swap arrays
            self.u_n[:] = self.u_np1[:]; self.v_n[:] = self.v_np1[:]; self.t_n[:,:] = self.t_np1[:,:]

#












#
    def do_ED(self):
        #==========================================================================
        # Compute Brunt-Vaisala frequency bvf for TKE production/destruction term
        #==========================================================================
        if self.lineos:
          rho,self.bvf = scm_oce.rho_eos_lin(
                                 self.t_np1[:,self.itemp], self.t_np1[:,self.isalt],
                                 self.z_r, self.alpha, self.beta, self.nz )
          #
          #Bsfc         = self.g*( self.alpha * (self.stflx[self.itemp]+self.srflx)
          #                        - self.beta  *  self.stflx[self.isalt]           )
        else:
          rho,self.bvf = scm_oce.rho_eos(self.t_np1[:,self.itemp],self.t_np1[:,self.isalt],self.z_r,self.z_w,self.nz)
        #  Bsfc,self.alpha,self.beta = scm_oce.compute_bsfc(self.t_np1[-1,self.itemp],self.t_np1[-1,self.isalt],self.stflx,self.srflx,self.ntra)
        #=======================================
        # Compute boundary conditions for TKE
        tke_sfc,tke_bot, flux_sfc = scm_tke.compute_tke_bdy(
                                      self.taux,  self.tauy, self.ED_tke_const  )
        #=======================================
        # Compute TKE production by shear
        self.shear = scm_tke.compute_shear(
                                  self.u_n  , self.v_n  ,
                                  self.u_np1, self.v_np1,
                                  self.akv  , self.z_r  , self.nz  )
        #================================================
        # Transfer from mean PE/KE to TKE by mass flux term
        if self.MF_tke and self.MF_tra: self.buoyMF[:] = -self.Fmass[:]*self.Bp[:]  ## corresponds to  - FM x Bp
        if not self.MF_tke : self.shearMF[:] = 0.                                   ## If not needed, cancel the transfer from mean KE to TKE by mass flux on momentum
        #================================================
        # Triple correlation term
        if self.MF_tke_trplCorr:
          self.triple_corr, self.tke_env = scm_mfc.compute_triplecorr(
                                        self.tke_n, self.tp[:,self.isalt+1], self.Fmass,
                                        self.up   , self.vp   ,self.wp  ,
                                        self.u_np1, self.v_np1, self.Hz ,
                                        self.z_r  , self.wtke, self.nz  )
        #=========================
        # Advance TKE
        self.tke_np1, self.Prdtl, self.eps, residual = scm_tke.advance_tke(
                                         self.tke_n, self.lupw, self.ldwn,
                                         self.akv  , self.akt , self.Hz  ,
                                         self.z_r  , self.bvf , self.buoyMF,
                                         self.shear, self.shearMF, self.triple_corr,
                                         self.wtke, self.dt, tke_sfc, tke_bot, flux_sfc,
                                         self.ED_tke_sfc_dirichlet, self.ED_tke_const, self.nz  )
        #===================================================
        # Finalize eddy-viscosity/diffusivity computation
        self.lupw,self.ldwn = scm_tke.compute_mxl(
                                         self.tke_np1, self.bvf , self.Hz,
                                         self.taux   , self.tauy, self.nz )
        self.akv,self.akt   = scm_tke.compute_ed (
                                         self.tke_np1, self.lupw, self.ldwn,
                                         self.Prdtl  , self.ED_extrap_sfc,
                                         self.ED_tke_const, self.nz  )
        #=========================
        # Apply EVD if necessary
        if self.ED_evd: scm_oce.compute_evd(
                                          self.bvf, self.akv, self.akt,
                                          10., self.nz  )
        # Vertically integrated energy budgets
        self.vint_TKE = sum( (self.z_r[1:self.nz]-self.z_r[0:self.nz-1])
                        *(self.tke_np1[1:self.nz]-self.tke_n[1:self.nz])/self.dt )
        #surface boundary contribution
        self.vint_TKE = self.vint_TKE - 0.5*self.inv_schmidt*(self.akv[-1]+self.akv[-2])*(self.tke_np1[-1]-self.tke_np1[-2])/self.Hz[-1]
        self.vint_TKE = self.vint_TKE - residual # remove the component artificially added
        self.tke_n[:] = self.tke_np1[:]
        #












#
    def do_MF(self):
        #=======================================================
        # Surface boundary conditions for wp,up,vp and tracers
        wp0,up0,vp0,tp0 = scm_mfc.compute_mf_bdy(
                                      self.u_np1[-2:]  , self.v_np1[-2:],
                                      self.t_np1[-2:,:], self.tke_n[-2:],
                                      self.Hz[-2:]     , self.ntra, 2 )
        #=================================================================
        # Compute the mean quantities used to constrain the mass flux eqns
        u_mean,v_mean,t_mean,dtke_m = scm_mfc.compute_mf_forcing(
                                      self.u_np1, self.v_np1, self.t_np1,
                                      self.tke_n, self.nz   , self.ntra )
        #=======================================================================================
        # Compute mass flux variables (Rio et al. 2010 closure or Pergaud et al. 2009 closure)
        if self.mass_flux_entr=='P09':
          self.ap,self.up,self.vp,self.wp,self.tp,self.Bp,self.ent,self.det  =  scm_mfc.mass_flux_p09(
                                    u_mean, v_mean, t_mean, self.z_w, self.Hz       ,
                                    tp0   , up0   , vp0   , wp0     , self.mf_params,
                                    self.alpha , self.beta, self.MF_small_ap,self.zinv ,
                                    self.nz , self.ntraMF , len(self.mf_params)  )
        if self.mass_flux_entr=='R10':
          self.ap,self.up,self.vp,self.wp,self.tp,self.Bp,self.ent,self.det, self.epsPlume = scm_mfc.mass_flux_r10(
                                    u_mean, v_mean, t_mean, self.tke_n, self.z_w, self.Hz,
                                    tp0   , up0   , vp0   , wp0     , self.mf_params,
                                    self.alpha , self.beta, self.MF_small_ap, self.lineos, self.zinv ,
                                    self.nz , self.ntraMF , len(self.mf_params)   )
#



#
    def output_init(self):
        fh01 = Dataset(self.output, mode='w',format="NETCDF4")
        fh01.createDimension('z_r',self.nz)
        fh01.createDimension('z_w',self.nz+1)
        fh01.createDimension('time',None)
        fh01.fcorSIN         = self.fcor
        fh01.fcorCOS         = self.ecor
        fh01.eddy_diff       = str(self.ED)
        fh01.evd             = str(self.ED_evd)
        fh01.mass_flux_tra   = str(self.MF_tra)
        fh01.mass_flux_dyn   = str(self.MF_dyn)
        fh01.mass_flux_tke   = str(self.MF_tke)
        fh01.mass_flux_tke_trplCorr = str(self.MF_tke_trplCorr)
        ocean_time = fh01.createVariable('ocean_time','f8',('time')); ocean_time[:] = 0.
        ocean_time.units = 'seconds'
        ocean_time.long_name = 'time since initialization'
        taux = fh01.createVariable('taux','f8',('time')); taux[:] = self.taux
        tauy = fh01.createVariable('tauy','f8',('time')); tauy[:] = self.tauy
        Qns  = fh01.createVariable('Qns','f8',('time'));   Qns[:] = self.stflx[self.itemp]*self.rho0*self.cp           #non-solar
        Qs   = fh01.createVariable( 'Qs','f8',('time'));    Qs[:] = self.srflx*self.rho0*self.cp                          #solar
        Fw   = fh01.createVariable( 'Fw','f8',('time'));    Fw[:] = self.stflx[self.isalt]            #freshwater
        zw   = fh01.createVariable('z_w','f8',('z_w')); zw[:] = self.z_w[:]
        zr   = fh01.createVariable('z_r','f8',('z_r')); zr[:] = self.z_r[:]
        var  = fh01.createVariable('u','f8',('time','z_r')); var[0,:] = self.u_n[:]; var.units = 'm s-1'; del var
        var  = fh01.createVariable('v','f8',('time','z_r')); var[0,:] = self.v_n[:]; var.units = 'm s-1'; del var
        var  = fh01.createVariable('temp','f8',('time','z_r')); var[0,:] = self.t_n[:,self.itemp]; var.units = 'Celsius'; del var
        var  = fh01.createVariable('salt','f8',('time','z_r')); var[0,:] = self.t_n[:,self.isalt]; var.units = 'psu'; del var
        var  = fh01.createVariable('Etot','f8',('time')); var[0] = self.vint_Etot
        var  = fh01.createVariable('Ekin','f8',('time')); var[0] = self.vint_Ekin
        var  = fh01.createVariable('Epot','f8',('time')); var[0] = self.vint_Epot
        var  = fh01.createVariable('Etke','f8',('time')); var[0] = self.vint_TKE
        var  = fh01.createVariable('Eeps','f8',('time')); var[0] = self.vint_Eps
        var  = fh01.createVariable('hmxl','f8',('time')); var[0] = self.hmxl
        var  = fh01.createVariable('WT','f8',('time','z_w')); var[0,:] = self.wted[:]+self.wtmf[:]; var.units = 's-2'; del var
        if self.ED:
            var  = fh01.createVariable('tke','f8',('time','z_w')); var[0,:] = self.tke_n[:]; var.units = 'm2 s-2'; del var
            var  = fh01.createVariable('Akv','f8',('time','z_w')); var[0,:] = self.akv[:]; var.units = 'm2 s-1'; del var
            var  = fh01.createVariable('Akt','f8',('time','z_w')); var[0,:] = self.akt[:]; var.units = 'm2 s-1'; del var
            var  = fh01.createVariable('bvf','f8',('time','z_w')); var[0,:] = self.bvf[:]; var.units = 's-2'; del var
            var  = fh01.createVariable('lup','f8',('time','z_w')); var[0,:] = self.lupw[:]; var.units = 'm'; del var
            var  = fh01.createVariable('ldw','f8',('time','z_w')); var[0,:] = self.ldwn[:]; var.units = 'm'; del var
            var  = fh01.createVariable('WT_ED','f8',('time','z_w')); var[0,:] = self.wted[:]; var.units = 's-2'; del var
            var  = fh01.createVariable('Prdtl','f8',('time','z_w')); var[0,:] = self.Prdtl[:]; var.units = 'none'; del var
        if self.MF_tra:
            var = fh01.createVariable('a_p','f8',('time','z_w')); var[0,:] = self.ap[:]; del var
            var = fh01.createVariable('zinv','f8',('time')); var[0] = self.zinv
            var = fh01.createVariable('w_p','f8',('time','z_w')); var[0,:] = self.wp[:]; del var
            var = fh01.createVariable('B_p','f8',('time','z_w')); var[0,:] = self.buoyMF[:]; del var
            var = fh01.createVariable('temp_p','f8',('time','z_w')); var[0,:] = self.tp[:,self.itemp]; del var
            var = fh01.createVariable('salt_p','f8',('time','z_w')); var[0,:] = self.tp[:,self.isalt]; del var
            var = fh01.createVariable('Ent','f8',('time','z_r')); var[0,:] = self.ent[:]; del var
            var = fh01.createVariable('Det','f8',('time','z_r')); var[0,:] = self.det[:]; del var
            var  = fh01.createVariable('WT_MF','f8',('time','z_w')); var[0,:] = self.wtmf[:]; del var
            var = fh01.createVariable('u_p','f8',('time','z_w')); var[0,:] = self.up[:]; del var
            var = fh01.createVariable('v_p','f8',('time','z_w')); var[0,:] = self.vp[:]; del var
            if self.MF_tke:
                var = fh01.createVariable('buoyMF','f8',('time','z_w')); var[0,:] = self.buoyMF[:]; del var
                var = fh01.createVariable('shearMF','f8',('time','z_w')); var[0,:] = self.shearMF[:]; del var
            if self.MF_tke_trplCorr:
                var = fh01.createVariable('we','f8',('time','z_w')); var[0,:] = self.triple_corr[:]; del var
                var = fh01.createVariable('tke_p','f8',('time','z_w')); var[0,:] = self.tp[:,self.isalt+1]; del var
        fh01.close()
#












#
    def output_state(self,TimeInSecs,kout):
        fh01 = Dataset(self.output, mode='a',format="NETCDF4")
        fh01.variables['ocean_time'][kout] = TimeInSecs
        fh01.variables['taux'][kout]       = self.taux
        fh01.variables['tauy'][kout]       = self.tauy
        fh01.variables['Qns'][kout]        = self.stflx[self.itemp]*self.cp*self.rho0
        fh01.variables['Qs'][kout]         = self.srflx
        fh01.variables['Fw'][kout]         = self.stflx[self.isalt]*self.cp*self.rho0
        fh01.variables['u'][kout,:]        = self.u_n[:]
        fh01.variables['v'][kout,:]        = self.v_n[:]
        fh01.variables['temp'][kout,:]     = self.t_n[:,self.itemp]
        fh01.variables['salt'][kout,:]     = self.t_n[:,self.isalt]
        fh01.variables['Etot'][kout]       = self.vint_Etot
        fh01.variables['Ekin'][kout]       = self.vint_Ekin
        fh01.variables['Epot'][kout]       = self.vint_Epot
        fh01.variables['Etke'][kout]       = self.vint_TKE
        fh01.variables['Eeps'][kout]       = self.vint_Eps
        fh01.variables['hmxl'][kout]       = self.hmxl
        fh01.variables['WT'][kout,:]       = self.wted[:]+self.wtmf[:]
        if self.ED:
            fh01.variables['tke'][kout,:] = self.tke_n[:]
            fh01.variables['Akv'][kout,:] = self.akv[:]
            fh01.variables['Akt'][kout,:] = self.akt[:]
            fh01.variables['bvf'][kout,:] = self.bvf[:]
            fh01.variables['lup'][kout,:] = self.lupw[:]
            fh01.variables['ldw'][kout,:] = self.ldwn[:]
            fh01.variables['WT_ED'][kout,:] = self.wted[:]
            fh01.variables['Prdtl'][kout,:] = self.Prdtl[:]
        if self.MF_tra:
            fh01.variables['zinv'][kout]  = self.zinv
            fh01.variables['a_p'][kout,:] = self.ap[:]
            fh01.variables['w_p'][kout,:] = self.wp[:]
            fh01.variables['B_p'][kout,:] = self.Bp[:]
            fh01.variables['temp_p'][kout,:] = self.tp[:,self.itemp]
            fh01.variables['salt_p'][kout,:] = self.tp[:,self.isalt]
            fh01.variables['Ent'][kout,:] = self.ent[:]
            fh01.variables['Det'][kout,:] = self.det[:]
            fh01.variables['WT_MF'][kout,:] = self.wtmf[:]
            fh01.variables['u_p'][kout,:] = self.up[:]
            fh01.variables['v_p'][kout,:] = self.vp[:]
            if self.MF_tke:
                fh01.variables['buoyMF'][kout,:] = self.buoyMF[:]
                fh01.variables['shearMF'][kout,:] = self.shearMF[:]
            if self.MF_tke_trplCorr:
                fh01.variables['we'][kout,:] = self.triple_corr[:]
                fh01.variables['tke_p'][kout,:] = self.tp[:,self.isalt+1]
        fh01.close()

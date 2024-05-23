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
from scmgls import scm_gls
#from scipy.io import netcdf
from netCDF4 import Dataset
###########################################
# Classes
###########################################
class SCM:
    """Single Column model class"""

    def __init__(self,  nz     = 100                    , dt       =  30.        , h0              =  1000.,
                        thetas = 6.5                    , hc       = 400.        , T0              =     2.,
                        N0     = 1.9620001275490499e-6  , gridType = 'analytic'  ,
                        mld_ini = -0.                   , mld_iniS = -0.,
                        cpoce  = 3985.                  , lin_eos = True        , rho0            =  1027.,
                        Tcoef  = 0.2048                 , Scoef   =  0.         , SaltCst         =    35.,
                        Tref   = 2.                     , Sref    = 35.         , lat0            =    45.,
                        sustr  = 0.                     , svstr   =  0.         , stflx           =  -500.,
                        srflx  = 0.                     , ssflx   =  0.         , diurnal_cycle   = False ,
                        btflx  = 'no flux'              ,
                        nbhours = 72                    , outfreq  =  1.        , output_filename = "scm_output.nc",
                        eddy_diff = True                , evd      = False      , eddy_diff_scheme = 'TKE',
                        mass_flux_tra   = False  ,
                        mass_flux_dyn = False           , mass_flux_tke = False , mass_flux_tke_trplCorr = False,
                        mass_flux_small_ap     = True   , extrap_ak_surf = True , tke_sfc_dirichlet = False,
                        eddy_diff_tke_const = 'NEMO'    , akvmin   = 1.e-4      , aktmin   = 1.e-5         ,
                        mxlmin  = 1.0,
                        Cent  = 0.55         , Cdet = -1         , wp_a =  1        , wp_b  = 1   ,
                        wp_bp = 0.0002       , up_c = 0.5        , vp_c = 0.5       , bc_ap = 0.1 ,
                        delta_bkg = 0.       , wp0=-1.e-08       ,
                        entr_scheme = 'R10'  , write_netcdf=False, avg_last_hour=False ):
        """[summary]
        Args:
            nz: Number of grid points. Defaults to 100.
            dt: Dynamical time-step. Defaults to 30.
            h0: Depth of the water column. Defaults to 1000.
            thetas: stretching parameter toward the surface for the vertical grid. Defaults to 6.5.
            hc: Stretching parameter for the vertical grid. Resolution is almost constant between 0 and hc. Defaults to 400.
            T0: Initial surface temperature. Defaults to 2 (Celsius).
            N0: Initial stratification. Defaults to 1.9620001275490499E-6.
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
        self.write_netcdf = write_netcdf
        self.do_avg    = avg_last_hour
        self.start_avg = (nbhours-1.)*3600.
        # physical constants
        self.rho0        = rho0; self.cp        = cpoce; self.g         = 9.81
        self.lineos      = True
        self.alpha       = Tcoef/self.rho0;      self.beta    = Scoef/self.rho0
        self.T0          = Tref           ;      self.S0      = Sref
        self.eos_params  = np.array([self.rho0,self.alpha,self.beta,self.T0,self.S0])
        ####################################
        # surface boundary conditions
        #-----------------------------------
        self.ustr_sfc           = sustr                      ## zonal wind stress       [m2/s2]
        self.vstr_sfc           = svstr                      ## meridional wind stress  [m2/s2]
        self.stflx              = np.zeros(2)
        cff                     = 1./(self.rho0*self.cp)
        self.stflx[self.itemp]  = stflx*cff                  ## non-penetrative heat flux  [K m / s]
        self.srflx              = srflx*cff                  ## solar radiation            [K m / s]
        self.stflx[self.isalt]  = ssflx                      ## freshwater flux         [psu m / s]
        self.btflx              = np.zeros(2)
        self.ustr_bot           = 0.
        self.vstr_bot           = 0.
        self.z0b                = 1.e-14
        self.DC                 = diurnal_cycle
        self.Qswmax             = srflx*cff
        ####################################
        # Eddy-diffusion parameters
        #-----------------------------------
        self.ED          = eddy_diff     ; self.ED_evd          = evd
        if self.ED:
          self.ED_scheme            = eddy_diff_scheme
          self.ED_tke_sfc_dirichlet = tke_sfc_dirichlet
          self.ED_extrap_sfc        = extrap_ak_surf
          epsmin = 1.e-12
          if self.ED_scheme=='TKE':
            self.ED_tke_const         = 0  ; self.inv_schmidt     = 1.; cm = 0.1
            if eddy_diff_tke_const=='MNH':
              self.ED_tke_const    = 1; cm = 0.126; self.inv_schmidt = 0.34/cm
            if eddy_diff_tke_const=='RS81':
              self.ED_tke_const    = 2; cm = 0.0667; self.inv_schmidt= 0.4/cm
            tkemin = np.power( akvmin/(cm*mxlmin), 2 )
            self.min_Threshold   = np.array([tkemin,akvmin,aktmin,mxlmin])
          elif self.ED_scheme=='Keps':
            self.pnm    = np.zeros(3)
            self.pnm[0] = 3.0; self.pnm[1] = -1.0; self.pnm[2] = 1.5
            self.betaCoef   = np.zeros(4)
            self.betaCoef[0] = 1.44 ; self.betaCoef[1] = 1.92
            #self.betaCoef[2] = -0.4 ; self.betaCoef[3] = 1.00
            self.betaCoef[2] = -1.83 ; self.betaCoef[3] = -1.83
            Sch_k = 1.; Sch_eps = 1.3  # Schmidt numbers
            self.OneOverSig_k   = 1./Sch_k; self.OneOverSig_psi = 1./Sch_eps
            c1 = 5.; c2 = 0.8; c3 = 1.968; c4 = 1.136
            a1  = 0.66666666667 - 0.5*c2; a2  = 1.            - 0.5*c3
            a3  = 1.0           - 0.5*c4; nn  = 0.5*c1
            self.cm0 =  np.power( (a2*a2 - 3.0*a3*a3 + 3.0*a1*nn)/(3.0*nn*nn), 0.25 )
            tkemin = np.power(mxlmin*epsmin/(self.cm0*self.cm0*self.cm0),2./3.)
            self.min_Threshold   = np.array([tkemin,akvmin,aktmin,mxlmin,epsmin])
            self.bdy_tke_sfc = np.zeros(2); self.bdy_tke_bot = np.zeros(2)
            self.bdy_eps_sfc = np.zeros(2); self.bdy_eps_bot = np.zeros(2)
            if  not tke_sfc_dirichlet: self.bdy_tke_sfc[0] = 1.; self.bdy_eps_sfc[0] = 1.
            # we assume Dirichlet condition at the bottom
            #if  bdy_tke_bot == 'Neu': self.bdy_tke_bot[0] = 1.
            #if  bdy_gls_bot == 'Neu': self.bdy_eps_bot[0] = 1.
          else:
            exit("Wrong entry for parameter eddy_diff_scheme (should be TKE or Keps)")
        ####################################
        # Mass flux options and parameters
        #-----------------------------------
        self.MF_tra      = mass_flux_tra  ; self.MF_dyn          = mass_flux_dyn
        self.MF_tke      = mass_flux_tke  ; self.MF_tke_trplCorr = mass_flux_tke_trplCorr
        self.mass_flux_entr  = entr_scheme ; self.MF_small_ap     = mass_flux_small_ap
        self.triple_corr_opt = 0
        if self.mass_flux_entr == 'R10_HB09': self.triple_corr_opt = 1
        if self.mass_flux_entr == 'R10_Wi11': self.triple_corr_opt = 2
        self.mf_params  = np.array([Cent,Cdet,wp_a,wp_b,wp_bp,up_c,vp_c,bc_ap,delta_bkg])
        self.wp0 = wp0                    ; self.bc_ap = bc_ap
        ####################################
        # define vertical grid
        #-----------------------------------
        if gridType!='ORCA75':
          Sc_r  = np.arange(-nz+0.5, 0.5, 1) / float(nz)
          Sc_w  = np.arange(-nz, 1, 1) / float(nz)
          Cs_r = (1.-np.cosh(thetas*Sc_r))/(np.cosh(thetas)-1.)
          Cs_w = (1.-np.cosh(thetas*Sc_w))/(np.cosh(thetas)-1.)
          self.z_w   = (hc*Sc_w+Cs_w*h0)*h0/(h0+hc)
          self.z_r   = (hc*Sc_r+Cs_r*h0)*h0/(h0+hc)
        else:
          # construct the ORCA 75 levels grid and extract the levels between 0 and h0
          nz   = 75
          zsur = -3958.95137127683; za2  = 100.760928500000
          za0  =  103.953009600000; za1  = 2.41595126900000
          zkth =  15.3510137000000; zkth2= 48.0298937200000
          zacr =  7.00000000000000; zacr2= 13.0000000000000
          Sc_r = np.arange( nz-0.5, 0.5 , -1.)
          Sc_w = np.arange( nz    , 0.  , -1.)
          z_w  = -( zsur + za0 * Sc_w + za1 * zacr * np.log(np.cosh((Sc_w-zkth )/zacr) ) + za2 * zacr2* np.log(np.cosh( (Sc_w-zkth2) / zacr2 ) )  )
          z_r  = -( zsur + za0 * Sc_r + za1 * zacr * np.log(np.cosh((Sc_r-zkth )/zacr) ) + za2 * zacr2* np.log(np.cosh( (Sc_r-zkth2) / zacr2 ) )  )
          Hz   = z_w[1:]-z_w[:-1]
          nbot = np.argmin(z_w <= -h0)
          self.z_w     = z_w[nbot-1:]
          self.nz      = len(self.z_w) - 1
          self.z_r     = z_r[nbot-1:]
          self.z_w[-1] = 0.
          #
        self.Hz    = self.z_w[1:]-self.z_w[:-1]
        zInvMin    = -self.Hz[-1]
        ####################################
        # define the initial state
        #-----------------------------------
        self.u_n      = np.zeros(self.nz); self.v_n      = np.zeros(self.nz)
        self.u_np1    = np.zeros(self.nz); self.v_np1    = np.zeros(self.nz)
        self.t_n      = np.zeros((self.nz,self.ntra), order='F')
        self.t_np1    = np.zeros((self.nz,self.ntra), order='F')
        self.t_n[:,self.isalt] = SaltCst                      # constant salinity
        self.t_np1[:,self.isalt] = self.t_n[:,self.isalt]     # constant salinity
        #
        self.T0      = T0;  self.N0      = N0;
        strat        = N0 / (self.g * self.alpha )
        self.t_n  [:,self.itemp] = T0
        self.t_np1[:,self.itemp] = T0
        for k in range(0, self.nz):
          if self.z_w[k] < mld_ini :
            self.t_n  [k,self.itemp] = T0 + 0.5*strat*( self.z_w[k+1]+self.z_w[k] ) + strat*(-mld_ini)
            self.t_np1[k,self.itemp] = T0 + 0.5*strat*( self.z_w[k+1]+self.z_w[k] ) + strat*(-mld_ini)
        self.btflx[self.itemp]  = 0.
        if btflx=='linear_continuation': self.btflx[self.itemp] = aktmin * N0 / (self.g * self.alpha )
        self.btflx[self.isalt]  = 0.
        ####################################
        # initialize arrays for EDDY-DIFFUSION (tke, bvf, lup,ldwn,tke,AkT,Akv)
        #-----------------------------------
        self.akt      = np.zeros(self.nz+1); self.akt[:]     = aktmin
        self.akv      = np.zeros(self.nz+1); self.akv[:]     = akvmin
        self.bvf      = np.zeros(self.nz+1); self.eps_n      = np.zeros(self.nz+1)
        if self.ED:
          self.tke_n    = np.zeros(self.nz+1); self.tke_n[:]   = tkemin
          self.tke_np1  = np.zeros(self.nz+1); self.tke_np1[:] = tkemin
          self.lupw     = np.zeros(self.nz+1); self.lupw[:]    = mxlmin
          self.ldwn     = np.zeros(self.nz+1); self.ldwn[:]    = mxlmin
          self.Prdtl      = np.zeros(self.nz+1)
          self.eps_n[:]   = epsmin  # TKE dissipation
          self.eps_np1  = np.zeros(self.nz+1); self.eps_np1[:] = epsmin
          self.shear    = np.zeros(self.nz+1); self.residual   = np.array(0.)
          if self.ED_scheme=='Keps':
            self.cmu      = np.zeros(self.nz+1); self.cmu     [:] = 0.1
            self.cmu_prim = np.zeros(self.nz+1); self.cmu_prim[:] = 0.1
            self.cmu_star = np.zeros(self.nz+1); self.cmu_star[:] = 0.1
            self.varT_n   = np.zeros(self.nz+1); self.varT_n  [:] = 1.e-14
            self.varT_np1 = np.zeros(self.nz+1); self.varT_np1[:] = 1.e-14
        else:
          self.ED_scheme = 'None'
        # MLD computation params
        self.rhoc   = 0.01    # density criterion for mixed layer depth (consistent with NEMO)
        self.hmxl   = np.array(zInvMin)
        ####################################
        # initialize arrays for MASS-FLUX
        #-----------------------------------
        if self.MF_tra or self.MF_dyn:
          self.zinvMin = np.array(zInvMin)
          self.zinv   = self.zinvMin
          self.ap     = np.zeros(self.nz+1); self.wp     = np.zeros(self.nz+1)
          self.tp     = np.zeros((self.nz+1, self.ntraMF), order='F');
          self.Bp     = np.zeros(self.nz+1); self.Fmass  = np.zeros(self.nz+1)
          self.ent    = np.zeros(self.nz  ); self.det    = np.zeros(self.nz  )
          self.up     = np.zeros(self.nz+1); self.vp     = np.zeros(self.nz+1)
          self.epsPlume = np.zeros(self.nz)  # TKE plume dissipation
          if self.mass_flux_entr=='R10corNT': self.vortp  = np.zeros(self.nz+1)
        ####################################
        # store MASS-FLUX solution (temperature, velocity)
        # in order to avoid writing a netcdf file
        #-----------------------------------
        dimt = int( nbhours // outfreq )
        self.t_history = np.zeros((self.nz,dimt), order='F')
        self.u_history = np.zeros((self.nz,dimt), order='F')
        self.v_history = np.zeros((self.nz,dimt), order='F')
        ####################################
        # edmf_diags
        ####################################
        if self.write_netcdf:
            self.wted   = np.zeros(self.nz+1); self.wtmf   = np.zeros(self.nz+1)
            self.wued   = np.zeros(self.nz+1); self.wumf   = np.zeros(self.nz+1)
            self.wved   = np.zeros(self.nz+1); self.wvmf   = np.zeros(self.nz+1)
        #
        self.buoyMF      = np.zeros(self.nz+1); self.shearMF     = np.zeros(self.nz+1)
        self.wtke        = np.zeros(self.nz  ); self.triple_corr = np.zeros(self.nz+1)
        ####################################
        # Energy diagnostics
        #-----------------------------------
        if self.write_netcdf:
            self.vint_Etot = np.array(1e-14); self.vint_Ekin = np.array(1e-14)
            self.vint_Epot = np.array(1e-14); self.vint_TKE  = np.array(1e-14)
            self.vint_Eps  = np.array(1e-14)
        ################################################
        # Average over the last hour of the simulation
        #-----------------------------------------------
        if self.do_avg:
            self.u_avg      = np.zeros(self.nz)
            self.v_avg      = np.zeros(self.nz)
            self.t_avg      = np.zeros((self.nz,self.ntra), order='F')
            self.tke_avg    = np.zeros(self.nz+1)
            self.wt_avg     = np.zeros(self.nz+1)
            self.akv_avg    = np.zeros(self.nz+1)
            self.akt_avg    = np.zeros(self.nz+1)
            self.navg       = 0.
################################################################################################################



# store MASS-FLUX solution (temperature)
#-----------------------------------
#
    def run_direct(self):
        #=============================
        kout = 0
        if self.write_netcdf:
            self.output_init( ); #to initialize .nc dataset
        #save initial temperature and velocity
        self.t_history[:,kout] = self.t_n[:,self.itemp]
        self.u_history[:,kout] = self.u_n[:]; self.v_history[:,kout] = self.v_n[:]
        #
        swr_frac = scm_oce.lmd_swfrac(self.Hz,self.nz)   ## Compute fraction of solar shortwave flux penetrating to specified depth
        if self.ED_scheme=='TKE': self.do_TKE( )                        ## Initialize eddy-diffusion scheme to compute eddy diffusivity/eddy viscosity for the first time-step
        if self.ED_scheme=='Keps': self.do_KEPS( )
        ####################################
        # Time LOOP
        #-----------------------------------
        for kt in range(self.nbsteps):
            # update current time
            time = self.dt*float(kt+1)                   ## time in seconds at the end of the current time-step
            modulation = 1.
            if self.DC: modulation = max(np.cos(2.*np.pi*(time/86400. - 0.5)),0.)
            self.srflx = modulation*self.Qswmax
            #===================================================
            # Advance tracers to n+1 (vertical diffusion only)
            #===================================================
            self.t_np1  = scm_oce.advance_tra_ed(
                                        self.t_n, self.stflx, self.srflx,
                                        swr_frac, self.btflx, self.Hz   , self.akt  ,
                                        self.z_w, self.eps_n, self.eos_params[1],
                                        self.dt , self.nz   , self.ntra  )
            #==============================================================
            # advance dynamics to n+1 (Coriolis + vertical viscosity only)
            #==============================================================
            self.u_np1, self.v_np1 = scm_oce.advance_dyn_cor_ed(
                                        self.u_n, self.v_n, self.ustr_sfc, self.vstr_sfc,
                                        self.ustr_bot, self.vstr_bot,
                                        self.Hz , self.akv, self.fcor, self.dt, self.nz  )
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
            #if self.ED:
            if self.ED_scheme=='TKE' : self.do_TKE( )
            if self.ED_scheme=='Keps': self.do_KEPS( )
            #==========================================================
            if time >= self.start_avg and self.do_avg:
              self.do_turb_fluxes (  )
              self.u_avg[:]   += self.u_np1[:]; self.v_avg[:]   += self.v_np1[:]
              self.t_avg[:]   += self.t_np1[:]; self.tke_avg[:] += self.tke_np1[:]
              self.wt_avg[:]  +=  - (self.wted[:]+self.wtmf[:])
              self.akv_avg[:] += self.akv[:]; self.akt_avg[:] += self.akt[:]
              self.navg       += 1.
            #==========================================================
            # Check for outputs & diagnostics
            #==========================================================
            if  time % self.outfreq == 0:
              kout+=1
              # save in np.array advanced temperature and velocity
              self.t_history[:,kout-1] = self.t_np1[:,self.itemp]
              self.u_history[:,kout-1] = self.u_np1[:]
              self.v_history[:,kout-1] = self.v_np1[:]

              if self.write_netcdf:
                # compute diagnostics
                if self.ED: self.do_diags_energy(  )
                self.do_turb_fluxes (  )
                #=======================================
                # Compute mixed layer depth from bvf
                #self.hmxl  = scm_oce.compute_mxl2(self.bvf,self.rhoc,self.z_r,-10.,self.rho0,self.nz)
                self.hmxl  = scm_oce.compute_mxl2(self.bvf,0.01,self.z_r,-0.5,self.rho0,self.nz)
                #=======================================
                # Write outputs in .nc file
                self.output_state(TimeInSecs=time,kout=kout)
            #===================================================================
            # swap arrays
            #===================================================================
            self.u_n[:] = self.u_np1[:]; self.v_n[:] = self.v_np1[:];
            self.t_n[:,:] = self.t_np1[:,:]; self.tke_n[:] = self.tke_np1[:]
            if self.ED_scheme=='Keps':
                self.eps_n [:] = self.eps_np1 [:]
                self.varT_n[:] = self.varT_np1[:]
            #===================================================================
            if kt<self.nbsteps-1: self.wtke[:]  = 0.                            ## diagnostics : reset the array containing w'e
        if self.write_netcdf:
          self.do_turb_fluxes (  )
        if self.do_avg:
          cff = 1./self.navg
          self.u_avg[:]   = cff*self.u_avg[:]
          self.v_avg[:]   = cff*self.v_avg[:]
          if self.ED: self.tke_avg[:] = cff*self.tke_avg[:]
          self.wt_avg[:]  = cff*self.wt_avg[:]
          self.t_avg[:]   = cff*self.t_avg[:]
          self.akv_avg[:] = cff*self.akv_avg[:]
          self.akt_avg[:] = cff*self.akt_avg[:]
#












#
    def do_turb_fluxes(self):
        # Turbulent fluxes diagnostics
        self.wted[  self.nz] = self.stflx[self.itemp]                            ## (w'theta')_ED
        self.wted[1:self.nz] = self.akt[1:self.nz] * (  self.t_np1[1:self.nz  ,self.itemp]
                                                      - self.t_np1[0:self.nz-1,self.itemp] ) / self.Hz[0:self.nz-1]
        self.wued[  self.nz] = self.ustr_sfc                                     ## (w'theta')_ED
        self.wved[  self.nz] = self.vstr_sfc
        self.wued[1:self.nz] = self.akv[1:self.nz] * (  self.u_np1[1:self.nz]
                                                      - self.u_np1[0:self.nz-1] ) / self.Hz[0:self.nz-1]
        self.wved[1:self.nz] = self.akv[1:self.nz] * (  self.v_np1[1:self.nz]
                                                      - self.v_np1[0:self.nz-1] ) / self.Hz[0:self.nz-1]
        if self.MF_tra:                                                         ## (w'theta')_MF
          self.wtmf[1:self.nz] =  self.Fmass[1:self.nz]*(
                                              self.tp   [1:self.nz,  self.itemp]
                                            - self.t_np1[0:self.nz-1,self.itemp] )
        if self.MF_dyn:                                                         ## (w'u')_MF
          self.wumf[1:self.nz] =  self.Fmass[1:self.nz]*(
                                               self.up   [1:self.nz]
                                             - self.u_np1[0:self.nz-1] )
          self.wvmf[1:self.nz] =  self.Fmass[1:self.nz]*(
                                               self.vp   [1:self.nz]
                                             - self.v_np1[0:self.nz-1] )
#












#
    def do_diags_energy(self):
        #Compute b_n and b_np1 for diagnostics purposes
        if self.lineos:
          rho_n,bvf = scm_oce.rho_eos_lin(
                            self.t_n[:,self.itemp], self.t_n[:,self.isalt],
                            self.z_r , self.eos_params , self.nz, len(self.eos_params)  )
          rho_np1,bvf = scm_oce.rho_eos_lin(
                            self.t_np1[:,self.itemp], self.t_np1[:,self.isalt],
                            self.z_r , self.eos_params , self.nz, len(self.eos_params)  )
        else:
          rho_n,bvf = scm_oce.rho_eos(
                            self.t_n[:,self.itemp], self.t_n[:,self.isalt],
                                             self.z_r, self.z_w, self.rho0, self.nz )
          rho_np1,bvf = scm_oce.rho_eos(
                            self.t_np1[:,self.itemp], self.t_np1[:,self.isalt],
                                             self.z_r, self.z_w, self.rho0, self.nz )
        # Compute buoyancy
        b_n   = -(self.g/self.rho0)*rho_n[:]
        b_np1 = -(self.g/self.rho0)*rho_np1[:]
        # Energy diagnostics
        self.vint_Ekin = 0.;self.vint_Epot = 0.
        self.vint_Eps  = 0.;self.vint_TKE  = 0.
        dissip         = 0.
        #
        for k in range(self.nz):
          self.vint_Ekin = self.vint_Ekin + 0.5*self.Hz[k] * (
                                             self.u_np1[k]*self.u_np1[k]
                                          -  self.u_n  [k]*self.u_n  [k]
                                          +  self.v_np1[k]*self.v_np1[k]
                                          -  self.v_n  [k]*self.v_n  [k] ) / self.dt   # m3/s3
          self.vint_Epot = self.vint_Epot + self.Hz[k]*(-self.z_r[k])*(b_np1[k]-b_n[k])/self.dt
          dissip = dissip + self.Hz[k]*self.cp*(self.t_np1[k,self.itemp]-self.t_n[k,self.itemp])/self.dt
        #
        for k in range(self.nz-1):
          self.vint_Eps = self.vint_Eps+(self.z_r[k+1]-self.z_r[k])*self.eps_n[k+1]
          self.vint_TKE = self.vint_TKE+(self.z_r[k+1]-self.z_r[k])*(self.tke_np1[k+1]-self.tke_n[k+1])/self.dt
        # remove surface contributions
        sfc_KE   = 0.5*( (self.u_np1[-1]+self.u_n[-1])*self.ustr_sfc
                       + (self.v_np1[-1]+self.v_n[-1])*self.vstr_sfc )
        sfc_Epot = -self.cp*self.stflx[self.itemp]
        sfc_Epot = sfc_Epot - 0.5*self.Hz[-1]*self.g*self.eos_params[1]*self.stflx[self.itemp]
        #
        sfc_TKE  = self.wtke[-1]
        self.vint_Ekin = self.vint_Ekin - sfc_KE
        self.vint_Epot = self.vint_Epot + dissip  + sfc_Epot
        self.vint_TKE  = self.vint_TKE  + sfc_TKE - self.residual
        self.vint_Etot = self.vint_Epot + self.vint_TKE + self.vint_Ekin
#












#
    def do_TKE(self):
        #==========================================================================
        # Compute Brunt-Vaisala frequency bvf for TKE production/destruction term
        #==========================================================================
        if self.lineos:
          rho,self.bvf = scm_oce.rho_eos_lin(
                                 self.t_np1[:,self.itemp], self.t_np1[:,self.isalt],
                                 self.z_r, self.eos_params,self.nz, len(self.eos_params) )
        else:
          rho,self.bvf = scm_oce.rho_eos(self.t_np1[:,self.itemp],self.t_np1[:,self.isalt],
                                 self.z_r,self.z_w,self.rho0,self.nz)
        #=======================================
        # Compute boundary conditions for TKE
        tkemin = self.min_Threshold[0]
        if self.MF_tke_trplCorr:
            tke_sfc,tke_bot, flux_sfc = scm_tke.compute_tke_bdy(
                                      self.ustr_sfc,  self.vstr_sfc, self.ED_tke_const,
                                      self.bc_ap, self.wp0, tkemin)
        else:
           tke_sfc,tke_bot, flux_sfc = scm_tke.compute_tke_bdy(
                                      self.ustr_sfc,  self.vstr_sfc, self.ED_tke_const,
                                      0.*self.bc_ap, 0.*self.wp0, tkemin)
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
          self.triple_corr = scm_mfc.compute_triplecorr(
                                        self.tke_n, self.tp[:,self.isalt+1], self.Fmass,
                                        self.up   , self.vp   ,self.wp  ,
                                        self.u_np1, self.v_np1, self.Hz ,
                                        self.z_r  , self.triple_corr_opt, self.wtke, self.nz  )
        #=========================
        # Advance TKE
        self.tke_np1, self.Prdtl, self.eps_n, self.residual = scm_tke.advance_tke(
                                         self.tke_n, self.lupw, self.ldwn,
                                         self.akv  , self.akt , self.Hz  ,
                                         self.z_r  , self.bvf , self.buoyMF,
                                         self.shear, self.shearMF, self.triple_corr,
                                         self.wtke, self.dt, tke_sfc, tke_bot, flux_sfc,
                                         self.ED_tke_sfc_dirichlet, self.ED_tke_const  ,
                                         tkemin, self.nz  )
        #===================================================
        # Finalize eddy-viscosity/diffusivity computation
        mxlmin = self.min_Threshold[3]
        self.lupw,self.ldwn = scm_tke.compute_mxl(
                                         self.tke_np1, self.bvf , self.Hz,
                                         self.ustr_sfc   , self.vstr_sfc, mxlmin, self.nz )
        Akvmin = self.min_Threshold[1]; Aktmin = self.min_Threshold[2]
        self.akv,self.akt   = scm_tke.compute_ed (
                                         self.tke_np1, self.lupw, self.ldwn,
                                         self.Prdtl  , self.ED_extrap_sfc,
                                         self.ED_tke_const, Akvmin, Aktmin, self.nz  )
        #=========================
        # Apply EVD if necessary
        if self.ED_evd: scm_oce.compute_evd(
                                          self.bvf, self.akv, self.akt,
                                          10., self.nz  )
        #












#
    def do_MF(self):
        #=======================================================
        # Surface boundary conditions for up,vp and tracers
        up0,vp0,tp0 = scm_mfc.compute_mf_bdy(
                                      self.u_np1[-2:]  , self.v_np1[-2:],
                                      self.t_np1[-2:,:], self.tke_n[-2:],
                                      self.Hz[-2:]     , self.ntra, 2 )
        wp0=self.wp0
        #=================================================================
        # Compute the mean quantities used to constrain the mass flux eqns
        u_mean,v_mean,t_mean,dtke_m = scm_mfc.compute_mf_forcing(
                                      self.u_np1, self.v_np1, self.t_np1,
                                      self.tke_n, self.nz   , self.ntra )
        #=======================================================================================
        # Compute mass flux variables (Rio et al. 2010 closure or Pergaud et al. 2009 closure)
        tkepmin = self.min_Threshold[0]
        mxlpmin = self.min_Threshold[3]
        #
        if self.mass_flux_entr=='R10' or self.mass_flux_entr=='R10_HB09' or self.mass_flux_entr=='R10_Wi11':
          self.ap,self.up,self.vp,self.wp,self.tp,self.Bp,self.ent,self.det, self.epsPlume = scm_mfc.mass_flux_r10(
                                    u_mean, v_mean, t_mean, self.tke_n, self.z_w, self.Hz,
                                    tp0   , up0   , vp0   , wp0     , self.mf_params,
                                    self.eos_params, tkepmin, mxlpmin, self.MF_small_ap, self.lineos, self.triple_corr_opt, self.zinv,
                                    self.nz , self.ntraMF , len(self.mf_params), len(self.eos_params)   )
        if self.mass_flux_entr=='R10corNT':
          self.ap,self.up,self.vp,self.wp,self.tp,self.Bp,self.ent,self.det, self.vortp, self.epsPlume = scm_mfc.mass_flux_r10_cor(
                                    u_mean, v_mean, t_mean, self.tke_n, self.z_w, self.Hz,
                                    tp0   , up0   , vp0   , wp0     , self.mf_params,
                                    self.eos_params, self.fcor, self.ecor, tkepmin, mxlpmin, self.MF_small_ap, self.lineos, self.zinv ,
                                    self.nz , self.ntraMF , len(self.mf_params), len(self.eos_params)   )
        self.zinv = min(self.zinvMin,self.zinv)
    #











    #
    def do_KEPS(self):
        #==========================================================================
        # Compute Brunt-Vaisala frequency bvf for TKE production/destruction term
        #==========================================================================
        if self.lineos:
          rho,self.bvf = scm_oce.rho_eos_lin(
                                 self.t_np1[:,self.itemp], self.t_np1[:,self.isalt],
                                 self.z_r, self.eos_params,self.nz, len(self.eos_params) )
        else:
          rho,self.bvf = scm_oce.rho_eos(self.t_np1[:,self.itemp],self.t_np1[:,self.isalt],
                                 self.z_r,self.z_w,self.rho0,self.nz)
        #=======================================================================
        # Compute TKE production by shear (Not multiplied by Akv !!)
        #=======================================================================
        self.shear2 = scm_gls.compute_shear(
                                  self.u_np1, self.v_np1, self.u_np1,
                                  self.v_np1, self.z_r  , self.nz  )
        #=======================================================================
        # Compute boundary conditions
        #=======================================================================
        tke1 = 0.5*(self.tke_n[self.nz]+self.tke_n[self.nz-1])
        tke2 = 0.5*(self.tke_n[      0]+self.tke_n[        1])
        tkemin = self.min_Threshold[0]
        epsmin = self.min_Threshold[4]
        #
        tke_sfc,tke_bot,ftke_sfc,ftke_bot,eps_sfc,eps_bot,feps_sfc,feps_bot = scm_gls.compute_tke_eps_bdy(
            self.ustr_sfc, self.vstr_sfc, self.ustr_bot, self.vstr_bot,
            self.z0b, tke1, self.Hz[-1], tke2, self.Hz[ 0], self.OneOverSig_psi, self.pnm, tkemin, epsmin)
        #
        if self.bdy_tke_sfc[0] < 0.5: self.bdy_tke_sfc[1] = tke_sfc
        else: self.bdy_tke_sfc[1] = ftke_sfc
        if self.bdy_tke_bot[0] < 0.5: self.bdy_tke_bot[1] = tke_bot
        else: self.bdy_tke_bot[1] = ftke_bot
        #
        if self.bdy_eps_sfc[0] < 0.5: self.bdy_eps_sfc[1] = eps_sfc
        else: self.bdy_eps_sfc[1] = feps_sfc
        if self.bdy_eps_bot[0] < 0.5: self.bdy_eps_bot[1] = eps_bot
        else: self.bdy_eps_bot[1] = feps_bot
        #
        #================================================
        #
        self.tke_np1[:], self.wtke[:] = scm_gls.advance_turb_tke(self.tke_n, self.bvf, self.shear2,
                            self.OneOverSig_k*self.akv,self.akv,self.akt,self.eps_n,self.Hz,
                                   self.dt,tkemin,self.bdy_tke_sfc,self.bdy_tke_bot,self.nz)
        #
        self.eps_np1[:] = scm_gls.advance_turb_eps(self.eps_n, self.bvf, self.shear2,
                                self.OneOverSig_psi*self.akv, self.cmu[:],
                                self.cmu_prim[:], self.tke_n, self.tke_np1,
                                self.Hz   ,self.dt         , self.betaCoef  ,
                                epsmin, self.bdy_eps_sfc,self.bdy_eps_bot,  self.nz)
        #=========================
        akvmin = self.min_Threshold[1]; aktmin = self.min_Threshold[2]
        self.akv,self.akt,self.cmu,self.cmu_prim, self.cmu_star = scm_gls.compute_ev_ed_filt_ng(
                        self.tke_np1, self.eps_np1, self.bvf,
                        self.shear2 , self.pnm, akvmin, aktmin, epsmin, self.nz )
        #=========================
        self.varT_np1[:] = scm_gls.advance_turb_vart(self.varT_n,self.akv,self.akt,self.t_np1[:,self.itemp],
                                                     self.eps_np1,self.tke_np1,self.cmu_star,
                                                     self.Hz, self.dt, self.alpha, 1.e-14, self.nz)
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
        fh01.rho0            = self.eos_params[0]
        fh01.alpha           = self.eos_params[1]
        fh01.beta            = self.eos_params[2]
        fh01.T0              = self.eos_params[3]
        fh01.S0              = self.eos_params[4]
        #self.min_Threshold
        fh01.tkemin          = self.min_Threshold[0]
        fh01.akvmin          = self.min_Threshold[1]
        fh01.aktmin          = self.min_Threshold[2]
        fh01.mxlmin          = self.min_Threshold[3]
        #
        ocean_time = fh01.createVariable('ocean_time','f8',('time')); ocean_time[:] = 0.
        ocean_time.units = 'seconds'
        ocean_time.long_name = 'time since initialization'
        taux = fh01.createVariable('taux','f8',('time')); taux[:] = self.ustr_sfc
        tauy = fh01.createVariable('tauy','f8',('time')); tauy[:] = self.vstr_sfc
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
        var  = fh01.createVariable('WT','f8',('time','z_w')); var[0,:] = self.wted[:]+self.wtmf[:]; var.units = 'K m s-1'; del var
        var  = fh01.createVariable('WU','f8',('time','z_w')); var[0,:] = self.wued[:]+self.wumf[:]; var.units = 'm2 s-2'; del var
        var  = fh01.createVariable('WV','f8',('time','z_w')); var[0,:] = self.wved[:]+self.wvmf[:]; var.units = 'm2 s-2'; del var

        if self.ED:
            var  = fh01.createVariable('tke','f8',('time','z_w')); var[0,:] = self.tke_n[:]; var.units = 'm2 s-2'; del var
            var  = fh01.createVariable('Akv','f8',('time','z_w')); var[0,:] = self.akv[:]; var.units = 'm2 s-1'; del var
            var  = fh01.createVariable('Akt','f8',('time','z_w')); var[0,:] = self.akt[:]; var.units = 'm2 s-1'; del var
            var  = fh01.createVariable('bvf','f8',('time','z_w')); var[0,:] = self.bvf[:]; var.units = 's-2'; del var
            var  = fh01.createVariable('lup','f8',('time','z_w')); var[0,:] = self.lupw[:]; var.units = 'm'; del var
            var  = fh01.createVariable('ldw','f8',('time','z_w')); var[0,:] = self.ldwn[:]; var.units = 'm'; del var
            var  = fh01.createVariable('WT_ED','f8',('time','z_w')); var[0,:] = self.wted[:]; var.units = 'K m s-1'; del var
            var  = fh01.createVariable('WU_ED','f8',('time','z_w')); var[0,:] = self.wued[:]; var.units = 'm2 s-2'; del var
            var  = fh01.createVariable('WV_ED','f8',('time','z_w')); var[0,:] = self.wved[:]; var.units = 'm2 s-2'; del var
            var  = fh01.createVariable('Prdtl','f8',('time','z_w')); var[0,:] = self.akv[:]/self.akt[:]; var.units = 'none'; del var
            var  = fh01.createVariable('epsil','f8',('time','z_w')); var[0,:] = self.eps_n[:]; var.units = 'm2 s-3'; del var
            if self.ED_scheme=='Keps':
                var = fh01.createVariable('varT','f8',('time','z_w')); var[0,:] = self.varT_n[:]; del var
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
            var  = fh01.createVariable('WU_MF','f8',('time','z_w')); var[0,:] = self.wumf[:]; del var
            var  = fh01.createVariable('WV_MF','f8',('time','z_w')); var[0,:] = self.wvmf[:]; del var
            var = fh01.createVariable('u_p','f8',('time','z_w')); var[0,:] = self.up[:]; del var
            var = fh01.createVariable('v_p','f8',('time','z_w')); var[0,:] = self.vp[:]; del var
            if self.MF_tke:
                var = fh01.createVariable('buoyMF','f8',('time','z_w')); var[0,:] = self.buoyMF[:]; del var
                var = fh01.createVariable('shearMF','f8',('time','z_w')); var[0,:] = self.shearMF[:]; del var
            if self.MF_tke_trplCorr:
                var = fh01.createVariable('we','f8',('time','z_w')); var[0,:] = self.triple_corr[:]; del var
                var = fh01.createVariable('tke_p','f8',('time','z_w')); var[0,:] = self.tp[:,self.isalt+1]; del var
            if self.mass_flux_entr=='R10corNT':
                var = fh01.createVariable('vort_p','f8',('time','z_w')); var[0,:] = self.vortp[:]; del var
        fh01.close()
#












#
    def output_state(self,TimeInSecs,kout):
        fh01 = Dataset(self.output, mode='a',format="NETCDF4")
        fh01.variables['ocean_time'][kout] = TimeInSecs
        fh01.variables['taux'][kout]       = self.ustr_sfc
        fh01.variables['tauy'][kout]       = self.vstr_sfc
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
        fh01.variables['WU'][kout,:]       = self.wued[:]+self.wumf[:]
        fh01.variables['WV'][kout,:]       = self.wved[:]+self.wvmf[:]


        if self.ED:
            fh01.variables['tke'][kout,:] = self.tke_n[:]
            fh01.variables['Akv'][kout,:] = self.akv[:]
            fh01.variables['Akt'][kout,:] = self.akt[:]
            fh01.variables['bvf'][kout,:] = self.bvf[:]
            fh01.variables['lup'][kout,:] = self.lupw[:]
            fh01.variables['ldw'][kout,:] = self.ldwn[:]
            fh01.variables['WT_ED'][kout,:] = self.wted[:]
            fh01.variables['WU_ED'][kout,:] = self.wued[:]
            fh01.variables['WV_ED'][kout,:] = self.wved[:]
            fh01.variables['Prdtl'][kout,:] = self.akv[:]/self.akt[:]
            fh01.variables['epsil'][kout,:] = self.eps_n[:]
            if self.ED_scheme=='Keps':
                fh01.variables['varT'][kout,:] = self.varT_n[:]
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
            fh01.variables['WU_MF'][kout,:] = self.wumf[:]
            fh01.variables['WV_MF'][kout,:] = self.wvmf[:]
            fh01.variables['u_p'][kout,:] = self.up[:]
            fh01.variables['v_p'][kout,:] = self.vp[:]
            if self.MF_tke:
                fh01.variables['buoyMF'][kout,:] = self.buoyMF[:]
                fh01.variables['shearMF'][kout,:] = self.shearMF[:]
            if self.MF_tke_trplCorr:
                fh01.variables['we'][kout,:] = self.triple_corr[:]
                fh01.variables['tke_p'][kout,:] = self.tp[:,self.isalt+1]
            if self.mass_flux_entr=='R10corNT':
                fh01.variables['vort_p'][kout,:] = self.vortp[:]
        fh01.close()

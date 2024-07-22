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
from datetime import datetime
from netCDF4 import Dataset
###########################################
# Classes
###########################################
class SCM:
    """Single Column model class"""

    def __init__(self,  nz     = 75  ,   dt =  1200. ,
                        initial_filename="../data/asicsmed/init_ASICS_m01d15.nc",
                        sfcforc_filename="../data/asicsmed/forc_ASICS_y2013.nc",
                        lat0 =   42.04, nbhours = 720., outfreq  = 1.,
                        output_filename="scm_asicsmed.nc" , eddy_diff = True ,
                        evd = False  , mass_flux_tra          = False  ,
                                       mass_flux_dyn          = False  ,
                                       mass_flux_tke          = False  ,
                                       mass_flux_tke_trplCorr = False  ,
                                       mass_flux_small_ap     = False  ,
                        lin_eos           = False , extrap_ak_surf      = True ,
                        tke_sfc_dirichlet = False , eddy_diff_scheme = 'TKE'   ,
                        eddy_diff_tke_const = 'NEMO',
                        akvmin   = 1.e-4      , aktmin   = 1.e-5 ,mxlmin  = 1.0 ,
                        Cent  = 0.99  , Cdet = 1.99 , wp_a =  1 , wp_b  = 1   ,
                        wp_bp = 0.0002, up_c = 0.5, vp_c = 0.5, bc_ap = 0.2 ,
                        delta_bkg = 0., wp0=-1.e-08, entr_scheme = 'R10'  ):
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
        dt1            = datetime(1850,1, 1,00,00)
        dt2            = datetime(2013,1,15,00,00)  # initial time for the simulation
        #dt3            = datetime(2013,3,26,00,00)  # final time for the simulation
        dt3            = datetime(2013,2,15,00,00)  # final time for the simulation
        timediff1      = dt2-dt1
        timediff2      = dt3-dt2
        self.inifrctime= int(timediff1.total_seconds()/3600.)
        self.nbsteps   = int( ( timediff2.total_seconds() )/self.dt )           ## number of time steps
        self.outfreq   = int( outfreq*3600. )                                   ## frequency of outputs in seconds
        #====================================
        self.frcname   = sfcforc_filename
        self.kfrcStr   = 336; self.kfrcEnd   = 1681
        #====================================
        # should check if self.outfreq is a multiple of the time-step
        self.output    = output_filename
        # physical constants
        self.rho0      = 1027; self.cp        = 3985.; self.g         = 9.81
        self.lineos    = False; self.eos_params  = np.array([self.rho0])                                                   ## no other option implemented yet
        ####################################
        # surface boundary conditions
        #-----------------------------------
        self.ustr_sfc  = 0.                                                     ## zonal wind stress       [m2/s2]
        self.vstr_sfc  = 0.                                                     ## meridional wind stress  [m2/s2]
        self.ustr_bot  = 0.
        self.vstr_bot  = 0.
        self.stflx     = np.zeros(2)
        self.stflx[self.itemp]  = 0.                                            ## non-penetrative heat flux  [K m / s]
        self.srflx     = 0.                                                     ## solar radiation            [K m / s]
        self.stflx[self.isalt]  = 0.
        self.btflx     = np.zeros(2)                                          ## freshwater flux         [psu m / s]
        self.btflx[self.itemp]  = 0.
        self.btflx[self.isalt]  = 0.
        self.z0b                = 1.e-14
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
            self.betaCoef[2] = -0.4; self.betaCoef[3] = 1.00
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
        self.mass_flux_entr = entr_scheme ; self.MF_small_ap     = mass_flux_small_ap
        self.mf_params  = np.array([Cent,Cdet,wp_a,wp_b,wp_bp,up_c,vp_c,bc_ap,delta_bkg])
        #self.mf_params  = np.array([Cent,2.*Cent,wp_a,wp_b,wp_bp,up_c,vp_c,bc_ap,0.5*wp_bp])
        self.wp0 = wp0                    ; self.bc_ap = bc_ap
        ####################################
        # define vertical grid
        #-----------------------------------
        zsur = -3958.95137127683; za2  = 100.760928500000
        za0  =  103.953009600000; za1  = 2.41595126900000
        zkth =  15.3510137000000; zkth2= 48.0298937200000
        zacr =  7.00000000000000; zacr2= 13.0000000000000
        Sc_r     = np.arange( nz-0.5, 0.5 , -1.)
        Sc_w     = np.arange( nz    , 0.  , -1.)    #/ float(nz)
        self.z_w      = -( zsur + za0 * Sc_w + za1 * zacr * np.log(np.cosh((Sc_w-zkth )/zacr) ) + za2 * zacr2* np.log(np.cosh( (Sc_w-zkth2) / zacr2 ) )  )
        self.z_r      = -( zsur + za0 * Sc_r + za1 * zacr * np.log(np.cosh((Sc_r-zkth )/zacr) ) + za2 * zacr2* np.log(np.cosh( (Sc_r-zkth2) / zacr2 ) )  )
        self.nz       = self.nz - 1  # Warning : this is inherent to the NEMO grid definition
        self.Hz       = self.z_w[1:]-self.z_w[:-1]
        zInvMin       = - self.Hz[-1]
        ####################################
        # define the initial state
        #-----------------------------------
        # read initial conditions in initial_filename
        self.u_n        = np.zeros(self.nz); self.v_n        = np.zeros(self.nz)
        self.u_np1      = np.zeros(self.nz); self.v_np1      = np.zeros(self.nz)
        self.u_np1_star = np.zeros(self.nz); self.v_np1_star = np.zeros(self.nz)
        self.uFlx       = np.zeros(self.nz+1); self.vFlx     = np.zeros(self.nz+1)
        fh01      = Dataset(initial_filename, mode='r',format="NETCDF4")
        tini      = np.squeeze(fh01.variables['votemper'][0,:-1,0,0]);
        sini      = np.squeeze(fh01.variables['vosaline'][0,:-1,0,0]);
        self.t_n        = np.zeros((self.nz,self.ntra), order='F')
        self.t_np1      = np.zeros((self.nz,self.ntra), order='F')
        self.t_np1_star = np.zeros((self.nz,self.ntra), order='F')
        self.tFlx       = np.zeros(self.nz+1)
        self.t_n  [:,self.isalt] = np.flip(sini)
        self.t_np1[:,self.isalt] = np.flip(sini)
        self.t_np1_star[:,self.isalt] = np.flip(sini)
        self.t_n  [:,self.itemp] = np.flip(tini)
        self.t_np1[:,self.itemp] = np.flip(tini)
        self.t_np1_star[:,self.itemp] = np.flip(tini)
        #------------------------------------
        self.alpha   = 0.; self.beta    = 0.
        ####################################
        # initialize arrays for EDDY-DIFFUSION (tke, bvf, lup,ldwn,tke,AkT,Akv)
        #-----------------------------------
        self.akt      = np.zeros(self.nz+1); self.akt[:]     = aktmin
        self.akv      = np.zeros(self.nz+1); self.akv[:]     = akvmin
        self.bvf      = np.zeros(self.nz+1); self.eps_n      = np.zeros(self.nz+1)
        if self.ED:
          self.tke_n    = np.zeros(self.nz+1); self.tke_n[:]   = tkemin
          self.tke_env  = np.zeros(self.nz+1); self.tke_env[:] = tkemin
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
        else:
          self.ED_scheme = 'None'
        # MLD computation params
        self.rhoc    = 0.01    # density criterion for mixed layer depth (consistent with NEMO)
        self.hmxl    = np.array(zInvMin)
        self.rhoc10  = 0.01            # density criterion for mixed layer depth (with reference level at 10m)
        self.rhoc300 = 0.01-0.00001    # density criterion for mixed layer depth (with reference level at 300m)
        self.zinv    = np.array(-10.)
        self.hmxl10  = np.array(-10.)
        self.hmxl300 = np.array(-300.)
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
        # edmf_diags
        ####################################
        self.wted   = np.zeros(self.nz+1); self.wtmf   = np.zeros(self.nz+1)
        self.wsed   = np.zeros(self.nz+1); self.wsmf   = np.zeros(self.nz+1)
        self.wued   = np.zeros(self.nz+1); self.wumf   = np.zeros(self.nz+1)
        self.wved   = np.zeros(self.nz+1); self.wvmf   = np.zeros(self.nz+1)
        self.buoyMF      = np.zeros(self.nz+1); self.shearMF     = np.zeros(self.nz+1)
        self.wtke        = np.zeros(self.nz  ); self.triple_corr = np.zeros(self.nz+1)
################################################################################################################
#-----------------------------------
#
    def run_direct(self):
        ####################################
        printMsg     = 3600.*24.
        write_netcdf = True
        # penalisation
        #-----------------------------------
        sigma    = np.zeros(self.nz)
        rr       = -(self.z_r[:]+2600.)/self.Hz[:]
        sigma    = 0.5*(1.+np.tanh(2.*(rr-0.64)))
        traLS    = np.zeros(2)
        #=============================================================================================================
        if write_netcdf: self.output_init( ); kout = 0;
        kfrc = 0
        swr_frac  = scm_oce.lmd_swfrac(self.Hz,self.nz)   ## Compute fraction of solar shortwave flux penetrating to specified depth
        irho0cp   = 1./(self.rho0*self.cp)
        #============================= read forcing for start of the simulation
        frc01     = Dataset(self.frcname, mode='r',format="NETCDF4")
        frc_time  = (frc01.variables['time'][self.kfrcStr:self.kfrcStr+self.kfrcEnd+1]-self.inifrctime)*3600.
        taux1     = np.squeeze(frc01.variables['TAUX2'][self.kfrcStr+kfrc,0,0]); taux2     = taux1
        tauy1     = np.squeeze(frc01.variables['TAUY2'][self.kfrcStr+kfrc,0,0]); tauy2     = tauy1
        Qs1       = np.squeeze(frc01.variables['FSOL2'][self.kfrcStr+kfrc,0,0]); Qs2       = Qs1
        Qnet1     = np.squeeze(frc01.variables['FNET2'][self.kfrcStr+kfrc,0,0]); Qnet2     = Qnet1
        Emp1      = np.squeeze(frc01.variables['EMP2'][self.kfrcStr+kfrc,0,0]) * 1.e+6; Emp2 = Emp1
        frc01.close()
        #=============================================================================================================
        salt = np.zeros(self.nz)  # array used to compute hmxl by ignoring salinity contribution to bvf
        if self.ED_scheme=='TKE': self.do_TKE( )     ## Initialize eddy-diffusion scheme to compute eddy diffusivity/eddy viscosity for the first time-step
        if self.ED_scheme=='Keps': self.do_KEPS( )
        ####################################
        # Time LOOP
        #-----------------------------------
        for kt in range(self.nbsteps):
            # update current time
            traLS[self.itemp] = self.t_n  [17,self.itemp]
            traLS[self.isalt] = self.t_n  [17,self.isalt]
            time = self.dt*float(kt+1)                   ## time in seconds at the end of the current time-step
            #====================================================================
            if time >= frc_time[kfrc]:   # we need to read a new time
              taux1 = taux2; tauy1 = tauy2; kfrc = kfrc + 1
              Qs1   = Qs2  ; Qnet1 = Qnet2; Emp1 = Emp2
              frc01      = Dataset(self.frcname, mode='r',format="NETCDF4")
              taux2 = np.squeeze(frc01.variables['TAUX2'][self.kfrcStr+kfrc,0,0])
              tauy2 = np.squeeze(frc01.variables['TAUY2'][self.kfrcStr+kfrc,0,0])
              Qs2   = np.squeeze(frc01.variables['FSOL2'][self.kfrcStr+kfrc,0,0])
              Qnet2 = np.squeeze(frc01.variables['FNET2'][self.kfrcStr+kfrc,0,0])
              Emp2  = np.squeeze(frc01.variables['EMP2'][self.kfrcStr+kfrc,0,0]) * 1.e+6
              frc01.close()
            #======= surface forcings
            if time < frc_time[0]:
              wght = 0.
            else:
              wght = (time-frc_time[kfrc-1])/(frc_time[kfrc]-frc_time[kfrc-1])
            self.ustr_sfc  = ((1.-wght)*taux1+wght*taux2) / self.rho0         # in m2/s2
            self.vstr_sfc  = ((1.-wght)*tauy1+wght*tauy2) / self.rho0         # in m2/s2
            self.srflx     = ((1.-wght)*Qs1 + wght*Qs2)*irho0cp
            self.stflx[self.itemp] = ((1.-wght)*(Qnet1-Qs1)+wght*(Qnet2-Qs2))*irho0cp
            coeffEMP               = self.t_n[-1,self.isalt] / self.rho0
            self.stflx[self.isalt] = coeffEMP*( (1.-wght)*Emp1 + wght*Emp2 )
            #===================================================
            # Advance tracers to n+1 (vertical diffusion only)
            #===================================================
            self.t_np1_star, self.tFlx = scm_oce.advance_tra_ed(
                                        self.t_n, self.stflx   , self.srflx,
                                        swr_frac, self.btflx   , self.Hz   , self.akt  ,
                                        self.z_r, 0.*self.eps_n, self.alpha,
                                        self.dt , self.nz      , self.ntra  )
            # penalization
            self.t_np1_star[:,self.itemp] = (self.t_np1_star[:,self.itemp]+sigma[:]*self.dt*traLS[self.itemp])/(1.+sigma[:]*self.dt)
            self.t_np1_star[:,self.isalt] = (self.t_np1_star[:,self.isalt]+sigma[:]*self.dt*traLS[self.isalt])/(1.+sigma[:]*self.dt)
            #==============================================================
            # advance dynamics to n+1 (Coriolis + vertical viscosity only)
            #==============================================================
            self.u_np1_star, self.v_np1_star, self.uFlx, self.vFlx = scm_oce.advance_dyn_cor_ed(
                                        self.u_n, self.v_n, self.ustr_sfc, self.vstr_sfc,
                                        self.ustr_bot, self.vstr_bot,
                                        self.Hz , self.akv, self.z_r, self.fcor, self.dt, self.nz  )
            # penalization
            self.u_np1_star[:] = np.exp(-self.dt*sigma[:])*self.u_np1_star[:]
            self.v_np1_star[:] = np.exp(-self.dt*sigma[:])*self.v_np1_star[:]
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
            if self.MF_tra:
                self.t_np1 = scm_oce.advance_tra_mf(
                                      self.t_np1_star, self.tp, self.Fmass, self.tFlx, self.Hz,
                                      self.dt, self.nz, self.ntra )
            else:
                self.t_np1[:,self.itemp] = self.t_np1_star[:,self.itemp]
                self.t_np1[:,self.isalt] = self.t_np1_star[:,self.isalt]
            # apply the MF term to the velocity equation + compute the transfer of KE to subgrid scales stored in shearMF
            if self.MF_dyn:
                self.u_np1, self.v_np1 = scm_oce.advance_dyn_mf(
                                      self.u_np1_star, self.v_np1_star, self.shearMF,
                                      self.u_n  , self.v_n  , self.up , self.vp ,
                                      self.Fmass, self.uFlx, self.vFlx, self.Hz   , self.dt , self.nz  )
            else:
                self.u_np1[:] = self.u_np1_star[:]
                self.v_np1[:] = self.v_np1_star[:]
            #==========================================================
            # Compute eddy-viscosity / eddy-diffusivity  (TKE scheme)
            #==========================================================
            if self.ED_scheme=='TKE' : self.do_TKE( )
            if self.ED_scheme=='Keps': self.do_KEPS( )
            #==========================================================
            # Check for outputs & diagnostics
            #==========================================================
            if  time % self.outfreq == 0 and write_netcdf:
                self.do_turb_fluxes (  )
                # Compute mixed layer depth from bvf
                salt[:] = 38.5     #salt[:] = self.t_np1[:,self.isalt]
                rho,bvf_loc = scm_oce.rho_eos( self.t_np1[:,self.itemp], salt,self.z_r, self.z_w, self.rho0 , self.nz )
                self.hmxl10,self.hmxl300  = scm_oce.compute_mxl(bvf_loc,self.rhoc10,self.rhoc300,self.z_r,self.rho0,self.nz)
                kout+=1
                #=======================================
                # Write outputs
                self.output_state(TimeInSecs=time,kout=kout)
            # swap arrays
            self.u_n[:]   = self.u_np1[:]; self.v_n  [:] = self.v_np1[:]; self.t_n[:,:] = self.t_np1[:,:]
            self.wtke[:]  = 0.           ; self.tke_n[:] = self.tke_np1[:]
            if self.ED_scheme=='Keps': self.eps_n[:] = self.eps_np1[:]
            #
            if time % printMsg == 0:
                print('Day : ',int(time/printMsg))
#












#
    def do_turb_fluxes(self):
        # Turbulent fluxes diagnostics
        self.wted[  self.nz] = self.stflx[self.itemp]                            ## (w'theta')_ED
        self.wted[1:self.nz] = self.akt[1:self.nz] * (  self.t_np1[1:self.nz  ,self.itemp]
                                                      - self.t_np1[0:self.nz-1,self.itemp] ) / self.Hz[0:self.nz-1]
        self.wsed[  self.nz] = self.stflx[self.isalt]                   ## (w'S')_ED
        self.wsed[1:self.nz] = self.akt[1:self.nz] * (  self.t_np1[1:self.nz  ,self.isalt]
                                                      - self.t_np1[0:self.nz-1,self.isalt] ) / self.Hz[0:self.nz-1]
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
          self.wsmf[1:self.nz] =  self.Fmass[1:self.nz]*(
                                              self.tp   [1:self.nz,  self.isalt]
                                            - self.t_np1[0:self.nz-1,self.isalt] )
        if self.MF_dyn:                                                         ## (w'u')_MF
          self.wumf[1:self.nz] =  self.Fmass[1:self.nz]*(
                                               self.up   [1:self.nz]
                                             - self.u_np1[0:self.nz-1] )
          self.wvmf[1:self.nz] =  self.Fmass[1:self.nz]*(
                                               self.vp   [1:self.nz]
                                             - self.v_np1[0:self.nz-1] )
#










#
    def do_TKE(self):
        #==========================================================================
        # Compute Brunt-Vaisala frequency bvf for TKE production/destruction term
        #==========================================================================
        rho,self.bvf = scm_oce.rho_eos(
                                 self.t_np1_star[:,self.itemp], self.t_np1_star[:,self.isalt],
                                 self.z_r, self.z_w,self.rho0, self.nz )
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
        #=======================================
        self.shear = scm_tke.compute_shear(
                                  self.u_n  , self.v_n  ,
                                  self.u_np1, self.v_np1,
                                  self.u_np1_star, self.v_np1_star,
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
                                        self.z_r  , 0, self.wtke, self.nz  )
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
    def do_MF(self):
        #=======================================================
        # Surface boundary conditions for wp,up,vp and tracers
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
        if self.mass_flux_entr=='R10':
          self.ap,self.up,self.vp,self.wp,self.tp,self.Bp,self.ent,self.det, self.epsPlume = scm_mfc.mass_flux_r10(
                                    u_mean, v_mean, t_mean, self.tke_n, self.z_w, self.Hz,
                                    tp0   , up0   , vp0   , wp0     , self.mf_params,
                                    self.eos_params, tkepmin, mxlpmin, self.MF_small_ap, self.lineos, 0, self.zinv ,
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
        self.akv,self.akt,self.cmu,self.cmu_prim = scm_gls.compute_ev_ed_filt(
                        self.tke_np1, self.eps_np1, self.bvf,
                        self.shear2 , self.pnm, akvmin, aktmin, epsmin, self.nz )
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
        fh01.rho0            = self.rho0
        #self.min_Threshold
        fh01.tkemin          = self.min_Threshold[0]
        fh01.akvmin          = self.min_Threshold[1]
        fh01.aktmin          = self.min_Threshold[2]
        fh01.mxlmin          = self.min_Threshold[3]
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
        var = fh01.createVariable('hmxl_10m','f8',('time')); var[0] = self.hmxl10
        var = fh01.createVariable('hmxl_300m','f8',('time')); var[0] = self.hmxl300
        var  = fh01.createVariable('WT','f8',('time','z_w')); var[0,:] = self.wted[:]; var.units = 's-2'; del var
        var  = fh01.createVariable('WS','f8',('time','z_w')); var[0,:] = self.wsed[:]; var.units = 's-2'; del var
        if self.ED:
            var  = fh01.createVariable('tke','f8',('time','z_w')); var[0,:] = self.tke_n[:]; var.units = 'm2 s-2'; del var
            var  = fh01.createVariable('Akv','f8',('time','z_w')); var[0,:] = self.akv[:]; var.units = 'm2 s-1'; del var
            var  = fh01.createVariable('Akt','f8',('time','z_w')); var[0,:] = self.akt[:]; var.units = 'm2 s-1'; del var
            var  = fh01.createVariable('bvf','f8',('time','z_w')); var[0,:] = self.bvf[:]; var.units = 's-2'; del var
            var  = fh01.createVariable('lup','f8',('time','z_w')); var[0,:] = self.lupw[:]; var.units = 'm'; del var
            var  = fh01.createVariable('ldw','f8',('time','z_w')); var[0,:] = self.ldwn[:]; var.units = 'm'; del var
            var  = fh01.createVariable('WT_ED','f8',('time','z_w')); var[0,:] = self.wted[:]; var.units = 's-2'; del var
            var  = fh01.createVariable('WS_ED','f8',('time','z_w')); var[0,:] = self.wsed[:]; var.units = 'none'; del var
        if self.MF_tra:
            var = fh01.createVariable('a_p','f8',('time','z_w')); var[0,:] = self.ap[:]; del var
            var = fh01.createVariable('zinv','f8',('time')); var[0] = self.zinv
            var = fh01.createVariable('w_p','f8',('time','z_w')); var[0,:] = self.wp[:]; del var
            var = fh01.createVariable('B_p','f8',('time','z_w')); var[0,:] = self.buoyMF[:]; del var
            var = fh01.createVariable('temp_p','f8',('time','z_w')); var[0,:] = self.tp[:,self.itemp]; del var
            var = fh01.createVariable('salt_p','f8',('time','z_w')); var[0,:] = self.tp[:,self.isalt]; del var
            var = fh01.createVariable('Ent','f8',('time','z_r')); var[0,:] = self.ent[:]; del var
            var = fh01.createVariable('Det','f8',('time','z_r')); var[0,:] = self.det[:]; del var
            var = fh01.createVariable('WT_MF','f8',('time','z_w')); var[0,:] = self.wtmf[:]; del var
            var = fh01.createVariable('WS_MF','f8',('time','z_w')); var[0,:] = self.wsmf[:]; del var
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
        fh01.variables['taux'][kout]       = self.ustr_sfc
        fh01.variables['tauy'][kout]       = self.vstr_sfc
        fh01.variables['Qns'][kout]        = self.stflx[self.itemp]*self.cp*self.rho0
        fh01.variables['Qs'][kout]         = self.srflx*self.cp*self.rho0
        fh01.variables['Fw'][kout]         = self.stflx[self.isalt]
        fh01.variables['u'][kout,:]        = self.u_n[:]
        fh01.variables['v'][kout,:]        = self.v_n[:]
        fh01.variables['temp'][kout,:]     = self.t_n[:,self.itemp]
        fh01.variables['salt'][kout,:]     = self.t_n[:,self.isalt]
        fh01.variables['hmxl_10m'][kout]   = self.hmxl10
        fh01.variables['hmxl_300m'][kout]  = self.hmxl300
        fh01.variables['WT'][kout,:]       = self.wted[:]+self.wtmf[:]
        if self.ED:
            fh01.variables['tke'][kout,:] = self.tke_n[:]
            fh01.variables['Akv'][kout,:] = self.akv[:]
            fh01.variables['Akt'][kout,:] = self.akt[:]
            fh01.variables['bvf'][kout,:] = self.bvf[:]
            fh01.variables['lup'][kout,:] = self.lupw[:]
            fh01.variables['ldw'][kout,:] = self.ldwn[:]
            fh01.variables['WT_ED'][kout,:] = self.wted[:]
            fh01.variables['WS_ED'][kout,:] = self.wsed[:]
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
            fh01.variables['WS_MF'][kout,:] = self.wsmf[:]
            fh01.variables['u_p'][kout,:] = self.up[:]
            fh01.variables['v_p'][kout,:] = self.vp[:]
            if self.MF_tke:
                fh01.variables['buoyMF'][kout,:] = self.buoyMF[:]
                fh01.variables['shearMF'][kout,:] = self.shearMF[:]
            if self.MF_tke_trplCorr:
                fh01.variables['we'][kout,:] = self.triple_corr[:]
                fh01.variables['tke_p'][kout,:] = self.tp[:,self.isalt+1]
        fh01.close()

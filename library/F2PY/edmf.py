from scmoce_new  import scm_oce_new
from scmtke_new  import scm_tke_new
from scmmfc_new  import scm_mfc_new
import numpy as np

def do_TKE(self):
    #==========================================================================
    # Compute Brunt-Vaisala frequency bvf for TKE production/destruction term
    #==========================================================================
    tkemin = self.min_Threshold[0]; mxlmin = self.min_Threshold[3]
    #
    if self.lineos:
      rho,self.bvf = scm_oce_new.rho_eos_lin(
                             self.t_np1_star[:,self.itemp], self.t_np1_star[:,self.isalt],
                             self.z_r, self.eos_params,self.nz, len(self.eos_params) )
    else:
      rho,self.bvf = scm_oce_new.rho_eos(self.t_np1_star[:,self.itemp],self.t_np1_star[:,self.isalt],
                             self.z_r,self.z_w,self.rho0,self.nz)
    #===========================================================================
    # Compute boundary conditions for TKE
    #===========================================================================
    coeffBC = 0.
    if self.MF_tke_trplCorr: coeffBC = 1.
    #
    tke_sfc,tke_bot, flux_sfc = scm_tke_new.compute_tke_bdy(
                                  self.ustr_sfc,  self.vstr_sfc, self.ED_tke_const,
                                  coeffBC*self.bc_ap, coeffBC*self.wp0, tkemin)
    #
    if self.bdy_tke_sfc[0] < 0.5: self.bdy_tke_sfc[1] = tke_sfc
    else: self.bdy_tke_sfc[1] = flux_sfc
    self.bdy_tke_bot[1] = tke_bot;  self.bdy_tke_sfc[2] = tke_sfc
    self.bdy_tke_bot[2] = tke_bot
    #============================================================================
    # Compute TKE production by shear (no multiplication by akv here)
    #============================================================================
    self.shear = scm_tke_new.compute_shear(
                              self.u_n  , self.v_n  , self.u_np1, self.v_np1,
                              self.u_np1_star, self.v_np1_star , self.z_r  , self.nz  )
    #============================================================================
    # Transfer from mean PE/KE to TKE by mass flux term
    #============================================================================
    if self.MF_Etransfer_to_tke and self.MF_tra: self.buoyMF[:] = -self.Fmass[:]*self.Bp[:]  ## corresponds to  - FM x Bp
    if not self.MF_Etransfer_to_tke : self.shearMF[:] = 0.                                   ## If not needed, cancel the transfer from mean KE to TKE by mass flux on momentum
    self.Bprod = -self.akt*self.bvf; self.Bprod[0]=0.; self.Bprod[-1]=0.        ## Store ED-related buoyancy production term
    #============================================================================
    # Witek 2011 triple correlation term (opt = 0 -> 3D vector norm, opt = 2 -> 1D vector norm, opt = 1 -> return 0)
    # self.NormVel contains the 0.5 \| u - up \|^2 term
    #============================================================================
    if self.MF_tke_trplCorr:
        self.dWKdz_MF, self.NormVel = scm_mfc_new.compute_modifiedwitek11(
                                        self.Fmass, self.up   , self.vp ,self.wp  ,
                                        self.u_np1, self.v_np1, self.z_r ,
                                        self.triple_corr_opt,  self.wtke,self.nz  )
    #============================================================================
    # Advance TKE (w/o mass flux term)
    #============================================================================
    self.tke_np1_star, self.Prdtl, self.eps_np1, self.residual = scm_tke_new.tkeeq(
                                       self.tke_n  , self.bvf , self.shear, self.akv, self.akt   ,
                                       self.lupw   , self.ldwn, self.Hz   , self.z_r, self.buoyMF,
                                       self.shearMF, self.dWKdz_MF, self.wtke, self.dt, tkemin,
                                       self.bdy_tke_sfc, self.bdy_tke_bot, self.ED_tke_const  , self.nz  )
    #============================================================================
    # Advance TKE (with mass flux term)
    #============================================================================
    if self.MF_tke_trplCorr:
        self.tke_np1, self.tkep, self.epsp = scm_tke_new.tkeeq_mf(self.tke_np1_star,
                                                               self.ap  , self.wp  , self.Fmass  , self.NormVel,
                                                               self.Hz  , self.z_r,  self.wtke   , self.mf_params,
                                                               self.zinv, tkemin   , mxlmin   , self.MF_small_ap,
                                                               self.triple_corr_opt, self.dt, len(self.mf_params), self.nz )
    else:
        self.tke_np1 = self.tke_np1_star
    #=============================================================================
    # Finalize eddy-viscosity/diffusivity computation
    #=============================================================================
    self.lupw,self.ldwn = scm_tke_new.compute_mxl(
                                     self.tke_np1 , self.bvf     , self.Hz ,
                                     self.ustr_sfc, self.vstr_sfc, mxlmin  , self.nz )

    self.akv,self.akt   = scm_tke_new.compute_ed (
                                     self.tke_np1 , self.lupw    , self.ldwn, self.Prdtl ,
                                     self.ED_extrap_sfc, self.ED_tke_const, self.nz  )
    #==============================================================================
    # Apply EVD if necessary
    #==============================================================================
    if self.ED_evd: scm_oce_new.compute_evd(
                                      self.bvf, self.akv, self.akt,
                                      10., self.nz  )
#
#
#
#
#
#
#
#
#
def do_MF(self):
    #===============================================================================
    # Surface boundary conditions for up,vp and tracers
    #===============================================================================
    up0,vp0,tp0 = scm_mfc_new.compute_mf_bdy(
                            self.u_np1_star[-2:]  , self.v_np1_star[-2:],
                            self.t_np1_star[-2:,:], self.Hz[-2:], self.ntra, 2 )
    wp0=self.wp0
    # Pergaud 2009 boundary conditions
    if self.bc_P09 != 'false':
        SqrTKE = np.sqrt(self.tke_n[-1]) # surface value of TKE
        self.Tmean_sfc  = tp0[self.itemp]
        self.Smean_sfc  = tp0[self.isalt]
        wp0    = -np.sqrt(0.66666667)*SqrTKE; self.wp0 = wp0
        #
        tp0[self.itemp]    = self.Tmean_sfc + self.beta_bc_P09*self.stflx[self.itemp]/SqrTKE
        tp0[self.isalt]    = self.Smean_sfc + self.beta_bc_P09*self.stflx[self.isalt]/SqrTKE
    #================================================================================
    # Compute mass flux variables (Rio et al. 2010 closure)
    #================================================================================
    if self.trad_coriolis_mod:
        # MODULATION due to traditional rotation terms
        B0 = self.g*self.alpha*self.stflx[self.itemp]
        Ro = np.minimum((np.abs(B0)/self.fcor)**(0.5)/(self.fcor*np.abs(self.zinv)) , (np.abs(B0)/self.fcor)**(0.5)/(self.fcor*1000.))
        # cff= np.tanh(Ro**0.3) #modulation coefficient, Wang 2006
        cff= np.tanh(Ro**0.37) #modulation coefficient, Wang 2006
        # cff=1
        # reduce integral lenth scale: L_int = cff*zinv (division by zinv is done in the fortran routine)
        #self.mf_params  = np.array([Cent,Cdet,wp_a,wp_b,wp_bp,up_c,vp_c,bc_ap,delta_bkg])
        self.mf_params[-1] = self.delta_bkg/cff
        self.mf_params[4] =      self.wp_bp/cff
        # # reduce lateral exchanges
        self.mf_params[0] = self.Cent/cff
        self.mf_params[1] = self.Cdet/cff
        # self.mf_params[2] = self.wp_a*cff
        # self.mf_params[3] = self.wp_b*cff
        ######################
    if self.mass_flux_entr=='R10' or self.mass_flux_entr=='R10_HB09' or self.mass_flux_entr=='R10_Wi11':
        self.ap,self.up,self.vp,self.wp,self.tp,self.Bp,self.ent,self.det, self.vortp = scm_mfc_new.mass_flux_r10(
                                    self.u_np1_star, self.v_np1_star, self.t_np1_star, self.z_w, self.Hz,
                                    tp0   , up0   , vp0   , wp0     , self.mf_params,
                                    self.eos_params, self.MF_small_ap, self.lineos, self.zinv, self.fcor,
                                    self.nz , self.ntra , len(self.mf_params), len(self.eos_params)   )
    #
    self.zinv = min(self.zinvMin,self.zinv)
    #
    if self.bc_P09 != 'false':
        self.tp[-1,self.itemp] = self.t_np1_star[-1,self.itemp] + self.beta_bc_P09*self.stflx[self.itemp]/SqrTKE
    else:
        self.tp[-1,self.itemp] = self.t_np1_star[-1,self.itemp]


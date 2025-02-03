from scmoce_new  import scm_oce_new
from scmkeps_new import scm_keps_new
#
def do_KEPS(self):
    #==========================================================================
    # Compute Brunt-Vaisala frequency bvf for TKE production/destruction term
    #==========================================================================
    if self.lineos:
        rho,self.bvf = scm_oce_new.rho_eos_lin(
                             self.t_np1[:,self.itemp], self.t_np1[:,self.isalt],
                             self.z_r, self.eos_params,self.nz, len(self.eos_params) )
    else:
        rho,self.bvf = scm_oce_new.rho_eos(self.t_np1[:,self.itemp],self.t_np1[:,self.isalt],
                             self.z_r,self.z_w,self.rho0,self.nz)
    #==========================================================================
    # Compute TKE production by shear (Not multiplied by Akv !!)
    #==========================================================================
    self.shear = scm_keps_new.compute_shear(
                              self.u_n, self.v_n, self.u_np1,
                              self.v_np1, self.z_r  , self.nz  )
    #==========================================================================
    # Compute boundary conditions
    #==========================================================================
    tke1 = self.tke_n[self.nz-1]
    tke2 = self.tke_n[        1]
    tke_sfc,tke_bot,ftke_sfc,ftke_bot,eps_sfc,eps_bot,feps_sfc,feps_bot = scm_keps_new.compute_tke_eps_bdy(
                self.ustr_sfc, self.vstr_sfc, self.ustr_bot, self.vstr_bot,
                self.z0b, self.z0s, tke1, self.Hz[-1], tke2, self.Hz[0], 
                self.OneOverSig_psi, self.tkemin, self.epsmin )
    #
    if self.bdy_tke_sfc[0] < 0.5: self.bdy_tke_sfc[1] = tke_sfc  ## Dirichlet BC for surface TKE
    else: self.bdy_tke_sfc[1] = ftke_sfc                         ## Neumann BC for surface TKE
    if self.bdy_tke_bot[0] < 0.5: self.bdy_tke_bot[1] = tke_bot  ## Dirichlet BC for bottom TKE
    else: self.bdy_tke_bot[1] = ftke_bot                         ## Neumann BC for bottom TKE
    #
    if self.bdy_eps_sfc[0] < 0.5: self.bdy_eps_sfc[1] = eps_sfc  ## Dirichlet BC for surface EPS
    else: self.bdy_eps_sfc[1] = feps_sfc                         ## Neumann BC for surface EPS
    if self.bdy_eps_bot[0] < 0.5: self.bdy_eps_bot[1] = eps_bot  ## Dirichlet BC for bottom EPS
    else: self.bdy_eps_bot[1] = feps_bot                         ## Neumann BC for bottom EPS
    #
    self.bdy_tke_sfc[2] = tke_sfc; self.bdy_tke_bot[2] = tke_bot ## Provide the Dirichlet values anyway 
    self.bdy_eps_sfc[2] = eps_sfc; self.bdy_eps_bot[2] = eps_bot ## Provide the Dirichlet values anyway
    #==========================================================================
    # Advance TKE equation 
    #==========================================================================
    self.tke_np1[:], self.wtke[:] = scm_keps_new.tkeeq(self.tke_n,self.bvf, self.shear,
                       self.akv,self.akt, self.eps_n,
                       self.Hz,self.dt,self.tkemin,
                       self.bdy_tke_sfc,self.bdy_tke_bot,self.nz)
    #==========================================================================
    # Advance dissipation equation
    #==========================================================================
    self.eps_np1[:] = scm_keps_new.dissipationeq(self.eps_n,self.bvf, self.shear,
                       self.akv,self.akt, self.tke_n,
                       self.Hz,self.dt,self.epsmin,
                       self.bdy_eps_sfc,self.bdy_eps_bot,self.nz)
    #=========================
    #self.akv,self.akt,self.cmu,self.cmu_prim,self.lscale = scm_keps_new.compute_ev_ed_eq(
    #                self.tke_np1, self.eps_np1, self.bvf,  self.nz )
    #==========================================================================
    # Compute stability functions and finalize the computation
    #==========================================================================
    self.akv,self.akt,self.cmu,self.cmu_prim,self.lscale = scm_keps_new.compute_ev_ed_weak(
                    self.tke_np1, self.eps_np1, self.shear, self.bvf,  self.nz )
    #

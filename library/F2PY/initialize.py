import numpy as np

def validate_config(self):
    """Validate and normalize the configuration."""
    # Check for eddy diff scheme
    ed  = ['TKE','Keps']
    cst = ['NEMO','MNH','RS81']
    self.eddy_diff = self.config["eddy_diff"]
    if self.eddy_diff:
        if self.config["eddy_diff_scheme"] in ed:
            print("Eddy diffusion scheme: ",self.config["eddy_diff_scheme"])
            self.ED_scheme = self.config["eddy_diff_scheme"]
        else:
            raise ValueError(f"Wrong eddy diffusion scheme (either TKE or Keps)")
        if self.ED_scheme=="TKE":
            if self.config["eddy_diff_tke_const"] in cst:
                print(f"TKE scheme constants: ",self.config["eddy_diff_tke_const"])
            else:
                raise ValueError(f"Wrong choice of tke constants (either NEMO, MNH, or RS81)")
    # Check for mass flux scheme
    ent = ['R10','R10_HB09','R10_Wi11']
    if self.config["mass_flux_tra"]:
        print("Mass flux scheme is activated")
        if self.config["entr_scheme"] in ent:
            print(f"Entrainement scheme: ",self.config["entr_scheme"])
        else:
            raise ValueError(f"Wrong choice of entrainment constants (either R10, R10_HB09, or R10_Wi11)")
    # Check for grid type
    grd = ['croco_new','ORCA75','croco_old']
    if self.config["gridType"] in grd:
        print(f"Vertical grid type: ",self.config["gridType"])
    else:
        raise ValueError(f"Wrong choice of vertical grid type (either croco_new, ORCA75, or croco_old)")

def initialize_physical_constants(self):
    """Initialize physical constants."""
    omega = 7.292116e-5; rad = np.pi / 180.0
    self.fcor   = 2.0 * omega * np.sin(rad * self.config["lat0"])
    self.ecor   = 2.0 * omega * np.cos(rad * self.config["lat0"])
    self.rho0   = self.config["rho0"]
    self.cp     = self.config["cpoce"]
    self.alpha  = self.config["alpha"]
    self.beta   = self.config["beta"]
    self.g      = 9.81
    self.lineos = self.config["lin_eos"]
    self.Tref   = self.config["Tref"]
    self.Sref   = self.config["Sref"]
    self.eos_params  = np.array([self.rho0,self.alpha,self.beta,self.Tref,self.Sref])



def initialize_boundary_conditions(self):
    """Initialize boundary boundary conditions."""
    self.ustr_sfc = self.config["sustr"]     ## zonal wind stress       [m2/s2]
    self.vstr_sfc = self.config["svstr"]     ## meridional wind stress  [m2/s2]
    self.stflx    = np.zeros(2)
    cff = 1.0 / (self.rho0 * self.cp)
    self.stflx[self.itemp] = self.config["stflx"]*cff  # Heat flux
    self.srflx    = self.config["srflx"]*cff  # Solar radiation
    self.stflx[self.isalt] = self.config["ssflx"]   # Freshwater flux
    self.btflx    = np.zeros(2)
    self.ustr_bot = 0.
    self.vstr_bot = 0.
    self.z0b      = 1.e-14
    self.DC       = self.config["diurnal_cycle"]
    self.Qswmax   = self.srflx
    self.bc_P09   = self.config["bc_P09"]
    self.beta_bc_P09 = self.config["beta_bc_P09"]



def initialize_eddy_diffusion(self):
    """Initialize eddy-diffusion parameters."""
    self.ED_evd        = self.config["evd"]
    self.akvmin        = self.config["akvmin"]
    self.aktmin        = self.config["aktmin"]
    self.ED_extrap_sfc = self.config["extrap_ak_surf"]
    self.tkemin        = 1.e-08
    self.epsmin        = 1.e-12
    self.min_Threshold = np.array([self.tkemin,self.akvmin,self.aktmin,self.config["mxlmin"]])
    if self.ED_scheme=='TKE':
        self.ED_tke_const      = 0  ; self.inv_schmidt     = 1.; cm = 0.1
        if self.config["eddy_diff_tke_const"]=='MNH':
            self.ED_tke_const    = 1; cm = 0.126; self.inv_schmidt = 0.34/cm
        if self.config["eddy_diff_tke_const"]=='RS81':
            self.ED_tke_const    = 2; cm = 0.0667; self.inv_schmidt= 0.4/cm
        self.bdy_tke_sfc    = np.zeros(3)
        self.bdy_tke_bot    = np.zeros(3)
        self.bdy_tke_sfc[0] = 1.               # Neumann BC for TKE
        if self.config["tke_sfc_dirichlet"]: self.bdy_tke_sfc[0] = 0. # Dirichlet BC for TKE
        self.bdy_tke_bot[0] = 0.
    elif self.ED_scheme=='Keps':
        Sch_k = 1.; Sch_eps = 1.3  # Schmidt numbers
        self.OneOverSig_k   = 1./Sch_k; self.OneOverSig_psi = 1./Sch_eps
        # boundary conditions for k-eps scheme
        self.bdy_tke_sfc = np.zeros(3); self.bdy_tke_bot = np.zeros(3)
        self.bdy_eps_sfc = np.zeros(3); self.bdy_eps_bot = np.zeros(3)
        # Neumann B.C. (default option in GOTM)
        self.bdy_tke_sfc[0] = 1.;self.bdy_eps_sfc[0] = 1.
        self.bdy_tke_bot[0] = 1.;self.bdy_eps_bot[0] = 1.



def initialize_mass_flux(self):
    """Initialize mass flux options."""
    self.MF_tra              = self.config["mass_flux_tra"]
    self.MF_dyn              = self.config["mass_flux_dyn"]
    self.MF_Etransfer_to_tke = self.config["mass_flux_tke"]
    self.MF_tke_trplCorr     = self.config["mass_flux_tke_trplCorr"]
    self.MF_small_ap         = self.config["mass_flux_small_ap"]
    self.mass_flux_entr      = self.config["entr_scheme"]
    if self.mass_flux_entr == 'R10': self.triple_corr_opt = 0
    if self.mass_flux_entr == 'R10_HB09': self.triple_corr_opt = 1
    if self.mass_flux_entr == 'R10_Wi11': self.triple_corr_opt = 2
    self.mf_params  =np.array([self.config["Cent"] ,self.config["Cdet"],
                               self.config["wp_a"] ,self.config["wp_b"],
                               self.config["wp_bp"],self.config["up_c"],
                               self.config["vp_c"] ,self.config["bc_ap"],
                               self.config["delta_bkg"]])
    self.wp0                 = self.config["wp0"]
    self.bc_ap               = self.config["bc_ap"]
    self.delta_bkg           = self.config["delta_bkg"]
    self.wp_bp               = self.config["wp_bp"]
    self.Cent                = self.config["Cent"]
    self.Cdet                = self.config["Cdet"]
    self.wp_a                = self.config["wp_a"]
    self.wp_b                = self.config["wp_b"]
    self.trad_coriolis_mod   = self.config["trad_coriolis_mod"]


def initialize_vertical_grid(self):
    """Define the vertical grid."""
    gtype = self.config["gridType"]
    if gtype=='croco_old':
        self.z_r    = np.zeros(self.nz); self.z_w    = np.zeros(self.nz+1)
        self.z_w[0] = -self.config["h0"]
        ds  = 1./self.nz
        cff = (self.config["h0"]-self.config["hc"])/np.sinh(self.config["thetas"])
        for k in range(self.nz,0,-1):
            sc_w     = ds * float(k-self.nz);     self.z_w[k  ] = self.config["hc"]*sc_w + cff*np.sinh(self.config["thetas"]*sc_w)
            sc_r     = ds*(float(k-self.nz)-0.5); self.z_r[k-1] = self.config["hc"]*sc_r + cff*np.sinh(self.config["thetas"]*sc_r)
    if gtype=='ORCA75':
        # construct the ORCA 75 levels grid and extract the levels between 0 and h0
        self.nz   = 75
        zsur = -3958.95137127683; za2  = 100.760928500000; za0  =  103.953009600000; za1  = 2.41595126900000
        zkth =  15.3510137000000; zkth2= 48.0298937200000; zacr =  7.00000000000000; zacr2= 13.0000000000000
        Sc_r = np.arange( self.nz-0.5, 0.5 , -1.); Sc_w = np.arange( self.nz    , 0.  , -1.)
        z_w  = -( zsur + za0 * Sc_w + za1 * zacr * np.log(np.cosh((Sc_w-zkth )/zacr) ) + za2 * zacr2* np.log(np.cosh( (Sc_w-zkth2) / zacr2 ) )  )
        z_r  = -( zsur + za0 * Sc_r + za1 * zacr * np.log(np.cosh((Sc_r-zkth )/zacr) ) + za2 * zacr2* np.log(np.cosh( (Sc_r-zkth2) / zacr2 ) )  )
        nbot = np.argmin(z_w <= -h0)
        self.z_w     = z_w[nbot-1:]; self.z_w[-1] = 0.
        self.nz      = len(self.z_w) - 1
        self.z_r     = z_r[nbot-1:]
        #
    if gtype=='croco_new':
        Sc_r  = np.arange(-self.nz+0.5, 0.5, 1) / float(self.nz)
        Sc_w  = np.arange(-self.nz, 1, 1) / float(self.nz)
        Cs_r = (1.-np.cosh(self.config["thetas"]*Sc_r))/(np.cosh(self.config["thetas"])-1.)
        Cs_w = (1.-np.cosh(self.config["thetas"]*Sc_w))/(np.cosh(self.config["thetas"])-1.)
        self.z_w   = (self.config["hc"]*Sc_w+Cs_w*self.config["h0"])*self.config["h0"]/(self.config["h0"]+self.config["hc"])
        self.z_r   = (self.config["hc"]*Sc_r+Cs_r*self.config["h0"])*self.config["h0"]/(self.config["h0"]+self.config["hc"])
        #
    self.Hz    = self.z_w[1:]-self.z_w[:-1]



def initialize_initial_state(self):
    """Define the initial state."""
    self.u_n = np.zeros(self.nz); self.u_np1 = np.zeros(self.nz)
    self.v_n = np.zeros(self.nz); self.v_np1 = np.zeros(self.nz)
    self.u_np1_star = np.zeros(self.nz)
    self.v_np1_star = np.zeros(self.nz)
    self.t_n        = np.zeros((self.nz,self.ntra), order='F')
    self.t_np1      = np.zeros((self.nz,self.ntra), order='F')
    self.t_np1_star = np.zeros((self.nz,self.ntra), order='F')
    # constant salinity
    self.t_n       [:,self.isalt] = self.config["SaltCst"]
    self.t_np1     [:,self.isalt] = self.t_n[:,self.isalt]
    self.t_np1_star[:,self.isalt] = self.t_n[:,self.isalt]
    # temperature
    self.T0      = self.config["T0"];  self.N0      = self.config["N0"];
    strat        = self.N0 / (self.g * self.alpha )
    self.t_n  [:,self.itemp] = self.T0
    self.t_np1[:,self.itemp] = self.T0
    for k in range(0, self.nz):
        if self.z_w[k] < self.config["mld_ini_temp"]:
            self.t_n  [k,self.itemp]      = self.T0 + 0.5*strat*( self.z_w[k+1]+self.z_w[k] ) + strat*(-self.config["mld_ini_temp"])
            self.t_np1[k,self.itemp]      = self.t_n  [k,self.itemp]
            self.t_np1_star[k,self.itemp] = self.t_np1[k,self.itemp]
    self.btflx[self.itemp]  = 0.
    if self.config["btflx"]=='linear_continuation': self.btflx[self.itemp] = self.aktmin * self.N0 / (self.g * self.alpha )
    self.btflx[self.isalt]  = 0.



def initialize_ed_arrays(self):
    self.akt      = np.zeros(self.nz+1); self.akt[:]     = self.aktmin
    self.akv      = np.zeros(self.nz+1); self.akv[:]     = self.akvmin
    self.bvf      = np.zeros(self.nz+1); self.eps_n      = np.zeros(self.nz+1)
    self.Bprod    = np.zeros(self.nz+1); self.eps_np1    = np.zeros(self.nz+1)
    if self.eddy_diff:
        self.tke_n   = np.zeros(self.nz+1); self.tke_n[:]   = self.tkemin
        self.tke_np1 = np.zeros(self.nz+1); self.tke_np1[:] = self.tkemin
        self.tke_np1_star = np.zeros(self.nz+1); self.tke_np1_star[:] = self.tkemin
        self.lupw    = np.zeros(self.nz+1); self.lupw[:]    = self.config["mxlmin"]
        self.ldwn    = np.zeros(self.nz+1); self.ldwn[:]    = self.config["mxlmin"]
        self.Prdtl   = np.zeros(self.nz+1)
        self.eps_n[:] = self.epsmin;
        self.eps_np1[:] = self.epsmin
        self.shear        = np.zeros(self.nz+1)
        self.residual   = np.array(0.)
        if self.ED_scheme=='Keps':
            self.cmu        = np.zeros(self.nz+1); self.cmu     [:] = 0.1
            self.cmu_prim   = np.zeros(self.nz+1); self.cmu_prim[:] = 0.1
            self.z0s        = 0.02
    self.buoyMF      = np.zeros(self.nz+1); self.shearMF     = np.zeros(self.nz+1)
    self.wtke        = np.zeros(self.nz  ) #ED+MF wtke
    self.dWKdz_MF    = np.zeros(self.nz+1)



def initialize_mf_arrays(self):
    self.zinvMin = np.array(-self.Hz[-1])
    self.zinv   = self.zinvMin
    self.ap     = np.zeros(self.nz+1); self.wp       = np.zeros(self.nz+1)
    self.tp     = np.zeros((self.nz+1, self.ntra), order='F');
    self.Bp     = np.zeros(self.nz+1); self.Fmass    = np.zeros(self.nz+1)
    self.ent    = np.zeros(self.nz  ); self.det      = np.zeros(self.nz  )
    self.up     = np.zeros(self.nz+1); self.vp       = np.zeros(self.nz+1)
    self.tkep   = np.zeros(self.nz+1); self.tkep[:]  = self.tkemin
    self.epsp   = np.zeros(self.nz);
    self.vortp  = np.zeros(self.nz+1); self.NormVel  = np.zeros(self.nz  )

def initialize_output_and_averaging(self):
    self.output       = self.config["output_filename"]
    self.write_netcdf = self.config["write_netcdf"]
    self.do_avg       = self.config["avg_last_hour"]
    self.start_avg    = (self.config["nbhours"]-1.)*3600.
    dimt = int( self.config["nbhours"] // self.config["outfreq"] )
    self.t_history = np.zeros((self.nz,dimt), order='F')
    self.u_history = np.zeros((self.nz,dimt), order='F')
    self.v_history = np.zeros((self.nz,dimt), order='F')
    self.tke_history = np.zeros((self.nz+1,dimt), order='F')
    if self.write_netcdf:
        self.wted   = np.zeros(self.nz+1); self.wtmf   = np.zeros(self.nz+1)
        self.wued   = np.zeros(self.nz+1); self.wumf   = np.zeros(self.nz+1)
        self.wved   = np.zeros(self.nz+1); self.wvmf   = np.zeros(self.nz+1)
        self.wbed   = np.zeros(self.nz+1); self.wbmf   = np.zeros(self.nz+1)
        self.wted[  self.nz] = self.stflx[self.itemp]
    if self.do_avg:
        self.u_avg      = np.zeros(self.nz)
        self.v_avg      = np.zeros(self.nz)
        self.t_avg      = np.zeros((self.nz,self.ntra), order='F')
        self.tke_avg    = np.zeros(self.nz+1)
        self.wt_avg     = np.zeros(self.nz+1)
        self.akv_avg    = np.zeros(self.nz+1)
        self.akt_avg    = np.zeros(self.nz+1)
        self.navg       = 0.

def initialize_diagnostics(self):
    # MLD computation params
    self.rhoc      = 0.01    # density criterion for mixed layer depth (consistent with NEMO)
    self.hmxl      = np.array(-self.Hz[-1])
    self.tFlx      = np.zeros(self.nz+1); self.tFlx[-1] = - self.stflx[self.itemp]
    self.uFlx      = np.zeros(self.nz+1)
    self.vFlx      = np.zeros(self.nz+1)
    self.vint_Etot = np.array(1e-14)
    self.vint_Ekin = np.array(1e-14)
    self.vint_Epot = np.array(1e-14)
    self.vint_TKE  = np.array(1e-14)
    self.vint_Eps  = np.array(1e-14)

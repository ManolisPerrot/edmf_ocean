MODULE scm_oce
  !!============================================================================<br />
  !!                       ***  MODULE  scm_oce  ***                            <br />
  !! Single Column Model (SCM) time-stepping and oceanic equation of state      <br />
  !!============================================================================<br />
  !!----------------------------------------------------------------------------<br />
  !!   advance_tra_ED     : integrate eddy-diffusion term for tracers           <br />
  !!   advance_tra_MF     : integrate mass-flux term for tracers                <br />
  !!   advance_dyn_MF     : integrate mass-flux term for dynamics               <br />
  !!   advance_dyn_Cor_ED : integrate eddy-viscosity and Coriolis terms for dynamics <br />
  !!   compute_evd        : compute Enhanced Vertical Diffusion (EVD)           <br />
  !!   lmd_swfrac         : compute fraction of solar penetration               <br />
  !!   rho_eos_lin        : linear equation of state                            <br />
  !!   tridiag_solve      : tridiagonal solver                                  <br />
  !!----------------------------------------------------------------------------<br />
    IMPLICIT NONE

CONTAINS
  !===================================================================================================
  SUBROUTINE advance_tra_ED(t_n,stflx,srflx,swr_frac,btflx,Hz,Akt,zw,eps,alpha,dt,N,ntra,t_np1)
  !---------------------------------------------------------------------------------------------------
    !!============================================================================<br />
    !!                  ***  ROUTINE advance_tra_ED  ***                          <br />
    !! ** Purposes : integrate vertical diffusion term for tracers                <br />
    !!============================================================================<br />
    !! \[ \overline{\phi}^{n+1,*} = \overline{\phi}^n + \Delta t \partial_z \left(  K_m \partial_z  \overline{\phi}^{n+1,*} \right) \]
    USE scm_par, ONLY: cp, grav
    IMPLICIT NONE
    INTEGER, INTENT(IN   )          :: N                !! number of vertical levels
    INTEGER, INTENT(IN   )          :: ntra             !! number of tracers to integrate
    REAL(8), INTENT(IN   )          :: dt               !! time-step [s]
    REAL(8), INTENT(IN   )          :: t_n  (1:N,ntra)  !! tracer at time step n
    REAL(8), INTENT(IN   )          :: stflx(  ntra)    !! surface tracer fluxes
    REAL(8), INTENT(IN   )          :: btflx(  ntra)    !! surface tracer fluxes
    REAL(8), INTENT(IN   )          :: srflx            !! surface radiative flux [W/m2]
    REAL(8), INTENT(IN   )          :: swr_frac(0:N)    !! fraction of solar penetration
    REAL(8), INTENT(IN   )          :: Hz      (1:N)    !! layer thickness [m]
    REAL(8), INTENT(IN   )          :: Akt     (0:N)    !! eddy-diffusivity [m2/s]
    REAL(8), INTENT(IN   )          :: eps     (0:N)    !! TKE dissipation [m2/s3]
    REAL(8), INTENT(IN   )          :: zw      (0:N)    !! depth at cell interfaces [m]
    REAL(8), INTENT(IN   )          :: alpha            !! thermal expension coefficient [C-1]
    REAL(8), INTENT(  OUT)          :: t_np1(1:N,ntra)  !! tracer at time step n+1
    ! local variables
    INTEGER                         :: k,itrc
    REAL(8)                         :: FC(0:N), ff(1:N),cffp,cffm
    !
    DO itrc=1,ntra
      FC(0) = 0.; FC(N) = 0.
      !=======================================================================
      !! 1 - Compute fluxes associated to solar penetration and surface bdy condition <br />
      !=======================================================================
      IF (itrc.eq.1) then
        FC(  N  ) = stflx(itrc)+srflx      !<-- surface heat flux (including latent and solar components)
        FC(1:N-1) = srflx*swr_frac(1:N-1)    !<-- penetration of solar heat flux
      ELSE
        FC(  N  ) = stflx(itrc)            !<-- salinity (fresh water flux)
        FC(1:N-1) = 0.D0
      ENDIF
      ! Apply flux divergence
      DO k=1,N
        t_np1(k,itrc)=Hz(k)*t_n(k,itrc)+dt*(FC(k  )-FC(k-1))
      ENDDO
      IF(itrc.eq.1) then
        DO k=1,N
          cffp  = eps(k  ) / ( cp-alpha*grav*zw(k  ) )
          cffm  = eps(k-1) / ( cp-alpha*grav*zw(k-1) )
          t_np1(k,itrc)=t_np1(k,itrc)+dt*0.5*Hz(k)*(cffp+cffm)
        ENDDO
      ENDIF
      !=======================================================================
      !! 2 - Implicit integration for vertical diffusion <br />
      !=======================================================================
      ! right hand side for the tridiagonal problem
      ff(1:N) = t_np1(1:N,itrc)
      ff(1  ) = ff(1) - dt*btflx(itrc)
      ! solve tridiagonal problem
      CALL tridiag_solve(N,Hz,Akt,ff,dt)
      ! update tracer
      t_np1(1:N,itrc) = ff(1:N)
      !-----------
    ENDDO
  !---------------------------------------------------------------------------------------------------
  END SUBROUTINE advance_tra_ED
  !===================================================================================================



  !===================================================================================================
  SUBROUTINE advance_tra_MF(t_np1,t_p,Fmass,Hz,dt,N,ntra)
  !---------------------------------------------------------------------------------------------------
    !!============================================================================<br />
    !!                  ***  ROUTINE advance_tra_MF  ***                          <br />
    !! ** Purposes : integrate mass flux term for tracers                         <br />
    !!============================================================================<br />
    USE scm_par
    IMPLICIT NONE
    INTEGER, INTENT(IN   )          :: N                  !! number of vertical levels
    INTEGER, INTENT(IN   )          :: ntra               !! number of tracers to integrate
    REAL(8), INTENT(IN   )          :: dt                 !! time-step [s]
    REAL(8), INTENT(IN   )          :: t_p  (0:N,ntra+1)  !! tracer properties in the plume
    REAL(8), INTENT(IN   )          :: Hz      (1:N)      !! layer thickness [m]
    REAL(8), INTENT(IN   )          :: Fmass (0:N       ) !! mass flux [m/s]
    REAL(8), INTENT(INOUT)          :: t_np1 (1:N,ntra  ) !! tracer at time step n+1
    ! local variables
    INTEGER                         :: k,itrc
    REAL(8)                         :: cff, FC(0:N)
    !
    DO itrc=1,ntra
      FC(0) = 0. !; FC(N) = 0.
      ! Compute fluxes associated to mass flux
      DO k = 1,N
        FC(k) = Fmass(k)*(t_p(k,itrc)-t_np1(k,itrc))  !! Compute fluxes associated to mass flux \( F_{k+1/2}^{\rm MF} = (a^{\rm p}w^{\rm p})_{k+1/2}\left( \phi^{\rm p}_{k+1/2} - \overline{\phi}_k^{n+1,\star} \right)  \) <br />
      ENDDO
      FC(N) = 0.
      ! Apply flux divergence
      DO k = 1,N
        t_np1(k,itrc)=t_np1(k,itrc)+dt*(FC(k)-FC(k-1))/Hz(k) !! \(   \overline{\phi}_k^{n+1} = \overline{\phi}_k^{n+1,\star} + \frac{\Delta t}{\Delta z_k} ( F_{k+1/2}^{\rm MF}-F_{k-1/2}^{\rm MF} )   \)
      ENDDO
    ENDDO
    !
  !---------------------------------------------------------------------------------------------------
  END SUBROUTINE advance_tra_MF
  !===================================================================================================

  !===================================================================================================
  SUBROUTINE advance_dyn_MF(u_np1,v_np1,shear2_MF,u_n,v_n,u_p,v_p,Fmass,Hz,dt,N)
  !---------------------------------------------------------------------------------------------------
    !!==========================================================================<br />
    !!                  ***  ROUTINE advance_dyn_MF  ***                        <br />
    !! ** Purposes : integrate mass flux term for dynamics and compute transfer
    !!                         of KE to TKE related to the mass flux            <br />
    !!==========================================================================<br />
    !! \[ \mathbf{u}^{n+1} = \mathbf{u}^{n+1,\star \star} - \Delta t \partial_z\left( a^{\rm p} w^{\rm p} \left( \mathbf{u}^{\rm p} - \mathbf{u}^{n+1,\star \star}  \right)   \right) \]
    !! Compute the shear term associated with mass flux
    !! \[ {\rm Sh}^{\rm p} = - a^{\rm p} w^{\rm p} \left( \mathbf{u}^{\rm p} - \mathbf{u}^{n+1,\star \star}  \right) \cdot \frac{1}{2} \left( \partial_z \mathbf{u}^n + \partial_z \mathbf{u}^{n+1}  \right) \]
    USE scm_par
    IMPLICIT NONE
    INTEGER, INTENT(IN   )          :: N                  !! number of vertical levels
    REAL(8), INTENT(IN   )          :: dt                 !! time-step [s]
    REAL(8), INTENT(INOUT)          :: u_np1 (1:N )       !! u-velocity component at time n+1 [m/s]
    REAL(8), INTENT(INOUT)          :: v_np1 (1:N )       !! v-velocity component at time n+1 [m/s]
    REAL(8), INTENT(INOUT)          :: shear2_MF(0:N)     !! TKE production term associated to mass flux [m2/s3]
    REAL(8), INTENT(IN   )          :: u_n   (1:N )       !! u-velocity component at time n [m/s]
    REAL(8), INTENT(IN   )          :: v_n   (1:N )       !! v-velocity component at time n [m/s]
    REAL(8), INTENT(IN   )          :: u_p   (0:N )       !! u-velocity component in the plume [m/s]
    REAL(8), INTENT(IN   )          :: v_p   (0:N )       !! v-velocity component in the plume [m/s]
    REAL(8), INTENT(IN   )          :: Hz      (1:N)      !! layer thickness [m]
    REAL(8), INTENT(IN   )          :: Fmass (0:N       ) !! mass flux [m/s]
    ! local variables
    INTEGER                         :: k
    REAL(8)                         :: cff, FCu(0:N), FCv(0:N)
    REAL(8)                         :: cffu,cffv
    !
    FCu(0) = 0. ; FCv(0) = 0.
    ! Compute fluxes associated to mass flux
    DO k = 1,N
      FCu(k) = Fmass(k)*(u_p(k)-u_np1(k))
      FCv(k) = Fmass(k)*(v_p(k)-v_np1(k))
    ENDDO
    FCu(N) = 0. ; FCv(N) = 0.
    ! Apply flux divergence
    DO k = 1,N
      u_np1(k)=u_np1(k)+dt*(FCu(k)-FCu(k-1))/Hz(k)
      v_np1(k)=v_np1(k)+dt*(FCv(k)-FCv(k-1))/Hz(k)
    ENDDO
    ! Compute shear term for TKE
    shear2_MF(0) = 0.
    shear2_MF(N) = 0.
    !
    DO k = 1,N-1   ! use upwind value for u_np1/v_np1
          cffu = Fcu(k)*0.5*( (u_np1(k+1)+u_n(k+1))-(u_np1(k)+u_n(k)) )
          cffv = FCv(k)*0.5*( (v_np1(k+1)+v_n(k+1))-(v_np1(k)+v_n(k)) )
          shear2_MF(k) = 2.*(cffu+cffv)/(Hz(k+1)+Hz(k))
    ENDDO
  !---------------------------------------------------------------------------------------------------
  END SUBROUTINE advance_dyn_MF
  !===================================================================================================


  !===================================================================================================
  SUBROUTINE advance_dyn_Cor_ED(u_n,v_n,ustr_sfc,vstr_sfc,ustr_bot,vstr_bot,Hz,Akv,fcor,dt,N,u_np1,v_np1)
  !---------------------------------------------------------------------------------------------------
    !!============================================================================<br />
    !!                  ***  ROUTINE advance_dyn_Cor_ED  ***                      <br />
    !! ** Purposes : integrate vertical viscosity and Coriolis terms for dynamics <br />
    !!============================================================================<br />
    !! 1- Compute Coriolis term <br />
    !! if n is even
    !! \begin{align*}
    !! u^{n+1,\star} &= u^n + \Delta t f v^n \\
    !! v^{n+1,\star} &= v^n - \Delta t f u^{n+1,\star}
    !! \end{align*}
    !! if n is odd
    !! \begin{align*}
    !! v^{n+1,\star} &= v^n - \Delta t f u^n \\
    !! u^{n+1,\star} &= u^n + \Delta t f v^{n+1,\star}
    !! \end{align*}
    IMPLICIT NONE
    INTEGER, INTENT(IN   )          :: N                  !! number of vertical levels
    REAL(8), INTENT(IN   )          :: dt                 !! time-step [s]
    REAL(8), INTENT(IN   )          :: fcor               !! Coriolis frequaency [s-1]
    REAL(8), INTENT(IN   )          :: u_n   (1:N )       !! u-velocity component at time n [m/s]
    REAL(8), INTENT(IN   )          :: v_n   (1:N )       !! v-velocity component at time n [m/s]
    REAL(8), INTENT(IN   )          :: Hz      (1:N)      !! layer thickness [m]
    REAL(8), INTENT(IN   )          :: Akv     (0:N)      !! eddy-viscosity [m2/s]
    REAL(8), INTENT(IN   )          :: ustr_sfc              !! zonal surface stress      [m2/s2]
    REAL(8), INTENT(IN   )          :: vstr_sfc              !! meridional surface stress [m2/s2]
    REAL(8), INTENT(IN   )          :: ustr_bot              !! zonal surface stress      [m2/s2]
    REAL(8), INTENT(IN   )          :: vstr_bot              !! meridional surface stress [m2/s2]
    REAL(8), INTENT(  OUT)          :: u_np1 (1:N )       !! u-velocity component at time n+1 [m/s]
    REAL(8), INTENT(  OUT)          :: v_np1 (1:N )       !! v-velocity component at time n+1 [m/s]
    ! local variables
    INTEGER                         :: k
    REAL(8)                         :: ff(1:N), cff, cff1
    REAL(8)                         :: gamma_Cor = 0.55
    !
    cff      = (dt*fcor)*(dt*fcor)
    cff1     = 1./(1+gamma_Cor*gamma_Cor*cff)
    DO k=1,N
      u_np1(k) = cff1*Hz(k)*( (1.-gamma_Cor*(1.-gamma_Cor)*cff)*u_n(k)  &
                                                      + dt*fcor*v_n(k) )
      v_np1(k) = cff1*Hz(k)*( (1.-gamma_Cor*(1.-gamma_Cor)*cff)*v_n(k)  &
                                                      - dt*fcor*u_n(k) )
    ENDDO
    !IF(nn==1) THEN             !! 1 - Coriolis term <br />
    !  DO k=1,N
    !    cff       = fcor * v_n(k)
    !    u_np1(k)  = u_n(k) + dt * cff
    !    cff       = fcor * u_np1(k)
    !    v_np1(k)  = Hz(k)*( v_n(k) - dt * cff  )
    !    u_np1(k)  = Hz(k)*( u_np1(k)           )
    !  ENDDO
    !ELSE
    !  DO k=1,N
    !    cff       = fcor * u_n(k)
    !    v_np1(k)  =         v_n(k)  - dt * cff
    !    cff       = fcor * v_np1(k)
    !    u_np1(k)  = Hz(k)*( u_n(k) + dt * cff  )
    !    v_np1(k)  = Hz(k)*( v_np1(k)           )
    !  ENDDO
    !ENDIF
    !! 2 - Apply surface and bottom forcing <br />
    u_np1(N)=u_np1(N) + dt*ustr_sfc    !<-- sustr is in m2/s2 here
    v_np1(N)=v_np1(N) + dt*vstr_sfc
    u_np1(1)=u_np1(1) - dt*ustr_bot    !<-- sustr is in m2/s2 here
    v_np1(1)=v_np1(1) - dt*vstr_bot
    !=======================================================================
    !! 3 - Implicit integration for vertical viscosity <br />
    !! \begin{align*}
    !! \mathbf{u}^{n+1,\star \star} &= \mathbf{u}^{n+1,\star} + \Delta t \partial_z \left(  K_m \partial_z  \mathbf{u}^{n+1,\star \star} \right)  \\
    !! \end{align*}
    !=======================================================================
    ! u-component
    !--------------------------------------------------------
    ff(1:N) = u_np1(1:N)
    CALL tridiag_solve(N,Hz,Akv,ff,dt) ! Invert tridiagonal matrix
    u_np1(1:N) = ff(1:N)    !EXP(-dt*sig(1:N))*ff(1:N)
    !
    ! v-component
    !--------------------------------------------------------
    ff(1:N) = v_np1(1:N)
    CALL tridiag_solve(N,Hz,Akv,ff,dt)    ! Invert tridiagonal matrix
    v_np1(1:N) = ff(1:N)
  !---------------------------------------------------------------------------------------------------
  END SUBROUTINE advance_dyn_Cor_ED
  !===================================================================================================


  !===================================================================================================
  SUBROUTINE compute_evd(bvf,Akv,Akt,AkEvd,N)
  !---------------------------------------------------------------------------------------------------
    !!==========================================================================<br />
    !!                  ***  ROUTINE compute_evd  ***                           <br />
    !! ** Purposes : compute enhanced vertical diffusion/viscosity
    !!                                   where the density profile is unstable  <br />
    !!==========================================================================<br />
    IMPLICIT NONE
    INTEGER, INTENT(IN   )         :: N                  !! number of vertical levels
    REAL(8), INTENT(IN   )         :: bvf(0:N)           !! Brunt Vaisala frequency [s-2]
    REAL(8), INTENT(INOUT)         :: Akv     (0:N)      !! eddy-viscosity [m2/s]
    REAL(8), INTENT(INOUT)         :: Akt     (0:N)      !! eddy-diffusivity [m2/s]
    REAL(8), INTENT(IN   )         :: AkEvd              !! value of enhanced diffusion [m2/s]
    ! local variables
    INTEGER                         :: k
    !
    DO k=1,N-1
      IF( bvf(k) <= -1.e-12 ) THEN
        Akv  (k  )  = AkEvd
        Akt  (k  )  = AkEvd
      ENDIF
    ENDDO
  !---------------------------------------------------------------------------------------------------
  END SUBROUTINE compute_evd
  !===================================================================================================


  !===================================================================================================
  SUBROUTINE lmd_swfrac(Hz,N,swr_frac)
  !---------------------------------------------------------------------------------------------------
    !!==========================================================================<br />
    !!                  ***  ROUTINE lmd_swfrac  ***                            <br />
    !! ** Purposes : Compute fraction of solar shortwave flux penetrating to specified
    !!                depth due to exponential decay in Jerlov water type.      <br />
    !!==========================================================================<br />
    IMPLICIT NONE
    INTEGER, INTENT(IN   )         :: N                  !! number of vertical levels
    REAL(8), INTENT(IN   )         :: Hz      (1:N)      !! layer thickness [m]
    REAL(8), INTENT(  OUT)         :: swr_frac(0:N)      !! fraction of solar penetration
    ! local variables
    INTEGER   ::  k
    REAL(8) swdk1,swdk2,xi1,xi2
    REAL(8) mu1,mu2, r1, attn1, attn2
    !----------
    mu1   = 0.35; mu2   = 23.0; r1    = 0.58
    attn1 = -1./mu1
    attn2 = -1./mu2
    swdk1 = r1               ! surface, then attenuate
    swdk2 = 1.-swdk1              ! them separately throughout
    swr_frac(N)=1.              ! the water column.
    !
    DO k=N,1,-1
      xi1=attn1*Hz(k)
      IF (xi1 .gt. -20.) then        ! this logic to avoid
        swdk1=swdk1*exp(xi1)        ! computing exponent for
      ELSE                           ! a very large argument
        swdk1=0.
      ENDIF
      xi2=attn2*Hz(k)
      IF (xi2 .gt. -20.) then
        swdk2=swdk2*exp(xi2)
      ELSE
        swdk2=0.
      ENDIF
      swr_frac(k-1)=swdk1+swdk2
    ENDDO
  !---------------------------------------------------------------------------------------------------
  END SUBROUTINE lmd_swfrac
  !===================================================================================================

  !===================================================================================================
  SUBROUTINE compute_mxl (bvf,rhoc10,rhoc300,zr,rhoRef,N,hmxl10,hmxl300)
  !---------------------------------------------------------------------------------------------------
    USE scm_par, ONLY : grav
    IMPLICIT NONE
    INTEGER, INTENT(IN   )       :: N                  !! number of vertical levels
    REAL(8), INTENT(IN   )       :: bvf(0:N)           !! Brunt Vaisala frequancy [s-2]
    REAL(8), INTENT(IN   )       :: zr (1:N)           !! depth at cell center [m]
    REAL(8), INTENT(IN   )       :: rhoc10             !! thermal expension coefficient [kg m-3]
    REAL(8), INTENT(IN   )       :: rhoc300            !! thermal expension coefficient [kg m-3]
    REAL(8), INTENT(IN   )       :: rhoRef
    REAL(8), INTENT(  OUT)       :: hmxl10             !! mixed layer depth [m]
    REAL(8), INTENT(  OUT)       :: hmxl300            !! mixed layer depth [m]
    ! local variables
    integer                        :: k
    real(8)                        :: cff_new,cff_old,bvf_c
    integer                        :: kstart
    ! find 10m depth
    kstart = N
    do while( zr(kstart) > -10.  )
      kstart = kstart - 1
    enddo
    !
    bvf_c = rhoc10*(grav/rhoRef)
    !
    hmxl10  = zr(kstart) ! initialize at the near bottom value
    cff_old = 0.
    !
    outer: DO k=kstart,2,-1
      cff_new = cff_old + MAX(bvf(k-1),0.)*(zr(k)-zr(k-1))
      !
      !IF(zr(k)<=-150.0 .AND. .not. found_150m) THEN
        !print*,'Int_10m_150m = ',cff_new
      !  found_150m = .true.
      !ENDIF
      !
      IF( cff_new >= bvf_c ) THEN
          !hmxl = zr(k-1); exit outer
          hmxl10 = ((cff_new-bvf_c)*zr(k)+(bvf_c-cff_old)*zr(k-1))/(cff_new-cff_old)
          exit outer
      ENDIF
      cff_old = cff_new
    ENDDO outer
    ! find 150m depth
    kstart = N
    do while( zr(kstart) > -295.  )
      kstart = kstart - 1.
    enddo
    !
    bvf_c = rhoc300*(grav/rhoRef)
    !
    hmxl300  = zr(kstart) ! initialize at the near bottom value
    cff_old = 0.
    outer2: DO k=kstart,2,-1
      cff_new = cff_old + MAX(bvf(k-1),0.)*(zr(k)-zr(k-1))
      IF( cff_new >= bvf_c ) THEN
          !hmxl = zr(k-1); exit outer
          hmxl300 = ((cff_new-bvf_c)*zr(k)+(bvf_c-cff_old)*zr(k-1))/(cff_new-cff_old)
          exit outer2
      ENDIF
      cff_old = cff_new
    ENDDO outer2
    !
  !---------------------------------------------------------------------------------------------------
  END SUBROUTINE compute_mxl
  !===================================================================================================



  !===================================================================================================
  SUBROUTINE compute_mxl2 (bvf,rhoc,zr,zref,rhoRef,N,hmxl)
  !---------------------------------------------------------------------------------------------------
    USE scm_par
    IMPLICIT NONE
    INTEGER, INTENT(IN   )       :: N                  !! number of vertical levels
    REAL(8), INTENT(IN   )       :: bvf(0:N)           !! Brunt Vaisala frequancy [s-2]
    REAL(8), INTENT(IN   )       :: zr (1:N)           !! depth at cell center [m]
    REAL(8), INTENT(IN   )       :: rhoc               !! thermal expension coefficient [kg m-3]
    REAL(8), INTENT(IN   )       :: zref
    REAL(8), INTENT(IN   )       :: rhoRef
    REAL(8), INTENT(  OUT)       :: hmxl              !! mixed layer depth [m]
    ! local variables
    integer                        :: k
    real(8)                        :: cff_new,cff_old,bvf_c
    integer                        :: kstart
    ! find 300m depth
    kstart = N
    do while( zr(kstart) > zref  )
      kstart = kstart - 1
    enddo
    !
    bvf_c = rhoc*(grav/rhoRef)
    !
    hmxl  = zr(kstart) ! initialize at the near bottom value
    cff_old = 0.
    !
    outer: DO k=kstart,2,-1
      cff_new = cff_old + MAX(bvf(k-1),0.)*(zr(k)-zr(k-1))
      !
      IF( cff_new >= bvf_c ) THEN
          hmxl = ((cff_new-bvf_c)*zr(k)+(bvf_c-cff_old)*zr(k-1))/(cff_new-cff_old)
          exit outer
      ENDIF
      cff_old = cff_new
    ENDDO outer
  !---------------------------------------------------------------------------------------------------
  END SUBROUTINE compute_mxl2
  !===================================================================================================





  !===================================================================================================
  SUBROUTINE rho_eos_lin (temp,salt,zr,eos_params,N,neos,rho,bvf)
  !---------------------------------------------------------------------------------------------------
    !!==========================================================================<br />
    !!                  ***  ROUTINE rho_eos_lin  ***                           <br />
    !! ** Purposes : Compute density anomaly and Brunt Vaisala frequency via linear
    !!                                                  Equation Of State (EOS) <br />
    !!==========================================================================<br />
    USE scm_par, ONLY : grav
    IMPLICIT NONE
    INTEGER, INTENT(IN   )         :: N,neos             !! number of vertical levels
    REAL(8), INTENT(IN   )         :: temp(1:N)          !! temperature [C]
    REAL(8), INTENT(IN   )         :: salt(1:N)          !! salinity [psu]
    REAL(8), INTENT(IN   )         :: zr (1:N  )         !! depth at cell centers [m]
    REAL(8), INTENT(IN   )         :: eos_params(neos)
    REAL(8), INTENT(  OUT)         :: bvf(0:N  )         !! Brunt Vaisala frequancy [s-2]
    REAL(8), INTENT(  OUT)         :: rho(1:N  )         !! density anomaly [kg/m3]
    ! local variables
    integer                        :: k
    real(8)                        :: cff
    real(8)                        :: alpha,beta,rhoRef,T0,S0
    !---------------------------------------------------------------------------
    rhoRef = eos_params(1); alpha  = eos_params(2); beta = eos_params(3)
    T0     = eos_params(4); S0     = eos_params(5)
    !---------------------------------------------------------------------------
    DO k=1,N
      rho(k)= rhoRef*( 1. - alpha*( temp(k) - T0) + beta*( salt(k) - S0) )    !! \(   \rho_{k} = \rho_0 \left( 1 - \alpha (\theta - 2) + \beta (S - 35)   \right)  \)  <br />
    ENDDO
    !----
    DO k=1,N-1
      cff    = 1./(zr(k+1)-zr(k))
      bvf(k) = -cff*(grav/rhoRef)*(rho(k+1)-rho(k))  !!  \(   (N^2)_{k+1/2} = - \frac{g}{\rho_0}  \frac{ \rho_{k+1}-\rho_{k} }{\Delta z_{k+1/2}} \)
    ENDDO
    bvf(0) = 0.
    bvf(N) = bvf(N-1)
  !---------------------------------------------------------------------------------------------------
  END SUBROUTINE rho_eos_lin
  !===================================================================================================


  !===================================================================================================
  SUBROUTINE rho_eos (temp,salt,zr,zw,rhoRef,N,rho,bvf)
  !---------------------------------------------------------------------------------------------------
    !!==========================================================================<br />
    !!                  ***  ROUTINE rho_eos  ***                           <br />
    !! ** Purposes : Compute density anomaly and Brunt Vaisala frequency via nonlinear
    !!                                                  Equation Of State (EOS) <br />
    !!==========================================================================<br />
    USE scm_par, ONLY : grav
    IMPLICIT NONE
    INTEGER, INTENT(IN   )         :: N                  !! number of vertical levels
    REAL(8), INTENT(IN   )         :: temp(1:N)          !! temperature [C]
    REAL(8), INTENT(IN   )         :: salt(1:N)          !! salinity [psu]
    REAL(8), INTENT(IN   )         :: zr (1:N  )         !! depth at cell centers [m]
    REAL(8), INTENT(IN   )         :: zw (0:N  )         !! depth at cell interfaces [m]
    REAL(8), INTENT(IN   )         :: rhoRef
    REAL(8), INTENT(  OUT)         :: bvf(0:N  )         !! Brunt Vaisala frequancy [s-2]
    REAL(8), INTENT(  OUT)         :: rho(1:N  )         !! density anomaly [kg/m3]
    ! local variables
    integer                        :: k
    real(8)                        :: cff,cff1,dr00,dpth
    real(8)                        :: Tt,Ts,sqrtTs,K0,K1,K2
    real(8)                        :: rho1(1:N),K_up(1:N),K_dw(1:N)
    ! EOS parameters
    real(8)                        :: r00,r01,r02,r03,r04,r05
    real(8)                        :: r10,r11,r12,r13,r14,r20
    real(8)                        :: rS0,rS1,rS2
    real(8)                        :: K00,K01,K02,K03,K04
    real(8)                        :: K10,K11,K12,K13
    real(8)                        :: KS0,KS1,KS2
    real(8)                        :: B00,B01,B02,B03,B10,B11,B12,BS1
    real(8)                        :: E00,E01,E02,E10,E11,E12
    !---------------------------------------------------------------------------
    parameter(r00=999.842594,   r01=6.793952E-2,  r02=-9.095290E-3,  &
              r03=1.001685E-4,  r04=-1.120083E-6, r05=6.536332E-9)
    !
    parameter(r10=0.824493,     r11=-4.08990E-3,  r12=7.64380E-5,    &
              r13=-8.24670E-7,  r14=5.38750E-9,   r20=4.8314E-4)

    parameter(rS0=-5.72466E-3,  rS1=1.02270E-4,   rS2=-1.65460E-6)
    !
    parameter(K00=19092.56,     K01=209.8925,     K02=-3.041638,    &
                                K03=-1.852732e-3, K04=-1.361629e-5)

    parameter(K10=104.4077,     K11=-6.500517,    K12=0.1553190,    &
                                                  K13=2.326469e-4)
    !
    parameter(KS0=-5.587545,    KS1=+0.7390729,   KS2=-1.909078e-2)
    !
    parameter(B00=0.4721788,    B01=0.01028859,   B02=-2.512549e-4, &
              B03=-5.939910e-7, B10=-0.01571896,  B11=-2.598241e-4, &
                                B12=7.267926e-6,  BS1=2.042967e-3)
    !
    parameter(E00=+1.045941e-5, E01=-5.782165e-10,E02=+1.296821e-7,  &
              E10=-2.595994e-7, E11=-1.248266e-9, E12=-3.508914e-9)
    !---------------------------------------------------------------------------
    dr00=r00-rhoRef
    ! Compute density anomaly via Equation Of State (EOS) for seawater
    DO k=1,N
      Tt       = temp(k)
      Ts       = salt(k)
      sqrtTs   = sqrt(Ts)
      !
      rho1(k)  = (dr00 +Tt*( r01+Tt*( r02+Tt*( r03+Tt*(r04+Tt*r05 ))))  &
                       +Ts*( r10+Tt*( r11+Tt*( r12+Tt*(r13+Tt*r14 )))   &
                          +sqrtTs*(rS0+Tt*(rS1+Tt*rS2 ))+Ts*r20 )       &
                 )
      !
      K0= Tt*( K01+Tt*( K02+Tt*( K03+Tt*K04 )))   &
         +Ts*( K10+Tt*( K11+Tt*( K12+Tt*K13 ))    &
               +sqrtTs*( KS0+Tt*( KS1+Tt*KS2 )) )
      !
      K1=B00+Tt*(B01+Tt*(B02+Tt*B03)) +Ts*( B10+Tt*( B11+Tt*B12 )+sqrtTs*BS1 )
      !
      K2=E00+Tt*(E01+Tt*E02) +Ts*(E10+Tt*(E11+Tt*E12))
      !
      dpth = -zr(k)
      cff  = K00-0.1*dpth
      cff1 = K0+dpth*(K1+K2*dpth)
      rho(k) = ( rho1(k)*cff*(K00+cff1)-0.1*dpth*rhoRef*cff1 )/(cff*(cff+cff1))
      ! For bvf computation
      dpth=0.-zw(k)
      K_up(k)=K0+dpth*(K1+K2*dpth)
      dpth=0.-zw(k-1)
      K_dw(k)=K0+dpth*(K1+K2*dpth)
    ENDDO
    cff=grav/rhoRef
    DO k=1,N-1
      cff1=0.1*(0.-zw(k))
      bvf(k)=-cff *(      (rho1(k+1)-rho1(k))*(K00+K_dw(k+1))*(K00+K_up(k))            &
             -cff1*( rhoRef*(K_dw(k+1)-K_up(k)) + K00*(rho1(k+1)-rho1(k))     &
                    + rho1(k+1)*K_dw(k+1) - rho1(k)*K_up(k)   &
        ) )/(  (K00+K_dw(k+1)-cff1)*(K00+K_up(k)-cff1)      &
                                        *(zr(k+1)-zr(k))  )
    ENDDO
    bvf(0) = 0.
    bvf(N) = 0.
  !---------------------------------------------------------------------------------------------------
  END SUBROUTINE rho_eos
  !===================================================================================================


  !===================================================================================================
  SUBROUTINE rho_eos2(temp,salt,zr,rhoRef,N,rho)
  !---------------------------------------------------------------------------------------------------
    !!==========================================================================<br />
    !!                  ***  ROUTINE rho_eos  ***                           <br />
    !! ** Purposes : Compute density anomaly and Brunt Vaisala frequency via nonlinear
    !!                                                  Equation Of State (EOS) <br />
    !!==========================================================================<br />
    USE scm_par, ONLY : grav
    IMPLICIT NONE
    INTEGER, INTENT(IN   )         :: N                  !! number of vertical levels
    REAL(8), INTENT(IN   )         :: temp(1:N)          !! temperature [C]
    REAL(8), INTENT(IN   )         :: salt(1:N)          !! salinity [psu]
    REAL(8), INTENT(IN   )         :: zr (1:N  )         !! depth at cell centers [m]
    REAL(8), INTENT(IN   )         :: rhoRef
    REAL(8), INTENT(  OUT)         :: rho(1:N  )         !! density anomaly [kg/m3]
    ! local variables
    integer                        :: k
    real(8)                        :: cff,cff1,dr00,dpth
    real(8)                        :: Tt,Ts,sqrtTs,K0,K1,K2
    real(8)                        :: rho1(1:N),K_up(1:N),K_dw(1:N)
    ! EOS parameters
    real(8)                        :: r00,r01,r02,r03,r04,r05
    real(8)                        :: r10,r11,r12,r13,r14,r20
    real(8)                        :: rS0,rS1,rS2
    real(8)                        :: K00,K01,K02,K03,K04
    real(8)                        :: K10,K11,K12,K13
    real(8)                        :: KS0,KS1,KS2
    real(8)                        :: B00,B01,B02,B03,B10,B11,B12,BS1
    real(8)                        :: E00,E01,E02,E10,E11,E12
    !---------------------------------------------------------------------------
    parameter(r00=999.842594,   r01=6.793952E-2,  r02=-9.095290E-3,  &
              r03=1.001685E-4,  r04=-1.120083E-6, r05=6.536332E-9)
    !
    parameter(r10=0.824493,     r11=-4.08990E-3,  r12=7.64380E-5,    &
              r13=-8.24670E-7,  r14=5.38750E-9,   r20=4.8314E-4)

    parameter(rS0=-5.72466E-3,  rS1=1.02270E-4,   rS2=-1.65460E-6)
    !
    parameter(K00=19092.56,     K01=209.8925,     K02=-3.041638,    &
                                K03=-1.852732e-3, K04=-1.361629e-5)

    parameter(K10=104.4077,     K11=-6.500517,    K12=0.1553190,    &
                                                  K13=2.326469e-4)
    !
    parameter(KS0=-5.587545,    KS1=+0.7390729,   KS2=-1.909078e-2)
    !
    parameter(B00=0.4721788,    B01=0.01028859,   B02=-2.512549e-4, &
              B03=-5.939910e-7, B10=-0.01571896,  B11=-2.598241e-4, &
                                B12=7.267926e-6,  BS1=2.042967e-3)
    !
    parameter(E00=+1.045941e-5, E01=-5.782165e-10,E02=+1.296821e-7,  &
              E10=-2.595994e-7, E11=-1.248266e-9, E12=-3.508914e-9)
    !---------------------------------------------------------------------------
    dr00=r00-rhoRef
    ! Compute density anomaly via Equation Of State (EOS) for seawater
    DO k=1,N
      Tt       = temp(k)
      Ts       = salt(k)
      sqrtTs   = sqrt(Ts)
      !
      rho1(k)  = (dr00 +Tt*( r01+Tt*( r02+Tt*( r03+Tt*(r04+Tt*r05 ))))  &
                       +Ts*( r10+Tt*( r11+Tt*( r12+Tt*(r13+Tt*r14 )))   &
                          +sqrtTs*(rS0+Tt*(rS1+Tt*rS2 ))+Ts*r20 )       &
                 )
      !
      K0= Tt*( K01+Tt*( K02+Tt*( K03+Tt*K04 )))   &
         +Ts*( K10+Tt*( K11+Tt*( K12+Tt*K13 ))    &
               +sqrtTs*( KS0+Tt*( KS1+Tt*KS2 )) )
      !
      K1=B00+Tt*(B01+Tt*(B02+Tt*B03)) +Ts*( B10+Tt*( B11+Tt*B12 )+sqrtTs*BS1 )
      !
      K2=E00+Tt*(E01+Tt*E02) +Ts*(E10+Tt*(E11+Tt*E12))
      !
      dpth = -zr(k)
      cff  = K00-0.1*dpth
      cff1 = K0+dpth*(K1+K2*dpth)
      rho(k) = ( rho1(k)*cff*(K00+cff1)-0.1*dpth*rhoRef*cff1 )/(cff*(cff+cff1))
    ENDDO
  !---------------------------------------------------------------------------------------------------
  END SUBROUTINE rho_eos2
  !===================================================================================================



  !===================================================================================================
  SUBROUTINE tridiag_solve(N,Hz,Ak,f,dt)
  !---------------------------------------------------------------------------------------------------
    !!============================================================================<br />
    !!                  ***  ROUTINE tridiag_solve  ***                           <br />
    !! ** Purposes : solve the tridiagonal problem associated with the implicit
    !!                         in time treatment of vertical diffusion/viscosity  <br />
    !!============================================================================<br />
    IMPLICIT NONE
    INTEGER, INTENT(IN   )    :: N                  !! number of vertical levels
    REAL(8), INTENT(IN   )    :: Hz(1:N)            !! layer thickness [m]
    REAL(8), INTENT(IN   )    :: dt                 !! time-step [s]
    REAL(8), INTENT(IN   )    :: Ak(0:N)            !! eddy diffusivity/viscosity [m2/s]
    REAL(8), INTENT(INOUT)    :: f(1:N)             !! (in: right-hand side) (out:solution of tridiagonal problem)
    ! local variables
    INTEGER                   :: k
    REAL(8)                   :: a(1:N),b(1:N),c(1:N),q(1:N)
    REAL(8)                   :: difA,difC,cff
    !
    DO k=2,N-1
      difA  =       -2.*dt* Ak  (k-1) / (Hz(k-1)+Hz(k))
      difC  =       -2.*dt* Ak  (k  ) / (Hz(k+1)+Hz(k))
      a (k) = difA
      c (k) = difC
      b (k) = Hz(k) - difA - difC
    ENDDO
    !++ Bottom BC
    a (1) = 0.
    difC  = - 2.*dt*Ak(1)/( Hz(2)+Hz(1) )
    c (1) = difC
    b (1) = Hz(1) - difC
    !++ Surface BC
    difA  = -2.*dt*Ak(N-1)/( Hz(N-1)+Hz(N))
    a (N) = difA
    c (N) = 0.
    b (N) = Hz(N) - difA
    !
    cff   = 1./b(1)
    q (1) = - c(1)*cff
    f (1) =   f(1)*cff
    DO k=2,N
      cff=1./(b(k)+a(k)*q(k-1))
      q(k)= -cff*c(k)
      f(k)=cff*( f(k)-a(k)*f(k-1) )
    ENDDO
    DO k=N-1,1,-1
      f(k)=f(k)+q(k)*f(k+1)
    ENDDO
  !---------------------------------------------------------------------------------------------------
  END SUBROUTINE tridiag_solve
  !===================================================================================================
  !
END MODULE scm_oce

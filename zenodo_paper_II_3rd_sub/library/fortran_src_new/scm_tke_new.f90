MODULE scm_tke_new
  !!============================================================================<br />
  !!                       ***  MODULE  scm_tke  ***                            <br />
  !! Eddy-diffusion closure: TKE scheme with diagnostic mixing lengths
  !!                         adapted from NEMO TKE turbulent closure model      <br />
  !!============================================================================<br />
  !!----------------------------------------------------------------------------<br />
  !!   compute_tke_bdy   : top and bottom boundary conditions for TKE           <br />
  !!   compute_shear     : compute shear production term                        <br />
  !!   advance_tke       : tke time stepping: advance tke at next time step     <br />
  !!   compute_mxl       : compute mixing length scale                          <br />
  !!   compute_ED        : compute avm and avt                                  <br />
  !!   tridiag_solve_tke : tridiagonal solver for TKE equation                  <br />
  !!----------------------------------------------------------------------------<br />
    IMPLICIT NONE

CONTAINS

  SUBROUTINE compute_tke_bdy(ustr_sfc, vstr_sfc, tke_const, bc_ap, wp0, tkemin, tke_sfc, tke_bot, flux_sfc)
    !!==========================================================================<br />
    !!                  ***  ROUTINE compute_tke_bdy  ***                       <br />
    !! ** Purposes : compute top and bottom boundary conditions for TKE equation<br />
    !!==========================================================================<br />
    USE scm_par, ONLY: cm_nemo,ceps_nemo,cm_mnh,ceps_mnh,cm_r81,ceps_r81, tke_min0
    IMPLICIT NONE
    REAL(8), INTENT(IN   )         :: ustr_sfc         !! zonal surface stress      [m2/s2]
    REAL(8), INTENT(IN   )         :: vstr_sfc         !! meridional surface stress [m2/s2]
    INTEGER, INTENT(IN   )         :: tke_const        !! choice of TKE constants
    REAL(8), INTENT(IN   )         :: bc_ap            !! choice of TKE constants
    REAL(8), INTENT(IN   )         :: wp0              !! surface value for plume vertical velocity [m/s]
    REAL(8), INTENT(IN   )         :: tkemin
    REAL(8), INTENT(  OUT)         :: tke_sfc          !! surface value for Dirichlet condition [m2/s2]
    REAL(8), INTENT(  OUT)         :: tke_bot          !! bottom value for Dirichlet condition [m2/s2]
    REAL(8), INTENT(  OUT)         :: flux_sfc         !! surface TKE ED flux for Neumann condition [m3/s3]
    !local variables
    REAL(8)                        :: ustar2_sfc, cm0inv2
    !
    ustar2_sfc = SQRT( ustr_sfc**2+vstr_sfc**2 )
    !
    IF( ustar2_sfc == 0. ) THEN     !! free convection case ( \( {\rm tke\_sfc\_dirichlet = True}  \) ) :
      tke_sfc  = tke_min0           !! \( k_{\rm sfc} = 10^{-4}\;{\rm m}^2\;{\rm s}^{-2} \)<br />
    ELSE
      IF(tke_const==0) THEN
        cm0inv2 = 1./SQRT(cm_nemo*ceps_nemo)
      ELSE IF(tke_const==1) THEN
        cm0inv2 = 1./SQRT( cm_mnh*ceps_mnh )
      ELSE
        cm0inv2 = 1./SQRT( cm_r81*ceps_r81 )
      ENDIF
      tke_sfc  = MAX( tkemin, cm0inv2*ustar2_sfc )
    ENDIF
    flux_sfc = 0.5*bc_ap*(wp0)**3  !! energetically consistent boundary condition \( F_{\rm sfc}^k = \left.  K_e \partial_z k \right|_{\rm sfc} \)
    tke_bot  = tkemin  !! bottom boundary condition : \( k_{\rm bot} = k_{\rm min}  \)
  !-----------------------------------------------------------------------------
  END SUBROUTINE compute_tke_bdy
  !=============================================================================

  SUBROUTINE compute_shear(u_n, v_n, u_np1, v_np1, u_np1_star, v_np1_star, zr, N, shear2)
    !!==========================================================================<br />
    !!                  ***  ROUTINE compute_shear  ***                         <br />
    !! ** Purposes : compute shear production term for TKE equation             <br />
    !!==========================================================================<br />
    IMPLICIT NONE
    INTEGER, INTENT(IN   )              :: N                               !! number of vertical levels
    REAL(8), INTENT(IN   )              :: u_n  (1:N),v_n  (1:N)           !! velocity components at time n    [m/s]
    REAL(8), INTENT(IN   )              :: u_np1(1:N),v_np1(1:N)           !! velocity components at time n+1  [m/s]
    REAL(8), INTENT(IN   )              :: u_np1_star(1:N),v_np1_star(1:N) !! velocity components at time n+1 (*)  [m/s]
    REAL(8), INTENT(IN   )              :: zr(1:N)                         !! depth at cell centers [m]
    REAL(8), INTENT(  OUT)              :: shear2(0:N)                     !! shear production term [m2/s3]
    ! local variables
    INTEGER                             :: k
    REAL(8)                             :: du,dv,cff
    shear2(0:N) = 0.
    DO k=1,N-1
      cff       = 1. / ( zr(k+1)-zr(k) )**2
      du        = cff*( u_np1_star(k+1)-u_np1_star(k) )*0.5*( u_n(k+1)+u_np1(k+1)-u_n(k)-u_np1(k) ) !! Shear production term using discretization from Perrot & Lemarie (2024) <br />
      !! \( {\rm Sh}_{k+1/2} = \frac{ 1 }{ \Delta z_{k+1/2}^2 } ( u_{k+1}^{n+1,\star} - u_{k}^{n+1,\star} ) ( u_{k+1}^{n+1/2} - u_{k}^{n+1/2} )  \)
      dv        = cff*( v_np1_star(k+1)-v_np1_star(k) )*0.5*( v_n(k+1)+v_np1(k+1)-v_n(k)-v_np1(k) )
      shear2(k) = du + dv
    ENDDO
  !---------------------------------------------------------------------------------------------------
  END SUBROUTINE compute_shear
  !===================================================================================================
  !
  !===================================================================================================
  SUBROUTINE tkeeq (tke_n,NN,SS,num,nuh,lup,ldwn,h,zr,NN_MF,SS_MF,dWKdz_MF,wtke,               &
                             dt,tke_min,bdy_sfc,bdy_bot,tke_const,N,tke_np1, pdlr, eps_np1, residual)
  !---------------------------------------------------------------------------------------------------
    !!==========================================================================<br />
    !!                  ***  ROUTINE advance_tke  ***                           <br />
    !! ** Purposes : tke time stepping, advance tke from time step n to n+1     <br />
    !!==========================================================================<br />
    USE scm_par, ONLY:ceps_nemo,Ric_nemo,ce_nemo,cm_nemo, &
                    ceps_mnh,Ric_mnh,ce_mnh,cm_mnh,ct_mnh, &
                    ceps_r81,Ric_r81,ce_r81,cm_r81,ct_r81,bshear,pdlrmin
    IMPLICIT NONE
    INTEGER, INTENT(IN   )  :: N,tke_const
    REAL(8), INTENT(IN   )  :: tke_n (0:N)
    REAL(8), INTENT(IN   )  :: NN    (0:N)
    REAL(8), INTENT(IN   )  :: SS    (0:N)
    REAL(8), INTENT(IN   )  :: num   (0:N)
    REAL(8), INTENT(IN   )  :: nuh   (0:N)
    REAL(8), INTENT(IN   )  :: lup   (0:N)
    REAL(8), INTENT(IN   )  :: ldwn  (0:N)
    REAL(8), INTENT(IN   )  :: h     (1:N)
    REAL(8), INTENT(IN   )  :: zr    (1:N)
    REAL(8), INTENT(IN   )  :: NN_MF (0:N)      !! TKE buoyancy forcing term associated with mass flux [m2/s3]
    REAL(8), INTENT(IN   )  :: SS_MF (0:N)      !! TKE shear forcing term associated with mass flux [m2/s3]
    REAL(8), INTENT(IN   )  :: dWKdz_MF(0:N) !! ap wp \| u - up \|^2 contribution to d(w'k)/dz term [m2/s3]
    REAL(8), INTENT(IN   )  :: dt
    REAL(8), INTENT(IN   )  :: tke_min
    REAL(8), INTENT(IN   )  :: bdy_sfc(3)       ! bdy_sfc(1) = 0. -> Dirichlet, bdy_sfc(1) = 1. -> Neumann
    REAL(8), INTENT(IN   )  :: bdy_bot(3)       ! bdy_bot(1) = 0. -> Dirichlet, bdy_top(1) = 1. -> Neumann
    REAL(8), INTENT(INOUT)  :: wtke(1:N)
    REAL(8), INTENT(  OUT)  :: tke_np1(0:N)
    REAL(8), INTENT(  OUT)  :: pdlr   (0:N)
    REAL(8), INTENT(  OUT)  :: eps_np1(0:N)
    REAL(8), INTENT(  OUT)  :: residual                     !! Diagnostics : TKE spuriously added to guarantee that tke >= tke_min [m3/s3]
    !local variables
    INTEGER                 :: k,k_ubc,k_lbc
    REAL(8)                 :: mxld(0:N),mxlm(0:N)
    REAL(8)                 :: rhs,rhsmin
    REAL(8)                 :: avh(0:N), prod, buoyan, diss, Ri
    REAL(8)                 :: DiffKup, DiffKdw, Lsour(0:N),Qsour(0:N)
    REAL(8)                 :: cm,ct,ceps,isch,Ric
    !=====================================================================
    Lsour(0:N) = 0.; Qsour(0:N) = 0.; avh(0:N) = 0.
    tke_np1(0:N) = tke_n(0:N)
    !=====================================================================
    ! Initialize constants
    IF(tke_const==0) THEN
      ceps = ceps_nemo; Ric  = Ric_nemo; isch = ce_nemo / cm_nemo
    ELSE IF(tke_const==1) THEN
      cm   = cm_mnh; ct   = ct_mnh; ceps = ceps_mnh
      Ric  = Ric_mnh; isch = ce_mnh / cm_mnh
    ELSE
      cm   = cm_r81; ct   = ct_r81; ceps = ceps_r81
      Ric  = Ric_r81; isch = ce_r81 / cm_r81
    ENDIF
    !=====================================================================
    DO k=1,N-1
      mxld(k) = SQRT( lup(k) * ldwn(k) ) !! Dissipative mixing length
      mxlm(k) = MIN ( lup(k),  ldwn(k) )
    ENDDO
    !=====================================================================
    pdlr(0:N) = 0. !! Inverse Prandtl number function of Richardson number
    DO k = 1,N-1
      prod    = SS(k)
      buoyan  = NN(k)
      Ri      = MAX( buoyan, 0. ) / ( prod + bshear )  !!
      pdlr(k) = MAX(  pdlrmin,  Ric / MAX( Ric , Ri ) )      !!
    END DO
    !=====================================================================
    DO k = 1,N-1
      avh(k)   =  isch*num(k)
      prod     =  num(k)*SS(k) + SS_MF(k)
      buoyan   = -nuh(k)*NN(k) + NN_MF(k)
      rhs      = prod + buoyan + dWKdz_MF(k)
      diss     = ceps*SQRT(tke_n(k))/mxld(k)
      rhsmin   = (tke_min-tke_n(k))/dt + diss*tke_min
      Qsour(k) = MAX(rhs,rhsmin)
      Lsour(k) = -diss
      IF(rhs<rhsmin) THEN
        residual = residual + (zr(k+1)-zr(k))*(rhsmin-rhs)
      ENDIF
    ENDDO
    !=====================================================================
    !  do diffusion step
    !=====================================================================
    k_ubc   = int(bdy_sfc(1))
    k_lbc   = int(bdy_bot(1))
    DiffKup = bdy_sfc(2)
    DiffKdw = bdy_bot(2)
    !
    call diff_face(N,dt,h,k_ubc,k_lbc,                          &
                  DiffKup,DiffKdw,avh,Lsour,Qsour,tke_np1)
    !  fill top and bottom value
    tke_np1(N   )  = tke_np1(N-1)
    tke_np1(0   )  = bdy_bot(3)
    !
    DO k = 1,N
      wtke(k) = wtke(k) - 0.5*(avh(k)+avh(k-1))*(tke_np1(k)-tke_np1(k-1))/h(k)
    ENDDO
    !=====================================================================
    !  Finalize computation with dissipation diagnostics
    !=====================================================================
    DO k = 1,N-1
      eps_np1(k) = ceps*tke_np1(k)*SQRT(tke_n(k))/mxld(k)
      tke_np1(k) = MAX( tke_np1(k   ), tke_min)
    ENDDO
    !  fill top and bottom value
    eps_np1(N   )  = 0.        ; eps_np1(0   )  = 0.
    !
  !---------------------------------------------------------------------------------------------------
  END SUBROUTINE tkeeq
  !===================================================================================================




  !===================================================================================================
  SUBROUTINE tkeeq_MF(tke_np1_star,a_p,w_p,Fmass,NormUmUp,h,zr,wtke,mf_params,zinv,  &
                        tkep_min,mxlp_min,small_ap,opt,dt,nparams,N,tke_np1,tke_p,eps_p)
  !---------------------------------------------------------------------------------------------------
    !!==========================================================================<br />
    !!                  ***  ROUTINE advance_tke  ***                           <br />
    !! ** Purposes : tke time stepping, advance tke from time step n to n+1     <br />
    !!==========================================================================<br />
    USE scm_par, ONLY:ceps_nemo,wpmin
    IMPLICIT NONE
    INTEGER, INTENT(IN   )  :: N,nparams
    REAL(8), INTENT(IN   )  :: tke_np1_star (0:N)
    REAL(8), INTENT(IN   )  :: a_p          (0:N)
    REAL(8), INTENT(IN   )  :: w_p          (0:N)
    REAL(8), INTENT(IN   )  :: Fmass        (0:N)
    REAL(8), INTENT(IN   )  :: NormUmUp     (1:N)
    REAL(8), INTENT(IN   )  :: h            (1:N)
    REAL(8), INTENT(IN   )  :: zr           (1:N)
    REAL(8), INTENT(INOUT)  :: wtke         (1:N)
    REAL(8), INTENT(IN   )  :: zinv
    REAL(8), INTENT(IN   )  :: tkep_min
    REAL(8), INTENT(IN   )  :: mxlp_min
    REAL(8), INTENT(IN   )  :: dt
    LOGICAL, INTENT(IN   )  :: small_ap
    INTEGER, INTENT(IN   )  :: opt          !! option for the computation of wtke (=0 PL24, =1 HB09, = 2 Wi11)
    REAL(8), INTENT(IN   )  :: mf_params(1:nparams)  !! parameters in the ODEs
    REAL(8), INTENT(  OUT)  :: tke_np1(0:N)
    REAL(8), INTENT(  OUT)  :: tke_p  (0:N)
    REAL(8), INTENT(  OUT)  :: eps_p  (0:N)
    !local variables
    INTEGER                 :: k
    REAL(8)                 :: Kpp,Kpm,k_env
    REAL(8)                 :: FC(1:N),imxld0(1:N)
    REAL(8)                 :: sigma,rhs,dissp,delta0,Cent,Cdet
    REAL(8)                 :: cff,lup,ldwn,ap,hEnt,hDet,dwp,D0,tke_min
    LOGICAL                 :: tke_comput_new = .true.
    !=====================================================================
    !  compute the tke of the plume using tke_np1_star
    !=====================================================================
    !=====================================
    tke_min = tkep_min
    Cent    = mf_params(1)
    Cdet    = mf_params(2)
    delta0  = mf_params(9)/ABS(zinv)
    !=====================================
    tke_p(N) = tke_np1_star(N)
    !
    IF( opt==0 .OR. opt==1 ) THEN
      !
      DO k = 1,N
        cff       = zr(k)
        lup       = MAX( -cff     , mxlp_min)
        ldwn      = MAX(  cff-zinv, mxlp_min)
        imxld0(k) = 1./SQRT(lup*ldwn)
      ENDDO
      !
      DO k = N,1,-1
        Kpp     = a_p(k)*w_p(k)*tke_p(k)
        ap      = 0.5*(a_p(k)+a_p(k-1))
        dissp   = ceps_nemo*imxld0(k)*tke_p(k)*SQRT(tke_p(k))
        dwp     = w_p(k)-w_p(k-1)
        ! Compute entrainment and detrainment rates multiplied by h
        D0      = MIN(0.5*h(k)*delta0*(w_p(k  )+w_p(k-1)),-2.*wpmin)
        IF( dwp >= 0. ) THEN
          hEnt    =  ap*Cent*dwp; hDet    = -ap*D0
        ELSE
          hEnt    =  0.         ; hDet    = -ap*( Cdet*dwp + D0 )
        ENDIF
        !
        sigma   = 0.; IF(.NOT. small_ap) sigma   = a_p(k)/( 1.+a_p(k))
        k_env   = tke_np1_star(k-1)+sigma*(tke_np1_star(k-1)-tke_p(k))
        rhs     = hEnt*k_env-hDet*tke_p(k)
        IF(opt==0) THEN  ! Perrot & Lemarie 2024 formulation
          rhs = rhs + hEnt*( 1.+sigma )*NormUmUp(k) - ap*dissp
        ENDIF
        !
        Kpm = Kpp - rhs
        !
        IF( a_p(k-1)>0. .AND. w_p(k-1) < -wpmin ) THEN
          tke_p(k-1)  = MAX( Kpm/(a_p(k-1)*w_p(k-1)), tkep_min )
        ELSE
          tke_p(k-1)  = tke_p(k  )
        ENDIF
        eps_p(k)      = dissp
        !
      ENDDO
    !
    ELSE  ! Witek et al. 2011 formulation (no mass flux on tke)
      DO k = 1,N-1
        tke_p(k) = tke_np1_star(k-1)
        eps_p(k) = 0.
      ENDDO
    ENDIF
    !=====================================================================
    !  do mass flux step
    !=====================================================================
    FC(1:N) = 0.
    DO k = 1,N-1
      FC(k)   = Fmass(k  )*( tke_p(k)-tke_np1_star(k-1) )
    ENDDO
    DO k = 1,N-1
      tke_np1(k) = MAX( tke_np1_star(k) + dt*(FC(k+1)-FC(k))/(zr(k+1)-zr(k)), tke_min )
      wtke   (k) = wtke(k) - FC(k)
    ENDDO
    !  fill top and bottom value
    tke_np1(N   )  = tke_np1(N-1); tke_np1(0   )  = tke_np1_star(0)
  !---------------------------------------------------------------------------------------------------
  END SUBROUTINE tkeeq_MF
  !===================================================================================================





  !===================================================================================================
  SUBROUTINE compute_mxl(tke, bvf, Hz, ustr_sfc, vstr_sfc, mxlmin, N, lup, ldwn)
  !---------------------------------------------------------------------------------------------------
    !!============================================================================<br />
    !!                  ***  ROUTINE compute_mxl  ***                             <br />
    !! ** Purposes : compute mixing length scales                                 <br />
    !!============================================================================<br />
    USE scm_par, ONLY:grav,vkarmn,rsmall,mxl_min0
    IMPLICIT NONE
    INTEGER, INTENT(IN   )  :: N              !! number of vertical levels
    REAL(8), INTENT(IN   )  :: tke(0:N)       !! turbulent kinetic energy [m2/s2]
    REAL(8), INTENT(IN   )  :: bvf(0:N)       !! Brunt Vaisala frequency [s-2]
    REAL(8), INTENT(IN   )  :: Hz(1:N)        !! layer thickness [m]
    REAL(8), INTENT(IN   )  :: ustr_sfc       !! surface stress [m2/s2]
    REAL(8), INTENT(IN   )  :: vstr_sfc       !! surface stress [m2/s2]
    REAL(8), INTENT(IN   )  :: mxlmin         !!
    REAL(8), INTENT(  OUT)  :: lup (0:N)      !! upward mixing length [m]
    REAL(8), INTENT(  OUT)  :: ldwn(0:N)      !! downward mixing length [m]
    ! local variables
    INTEGER                 :: k
    REAL(8)                 :: rn2,ld80(0:N),raug,ustar2_sfc
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !  Mixing length
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ustar2_sfc =  SQRT( ustr_sfc**2+vstr_sfc**2 )  !! \( u_{\star}^2 = \sqrt{\tau_x^2 + \tau_y^2} \)<br />
    !
    lup (0:N)  = mxlmin
    ldwn(0:N)  = mxlmin
    !
    DO k = 0, N              ! interior value : l=sqrt(2*e/n^2)
      rn2     = MAX( bvf(k) , rsmall )
      ld80(k) = MAX( mxlmin,  SQRT(2.*tke(k) / rn2 ) ) !! Buoyancy length scale : \(  (l_{\rm up})_{k+1/2}=(l_{\rm dwn})_{k+1/2}=(l_{D80})_{k+1/2} = \sqrt{\frac{2 k_{k+1/2}^{n+1}}{ \max( (N^2)_{k+1/2}, (N^2)_{\min} ) }}     \)<br />
    END DO
    !! Physical limits for the mixing lengths <br />
    ldwn(0  ) = 0.
    DO k = 1, N
      ldwn(k) = MIN( ldwn(k-1) + Hz(k  ) , ld80(k) )   !! Limit \( (l_{\rm dwn})_{k+1/2} \) such that \( \partial_z (l_{\rm dwn})_{k} \le 1 \) the bottom boundary condition is \( (l_{\rm dwn})_{1/2} = l_{\min} \) <br />
    END DO
    ! surface mixing length = F(stress)=vkarmn*2.e5*taum/(rho0*g)
    raug      = vkarmn * 2.e5 / grav
    lup(N)   = MAX( mxl_min0, raug * ustar2_sfc ) ! surface boundary condition
    !
    DO k = N-1,0,-1
      lup(k) = MIN( lup(k+1) + Hz(k+1) , ld80(k) )   !! Limit \( (l_{\rm up})_{k-1/2} \) such that \( - \partial_z (l_{\rm up})_{k} \le 1 \) the surface boundary condition is \( (l_{\rm up})_{N+1/2} = \frac{\kappa}{g} (2 \times 10^{-5}) u_{\star}^2 \) <br />
    END DO
    lup(N) = 0.   !<-- ensures that avm(jpk) = 0.
    !
  !---------------------------------------------------------------------------------------------------
  END SUBROUTINE compute_mxl
  !===================================================================================================

  !===================================================================================================
  SUBROUTINE compute_ED(tke,lup,ldwn,pdlr,extrap_sfc,tke_const,N,Akv,Akt)
  !---------------------------------------------------------------------------------------------------
    !!============================================================================<br />
    !!                  ***  ROUTINE compute_ED  ***                              <br />
    !! ** Purposes : compute the vertical eddy viscosity and diffusivity          <br />
    !!============================================================================<br />
    USE scm_par, ONLY: cm_nemo,cm_mnh,cm_r81,avmolU,avmolT
    IMPLICIT NONE
    INTEGER, INTENT(IN   )  :: N                   !! number of vertical levels
    REAL(8), INTENT(IN   )  :: tke(0:N)            !! turbulent kinetic energy [m2/s2]
    REAL(8), INTENT(IN   )  :: pdlr(0:N)           !! inverse turbulent Prandtl number
    REAL(8), INTENT(IN   )  :: lup(0:N)            !! upward mixing length [m]
    REAL(8), INTENT(IN   )  :: ldwn(0:N)           !! downward mixing length [m]
    LOGICAL, INTENT(IN   )  :: extrap_sfc          !! (T) extrapolate eddy coefficients to the surface
    INTEGER, INTENT(IN   )  :: tke_const           !! choice of TKE constants
    REAL(8), INTENT(  OUT)  :: Akv(0:N)            !! eddy-viscosity [m2/s]
    REAL(8), INTENT(  OUT)  :: Akt(0:N)            !! eddy-diffusivity [m2/s]
    ! local variables
    INTEGER                 :: k
    REAL(8)                 :: mxlm,av,cm
    !* Initialize constants
    IF(tke_const==0) THEN
      cm = cm_nemo
    ELSE IF(tke_const==1) THEN
      cm = cm_mnh
    ELSE
      cm = cm_r81
    ENDIF
    !
    DO k = 0, N
      mxlm          = MIN ( lup(k),  ldwn(k) )      !! Compute "master" mixing length \( (l_m)_{k+1/2} = \min( (l_{\rm up})_{k+1/2}, (l_{\rm dwn})_{k+1/2} ) \)<br />
      av            = cm * mxlm * SQRT(tke(k))      !! Compute eddy-viscosity \( (K_m)_{k+1/2} = C_m l_m \sqrt{k_{k+1/2}^{n+1}}  \)<br />
      Akv  (k  )    = av + avmolU
      Akt  (k  )    = pdlr(k) * av + avmolT          !! Compute eddy-diffusivity \( (K_s)_{k+1/2} = ({\rm Pr}_t)^{-1}_{k+1/2}   (K_m)_{k+1/2} \)<br />
    END DO
    Akv(N) = 0.; Akt(N) = 0.
    !Warning : extrapolation ignores the variations of Hz with depth
    IF(extrap_sfc) THEN
      Akv(N) = 1.5*Akv(N-1)-0.5*Akv(N-2) !! if \( {\rm extrap\_sfc = True} \qquad \rightarrow \qquad (K_m)_{N+1/2} = \frac{3}{2} (K_m)_{N-1/2} - \frac{1}{2} (K_m)_{N-3/2} \)
      Akt(N) = 1.5*Akt(N-1)-0.5*Akt(N-2)
    ENDIF
  !---------------------------------------------------------------------------------------------------
  END SUBROUTINE compute_ED
  !===================================================================================================
  !
  !===================================================================================================
  SUBROUTINE diff_face(N,dt,h,Bcup,Bcdw,Yup,Ydw,nuY,Lsour,Qsour,Y)
  !---------------------------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN   )  :: N
    REAL(8), INTENT(IN   )  :: dt
    REAL(8), INTENT(IN   )  :: h         (1:N)      ! layer thickness
    INTEGER, INTENT(IN   )  :: Bcup                 !  type of upper BC
    INTEGER, INTENT(IN   )  :: Bcdw                 !  type of lower BC
    REAL(8), INTENT(IN   )  :: Yup                  !  value of upper BC
    REAL(8), INTENT(IN   )  :: Ydw                  !  value of lower BC
    REAL(8), INTENT(IN   )  :: nuY(0:N)
    REAL(8), INTENT(IN   )  :: Lsour(0:N)
    REAL(8), INTENT(IN   )  :: Qsour(0:N)
    REAL(8), INTENT(INOUT)  :: Y(0:N)
    INTEGER,parameter       :: Dirichlet      = 0
    INTEGER,parameter       :: Neumann        = 1
    ! local variables
    INTEGER  :: k
    REAL(8)  :: a,c,l
    REAL(8)  :: au(1:N-1),bu(1:N-1),cu(1:N-1),du(1:N-1)
    !
    au(1:N-1) = 0.;bu(1:N-1) = 0.
    cu(1:N-1) = 0.;du(1:N-1) = 0.
    !  set up matrix
    DO k=2,N-2
      c  = dt*( nuY(k+1) + nuY(k  ) )  / ( h(k)+h(k+1) ) / h(k+1)
      a  = dt*( nuY(k  ) + nuY(k-1) )  / ( h(k)+h(k+1) ) / h(k  )
      l  = dt*Lsour(k)

      cu(k) =-c
      au(k) =-a
      bu(k) = 1. + (a + c) - l
      du(k) = Y(k) + dt*Qsour(k)
    ENDDO
    !   set up upper boundary condition
    SELECT CASE(Bcup)
      CASE(Neumann)
        a = dt*( nuY(N-1) + nuY(N-2) )  / ( h(N-1)+h(N) ) / h(N-1)
        l = dt*Lsour(N-1)
        au(N-1) =-a
        bu(N-1) = 1. + a - l
        du(N-1) = Y(N-1) + dt*Qsour(N-1)       &
                + 2.0*dt*Yup/( h(N-1)+h(N) )
      CASE(Dirichlet)
        au(N-1) = 0.
        bu(N-1) = 1.
        du(N-1) = Yup
    END SELECT
    !   set up lower boundary condition
    SELECT CASE(Bcdw)
      CASE(Neumann)
        c = dt*( nuY(2) + nuY(1) )  / ( h(1)+h(2) ) / h(2)
        l = dt*Lsour(1)
        cu(1) =-c
        bu(1) = 1. + c - l
        du(1) = Y(1) + dt*Qsour(1)          &
                + 2.0*dt*Ydw/( h(1)+h(2) )
      CASE(Dirichlet)
          bu(1) = 1.
          cu(1) = 0.
          du(1) = Ydw
    END SELECT
    !  solve linear system
    CALL tridiagonal(au,bu,cu,du,Y,N)
    return
  !---------------------------------------------------------------------------------------------------
  END SUBROUTINE diff_face
  !===================================================================================================


  !===================================================================================================
  SUBROUTINE tridiagonal(au,bu,cu,du,f,N)
  !---------------------------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN   )  :: N
    REAL(8), INTENT(IN   )  :: au(1:N-1)
    REAL(8), INTENT(IN   )  :: bu(1:N-1)
    REAL(8), INTENT(IN   )  :: cu(1:N-1)
    REAL(8), INTENT(IN   )  :: du(1:N-1)
    REAL(8), INTENT(INOUT)  :: f(0:N)
    ! local variables
    INTEGER                 :: k
    REAL(8)                 :: ru(1:N-1),cff
    REAL(8)                 :: qu(1:N-1)
    !-------------------------------------
    ru(1:N-1) = 0.; qu(1:N-1) = 0.
    !
    cff = 1./bu(N-1)
    ru(N-1)=cff*au(N-1)
    qu(N-1)=cff*du(N-1)
    DO k=N-2,2,-1
      cff  = 1./(bu(k)-cu(k)*ru(k+1))
      ru(k)=cff*au(k)
      qu(k)=cff*(du(k)-cu(k)*qu(k+1))
    ENDDO
    qu(1)=(du(1)-cu(1)*qu(2))/(bu(1)-cu(1)*ru(2))
    f(1)=qu(1)
    DO k=2,N-1
       f(k)=qu(k)-ru(k)*f(k-1)
    ENDDO
    return
  !---------------------------------------------------------------------------------------------------
  END SUBROUTINE tridiagonal
  !===================================================================================================
!
END MODULE scm_tke_new

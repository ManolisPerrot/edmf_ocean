      MODULE scm_gls
      !
        IMPLICIT NONE
      !
      CONTAINS
      !++
      SUBROUTINE compute_tke_eps_bdy(ustr_sfc, vstr_sfc, ustr_bot, vstr_bot, z0b,       &
                                     tke_n_sfc, Hz_sfc, tke_n_bot, Hz_bot, OneOverSig,  &
                                     pnm, tke_min , eps_min,                            &
                                     tke_sfc ,  tke_bot, ftke_sfc, ftke_bot,            &
                                     eps_sfc ,  eps_bot, feps_sfc, feps_bot )
          !!==========================================================================<br />
          !!                  ***  ROUTINE compute_tke_gls_bdy  ***                       <br />
          !! ** Purposes : compute top and bottom boundary conditions for TKE and GLS equation<br />
          !!==========================================================================<br />
          USE scm_par
          IMPLICIT NONE
          REAL(8), INTENT(IN   )         :: ustr_sfc     !! zonal surface stress      [m2/s2]
          REAL(8), INTENT(IN   )         :: vstr_sfc     !! meridional surface stress [m2/s2]
          REAL(8), INTENT(IN   )         :: ustr_bot     !! zonal bottom stress      [m2/s2]
          REAL(8), INTENT(IN   )         :: vstr_bot     !! meridional bottom stress [m2/s2]
          REAL(8), INTENT(IN   )         :: z0b          !! bottom roughness length [m]
          REAL(8), INTENT(IN   )         :: tke_n_sfc    !!  [m2/s2]
          REAL(8), INTENT(IN   )         :: Hz_sfc       !!  [m]
          REAL(8), INTENT(IN   )         :: tke_n_bot    !!  [m2/s2]
          REAL(8), INTENT(IN   )         :: Hz_bot       !!  [m]
          REAL(8), INTENT(IN   )         :: OneOverSig   !!  [-]
          REAL(8), INTENT(IN   )         :: tke_min      !!  [m2/s2]
          REAL(8), INTENT(IN   )         :: eps_min      !!  [m2/s3]
          REAL(8), INTENT(IN   )         :: pnm(3)
          REAL(8), INTENT(  OUT)         :: tke_sfc      !! surface value for Dirichlet condition [m2/s2]
          REAL(8), INTENT(  OUT)         :: tke_bot      !! bottom value for Dirichlet condition [m2/s2]
          REAL(8), INTENT(  OUT)         :: ftke_sfc     !! surface TKE flux for Neumann condition [m3/s3]
          REAL(8), INTENT(  OUT)         :: ftke_bot     !! bottom  TKE flux for Neumann condition [m3/s3]
          REAL(8), INTENT(  OUT)         :: eps_sfc      !! surface EPS value for Dirichlet condition [m2/s3]
          REAL(8), INTENT(  OUT)         :: eps_bot      !! bottom  EPS value for Dirichlet condition [m2/s3]
          REAL(8), INTENT(  OUT)         :: feps_sfc     !! surface EPS flux for Neumann condition [m3/s4]
          REAL(8), INTENT(  OUT)         :: feps_bot     !! bottom  EPS flux for Neumann condition [m3/s4]
          !local variables
          REAL(8)                        :: ustar2_sfc, ustar2_bot
          REAL(8)                        :: cm0, cm0inv2, rp, rn, rm
          REAL(8)                        :: z0_s, z0_b, lgthsc
          !
          cm0 =  ( (a2*a2 - 3.0*a3*a3 + 3.0*a1*nn)/(3.0*nn*nn) )**0.25
          rp = pnm(1); rn = pnm(2); rm = pnm(3)
          cm0inv2 = 1./cm0**2
          ! velocity scales
          ustar2_sfc   =  SQRT( ustr_sfc**2+vstr_sfc**2 )
          ustar2_bot   =  SQRT( ustr_bot**2+vstr_bot**2 )
          ! TKE Dirichlet bdy condition
          tke_sfc  = MAX( tke_min, cm0inv2*ustar2_sfc )
          tke_bot  = MAX( tke_min, cm0inv2*ustar2_bot )
          ! TKE Neumann bdy condition
          ftke_sfc = 0.
          ftke_bot = 0.
          !
          z0_s    = MAX( Zosmin , chk*ustar2_sfc )   !<-- Charnock
          lgthsc  = vkarmn*(0.5*Hz_sfc+z0_s)
          !
          eps_sfc = MAX(eps_min,(cm0**rp)*(lgthsc**rn)*(tke_n_sfc**rm))
          !
          feps_sfc = -rn*cm0**(rp+1.)*vkarmn*OneOverSig*(tke_n_sfc**(rm+0.5))*(lgthsc**rn)
          ! bottom
          z0_b = MAX( z0b, Zobmin )
          !
          lgthsc = vkarmn*(0.5*Hz_bot+z0_b)
          eps_bot = MAX(eps_min,(cm0**rp)*(lgthsc**rn)*(tke_n_bot**rm))
          !
          feps_bot = -rn*cm0**(rp+1.)*vkarmn*OneOverSig*(tke_n_bot**(rm+0.5))*(lgthsc**rn)
          !
      !-----------------------------------------------------------------------------
      END SUBROUTINE compute_tke_eps_bdy
      !=============================================================================

      SUBROUTINE compute_shear(u_n, v_n, u_np1, v_np1, zr, N, shear2)
      !!==========================================================================<br />
          !!                  ***  ROUTINE compute_shear  ***                         <br />
          !! ** Purposes : compute shear production term for TKE equation             <br />
          !!==========================================================================<br />
          IMPLICIT NONE
          INTEGER, INTENT(IN   )              :: N                      !! number of vertical levels
          REAL(8), INTENT(IN   )              :: u_n  (1:N),v_n  (1:N)  !! velocity components at time n    [m/s]
          REAL(8), INTENT(IN   )              :: u_np1(1:N),v_np1(1:N)  !! velocity components at time n+1  [m/s]
          REAL(8), INTENT(IN   )              :: zr(1:N)                !! depth at cell centers [m]
          REAL(8), INTENT(  OUT)              :: shear2(0:N)            !! shear production term [m2/s3]
          ! local variables
          INTEGER                             :: k
          REAL(8)                             :: du,dv,cff
          shear2(0:N) = 0.
          DO k=1,N-1
            cff       = 1. / ( zr(k+1)-zr(k) )**2
            du        = cff*( u_np1(k+1)-u_np1(k) )*0.5*( u_n(k+1)+u_np1(k+1)-u_n(k)-u_np1(k) ) !! Shear production term using discretization from Burchard (2002) <br />
            !! \( {\rm Sh}_{k+1/2} = \frac{ 1 }{ \Delta z_{k+1/2}^2 } ( u_{k+1}^n - u_{k}^n ) ( u_{k+1}^{n+1/2} - u_{k}^{n+1/2} )  \)
            dv        = cff*( v_np1(k+1)-v_np1(k) )*0.5*( v_n(k+1)+v_np1(k+1)-v_n(k)-v_np1(k) )
            shear2(k) = du + dv
          ENDDO
      !---------------------------------------------------------------------------------------------------
      END SUBROUTINE compute_shear
      !===================================================================================================
      !
      !===================================================================================================
      SUBROUTINE compute_EV_ED (tke,gls,bvf,shear2,pnm,nuwm,nuws,eps_min,N,Akv,Akt,c_mu,c_mu_prim)
      !---------------------------------------------------------------------------------------------------
          USE scm_par
          IMPLICIT NONE
          INTEGER, INTENT(IN   )  :: N
          REAL(8), INTENT(IN   )  :: tke(0:N)
          REAL(8), INTENT(INOUT)  :: gls(0:N)
          REAL(8), INTENT(IN   )  :: bvf(0:N)
          REAL(8), INTENT(IN   )  :: shear2(0:N)
          REAL(8), INTENT(IN   )  :: pnm(3)
          REAL(8), INTENT(IN   )  :: nuwm
          REAL(8), INTENT(IN   )  :: nuws
          REAL(8), INTENT(IN   )  :: eps_min
          REAL(8), INTENT(  OUT)  :: Akv(0:N)            !! eddy-viscosity [m2/s]
          REAL(8), INTENT(  OUT)  :: Akt(0:N)            !! eddy-diffusivity [m2/s]
          REAL(8), INTENT(  OUT)  :: c_mu(0:N)
          REAL(8), INTENT(  OUT)  :: c_mu_prim(0:N)
          ! local variables
          INTEGER                 :: k
          REAL(8)                 :: L_lim, rp,rn,rm,cm0,cff
          REAL(8)                 :: e1,e2,e3,alpha_n,alpha_m
          REAL(8)                 :: alpha_m_max,alpha_n_min
          REAL(8)                 :: Denom, epsilon
          !
          cm0 =  ( (a2*a2 - 3.0*a3*a3 + 3.0*a1*nn)/(3.0*nn*nn) )**0.25
          rp = pnm(1); rn = pnm(2); rm = pnm(3)
          e1 =  3.0 + 1.*rp / rn; e2 =  1.5 + 1.*rm / rn; e3 = -1.0 / rn
          !
          ! Minimum value of alpha_n to ensure that alpha_m is positive
          !
          alpha_n_min = 0.5*( - (sf_d1+sf_nb0)+sqrt( (sf_d1+sf_nb0)**2   &
                            - 4. * sf_d0 *( sf_d4 + sf_nb1 ) ) )         &
                                            / ( sf_d4 + sf_nb1 )
          !--------------------------------------------------------------------
          !    Compute Akv & Lscale
          !--------------------------------------------------------------------
          c_mu(0:N) = 0.; c_mu_prim(0:N) = 0.
          DO k=1,N-1
            !
            ! Galperin limitation : l <= l_lim
            L_lim = galp * sqrt( 2.* tke(k)) / ( sqrt(max(rsmall, bvf(k)))  )
            !
            ! Limitation on psi (use MAX because rn is negative)
            cff    = (cm0**rp)*(L_lim**rn)*(tke(k)**rm)
            gls(k) = MAX( gls( k ),cff )
            !
            ! Dissipation rate epsilon
            epsilon = (cm0**e1) * ( tke(k)**e2 ) * ( gls(k)**e3 )
            ! Limitation of epsilon
            epsilon = MAX(epsilon, eps_min)
            !
            ! Compute alpha_n and alpha_m
            cff       = ( tke(k)/epsilon )**2
            alpha_m   = cff*  shear2(k)
            alpha_n   = cff*     bvf(k)
            !
            ! Limitation of alpha_n and alpha_m
            alpha_n   = MIN( MAX(0.73*alpha_n_min,alpha_n) , 1.0e10 )
            alpha_m_max = ( lim_am0 + lim_am1 * alpha_n             &
                     + lim_am2 * alpha_n**2+ lim_am3 * alpha_n**3) /  &
               ( lim_am4 + lim_am5 * alpha_n + lim_am6 * alpha_n**2 )
            alpha_m = MIN(alpha_m , alpha_m_max)
            !
            ! Compute stability functions
            Denom = sf_d0  + sf_d1*alpha_n +  sf_d2*alpha_m        &
                             + sf_d3*alpha_n*alpha_m                 &
                             + sf_d4*alpha_n**2 + sf_d5*alpha_m**2
            cff   = 1./Denom
            c_mu(k)      = cff*(sf_n0  + sf_n1 *alpha_n + sf_n2 *alpha_m)
            c_mu_prim(k) = cff*(sf_nb0 + sf_nb1*alpha_n + sf_nb2*alpha_m)
            !
            ! Finalize the computation of Akv and Akt
            cff = tke(k)**2 / epsilon
            Akv(k) = MAX( cff*c_mu(k)      , nuwm )
            Akt(k) = MAX( cff*c_mu_prim(k) , nuws )
          ENDDO
          !
          Akv(N) = 0.; Akv(0) = 0.
          Akt(N) = 0.; Akt(0) = 0.
          !Akv(N) = MAX( 1.5*Akv(N-1)-0.5*Akv(N-2), nuwm)
          !Akv(0) = MAX( 1.5*Akv(  1)-0.5*Akv(  2), nuwm)
          !Akt(N) = MAX( 1.5*Akt(N-1)-0.5*Akt(N-2), nuws)
          !Akt(0) = MAX( 1.5*Akt(  1)-0.5*Akt(  2), nuws)
      !---------------------------------------------------------------------------------------------------
      END SUBROUTINE compute_EV_ED
      !===================================================================================================





      !===================================================================================================
      SUBROUTINE compute_EV_ED_filt(tke,gls,bvf,shear2,pnm,nuwm,nuws,eps_min,N,Akv,Akt,c_mu,c_mu_prim)
      !---------------------------------------------------------------------------------------------------
          USE scm_par
          IMPLICIT NONE
          INTEGER, INTENT(IN   )  :: N
          REAL(8), INTENT(IN   )  :: tke(0:N)
          REAL(8), INTENT(INOUT)  :: gls(0:N)
          REAL(8), INTENT(IN   )  :: bvf(0:N)
          REAL(8), INTENT(IN   )  :: shear2(0:N)
          REAL(8), INTENT(IN   )  :: pnm(3)
          REAL(8), INTENT(IN   )  :: nuwm
          REAL(8), INTENT(IN   )  :: nuws
          REAL(8), INTENT(IN   )  :: eps_min
          REAL(8), INTENT(  OUT)  :: Akv(0:N)            !! eddy-viscosity [m2/s]
          REAL(8), INTENT(  OUT)  :: Akt(0:N)            !! eddy-diffusivity [m2/s]
          REAL(8), INTENT(  OUT)  :: c_mu(0:N)
          REAL(8), INTENT(  OUT)  :: c_mu_prim(0:N)
          ! local variables
          INTEGER                 :: k
          REAL(8)                 :: L_lim, rp,rn,rm,cm0,cff,c_filt(2:N-2)
          REAL(8)                 :: e1,e2,e3,alpha_n,alpha_m,aN(0:N),aM(0:N)
          REAL(8)                 :: alpha_m_max,alpha_n_min,epsil(1:N-1)
          REAL(8)                 :: Denom, epsilon
          REAL(8), PARAMETER      :: filter_cof = 0.5
          !
          cm0 =  ( (a2*a2 - 3.0*a3*a3 + 3.0*a1*nn)/(3.0*nn*nn) )**0.25
          rp = pnm(1); rn = pnm(2); rm = pnm(3)
          e1 =  3.0 + 1.*rp / rn; e2 =  1.5 + 1.*rm / rn; e3 = -1.0 / rn
          !
          ! Minimum value of alpha_n to ensure that alpha_m is positive
          !
          alpha_n_min = 0.5*( - (sf_d1+sf_nb0)+sqrt( (sf_d1+sf_nb0)**2   &
                            - 4. * sf_d0 *( sf_d4 + sf_nb1 ) ) )         &
                                            / ( sf_d4 + sf_nb1 )
          !--------------------------------------------------------------------
          !    Compute Akv & Lscale
          !--------------------------------------------------------------------
          c_mu(0:N) = 0.; c_mu_prim(0:N) = 0.
          DO k=1,N-1
            !
            ! Galperin limitation : l <= l_lim
            L_lim = galp * sqrt( 2.* tke(k)) / ( sqrt(max(1.e-14, bvf(k)))  )
            !
            ! Limitation on psi (use MAX because rn is negative)
            cff    = (cm0**rp)*(L_lim**rn)*(tke(k)**rm)
            gls(k) = MAX( gls( k ),cff )
            !
            epsilon = (cm0**e1) * ( tke(k)**e2 ) * ( gls(k)**e3 )
            epsilon = MAX(epsilon, eps_min)
            ! Compute alpha_n and alpha_m
            cff       = ( tke(k)/epsilon )**2
            aM(k)     = cff*  shear2(k)
            aN(k)     = cff*     bvf(k)
          ENDDO
          aM(N) = 0.; aN(N) = 0.
          aM(0) = 0.; aN(0) = 0.
          !
          DO k=1,N-1
            !
            alpha_n = aN(k)+filter_cof*(0.5*aN(k+1)-aN(k)+0.5*aN(k-1))
            alpha_m = aM(k)+filter_cof*(0.5*aM(k+1)-aM(k)+0.5*aM(k-1))
            ! Limitation of alpha_n and alpha_m
            alpha_n   = MIN( MAX(0.73*alpha_n_min,alpha_n) , 1.0e10 )
            alpha_m_max = ( lim_am0 + lim_am1 * alpha_n             &
                     + lim_am2 * alpha_n**2+ lim_am3 * alpha_n**3 ) /  &
               ( lim_am4 + lim_am5 * alpha_n + lim_am6 * alpha_n**2 )
            alpha_m = MIN(alpha_m , alpha_m_max)
            !
            ! Compute stability functions
            Denom = sf_d0  + sf_d1*alpha_n +  sf_d2*alpha_m        &
                             + sf_d3*alpha_n*alpha_m                 &
                             + sf_d4*alpha_n**2 + sf_d5*alpha_m**2
            cff   = 1./Denom
            c_mu(k)      = cff*(sf_n0  + sf_n1 *alpha_n + sf_n2 *alpha_m)
            c_mu_prim(k) = cff*(sf_nb0 + sf_nb1*alpha_n + sf_nb2*alpha_m)
          ENDDO
          !
          DO k = 1,N-1
            epsilon = (cm0**e1) * ( tke(k)**e2 ) * ( gls(k)**e3 )
            epsilon = MAX(epsilon, eps_min)
            ! Finalize the computation of Akv and Akt
            cff = tke(k)**2 / epsilon
            Akv(k) = MAX( cff*c_mu(k)      , nuwm )
            Akt(k) = MAX( cff*c_mu_prim(k) , nuws )
          ENDDO
          !
          !Akv(N) = 0.; Akv(0) = 0.
          !Akt(N) = 0.; Akt(0) = 0.
          Akv(N) = MAX( 1.5*Akv(N-1)-0.5*Akv(N-2), nuwm)
          Akv(0) = MAX( 1.5*Akv(  1)-0.5*Akv(  2), nuwm)
          Akt(N) = MAX( 1.5*Akt(N-1)-0.5*Akt(N-2), nuws)
          Akt(0) = MAX( 1.5*Akt(  1)-0.5*Akt(  2), nuws)
      !---------------------------------------------------------------------------------------------------
      END SUBROUTINE compute_EV_ED_filt
      !===================================================================================================






      !===================================================================================================
      SUBROUTINE advance_turb_tke (tke_n,bvf,shear2,Ak_tke,Akv,Akt,diss,Hz,dt,tke_min,bdy_sfc,bdy_bot,N,tke_np1,wtke)
      !---------------------------------------------------------------------------------------------------
          USE scm_par
          IMPLICIT NONE
          INTEGER, INTENT(IN   )  :: N
          REAL(8), INTENT(IN   )  :: tke_n (0:N)
          REAL(8), INTENT(IN   )  :: bvf   (0:N)
          REAL(8), INTENT(IN   )  :: shear2(0:N)
          REAL(8), INTENT(IN   )  :: Ak_tke(0:N)
          REAL(8), INTENT(IN   )  :: Akv   (0:N)
          REAL(8), INTENT(IN   )  :: Akt   (0:N)
          REAL(8), INTENT(IN   )  :: diss  (0:N)
          REAL(8), INTENT(IN   )  :: Hz    (1:N)
          REAL(8), INTENT(IN   )  :: dt
          REAL(8), INTENT(IN   )  :: tke_min
          REAL(8), INTENT(IN   )  :: bdy_sfc(2)       ! bdy_sfc(1) = 0. -> Dirichlet, bdy_sfc(1) = 1. -> Neumann
          REAL(8), INTENT(IN   )  :: bdy_bot(2)       ! bdy_bot(1) = 0. -> Dirichlet, bdy_top(1) = 1. -> Neumann
          REAL(8), INTENT(  OUT)  :: tke_np1(0:N)
          REAL(8), INTENT(  OUT)  :: wtke   (1:N)
          !local variables
          INTEGER                 :: k
          REAL(8)                 :: cff,Sprod,Bprod
          REAL(8)                 :: aa(0:N),bb(0:N),cc(0:N),rhs(0:N),q(0:N),invG
          !==================================
          ! fill in the tridiagonal matrix
          !==================================
          ! Off-diagonal terms for the tridiagonal problem
          cff=-0.5*dt
          DO k = 1,N-1
            aa(k) = cff*( Ak_tke(k)+Ak_tke(k-1) ) / Hz(k  )
            cc(k) = cff*( Ak_tke(k)+Ak_tke(k+1) ) / Hz(k+1)
          ENDDO
          !
          DO k = 1, N-1
            ! Shear and buoyancy production
            Sprod =  Akv(k)*shear2(k)
            Bprod = -Akt(k)*bvf(k)
            ! Diagonal and rhs term
            cff   = 0.5*( Hz(k) + Hz(k+1) )
            invG  = 1./tke_n(k)
            IF( (Bprod + Sprod) .gt. 0.) THEN
              rhs(k) = cff*( tke_n(k) + dt*(Bprod+Sprod) )
              bb (k) = cff*( 1.+dt*diss(k)*invG )-aa(k)-cc(k)
            ELSE
              rhs(k) = cff*( tke_n(k) + dt*Sprod  )
              bb (k) = cff*(1.+dt*(diss(k)-Bprod)*invG)-aa(k)-cc(k)
            ENDIF
          ENDDO
          IF(bdy_sfc(1)<0.5) THEN
            ! Dirichlet boundary condition (surface)
            aa (N) = 0.; bb(N) = 1.; cc(N) = 0.
            rhs(N) = bdy_sfc(2)
          ELSE
            ! Neumann boundary condition (surface)
            aa(N)  = -0.5*(Ak_tke(N)+Ak_tke(N-1))
            bb(N)  = +0.5*(Ak_tke(N)+Ak_tke(N-1))
            cc(N)  = 0.
            rhs(N) = Hz(N)*bdy_sfc(2)
          ENDIF
          IF(bdy_bot(1)<0.5) THEN
            ! Dirichlet boundary condition (bottom)
            aa(0) = 0.; bb(0) = 1.; cc(0) = 0.
            rhs(0) = bdy_bot(2)
          ELSE
            ! Neumann boundary condition (bottom)
            aa(0)  = 0.
            bb(0)  = -0.5*(Ak_tke(0)+Ak_tke(1))
            cc(0)  =  0.5*(Ak_tke(0)+Ak_tke(1))
            rhs(0) = Hz(1)*bdy_bot(2)
          ENDIF
          ! Solve tridiagonal problem
          cff    = 1./bb(0)
          q  (0) = - cc(0)*cff
          rhs(0) =  rhs(0)*cff
          DO k=1,N
            cff=1./(bb(k)+aa(k)*q(k-1))
            q(k)= -cff*cc(k)
            rhs(k)=cff*( rhs(k)-aa(k)*rhs(k-1) )
          ENDDO
          DO k=N-1,0,-1
            rhs(k)=rhs(k)+q(k)*rhs(k+1)
          ENDDO
          tke_np1(0:N) = MAX( rhs(0:N), tke_min )
          !
          DO k = 1,N
            wtke(k) = - 0.5*(Ak_tke(k)+Ak_tke(k-1))*(tke_np1(k)-tke_np1(k-1))/Hz(k)
          ENDDO
      !---------------------------------------------------------------------------------------------------
      END SUBROUTINE advance_turb_tke
      !===================================================================================================



      !===================================================================================================
      SUBROUTINE advance_turb_eps (eps_n,bvf,shear2,Ak_eps,c_mu,c_mu_prim,tke_n,tke_np1,Hz,dt,  &
                                            beta,eps_min,bdy_sfc,bdy_bot,N,eps_np1)
      !---------------------------------------------------------------------------------------------------
          USE scm_par
          IMPLICIT NONE
          INTEGER, INTENT(IN   )  :: N
          REAL(8), INTENT(IN   )  :: eps_n     (0:N)
          REAL(8), INTENT(IN   )  :: bvf       (0:N)
          REAL(8), INTENT(IN   )  :: shear2    (0:N)
          REAL(8), INTENT(IN   )  :: Ak_eps    (0:N)
          REAL(8), INTENT(IN   )  :: c_mu      (0:N)
          REAL(8), INTENT(IN   )  :: c_mu_prim (0:N)
          REAL(8), INTENT(IN   )  :: tke_n     (0:N)
          REAL(8), INTENT(IN   )  :: tke_np1   (0:N)
          REAL(8), INTENT(IN   )  :: Hz    (1:N)
          REAL(8), INTENT(IN   )  :: beta(4)
          REAL(8), INTENT(IN   )  :: dt
          REAL(8), INTENT(IN   )  :: eps_min
          REAL(8), INTENT(IN   )  :: bdy_sfc(2)       ! bdy_sfc(1) = 0. -> Dirichlet, bdy_sfc(1) = 1. -> Neumann
          REAL(8), INTENT(IN   )  :: bdy_bot(2)       ! bdy_bot(1) = 0. -> Dirichlet, bdy_top(1) = 1. -> Neumann
          REAL(8), INTENT(  OUT)  :: eps_np1(0:N)
          !local variables
          INTEGER                 :: k
          REAL(8)                 :: cff,Sprod,Bprod,beta1,beta2,beta3m,beta3p
          REAL(8)                 :: aa(0:N),bb(0:N),cc(0:N),rhs(0:N),q(0:N),invG
          !
          beta1  = beta(1); beta2  = beta(2)
          beta3m = beta(3); beta3p = beta(4)
          !==================================
          ! fill in the tridiagonal matrix
          !==================================
          ! Off-diagonal terms for the tridiagonal problem
          cff=-0.5*dt
          DO k = 1,N-1
            aa(k) = cff*( Ak_eps(k)+Ak_eps(k-1) ) / Hz(k  )
            cc(k) = cff*( Ak_eps(k)+Ak_eps(k+1) ) / Hz(k+1)
          ENDDO
          !
          DO k = 1, N-1
            ! Shear and buoyancy production
            Sprod =  beta1*c_mu(k)*tke_n(k)*shear2(k)
            Bprod = -c_mu_prim(k)*tke_n(k)*( beta3m*MAX(bvf(k),0.) + beta3p*MIN(bvf(k),0.) )
            ! Diagonal and rhs term
            cff   = 0.5*( Hz(k) + Hz(k+1) )
            IF( (Bprod + Sprod) .gt. 0.) THEN
              rhs(k) = cff*( eps_n(k) + dt*(Bprod+Sprod) )
              bb (k) = cff*( 1.+dt*beta2*eps_n(k)/tke_np1(k) )-aa(k)-cc(k)
            ELSE
              invG   = 1./eps_n(k)
              rhs(k) = cff*( eps_n(k) + dt*Sprod )
              bb (k) = cff*(1.+dt*beta2*eps_n(k)/tke_np1(k)-dt*invG*Bprod)-aa(k)-cc(k)
            ENDIF
          ENDDO
          IF(bdy_sfc(1)<0.5) THEN
            ! Dirichlet boundary condition (surface)
            aa (N) = 0.; bb(N) = 1.; cc(N) = 0.
            rhs(N) = bdy_sfc(2)
          ELSE
            ! Neumann boundary condition (surface)
            aa(N)  = -0.5*(Ak_eps(N)+Ak_eps(N-1))
            bb(N)  = +0.5*(Ak_eps(N)+Ak_eps(N-1))
            cc(N)  = 0.
            rhs(N) = Hz(N)*bdy_sfc(2)
          ENDIF
          IF(bdy_bot(1)<0.5) THEN
            ! Dirichlet boundary condition (bottom)
            aa(0) = 0.; bb(0) = 1.; cc(0) = 0.
            rhs(0) = bdy_bot(2)
          ELSE
            ! Neumann boundary condition (bottom)
            aa(0)  = 0.
            bb(0)  = -0.5*(Ak_eps(0)+Ak_eps(1))
            cc(0)  =  0.5*(Ak_eps(0)+Ak_eps(1))
            rhs(0) = Hz(1)*bdy_bot(2)
          ENDIF
          ! Solve tridiagonal problem
          cff    = 1./bb(0)
          q  (0) = - cc(0)*cff
          rhs(0) =  rhs(0)*cff
          DO k=1,N
            cff=1./(bb(k)+aa(k)*q(k-1))
            q(k)= -cff*cc(k)
            rhs(k)=cff*( rhs(k)-aa(k)*rhs(k-1) )
          ENDDO
          DO k=N-1,0,-1
            rhs(k)=rhs(k)+q(k)*rhs(k+1)
          ENDDO
          eps_np1(0:N) = MAX( rhs(0:N), eps_min )
      !---------------------------------------------------------------------------------------------------
    END SUBROUTINE advance_turb_eps
      !===================================================================================================

      END MODULE scm_gls

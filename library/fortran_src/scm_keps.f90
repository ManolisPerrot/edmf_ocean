      MODULE scm_keps
      !
        IMPLICIT NONE
      !
      CONTAINS
      !++
      SUBROUTINE compute_tke_eps_bdy(ustr_sfc, vstr_sfc, ustr_bot, vstr_bot, z0b, z0_s, &
                                       tke_n_sfc, Hz_sfc, tke_n_bot, Hz_bot, OneOverSig,  &
                                       tke_min,   eps_min,                                &
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
            REAL(8), INTENT(IN   )         :: z0_s         !! surface roughness length [m]
            REAL(8), INTENT(IN   )         :: tke_n_sfc    !!  [m2/s2]
            REAL(8), INTENT(IN   )         :: Hz_sfc       !!  [m]
            REAL(8), INTENT(IN   )         :: tke_n_bot    !!  [m2/s2]
            REAL(8), INTENT(IN   )         :: Hz_bot       !!  [m]
            REAL(8), INTENT(IN   )         :: OneOverSig   !!  [-]
            REAL(8), INTENT(IN   )         :: tke_min      !!  [m2/s2]
            REAL(8), INTENT(IN   )         :: eps_min
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
            REAL(8)                        :: z0_b, lgthsc
            REAL(8), PARAMETER             :: rm =  1.5
            REAL(8), PARAMETER             :: rn = -1.0
            REAL(8), PARAMETER             :: rp =  3.0
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
            lgthsc  = vkarmn*(Hz_sfc+z0_s)
            eps_sfc = MAX(eps_min,(cm0**3)*tke_n_sfc*SQRT(tke_n_sfc)/lgthsc)
            !
            lgthsc  = vkarmn*(0.5*Hz_sfc+z0_s)
            feps_sfc = (cm0**4)*tke_n_sfc*tke_n_sfc*vkarmn*OneOverSig/lgthsc
            ! bottom
            z0_b = MAX( z0b, Zobmin )
            !
            lgthsc = vkarmn*(Hz_bot+z0_b)
            eps_bot = MAX(eps_min,(cm0**3)*tke_n_bot*SQRT(tke_n_bot)/lgthsc)
            !
            feps_bot = (cm0**4)*tke_n_bot*tke_n_bot*vkarmn*OneOverSig/lgthsc
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
            !! \( {\rm Sh}_{k+1/2} = \frac{ (K_m)_{k+1/2} }{ \Delta z_{k+1/2}^2 } ( u_{k+1}^n - u_{k}^n ) ( u_{k+1}^{n+1/2} - u_{k}^{n+1/2} )  \)
            dv        = cff*( v_np1(k+1)-v_np1(k) )*0.5*( v_n(k+1)+v_np1(k+1)-v_n(k)-v_np1(k) )
            shear2(k) = du + dv
          ENDDO
          !
          shear2(0) = shear2(  1)
          shear2(N) = shear2(N-1)
          !
      !---------------------------------------------------------------------------------------------------
      END SUBROUTINE compute_shear
      !===================================================================================================



      !===================================================================================================
      SUBROUTINE compute_EV_ED_eq (tke,epsil,bvf,nuwm,nuws,N,Akv,Akt,c_mu,c_mu_prim,Lscale)
      !---------------------------------------------------------------------------------------------------
            USE scm_par
            IMPLICIT NONE
            INTEGER, INTENT(IN   )  :: N
            REAL(8), INTENT(IN   )  :: tke(0:N)
            REAL(8), INTENT(INOUT)  :: epsil(0:N)
            REAL(8), INTENT(IN   )  :: bvf(0:N)
            REAL(8), INTENT(IN   )  :: nuwm
            REAL(8), INTENT(IN   )  :: nuws
            REAL(8), INTENT(  OUT)  :: Akv(0:N)            !! eddy-viscosity [m2/s]
            REAL(8), INTENT(  OUT)  :: Akt(0:N)            !! eddy-diffusivity [m2/s]
            REAL(8), INTENT(  OUT)  :: c_mu(0:N)
            REAL(8), INTENT(  OUT)  :: c_mu_prim(0:N)
            REAL(8), INTENT(  OUT)  :: Lscale(0:N)
            ! local variables
            INTEGER                 :: k
            REAL(8)                 :: dCm,nCm,nCmp
            REAL(8)                 :: cff,cff0,cff1,cff2,tau2
            REAL(8)                 :: Denom, epsilon, aN(0:N),aS(0:N)
            !
            c_mu(0:N) = 0.; c_mu_prim(0:N) = 0.
            !
            DO k=0,N
               tau2   = tke(k)*tke(k) / ( epsil(k)*epsil(k) )
               aN(k)  = tau2 * bvf(k)
            ENDDO
            !
            IF ( ABS(n2-d5) .lt. small ) THEN
              DO k=1,N-1
                aN(k) = max(aN(k),anLimitFact*anMin)
                cff0  = -d0 - (d1 + nt0)*an(k) - (d4 + nt1)*an(k)*an(k)
                cff1  = -d2 + n0 +  (n1-d3-nt2)*an(k)
                aS(k) =  - cff0/cff1
                dCm  = d0  +  d1*an(k) +  d2*aS(k) + d3*an(k)*as(k) + d4*an(k)*an(k) + d5*aS(k)*aS(k)
                nCm  = n0  +  n1*an(k) +  n2*aS(k)
                nCmp = nt0 + nt1*an(k) + nt2*aS(k)
                c_mu     (k) =  cm3_inv*nCm /dCm
                c_mu_prim(k) =  cm3_inv*nCmp/dCm
              ENDDO
            ELSE
              DO k=1,N-1
                an(k) = max(aN(k),anLimitFact*anMin)
  !          compute the equilibrium value of as
                cff0  = -d0 - (d1 + nt0)*aN(k) - (d4 + nt1)*aN(k)*aN(k)
                cff1  = -d2 + n0 + (n1-d3-nt2)*aN(k)
                cff2  =  n2-d5
                aS(k) =  (-cff1 + sqrt(cff1*cff1-4.*cff0*cff2) ) / (2.*cff2)
   !          compute stability function
                dCm  = d0  +  d1*aN(k) +  d2*aS(k) + d3*an(k)*as(k) + d4*an(k)*an(k) + d5*aS(k)*aS(k)
                nCm  = n0  +  n1*aN(k) +  n2*aS(k)
                nCmp = nt0 + nt1*aN(k) + nt2*aS(k)
                c_mu     (k) =  cm3_inv*nCm /dCm
                c_mu_prim(k) =  cm3_inv*nCmp/dCm
              ENDDO
            ENDIF
            !
            !c_mu     (1:N-1) = 0.25
            !c_mu_prim(1:N-1) = 0.25
            !
            DO k=1,N-1
            ! compute dissipative scale
               Lscale(k) = cm3*SQRT(tke(k)*tke(k)*tke(k))/epsil(k)
               cff       =  SQRT(tke(k))*Lscale(k)
               Akv(k)    =  cff*c_mu(k) + avmolU
               Akt(k)    =  cff*c_mu_prim(k) + avmolT
            ENDDO
            !
            Akv(N) = 0.; Akv(0) = 0.
            Akt(N) = 0.; Akt(0) = 0.
            Lscale(N) = 0.; Lscale(0) = 0.
            !
        !---------------------------------------------------------------------------------------------------
        END SUBROUTINE compute_EV_ED_eq
        !===================================================================================================


        !===================================================================================================
        SUBROUTINE compute_EV_ED_weak (tke,epsil,shear2,bvf,nuwm,nuws,N,Akv,Akt,c_mu,c_mu_prim,Lscale)
        !---------------------------------------------------------------------------------------------------
            USE scm_par
            IMPLICIT NONE
            INTEGER, INTENT(IN   )  :: N
            REAL(8), INTENT(IN   )  :: tke(0:N)
            REAL(8), INTENT(INOUT)  :: epsil(0:N)
            REAL(8), INTENT(IN   )  :: shear2(0:N)
            REAL(8), INTENT(IN   )  :: bvf(0:N)
            REAL(8), INTENT(IN   )  :: nuwm
            REAL(8), INTENT(IN   )  :: nuws
            REAL(8), INTENT(  OUT)  :: Akv(0:N)            !! eddy-viscosity [m2/s]
            REAL(8), INTENT(  OUT)  :: Akt(0:N)            !! eddy-diffusivity [m2/s]
            REAL(8), INTENT(  OUT)  :: c_mu(0:N)
            REAL(8), INTENT(  OUT)  :: c_mu_prim(0:N)
            REAL(8), INTENT(  OUT)  :: Lscale(0:N)
            ! local variables
            INTEGER                 :: k
            REAL(8)                 :: dCm,nCm,nCmp,asMax
            REAL(8)                 :: cff,cff0,cff1,cff2,tau2
            REAL(8)                 :: asMaxNum, asMaxDen, aN(0:N),aS(0:N)
            !
            c_mu(0:N) = 0.; c_mu_prim(0:N) = 0.
            !
            DO k=0,N
               tau2   = tke(k)*tke(k) / ( epsil(k)*epsil(k) )
               aN(k)  = tau2 * bvf(k)
               aS(k)  = tau2 * shear2(k)
            ENDDO
            !
            DO k=1,N-1
              aN(k) = MAX(aN(k),anLimitFact*anMin)
              ! compute maximum as for shear instability
              asMaxNum  = d0*n0 + (d0*n1+d1*n0)*aN(k) + (d1*n1+d4*n0)*aN(k)*aN(k) + d4*n1*aN(k)*aN(k)*aN(k)
              asMaxDen  = d2*n0 + (d2*n1+d3*n0)*aN(k) + d3*n1*aN(k)*aN(k)
              asMax     = asMaxNum / asMaxDen
              !       clip as at miximum value
              aS(k) = MIN(aS(k),asLimitFact*asMax)
              !
              dCm  = d0  +  d1*aN(k) +  d2*aS(k) + d3*aN(k)*aS(k) + d4*aN(k)*aN(k) + d5*aS(k)*aS(k)
              nCm  = n0  +  n1*aN(k) +  n2*aS(k)
              nCmp = nt0 + nt1*aN(k) + nt2*aS(k)
              !
              c_mu     (k) =  cm3_inv*nCm /dCm
              c_mu_prim(k) =  cm3_inv*nCmp/dCm
            ENDDO
            !
            DO k=1,N-1
            ! compute dissipative scale
               Lscale(k) = cm3*SQRT(tke(k)*tke(k)*tke(k))/epsil(k)
               cff       =  SQRT(tke(k))*Lscale(k)
               Akv(k)    =  cff*c_mu(k) + avmolU
               Akt(k)    =  cff*c_mu_prim(k) + avmolT
            ENDDO
            !
            Akv(N) = 0.; Akv(0) = 0.
            Akt(N) = 0.; Akt(0) = 0.
            Lscale(N) = 0.; Lscale(0) = 0.
        !---------------------------------------------------------------------------------------------------
        END SUBROUTINE compute_EV_ED_weak
        !===================================================================================================




        !===================================================================================================
          SUBROUTINE dissipationeq (eps_n,NN,SS,num,nuh,tke_n,h,dt,eps_min,bdy_sfc,bdy_bot,N,eps_np1)
          !---------------------------------------------------------------------------------------------------
              USE scm_par
              IMPLICIT NONE
              INTEGER, INTENT(IN   )  :: N
              REAL(8), INTENT(IN   )  :: eps_n (0:N)
              REAL(8), INTENT(IN   )  :: NN    (0:N)
              REAL(8), INTENT(IN   )  :: SS    (0:N)
              REAL(8), INTENT(IN   )  :: num   (0:N)
              REAL(8), INTENT(IN   )  :: nuh   (0:N)
              REAL(8), INTENT(IN   )  :: tke_n (0:N)
              REAL(8), INTENT(IN   )  :: h     (1:N)
              REAL(8), INTENT(IN   )  :: dt
              REAL(8), INTENT(IN   )  :: eps_min
              REAL(8), INTENT(IN   )  :: bdy_sfc(3)       ! bdy_sfc(1) = 0. -> Dirichlet, bdy_sfc(1) = 1. -> Neumann
              REAL(8), INTENT(IN   )  :: bdy_bot(3)       ! bdy_bot(1) = 0. -> Dirichlet, bdy_top(1) = 1. -> Neumann
              REAL(8), INTENT(  OUT)  :: eps_np1(0:N)
              !local variables
              INTEGER                 :: k,eps_ubc,eps_lbc
              REAL(8)                 :: avh(0:N), prod, buoyan, diss, EpsOverTke,ce3
              REAL(8)                 :: DiffEpsup, DiffEpsdw, Lsour(0:N),Qsour(0:N), cnpar
              !==================================
              Lsour(0:N) = 0.; Qsour(0:N) = 0.; avh(0:N) = 0.
              eps_np1(0:N) = eps_n(0:N)
              !==================================
              DO k = 1,N-1
                avh(k)     = num(k)/sig_e
                EpsOverTke = eps_n(k)/tke_n(k)
                prod       = ce1*EpsOverTke*num(k)*SS(k)
                buoyan     = -nuh(k)*EpsOverTke*NN(k)
                ce3        = ce3minus
                IF (buoyan .gt. 0.) ce3=ce3plus
                diss       = ce2*EpsOverTke*eps_n(k)
                IF (prod+buoyan.gt.0) then
                  Qsour(k) = prod+ce3*buoyan
                  Lsour(k) = -diss/eps_n(k)
                ELSE
                  Qsour(k) = prod
                  Lsour(k) = -(diss-ce3*buoyan)/eps_n(k)
                ENDIF
              ENDDO
              !  do diffusion step
              eps_ubc   = int(bdy_sfc(1))
              eps_lbc   = int(bdy_bot(1))
              DiffEpsup = bdy_sfc(2)
              DiffEpsdw = bdy_bot(2)
              !
              cnpar = 1.0
              CALL diff_face(N,dt,cnpar,h,eps_ubc,eps_lbc,                          &
                                DiffEpsup,DiffEpsdw,avh,Lsour,Qsour,eps_np1)

              !  fill top and bottom value with something nice
              !  (only for output)
              eps_np1(N   )  = bdy_sfc(3)
              eps_np1(0   )  = bdy_bot(3)
              !  clip at k_min
              DO k=0,N
                eps_np1(k) = MAX(eps_np1(k),eps_min)
              ENDDO
          !---------------------------------------------------------------------------------------------------
          END SUBROUTINE dissipationeq
          !===================================================================================================







          !===================================================================================================
          SUBROUTINE tkeeq (tke_n,NN,SS,num,nuh,eps_n,h,dt,tke_min,bdy_sfc,bdy_bot,N,tke_np1,wtke)
          !---------------------------------------------------------------------------------------------------
              USE scm_par
              IMPLICIT NONE
              INTEGER, INTENT(IN   )  :: N
              REAL(8), INTENT(IN   )  :: tke_n (0:N)
              REAL(8), INTENT(IN   )  :: NN    (0:N)
              REAL(8), INTENT(IN   )  :: SS    (0:N)
              REAL(8), INTENT(IN   )  :: num   (0:N)
              REAL(8), INTENT(IN   )  :: nuh   (0:N)
              REAL(8), INTENT(IN   )  :: eps_n (0:N)
              REAL(8), INTENT(IN   )  :: h     (1:N)
              REAL(8), INTENT(IN   )  :: dt
              REAL(8), INTENT(IN   )  :: tke_min
              REAL(8), INTENT(IN   )  :: bdy_sfc(3)       ! bdy_sfc(1) = 0. -> Dirichlet, bdy_sfc(1) = 1. -> Neumann
              REAL(8), INTENT(IN   )  :: bdy_bot(3)       ! bdy_bot(1) = 0. -> Dirichlet, bdy_top(1) = 1. -> Neumann
              REAL(8), INTENT(  OUT)  :: tke_np1(0:N)
              REAL(8), INTENT(  OUT)  :: wtke(1:N)
              !local variables
              INTEGER                 :: k,k_ubc,k_lbc
              REAL(8)                 :: avh(0:N), prod, buoyan, diss
              REAL(8)                 :: DiffKup, DiffKdw, Lsour(0:N),Qsour(0:N), cnpar
              !==================================
              Lsour(0:N) = 0.; Qsour(0:N) = 0.; avh(0:N) = 0.
              tke_np1(0:N) = tke_n(0:N)
              !==================================
              DO k = 1,N-1
                avh(k)  =  num(k)/sig_k
                prod    =  num(k)*SS(k)
                buoyan  = -nuh(k)*NN(k)
                diss    =  eps_n(k)
                IF (prod+buoyan.gt.0) then
                  Qsour(k) = prod+buoyan
                  Lsour(k) = -diss/tke_n(k)
                ELSE
                  Qsour(k) = prod
                  Lsour(k) = -(diss-buoyan)/tke_n(k)
                ENDIF
              ENDDO
              !  do diffusion step
              k_ubc   = int(bdy_sfc(1))
              k_lbc   = int(bdy_bot(1))
              DiffKup = bdy_sfc(2)
              DiffKdw = bdy_bot(2)
              !
              cnpar = 1.0
              call diff_face(N,dt,cnpar,h,k_ubc,k_lbc,                          &
                                DiffKup,DiffKdw,avh,Lsour,Qsour,tke_np1)

              !  fill top and bottom value with something nice
              !  (only for output)
              tke_np1(N   )  = bdy_sfc(3)
              tke_np1(0   )  = bdy_bot(3)
              !  clip at k_min
              DO k=0,N
                tke_np1(k) = MAX(tke_np1(k),tke_min)
              ENDDO
              !
              DO k=1,N
                wtke(k) = - 0.5*(avh(k)+avh(k-1))*( tke_np1(k)-tke_np1(k-1) )/h(k)
              ENDDO
          !---------------------------------------------------------------------------------------------------
          END SUBROUTINE tkeeq
          !===================================================================================================







          !===================================================================================================
          SUBROUTINE diff_face(N,dt,cnpar,h,Bcup,Bcdw,Yup,Ydw,nuY,Lsour,Qsour,Y)
          !---------------------------------------------------------------------------------------------------
              IMPLICIT NONE
              INTEGER, INTENT(IN   )  :: N
              REAL(8), INTENT(IN   )  :: dt
              REAL(8), INTENT(IN   )  :: cnpar                ! level of implicitness
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

                cu(k) =-cnpar*c
                au(k) =-cnpar*a
                bu(k) = 1. + cnpar*(a + c) - l
                du(k) = (1. - (1.-cnpar)*(a        + c        ) )*Y(k)         &
                            + (1.-cnpar)*(a*Y(k-1) + c*Y(k+1) ) + dt*Qsour(k)
              ENDDO
              !   set up upper boundary condition
              SELECT CASE(Bcup)
                CASE(Neumann)
                  a = dt*( nuY(N-1) + nuY(N-2) )  / ( h(N-1)+h(N) ) / h(N-1)
                  l = dt*Lsour(N-1)
                  au(N-1) =-cnpar*a
                  bu(N-1) = 1. + cnpar*a - l
                  du(N-1) = (1. - (1.-cnpar)*a)*Y(N-1)                  &
                          + (1. - cnpar)*a*Y(N-2) + dt*Qsour(N-1)       &
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
                  cu(1) =-cnpar*c
                  bu(1) = 1. + cnpar*c - l
                  du(1) = (1. - (1.-cnpar)*c)*Y(1)                      &
                          + (1. - cnpar)*c*Y(2)  + dt*Qsour(1)          &
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



    END MODULE scm_keps

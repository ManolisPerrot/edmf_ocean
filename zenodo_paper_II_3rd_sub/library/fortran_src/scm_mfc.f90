MODULE scm_mfc
  !!============================================================================<br />
  !!                       ***  MODULE  scm_mfc  ***                            <br />
  !!        Mass-flux (MF) convection scheme for oceanic deep convection        <br />
  !!============================================================================<br />
  !!--------------------------------------------------------------------------<br />
  !!   compute_MF_bdy    : top boundary conditions for plume properties       <br />
  !!   compute_MF_forcing: compute mean quantities used to force the MF equations <br />
  !!   compute_tripleCorr: compute contribution of MF terms to turbulent flux w'e <br />
  !!   eos_val_lin   : linear equation of state                           <br />
  !!   get_w_p_R10   : compute plume vertical velocity for Rio et al. 2010 closure <br />
  !!   get_a_p_R10   : compute plume fractional area for Rio et al. 2010 closure   <br />
  !!   get_t_p_R10   : compute plume tracer properties for Rio et al. 2010 closure <br />
  !!   get_w_p_P09   : compute plume vertical velocity for Pergaud et al. 2009 closure <br />
  !!   get_a_p_P09   : compute plume fractional area for Pergaud et al. 2009 closure   <br />
  !!   get_t_p_P09   : compute plume tracer properties for Pergaud et al. 2009 closure <br />
  !!   mass_flux_R10 : computation of mass-flux and plume properties for Rio et al. 2010 closure<br />
  !!   mass_flux_P09 : computation of mass-flux and plume properties for Pergaud et al. 2009 closure<br />
  !!----------------------------------------------------------------------------<br />
  IMPLICIT NONE

CONTAINS
  !===================================================================================================
  SUBROUTINE compute_MF_bdy(u_m,v_m,t_m,tke,Hz,ntra,npts,up0,vp0,tp0)
  !---------------------------------------------------------------------------------------------------
    !!==========================================================================<br />
    !!                  ***  ROUTINE compute_MF_bdy  ***                        <br />
    !! ** Purposes : compute top initial condition for mass flux equation       <br />
    !!==========================================================================<br />
    USE scm_par
    IMPLICIT NONE
    INTEGER, INTENT(IN   )         :: ntra                 !! number of tracers
    INTEGER, INTENT(IN   )         :: npts                 !! number of points used for extrapolation
    REAL(8), INTENT(IN   )         :: tke  (1:npts)        !! turbulent kinetic energy [m2/s2]
    REAL(8), INTENT(IN   )         :: u_m  (1:npts)        !! mean zonal velocity [m/s]
    REAL(8), INTENT(IN   )         :: v_m  (1:npts)        !! mean meridional velocity [m/s]
    REAL(8), INTENT(IN   )         :: t_m  (1:npts,1:ntra) !! mean tracers
    REAL(8), INTENT(IN   )         :: Hz   (1:npts)        !! layer thickness [m]
    !REAL(8), INTENT(  OUT)         :: wp0                  !! vertical plume velocity at the surface [m/s]
    REAL(8), INTENT(  OUT)         :: up0                  !! zonal plume velocity at the surface [m/s]
    REAL(8), INTENT(  OUT)         :: vp0                  !! meridional plume velocity at the surface [m/s]
    REAL(8), INTENT(  OUT)         :: tp0(1:ntra+1)        !! tracer plume properties
    ! local variables
    INTEGER                        :: itrc
    REAL(8)                        :: cff
    !
    !wp0   = - wpmin  !! \( w^{\rm p}_{\rm sfc} = - w^{\rm p}_{\min} \)<br />
    !
    cff = 1./(Hz(npts-1)+Hz(npts))
    DO itrc = 1,ntra   ! temperature & salinity
      tp0 (itrc) = cff*((2.*Hz(npts)+Hz(npts-1))*t_m(npts  ,itrc)  &
                                           -Hz(npts)*t_m(npts-1,itrc)  ) !! \( \phi^{\rm p}_{\rm sfc} = \frac{ (2 h_N+h_{N-1}) \overline{\phi}_N - h_N \overline{\phi}_{N-1} }{h_N + h_{N-1}} \) <br />
      !tp0 (itrc) = (1.+0.75)*t_m(npts  ,itrc) - 0.75*t_m(npts-1,itrc)
    ENDDO
    ! TKE
    tp0 ( ntra+1 ) = tke(npts) !! \( k^{\rm p}_{\rm sfc} = k_{N+1/2} \) <br />
    !
    up0     = cff*((2.*Hz(npts)+Hz(npts-1))*u_m(npts)  &
                                         -Hz(npts)*u_m(npts-1)  ) !! \( u^{\rm p}_{\rm sfc} = \frac{ (2 h_N+h_{N-1}) \overline{u}_N - h_N \overline{u}_{N-1} }{h_N + h_{N-1}} \) <br />
    vp0     = cff*((2.*Hz(npts)+Hz(npts-1))*v_m(npts  )  &
                                         -Hz(npts)*v_m(npts-1)  ) !! \( v^{\rm p}_{\rm sfc} = \frac{ (2 h_N+h_{N-1}) \overline{v}_N - h_N \overline{v}_{N-1} }{h_N + h_{N-1}} \) <br />
  !---------------------------------------------------------------------------------------------------
  END SUBROUTINE compute_MF_bdy
  !===================================================================================================


  !===================================================================================================
  SUBROUTINE compute_MF_forcing(u_np1,v_np1,t_np1,tke,N,ntra,u_m,v_m,t_m,dtke_m)
  !---------------------------------------------------------------------------------------------------
    !!==========================================================================<br />
    !!                  ***  ROUTINE compute_MF_forcing  ***                    <br />
    !! ** Purposes : compute mean quantities used to force the MF equations     <br />
    !!==========================================================================<br />
    IMPLICIT NONE
    INTEGER, INTENT(IN   )         :: N                    !! number of vertical levels
    INTEGER, INTENT(IN   )         :: ntra                 !! number of tracers
    REAL(8), INTENT(IN   )         :: u_np1(1:N)           !! zonal velocity [m/s]
    REAL(8), INTENT(IN   )         :: v_np1(1:N)           !! meridional velocity [m/s]
    REAL(8), INTENT(IN   )         :: t_np1(1:N,ntra)      !! tracers
    REAL(8), INTENT(IN   )         :: tke  (0:N)           !! turbulent kinetic energy [m/s]
    REAL(8), INTENT(  OUT)         :: u_m(1:N)             !! mean zonal velocity for mass-flux equations [m/s]
    REAL(8), INTENT(  OUT)         :: v_m(1:N)             !! mean meridional velocity for mass-flux equations [m/s]
    REAL(8), INTENT(  OUT)         :: t_m(1:N,ntra+1)      !! mean tracer for mass-flux equations [m/s]
    REAL(8), INTENT(  OUT)         :: dtke_m(1:N)      !! mean tracer for mass-flux equations [m/s]
    ! local variables
    INTEGER                        :: k,itrc
    !
    DO itrc = 1,ntra  ! temperature & salinity
      DO k = 1,N
        t_m(k,itrc) = t_np1(k,itrc) !! \( \overline{\phi}_{k} = \phi^{\rm scm}_k \)  <br />
      ENDDO
    ENDDO
    ! TKE
    DO k = 1,N
      t_m(k,ntra+1) = 0.5*(tke(k)+tke(k-1))  !! \( k_{k+1/2} = \frac{1}{2} (k^{\rm scm}_{k+1/2}+k^{\rm scm}_{k-1/2}) \)  <br />
      dtke_m( k )   = tke(k)-tke(k-1)
    ENDDO
    !
    DO k = 1,N
      u_m(k) = u_np1(k) !! \( \overline{u}_{k} = u^{\rm scm}_k \)  <br />
      v_m(k) = v_np1(k) !! \( \overline{v}_{k} = v^{\rm scm}_k \)  <br />
    ENDDO
  !---------------------------------------------------------------------------------------------------
  END SUBROUTINE compute_MF_forcing
  !===================================================================================================

  !===================================================================================================
  SUBROUTINE compute_tripleCorr(tke,tke_p,Fmass,u_p,v_p,w_p,u_np1,v_np1,Hz,zr,opt,wtke,N,trplCorr)
  !---------------------------------------------------------------------------------------------------
    !!==========================================================================<br />
    !!                  ***  ROUTINE compute_tripleCorr  ***                    <br />
    !! ** Purposes : compute contribution of MF terms to turbulent flux w'e     <br />
    !!   \[ \overline{w'e}_{\rm MF} = a^{\rm p} w^{\rm p} \underbrace{\left[ \frac{1}{1-a^{\rm p}} \right]}_{\rm small\_ap = False} \left\{ \left(  k^{\rm p} - k^n \right)
    !!    + \frac{1}{2} \| \mathbf{v}^p - \mathbf{v}^{n+1} \|^2 \right\}  \]     <br />
    !!==========================================================================<br />
    USE scm_par
    IMPLICIT NONE
    INTEGER, INTENT(IN   )         :: N                    !! number of vertical levels
    REAL(8), INTENT(IN   )         ::  tke  (0:N)          !! mean turbulent kinetic energy [m/s]
    REAL(8), INTENT(IN   )         ::  tke_p(0:N)          !! plume turbulent kinetic energy [m/s]
    REAL(8), INTENT(IN   )         ::  w_p(0:N)            !! plume vertical velocity [m/s]
    REAL(8), INTENT(IN   )         ::  Fmass(0:N)          !! mass flux [m/s]
    REAL(8), INTENT(IN   )         ::  u_p  (0:N)          !! plume zonal velocity [m/s]
    REAL(8), INTENT(IN   )         ::  v_p  (0:N)          !! plume meridional velocity [m/s]
    REAL(8), INTENT(IN   )         ::  u_np1(1:N)          !! mean zonal velocity [m/s]
    REAL(8), INTENT(IN   )         ::  v_np1(1:N)          !! mean meridional velocity [m/s]
    REAL(8), INTENT(IN   )         ::  Hz  (1:N)           !! layer thickness [m]
    REAL(8), INTENT(IN   )         ::  zr  (1:N)           !! depth of cell centers [m]
    INTEGER, INTENT(IN   )         ::  opt                 !! option for the computation of wtke (=0 PL24, =1 HB09, = 2 Wi11)
    REAL(8), INTENT(INOUT)         ::  wtke(1:N)           !! turbulent w'e flux [m3/s3]
    REAL(8), INTENT(  OUT)         ::  trplCorr(0:N)       !! contribution of mass flux to divergence of w'e flux [m2/s3]
    ! local variables
    INTEGER                        :: k
    REAL(8)                        :: cff, cff1, cff2, cff3
    REAL(8)                        :: FC(1:N)
    REAL(8)                        :: Fm,wpr,upr,vpr
    !===============================================================
    ! Compute fluxes at rho-points and flux divergence at w-points
    FC      (1:N) = 0.
    FC(N) = 0.5*Fmass(N)*(w_p(N))**2 ! attention Fmass=-apwp
    wtke(N) = wtke(N) - FC(N)
    trplCorr(0:N) = 0.
    !
    DO k = 1,N-1
      Fm    = Fmass(k  )
      wpr   = w_p(k  )
      upr   = u_p(k  )
      vpr   = v_p(k  )
      !
      SELECT CASE( opt )
      CASE( 0  ) ! all terms
        FC(k) = Fm*( tke_p(k)-tke(k-1) ) !! \( F_{k} = (a^{\rm p} w^{\rm p})_{k+1/2} ( k^{\rm p}_{k+1/2} - k_{k-1/2}^n ) \)<br />
        FC(k) = FC(k) + 0.5*Fm*( (upr-u_np1(k))**2 +  (vpr-v_np1(k))**2   +   wpr*wpr ) !! \( F_{k} = F_{k} + \frac{(a^{\rm p} w^{\rm p})_{k+1/2}}{2}\left(  \mathbf{v}^{\rm p}_{k+1/2} - \mathbf{v}^{n+1}_{k} \right)^2  \) <br />
      CASE( 1  ) ! only MF on tke
        FC(k) = Fm*( tke_p(k)-tke(k-1) ) !! \( F_{k} = (a^{\rm p} w^{\rm p})_{k+1/2} ( k^{\rm p}_{k+1/2} - k_{k-1/2}^n ) \)<br />
      CASE( 2  ) ! only ap wp^3 term
        FC(k) = 0.5*Fm*( wpr*wpr )
      END SELECT
      !
      wtke(k) = wtke(k) - FC(k)
    ENDDO
    !
    DO k = 1,N-1
      trplCorr(k) = (FC(k+1)-FC(k))/(zr(k+1)-zr(k)) !! \( {\rm TOM}_{k+1/2} = \frac{ F_{k+1}-F_{k} }{ \Delta z_{k+1/2} }   \) <br />
    ENDDO
    !!@note Computation of environmental TKE for diagnostics is missing @endnote
    !DO k = 1,N
    !  upr        = u_p(k) - u_np1(k)
    !  vpr        = v_p(k) - v_np1(k)
    !  wpr        = w_p(k)
    !  tke_env(k) = tke(k) - a_p(k)*( tke_p(k) + wpr*wpr + upr*upr + vpr*vpr )
    !ENDDO
    !tke_env(0:N) = 0.
  !---------------------------------------------------------------------------------------------------
  END SUBROUTINE compute_tripleCorr
  !===================================================================================================


  !===================================================================================================
  SUBROUTINE eos_val_lin(temp,salt,eos_params,neos,rho)
  !---------------------------------------------------------------------------------------------------
    !!==========================================================================<br />
    !!                  ***  ROUTINE eos_val_lin  ***                           <br />
    !! ** Purposes : compute density anomaly from linear equation of state      <br />
    !!==========================================================================<br />
    IMPLICIT NONE
    INTEGER, INTENT(IN   )              :: neos
    REAL(8), INTENT(IN   )              :: temp    !! temperature [C]
    REAL(8), INTENT(IN   )              :: salt    !! salinity [psu]
    REAL(8), INTENT(IN   )              :: eos_params(neos)
    REAL(8), INTENT(INOUT)              :: rho     !! density anomaly
    !
    REAL(8)                             :: alpha,beta,rhoRef,T0,S0    !! haline expension coefficient [psu-1]
    !---------------------------------------------------------------------------
    rhoRef = eos_params(1); alpha  = eos_params(2); beta = eos_params(3)
    T0     = eos_params(4); S0     = eos_params(5)
    !---------------------------------------------------------------------------
    ! Compute density anomaly via linear Equation Of State (EOS)
    !-------
    rho = rhoRef*( 1. - alpha*( temp - T0 ) + beta*( salt - S0 ) )  !! \(   \rho_{k} = \rho_0 \left( 1 - \alpha (\theta - 2) + \beta (S - 35)   \right)  \)  <br />
    !!@note Reference values \( \theta_0 = 2^{o}C \) and \( S_0 = 35\;{\rm psu}\) are hardcoded !!  @endnote <br />
    return
  !---------------------------------------------------------------------------------------------------
  END SUBROUTINE eos_val_lin
  !===================================================================================================


  !===================================================================================================
  SUBROUTINE eos_val (temp,salt,zr,rhoRef,rho)
  !---------------------------------------------------------------------------------------------------
    !!==========================================================================<br />
    !!                  ***  ROUTINE rho_eos  ***                           <br />
    !! ** Purposes : Compute density anomaly and Brunt Vaisala frequency via nonlinear
    !!                                                  Equation Of State (EOS) <br />
    !!==========================================================================<br />
    USE scm_par
    IMPLICIT NONE
    REAL(8), INTENT(IN   )         :: temp               !! temperature [C]
    REAL(8), INTENT(IN   )         :: salt               !! salinity [psu]
    REAL(8), INTENT(IN   )         :: zr                 !! depth [m]
    REAL(8), INTENT(IN   )         :: rhoRef             !! reference density
    REAL(8), INTENT(INOUT)         :: rho                !! density anomaly [kg/m3]
    ! local variables
    real(8)                        :: cff,cff1,dr00,dpth,rho1
    real(8)                        :: Tt,Ts,sqrtTs,K0,K1,K2
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
    Tt       = temp
    Ts       = salt
    sqrtTs   = sqrt(Ts)
    !
    rho1  = (dr00 +Tt*( r01+Tt*( r02+Tt*( r03+Tt*(r04+Tt*r05 ))))  &
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
    dpth = -zr
    cff  = K00-0.1*dpth
    cff1 = K0+dpth*(K1+K2*dpth)
    rho  = ( rho1*cff*(K00+cff1)-0.1*dpth*rhoRef*cff1 )/(cff*(cff+cff1))
  !---------------------------------------------------------------------------------------------------
  END SUBROUTINE eos_val
  !===================================================================================================


  !===================================================================================================
  SUBROUTINE get_w_p_R10(wpm,wpp,aa,bb,bpr,beta1,Hz,Bp,wp_min,h,found)
  !---------------------------------------------------------------------------------------------------
    !!==========================================================================<br />
    !!                  ***  ROUTINE get_w_p_R10  ***                           <br />
    !! ** Purposes : compute plume vertical velocity for Rio et al. (2010)
    !!                        entrainment/detrainement closure                  <br />
    !!==========================================================================<br />
    IMPLICIT NONE
    REAL(8), INTENT(INOUT)      :: wpm          !! vertical velocity at the bottom of the grid cell [m/s]
    REAL(8), INTENT(INOUT)      :: h            !! distance from the top of the grid cell where w_p = w_p_min [m]
    LOGICAL, INTENT(INOUT)      :: found        !! (=T) the bottom of the plume is reached
    REAL(8), INTENT(IN   )      :: wpp          !! vertical velocity at the top of the grid cell [m/s]
    REAL(8), INTENT(IN   )      :: Bp           !! plume buoyancy forcing term [m/s2]
    REAL(8), INTENT(IN   )      :: wp_min
    REAL(8), INTENT(IN   )      :: Hz           !! thickness of the grid cell [m]
    REAL(8), INTENT(IN   )      :: aa           !! parameter of the MF scheme
    REAL(8), INTENT(IN   )      :: bb           !! parameter of the MF scheme
    REAL(8), INTENT(IN   )      :: bpr          !! parameter of the MF scheme
    REAL(8), INTENT(IN   )      :: beta1        !! parameter of the MF scheme
    !REAL(8), INTENT(IN   )      :: wpmin        !! minimum vertical velociy

    ! local variables
    REAL(8)                     :: cff1,cff,rhsw,wpm2
    !
    rhsw = bpr*(wpp*wpp+wpp*wpp)+2.*aa*Bp       !! \(   {\rm rhs}_{k}^{\star} = 2 b'( w^{\rm p}_{k+1/2} )^2 + 2 a B_{k}^{\rm p}  \)  <br />
    cff  = 1.; IF(rhsw < 0.) cff = 1. + bb*beta1 !! if \( {\rm rhs}_{k}^{\star} < 0 \Rightarrow \alpha_w = 1 + b \beta_1 \), \( \alpha_w = 1 \) otherwise <br />
    cff1 = 1./(cff+Hz*bpr)
    wpm2 = cff1*( (cff-Hz*bpr)*wpp*wpp - aa*Hz*2.*Bp )  !! \(  ( w^{\rm p}_{k-1/2} )^2 = \frac{ (\alpha_w - h_k b') ( w^{\rm p}_{k+1/2} )^2 - 2 a h_k B_{k}^{\rm p} }{\alpha_w + h_k b'}\) <br />
! in case the sign of the rhs has changed
    rhsw = bpr*(wpp*wpp+wpm2)+2.*aa*Bp !! \(   {\rm rhs}_{k} = b'\left( ( w^{\rm p}_{k-1/2} )^2 + ( w^{\rm p}_{k+1/2} )^2 \right) + 2 a B_{k}^{\rm p}  \)  <br />
    cff  = 1.; IF(rhsw < 0.) cff = 1. + bb*beta1 !! if \( {\rm rhs}_{k} < 0 \Rightarrow \alpha_w = 1 + b \beta_1 \), \( \alpha_w = 1 \) otherwise <br />
    cff1 = 1./(cff+Hz*bpr)
    wpm2 = cff1*( (cff-Hz*bpr)*wpp*wpp - aa*Hz*2.*Bp ) !! \(  ( w^{\rm p}_{k-1/2} )^2 = \frac{ (\alpha_w - h_k b') ( w^{\rm p}_{k+1/2} )^2 - 2 a h_k B_{k}^{\rm p} }{\alpha_w + h_k b'}\) <br />
    ! finalize computation
    wpm  = -SQRT( MAX(wpm2,wp_min*wp_min) )
    !
    IF(cff==1.) THEN                    ! the bottom of the plume is reached only if we are in the detrainment zone
      h = (wpp*wpp-wp_min*wp_min) /   &
                    (2.*aa*Bp+bpr*(wpp*wpp+wp_min*wp_min))
      IF(h>0. .AND. h<Hz) found = .true.
    ENDIF
    return
  !---------------------------------------------------------------------------------------------------
  END SUBROUTINE get_w_p_R10
  !===================================================================================================



  !===================================================================================================
  SUBROUTINE get_a_p_R10(apm,app,wpm,wpp,beta1,beta2,hk,delta0,wp_min)
  !---------------------------------------------------------------------------------------------------
    !!==========================================================================<br />
    !!                  ***  ROUTINE get_a_p_R10  ***                           <br />
    !! ** Purposes : compute plume fractional area for Rio et al. (2010)
    !!                        entrainment/detrainement closure                  <br />
    !!==========================================================================<br />
    IMPLICIT NONE
    REAL(8), INTENT(INOUT)      :: apm    !! fractional area at the bottom of the grid cell
    REAL(8), INTENT(IN   )      :: app    !! fractional area at the top   of the grid cell
    REAL(8), INTENT(IN   )      :: wpm    !! vertical velocity at the bottom of the grid cell [m/s]
    REAL(8), INTENT(IN   )      :: wpp    !! vertical velocity at the top of the grid cell [m/s]
    REAL(8), INTENT(IN   )      :: beta1  !! parameter of the MF scheme for the entrainment zone
    REAL(8), INTENT(IN   )      :: hk     !! thickness of the grid cell [m]
    REAL(8), INTENT(IN   )      :: delta0 !! background detrainement in the entrainment zone [m-1]
    REAL(8), INTENT(IN   )      :: wp_min
    REAL(8), INTENT(IN   )      :: beta2  !! parameter of the MF scheme for the detrainment zone
    !local variables
    REAL(8)                 :: cff1,cff,EmD
    !
    EmD = Ent_R10(beta1,wpp,wpm)  + &
          Det_R10(beta1,beta2,wpp,wpm,delta0,hk,wp_min) !! \(   {\rm EmD}_k = \frac{h_k}{a^{\rm p}} \left( E_k - D_k \right) = {\rm Ent\_R10} + {\rm Det\_R10}  \)  <br />
    cff  = 1. / (2.*wpm+EmD)
    cff1 = app *(2.*wpp-EmD)
    apm  = MAX(cff*cff1,0.)   !! \(  a^{\rm p}_{k-1/2} = a^{\rm p}_{k+1/2} \left(  \frac{2 w^{\rm p}_{k+1/2} - {\rm EmD}_k}{2 w^{\rm p}_{k-1/2} + {\rm EmD}_k}   \right)   \) <br />
    return
  !---------------------------------------------------------------------------------------------------
  END SUBROUTINE get_a_p_R10
  !===================================================================================================

  REAL(8) FUNCTION Ent_R10(beta1,wpp,wpm)   !! Entrainment \( {\rm Ent\_R10} = \frac{E_k h_k}{a^p} = \max\left( \beta_1(w_{k+1/2}^p-w_{k-1/2}^p), 0 \right)  \)
    IMPLICIT NONE
    REAL(8), INTENT(IN)    :: beta1  !! parameter of the MF scheme for the entrainment zone
    REAL(8), INTENT(IN)    :: wpp    !! \( w^{\rm p}_{k+1/2} \)
    REAL(8), INTENT(IN)    :: wpm    !! \( w^{\rm p}_{k-1/2} \)
    Ent_R10 = MAX( beta1*(wpp-wpm), 0. )
  END FUNCTION Ent_R10

  REAL(8) FUNCTION Det_R10(beta1,beta2,wpp,wpm,delta0,hk,wp_min)
    !! Detrainment \( {\rm Det\_R10} = -\frac{D_k h_k}{a^p} = \min\left( \beta_2(w_{k+1/2}^p-w_{k-1/2}^p), 0 \right) + \min\left(-2 w^{\rm p}_{\min}, h_k \delta_0 \frac{w_{k+1/2}^p+w_{k-1/2}^p}{2} \right) \)  <br />
    !! The minimum detrainment \( -2 w^{\rm p}_{\min} \) ensures that \( a^{\rm p}_{k-1/2} = 0 \) as soon as \( w^{\rm p}_{k+1/2} = w^{\rm p}_{k-1/2} = -w^{\rm p}_{\min} \)
    USE scm_par, ONLY: wpmin
    IMPLICIT NONE
    REAL(8), INTENT(IN)    :: beta1   !! parameter of the MF scheme for the entrainment zone
    REAL(8), INTENT(IN)    :: wpp     !! \( w^{\rm p}_{k+1/2} \)
    REAL(8), INTENT(IN)    :: wpm     !! \( w^{\rm p}_{k-1/2} \)
    REAL(8), INTENT(IN)    :: delta0  !! background detrainment [m-1]
    REAL(8), INTENT(IN)    :: hk      !! Thickness \( h_k /) of layer k
    REAL(8), INTENT(IN)    :: wp_min
    REAL(8), INTENT(IN)    :: beta2   !! parameter of the MF scheme for the detrainment zone
    REAL(8)                :: D0,D1
    !REAL(8), INTENT(IN)    :: wpmin   !! minimum vertical velociy
    D0 = 0.5*hk*delta0*(wpp+wpm)
    Det_R10 = MIN( beta2*(wpp-wpm), 0. ) + MIN(D0, -2.*wp_min)
  END FUNCTION Det_R10

  !===================================================================================================
  SUBROUTINE get_t_p_R10(tpm,tpp,te,apm,app,wpm,wpp,beta1,beta2,hk,delta0,wp_min)
  !---------------------------------------------------------------------------------------------------
    !!==========================================================================<br />
    !!                  ***  ROUTINE get_t_p_R10  ***                           <br />
    !! ** Purposes : compute plume tracer properties for Rio et al. (2010)
    !!                        entrainment/detrainement closure                  <br />
    !!==========================================================================<br />
    IMPLICIT NONE
    REAL(8), INTENT(INOUT)      :: tpm    !! tracer value at the bottom of the grid cell
    REAL(8), INTENT(IN   )      :: tpp    !! tracer value at the top   of the grid cell
    REAL(8), INTENT(IN   )      :: te     !! environmental value for tracer in the grid cell
    REAL(8), INTENT(IN   )      :: apm    !! fractional area at the bottom of the grid cell
    REAL(8), INTENT(IN   )      :: app    !! fractional area at the top   of the grid cell
    REAL(8), INTENT(IN   )      :: wpm    !! vertical velocity at the bottom of the grid cell [m/s]
    REAL(8), INTENT(IN   )      :: wpp    !! vertical velocity at the top of the grid cell [m/s]
    REAL(8), INTENT(IN   )      :: beta1  !! parameter of the MF scheme
    REAL(8), INTENT(IN   )      :: hk     !! thickness of the grid cell [m]
    REAL(8), INTENT(IN   )      :: delta0 !! background detrainement in the entrainment zone [m-1]
    REAL(8), INTENT(IN   )      :: beta2  !! increase the detrainement coefficient beta1 [m-1]
    REAL(8), INTENT(IN   )      :: wp_min
    ! local variables
    REAL(8)                 :: cffm,cffp,dwpm,dwpp,ap,cff
    !
    cffp = app*wpp*tpp         !! \(  \Phi^{\rm p}_{k+1/2} = a^{\rm p}_{k+1/2} w^{\rm p}_{k+1/2} \phi^{\rm p}_{k+1/2}  \)  <br />
    ap   = 0.5*(app+apm)       !! \(  a^{\rm p}_{k} = \frac{a^{\rm p}_{k+1/2}+a^{\rm p}_{k-1/2}}{2}  \)  <br />
    dwpp = Ent_R10(beta1,wpp,wpm)
    dwpm = Det_R10(beta1,beta2,wpp,wpm,delta0,hk,wp_min)
    cffm = cffp - ap*(dwpp*te+dwpm*tpp) !! \( \Phi^{\rm p}_{k-1/2} = \Phi^{\rm p}_{k+1/2} - a^{\rm p}_{k} \left( \underbrace{\frac{E_k h_k}{a^p}}_{\rm Ent\_R10} \phi_{k}^{\rm e} \underbrace{- \frac{D_k h_k}{a^p}}_{\rm Det\_R10} \phi_{k+1/2}^{\rm p}     \right)    \) <br />
!    IF(apm>0.) THEN
    IF( apm>0. .AND. wpm < -wp_min ) THEN
      tpm  = cffm/(apm*wpm)  !! \( \phi^{\rm p}_{k-1/2} = \frac{\Phi^{\rm p}_{k-1/2}}{( a^{\rm p} w^{\rm p} )_{k-1/2}} \)  <br />
    ELSE
      tpm  = tpp
    ENDIF
    !
    return
  !---------------------------------------------------------------------------------------------------
  END SUBROUTINE get_t_p_R10
  !===================================================================================================

  !===================================================================================================
  SUBROUTINE get_vort_p_R10(Vortpm,Vortpp,fcor,apm,app,wpm,wpp,beta1,beta2,hk,delta0,wp_min)
  !---------------------------------------------------------------------------------------------------
    !!==========================================================================<br />
    !!                  ***  ROUTINE get_t_p_R10  ***                           <br />
    !! ** Purposes : compute plume tracer properties for Rio et al. (2010)
    !!                        entrainment/detrainement closure                  <br />
    !!==========================================================================<br />
    IMPLICIT NONE
    REAL(8), INTENT(INOUT)      :: Vortpm    !! tracer value at the bottom of the grid cell
    REAL(8), INTENT(IN   )      :: Vortpp    !! tracer value at the top   of the grid cell
    REAL(8), INTENT(IN   )      :: fcor      !! environmental value for tracer in the grid cell
    REAL(8), INTENT(IN   )      :: apm       !! fractional area at the bottom of the grid cell
    REAL(8), INTENT(IN   )      :: app       !! fractional area at the top   of the grid cell
    REAL(8), INTENT(IN   )      :: wpm       !! vertical velocity at the bottom of the grid cell [m/s]
    REAL(8), INTENT(IN   )      :: wpp       !! vertical velocity at the top of the grid cell [m/s]
    REAL(8), INTENT(IN   )      :: beta1     !! parameter of the MF scheme
    REAL(8), INTENT(IN   )      :: hk        !! thickness of the grid cell [m]
    REAL(8), INTENT(IN   )      :: delta0    !! background detrainement in the entrainment zone [m-1]
    REAL(8), INTENT(IN   )      :: beta2     !! increase the detrainement coefficient beta1 [m-1]
    REAL(8), INTENT(IN   )      :: wp_min
    ! local variables
    REAL(8)                 :: cffm,num,hDet,dwpp,ap,Denom
    !
    ap     = 0.5*(app+apm)       !! \(  a^{\rm p}_{k} = \frac{a^{\rm p}_{k+1/2}+a^{\rm p}_{k-1/2}}{2}  \)  <br />
    hDet   = ap*Det_R10(beta1,beta2,wpp,wpm,delta0,hk,wp_min) != - Dj x hj
    num    = apm*wpm+app*wpp+hDet
    IF(ABS(num) > 1.e-20) THEN
      Denom  = Vortpp*( apm*wpm+app*wpp-hDet ) - 2.*fcor*(app*wpp-apm*wpm)
      Vortpm = Denom / num
    ELSE
      Vortpm = 0.
    ENDIF
    !print*,'Vortpm = ',Vortpm
    !
    return
  !---------------------------------------------------------------------------------------------------
  END SUBROUTINE get_vort_p_R10
  !===================================================================================================


  !===================================================================================================
  SUBROUTINE get_tke_p_R10( tkep_m,tkep_p, tke_m, tke_p, wpm,wpp,normvel,epsilon,hk,beta1 )
  !---------------------------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(8), INTENT(INOUT)      :: tkep_m    !! plume TKE value at the bottom of the grid cell
    REAL(8), INTENT(IN   )      :: tkep_p    !! plume TKE value at the top   of the grid cell
    REAL(8), INTENT(IN   )      :: tke_m     !! external forcing for the plume TKE
    REAL(8), INTENT(IN   )      :: tke_p     !! external forcing for the plume TKE
    REAL(8), INTENT(IN   )      :: wpm       !! vertical velocity at the bottom of the grid cell [m/s]
    REAL(8), INTENT(IN   )      :: wpp       !! vertical velocity at the top of the grid cell [m/s]
    REAL(8), INTENT(IN   )      :: normvel   !! \( \| \mathbf{u}^{\rm p} - \mathbf{u} \|^2  \) [m2/s-2]
    REAL(8), INTENT(IN   )      :: hk        !! thickness of the grid cell [m]
    REAL(8), INTENT(IN   )      :: epsilon   !! plume TKE dissipation term [m2/s-3]
    REAL(8), INTENT(IN   )      :: beta1     !! parameter of the MF scheme
    REAL(8)                     :: ent,cff
    !---------------------------------------------------------------------------------------------------
    ent    = beta1*MAX(wpp-wpm,0.)
    cff    = 1./( wpp + wpm - ent )
    tkep_m = cff*( 2.*hk*epsilon+tkep_p*(wpp + wpm)-ent*( tke_m+tke_p-tkep_p+normvel )  )
    return
  END SUBROUTINE get_tke_p_R10
  !===================================================================================================

  !===================================================================================================
  SUBROUTINE get_dtke_p_R10(tkep_m,tkep_p,dtke_m,normvel,wpm,wpp,hk,epsilon,beta1)
  !---------------------------------------------------------------------------------------------------
    !!==========================================================================<br />
    !!                  ***  ROUTINE get_dtke_p_R10  ***                           <br />
    !! ** Purposes : compute plume TKE for Rio et al. (2010)
    !!                        entrainment/detrainement closure                  <br />
    !!==========================================================================<br />
    IMPLICIT NONE
    REAL(8), INTENT(INOUT)      :: tkep_m    !! plume TKE value at the bottom of the grid cell
    REAL(8), INTENT(IN   )      :: tkep_p    !! plume TKE value at the top   of the grid cell
    REAL(8), INTENT(IN   )      :: dtke_m    !! external forcing for the plume TKE
    REAL(8), INTENT(IN   )      :: wpm       !! vertical velocity at the bottom of the grid cell [m/s]
    REAL(8), INTENT(IN   )      :: wpp       !! vertical velocity at the top of the grid cell [m/s]
    REAL(8), INTENT(IN   )      :: normvel   !! \( \| \mathbf{u}^{\rm p} - \mathbf{u} \|^2  \) [m2/s-2]
    REAL(8), INTENT(IN   )      :: hk        !! thickness of the grid cell [m]
    REAL(8), INTENT(IN   )      :: epsilon   !! plume TKE dissipation term [m2/s-3]
    REAL(8), INTENT(IN   )      :: beta1     !! parameter of the MF scheme
    ! local variables
    REAL(8)                 :: cffm,cffp,dwp,cff,chi
    !
    dwp     = MAX(wpp-wpm,0.)    !! \(  (\delta w^{\rm p})_{k}^{+} = \max( w^{\rm p}_{k+1/2}-w^{\rm p}_{k-1/2}  ,0) \)  <br />
    cff     = 1./(wpp+wpm)
    chi     = cff*beta1*dwp      !! \(  \chi_k = \frac{\beta_1}{w^{\rm p}_{k+1/2}+w^{\rm p}_{k-1/2}} (\delta w^{\rm p})_{k}^{+} \equiv \frac{E_k h_k}{2 a^{\rm p} w^{\rm p}}     \) <br />
    tkep_m  = ( (1.+chi)*tkep_p - chi*normvel + dtke_m + 2.*hk*cff*epsilon)/(1.-chi) !! \( (e^{\rm p}-e)_{k-1/2} = \frac{(1+\chi_k)(e^{\rm p}-e)_{k+1/2} - \chi_k \| \mathbf{u}^{\rm p} - \mathbf{u} \|_k^2 + (e_{k+1/2}-e_{k-1/2}) + \frac{2 h_k}{w^{\rm p}_{k+1/2}+w^{\rm p}_{k-1/2}} \epsilon_{k+1/2}^{\rm p}   }{1-\chi_k}    \)
    !
    return
  !---------------------------------------------------------------------------------------------------
  END SUBROUTINE get_dtke_p_R10
  !===================================================================================================

  !===================================================================================================
  SUBROUTINE get_u_p_R10(tpm,tpp,te,apm,app,wpm,wpp,beta1,beta2,hk,delta0,frc,cor,wp_min)
  !---------------------------------------------------------------------------------------------------
    !!==========================================================================<br />
    !!                  ***  ROUTINE get_t_p_R10  ***                           <br />
    !! ** Purposes : compute plume tracer properties for Rio et al. (2010)
    !!                        entrainment/detrainement closure                  <br />
    !!==========================================================================<br />
    IMPLICIT NONE
    REAL(8), INTENT(INOUT)      :: tpm    !! tracer value at the bottom of the grid cell
    REAL(8), INTENT(IN   )      :: tpp    !! tracer value at the top   of the grid cell
    REAL(8), INTENT(IN   )      :: te     !! environmental value for tracer in the grid cell
    REAL(8), INTENT(IN   )      :: apm    !! fractional area at the bottom of the grid cell
    REAL(8), INTENT(IN   )      :: app    !! fractional area at the top   of the grid cell
    REAL(8), INTENT(IN   )      :: wpm    !! vertical velocity at the bottom of the grid cell [m/s]
    REAL(8), INTENT(IN   )      :: wpp    !! vertical velocity at the top of the grid cell [m/s]
    REAL(8), INTENT(IN   )      :: beta1  !! parameter of the MF scheme
    REAL(8), INTENT(IN   )      :: hk     !! thickness of the grid cell [m]
    REAL(8), INTENT(IN   )      :: delta0 !! background detrainement in the entrainment zone [m-1]
    REAL(8), INTENT(IN   )      :: beta2  !! increase the detrainement coefficient beta1 [m-1]
    REAL(8), INTENT(IN   )      :: cor    !! traditional Coriolis term [m2/s2]
    REAL(8), INTENT(IN   )      :: wp_min
    REAL(8), INTENT(IN   )      :: frc    !! forcing term for up/vp equation [m/s]
    ! local variables
    REAL(8)                 :: cffm,cffp,dwpm,dwpp,ap,apwp,cff
    !
    cffp = app*wpp*tpp         !! \(  \Phi^{\rm p}_{k+1/2} = a^{\rm p}_{k+1/2} w^{\rm p}_{k+1/2} \phi^{\rm p}_{k+1/2}  \)  <br />
    ap   = 0.5*(app+apm)       !! \(  a^{\rm p}_{k} = \frac{a^{\rm p}_{k+1/2}+a^{\rm p}_{k-1/2}}{2}  \)  <br />
    apwp = 0.5*(app*wpp+apm*wpm)
    !
    ap   = app
    apwp = app*wpp
    !
    dwpp = Ent_R10(beta1,wpp,wpm)
    dwpm = Det_R10(beta1,beta2,wpp,wpm,delta0,hk,wp_min)
    cffm = cffp - ap*(dwpp*te+dwpm*tpp+cor) - apwp*frc !! \( \Phi^{\rm p}_{k-1/2} = \Phi^{\rm p}_{k+1/2} - a^{\rm p}_{k} \left( \underbrace{\frac{E_k h_k}{a^p}}_{\rm Ent\_R10} \phi_{k}^{\rm e} \underbrace{- \frac{D_k h_k}{a^p}}_{\rm Det\_R10} \phi_{k+1/2}^{\rm p}     \right)    \) <br />
    IF( apm>0. .AND. wpm < -wp_min ) THEN
    !IF( apm>0. ) THEN
      tpm  = cffm/(apm*wpm)  !! \( \phi^{\rm p}_{k-1/2} = \frac{\Phi^{\rm p}_{k-1/2}}{( a^{\rm p} w^{\rm p} )_{k-1/2}} \)  <br />
    ELSE
      tpm  = tpp
    ENDIF
    !
    !IF(ABS(tpm)>100.) then
    !  print*,'Plume vel = ',tpm,apm,wpm,app,wpp,frc,cor
    !ENDIF
    !
    return
  !---------------------------------------------------------------------------------------------------
  END SUBROUTINE get_u_p_R10
  !===================================================================================================





  !===================================================================================================
  SUBROUTINE mass_flux_R10(u_m,v_m,t_m,tke_m,z_w,Hz,tp0,up0,vp0,wp0,mf_params,eos_params,  &
                           tkep_min,mxlp_min,small_ap,lin_eos,opt,fcor,zinv, N,ntra,nparams,neos,   &
                            a_p,u_p,v_p,w_p,t_p,B_p,ent,det,eps,vort_p)
  !---------------------------------------------------------------------------------------------------
    !!==========================================================================<br />
    !!                  ***  ROUTINE mass_flux_R10  ***                         <br />
    !! ** Purposes : solve mass-flux equations for Rio et al. (2010)
    !!                        entrainment/detrainement closure                  <br />
    !! \begin{align*} E &= a^{\rm p} \beta_1 \max(0, \partial_z w^{\rm p}) \\
    !!                D &= - a^{\rm p} \beta_2 \min\left( 0 , \partial_z w^{\rm p} \right) -  a^{\rm p} w^{\rm p} \delta_0  \end{align*} <br />
    !!@note we must have \( \beta_1 < 1 \) and \( 1 < \beta_2 < 2 \) otherwise unphysical values of \( a^{\rm p} \) are obtained  @endnote <br />
    !! \begin{align*}
    !! \partial_z(a^{\rm p} w^{\rm p}) &= E-D  \\
    !! \partial_z(a^{\rm p} w^{\rm p} \phi^{\rm p}) &= E \phi^{\rm e} - D \phi^{\rm p}
    !! =   E \left\{ \overline{\phi} + \underbrace{\frac{a^{\rm p}}{1-a^{\rm p}} (\overline{\phi}-\phi^{\rm p})}_{\rm mass\_flux\_small\_ap=False} \right\} - D \phi^{\rm p} \\
    !! w^{\rm p} \partial_z w^{\rm p} &= \left( 1 + \underbrace{\frac{a^{\rm p}}{1-a^{\rm p}}}_{\rm mass\_flux\_small\_ap=False} \right) \left[ b' + b \epsilon \right] (w^{\rm p})^2 + a B^{\rm p} \\
    !! \partial_z(e^{\rm p}-e) &= \left( \frac{E}{-a^{\rm p} w^{\rm p}} \right) \left[ \left( 1 - \underbrace{\frac{a^{\rm p}}{1-a^{\rm p}}}_{\rm mass\_flux\_small\_ap=False}  \right) (e^{\rm p}-e) - \frac{1}{2} \left( 1 + \underbrace{\frac{a^{\rm p}}{1-a^{\rm p}}}_{\rm mass\_flux\_small\_ap=False}  \right) \| \mathbf{v}^{\rm p} - \mathbf{v} \|^2    \right] - \partial_z e - \frac{c_\epsilon}{w^{\rm p} l_{\epsilon}} e^{\rm p} \sqrt{e^{\rm p}} \\
    !! \partial_z(a^{\rm p} w^{\rm p} [u^{\rm p}-C_u \overline{u}]) &=  E \left\{  (1-C_u) \overline{u} + \underbrace{\frac{a^{\rm p}}{1-a^{\rm p}}\left( (1-C_u) \overline{u} -[u^{\rm p}-C_u \overline{u}]   \right)}_{\rm mass\_flux\_small\_ap=False}  \right\} - D [u^{\rm p}-C_u \overline{u}]
    !! \end{align*}<br />
    !
    ! \partial_z(a^{\rm p} w^{\rm p} e^{\rm p}) &=  E \left\{ e + \frac{1}{2} \| \mathbf{v}^{\rm p} - \mathbf{v} \|^2 + \underbrace{\frac{a^{\rm p}}{1-a^{\rm p}} \left( (e-e^{\rm p}) + \frac{1}{2} \| \mathbf{v}^{\rm p} - \mathbf{v} \|^2 \right)}_{\rm mass\_flux\_small\_ap=False} \right\} - D e^{\rm p}  \\
    !!==========================================================================<br />
    USE scm_par, ONLY : grav,ceps_nemo,wpmin
    IMPLICIT NONE
    INTEGER, INTENT(IN   )                 :: N                     !! number of vertical levels
    INTEGER, INTENT(IN   )                 :: ntra                  !! number of tracers
    INTEGER, INTENT(IN   )                 :: nparams               !! number of parameters in the EDOs
    INTEGER, INTENT(IN   )                 :: neos                  !! number of parameters in the EDOs
    REAL(8), INTENT(IN   )                 :: z_w(0:N)              !! depth at cell interfaces [m]
    REAL(8), INTENT(IN   )                 :: Hz(1:N)               !! layer thickness [m]
    REAL(8), INTENT(IN   )                 :: u_m(1:N     )         !! mean zonal velocity [m/s]
    REAL(8), INTENT(IN   )                 :: v_m(1:N     )         !! mean meridional velocity [m/s]
    REAL(8), INTENT(IN   )                 :: t_m(1:N,ntra)         !! mean tracer
    REAL(8), INTENT(IN   )                 :: tke_m(0:N)            !! mean TKE [m2/s2]
    REAL(8), INTENT(IN   )                 :: up0                   !! surface value for plume zonal velocity [m/s]
    REAL(8), INTENT(IN   )                 :: vp0                   !! surface value for plume meridional velocity [m/s]
    REAL(8), INTENT(IN   )                 :: wp0                   !! surface value for plume vertical velocity [m/s]
    REAL(8), INTENT(IN   )                 :: tkep_min
    REAL(8), INTENT(IN   )                 :: mxlp_min
    REAL(8), INTENT(IN   )                 :: tp0(1:ntra)           !! surface value for plume tracers
    REAL(8), INTENT(IN   )                 :: mf_params(1:nparams)  !! parameters in the ODEs
    REAL(8), INTENT(IN   )                 :: eos_params(1:neos)    !! parameters in the EOS (for lin eos only)
    LOGICAL, INTENT(IN   )                 :: small_ap              !! (T) small area approximation (F) no approximation
    LOGICAL, INTENT(IN   )                 :: lin_eos               !!
    INTEGER, INTENT(IN   )                 ::  opt                 !! option for the computation of wtke (=0 PL24, =1 HB09, = 2 Wi11)
    REAL(8), INTENT(IN   )                 :: fcor
    REAL(8), INTENT(  OUT)                 :: a_p(0:N)              !! fractional area occupied by the plume
    REAL(8), INTENT(  OUT)                 :: w_p(0:N)              !! vertical velocity in the plume [m/s]
    REAL(8), INTENT(  OUT)                 :: u_p(0:N)              !! zonal velocity in the plume [m/s]
    REAL(8), INTENT(  OUT)                 :: v_p(0:N)              !! meridional velocity in the plume [m/s]
    REAL(8), INTENT(  OUT)                 :: t_p(0:N,ntra)         !! tracer values in the plume
    REAL(8), INTENT(  OUT)                 :: B_p(0:N)              !! buoyancy forcing  [m/s2]
    REAL(8), INTENT(  OUT)                 :: ent(1:N)              !! diagnostics : entrainment [m-1]
    REAL(8), INTENT(  OUT)                 :: det(1:N)              !! diagnostics : detrainment [m-1]
    REAL(8), INTENT(  OUT)                 :: eps(1:N)          !! diagnostics : TKE dissipation [m2 s-3]
    REAL(8), INTENT(  OUT)                 :: vort_p(0:N)
    REAL(8), INTENT(INOUT)                 :: zinv                  !! depth at which w_p = wmin  [m]
    ! local variables
    REAL(8)                                :: delta0
    INTEGER                                :: k,itrc,iter,tke_comput
    REAL(8)                                :: cff,cff1, rho_m(1:N), N2sfc
    REAL(8)                                :: tke_env, u_env,v_env,t_env
    REAL(8)                                :: aa,bb,bp,beta1,zbb,zbp
    REAL(8)                                :: hinv,rho_p,rho_w,idwdz0
    REAL(8)                                :: cffu,cffv,cffw,apr,normVel
    REAL(8)                                :: temp_p,salt_p, Cu, Cv, tpm, apm, beta2, cffmax
    LOGICAL                                :: found    = .false.
    REAL(8)                                :: mxld(0:N), imxld0(1:N)
    REAL(8)                                :: lup, ldwn, epsilon, dtke, rn2,cffp
    !REAL(8)                                :: wpmin
    !=======================================================================
    tke_comput = 1
    !=======================================================================
    beta1 = mf_params(1); aa     = mf_params(3)
    bb    = mf_params(4); bp     = mf_params(5)/ABS(zinv)
    Cu    = mf_params(6); Cv     = mf_params(7)
    beta2 = mf_params(2); delta0 = mf_params(9)/ABS(zinv)
    !=======================================================================
    ! initialize plume properties with surface values
    a_p(N) = mf_params(8)         ; a_p(0:N-1       ) = 0.
    w_p(N) = wp0                  ; w_p(0:N-1       ) = 0.
    u_p(N) = (1.-Cu)*up0          ; u_p(0:N-1       ) = 0.
    v_p(N) = (1.-Cv)*vp0          ; v_p(0:N-1       ) = 0.
    t_p(N,1:ntra) = tp0(1:ntra)   ; t_p(0:N-1,1:ntra) = 0.
    B_p(0:N) = 0.; ent(1:N) = 0.  ; det(1:N) = 0.  ; eps(1:N) = 0.
    !! make a motivated CHOICE for the b.c. of vorticity... 
    ! vort_p(N) = fcor
    vort_p(N)=0.
    !
    IF(tke_comput == 0) THEN
        t_p(N,ntra) = 0.  ! in this case t_p(ntra) contains t_p - tke
    ELSE
        t_p(N,ntra) = tke_m(N)
    ENDIF
    !
    DO k = 1,N
      cff       = 0.5*(z_w(k)+z_w(k-1))
      lup       = MAX( -cff     , mxlp_min)
      ldwn      = MAX(  cff-zinv, mxlp_min)
      imxld0(k) = 1./SQRT(lup*ldwn)
    ENDDO
    imxld0(N) = 0.5*(3.*imxld0(N-1)-imxld0(N-2))
    !imxld0(1:N) = 0.
    !
    !=======================================================================
    ! Linear eos version
    !=======================================================================
    IF(lin_eos) THEN
      DO k = 1,N
        CALL eos_val_lin(t_m(k,1),t_m(k,2),eos_params,neos,rho_m(k))
      ENDDO
    ELSE
      DO k = 1,N
        CALL eos_val(t_m(k,1),t_m(k,2),0.5*(z_w(k)+z_w(k-1)),eos_params(1),rho_m(k))
      ENDDO
    ENDIF
    !=======================================================================
    N2sfc  = -(grav/eos_params(1))*(rho_m(N)-rho_m(N-1))/Hz(N)
    hinv   = z_w(N)
    !=======================================================================
    DO k=N,1,-1
      ! Compute B_p
      temp_p = t_p(k,1); salt_p = t_p(k,2)
      !
      IF(lin_eos) THEN; CALL eos_val_lin(temp_p,salt_p,eos_params,neos,rho_p)
      ELSE; CALL eos_val(temp_p,salt_p,0.5*(z_w(k)+z_w(k-1)),eos_params(1),rho_p); ENDIF
      !! Compute \( B^{\rm p}_{k} \) \[ B^{\rm p}_{k} = - \frac{g}{\rho_0} \left( \rho^{\rm p}_{k+1/2} - \overline{\rho}_k \right) \]
      B_p(k) = - grav * ( rho_p - rho_m(k) ) / eos_params(1)
      ! Wp equation first
      !! If \( {\rm small\_ap = False} : (b',b) \rightarrow \frac{(b',b)}{1-a^{\rm p}_{k+1/2}} \) <br />
      cff    = 1.
      IF(.not.small_ap) cff = 1./(1.-a_p(k))
      zbb = cff*bb; zbp = cff*bp
      !! Compute \( w^{\rm p}_{k-1/2} \) :: call \(  {\rm get\_w\_p\_R10}  \)
      !  \[ (w^{\rm p})^{2}_{k+1/2} - (w^{\rm p})^{2}_{k-1/2} =
      !  h_k (b' + b \epsilon_k) \left((w^{\rm p})^{2}_{k+1/2} + (w^{\rm p})^{2}_{k-1/2})\right)
      ! + 2 a h_k B^{\rm p}_{k}\]
      CALL get_w_p_R10(w_p(k-1),w_p(k),aa,zbb,zbp,beta1,Hz(k),B_p(k),wpmin,hinv,found)
      ! diagnostics
      cff       = (w_p(k)-w_p(k-1))/(0.5*Hz(k)*(w_p(k)+w_p(k-1)))
      ent  (k)  = MAX(0., -beta1*cff)
      det  (k)  = MAX(0.,  beta2*cff) + delta0
      !IF( det(k)-1. == det(k) ) THEN
      !  print*,'det(k) = ',det(k),cff,w_p(k),beta2
      !ENDIF
      !
      IF(found) THEN
        zinv     = z_w(k)-hinv; found = .false.
      ENDIF
      !! Compute \( a^{\rm p}_{k-1/2} \) :: call \(  {\rm get\_a\_p\_R10}  \)
      ! \[  (a^{\rm p} w^{\rm p})_{k-\frac{1}{2}}  =   (a^{\rm p} w^{\rm p})_{k+\frac{1}{2}}
      ! - \beta_1 \left(\frac{a^{\rm p}_{k+\frac{1}{2}} + a^{\rm p}_{k-\frac{1}{2}} }{2}\right)
      ! \left(  \max(0, (\delta_z w^{\rm p})_k) + \min\left( \frac{w^{\rm p} h_k}{\beta_1} \delta_0 ,
      ! (\delta_z w^{\rm p})_k + \frac{w^{\rm p} h_k}{\beta_1} (\delta_1)_k \right) \right) \\  \]
      CALL get_a_p_R10(a_p(k-1),a_p(k),w_p(k-1),w_p(k),beta1,beta2,Hz(k),delta0,wpmin)
      cff = a_p(k)/(1.-a_p(k))
      IF(small_ap) cff = 0.
      ! Compute tracers (except TKE_p)
      DO itrc = 1,ntra-1
        t_env = t_m(k,itrc) + cff*(t_m(k,itrc)-t_p(k,itrc)) !! Compute environment \( \phi^{\rm e}_k\)
              !!\begin{align*}
              !! \phi^{\rm e}_k &= \overline{\phi}_k  \hspace{7cm} \mbox{small_ap = True} \\
              !! \phi^{\rm e}_k &= \overline{\phi}_k + \left( \frac{a^{\rm p}_{k+1/2}}{1-a^{\rm p}_{k+1/2}} \right) ( \overline{\phi}_k - \phi^{\rm p}_{k+1/2} ) \hspace{1cm} \mbox{small_ap = False}
              !! \end{align*}
        CALL get_t_p_R10(t_p(k-1,itrc),t_p(k,itrc),t_env,  &
                     a_p(k-1),a_p(k),w_p(k-1),w_p(k),beta1,beta2,Hz(k),delta0,wpmin) !! Compute \( \phi^{\rm p}_{k-1/2}\) :: call \(  {\rm get\_t\_p\_R10}  \)
      ENDDO
      ! Compute up - Cu umean
      u_env = (1.-Cu)*u_m(k) + cff*( (1.-Cu)*u_m(k) - u_p(k) )
      CALL get_t_p_R10(u_p(k-1),u_p(k),u_env,a_p(k-1),a_p(k),w_p(k-1),w_p(k),   &
                                                      beta1,beta2,Hz(k),delta0,wpmin) !! Compute \( u^{\rm p}_{k-1/2}\) :: call \(  {\rm get\_t\_p\_R10}  \)
      ! Compute vp - Cv vmean
      v_env = (1.-Cv)*v_m(k) + cff*( (1.-Cv)*v_m(k) - v_p(k) )
      CALL get_t_p_R10(v_p(k-1),v_p(k),v_env,a_p(k-1),a_p(k),w_p(k-1),w_p(k),   &
                                                      beta1,beta2,Hz(k),delta0,wpmin) !! Compute \( v^{\rm p}_{k-1/2}\) :: call \(  {\rm get\_t\_p\_R10}  \)
      ! Compute TKEplume - TKEmean
      IF(tke_comput == 0) THEN
        cff       = a_p(k)/(1.-a_p(k))
        IF(small_ap) cff = 0.
        cffw      = 0.5*(w_p(k)+w_p(k-1))                    ! corresponds to wp
        cffu      = 0.5*(u_p(k)+u_p(k-1)) + (Cu-1.)*u_m(k)   ! corresponds to up-um
        cffv      = 0.5*(v_p(k)+v_p(k-1)) + (Cv-1.)*v_m(k)   ! corresponds to vp-vm
        normVel   = (1.+cff)*(cffu*cffu+cffv*cffv+cffw*cffw)/(1.-cff)  !
        dtke      = tke_m(k)-tke_m(k-1)
        epsilon   = ceps_nemo * (t_p(k,ntra)+tke_m(k))   &
                                        * SQRT(t_p(k,ntra)+tke_m(k)) * imxld0(k)
        IF(opt==1) THEN  ! HB09 tke_p equation is simply   d tke_p / dz = E tke - D tke_p
          epsilon = 0.; normVel = 0.
        ENDIF
        CALL get_dtke_p_R10(t_p(k-1,ntra),t_p(k,ntra),dtke,normVel,w_p(k-1),w_p(k),  &
                                                   Hz(k),epsilon,(1.-cff)*beta1) !! Compute \( e^{\rm p}_{k-1/2}\) :: call \(  {\rm get\_dtke\_p\_R10}  \)
        eps(k)    = epsilon
      ELSE
      !=====================================================================
        cff = 1./(1.-a_p(k))
        IF(small_ap) cff = 1.  ! small ap limit
        cffw      = 0.5*(w_p(k)+w_p(k-1))                    ! corresponds to wp
        cffu      = 0.5*(u_p(k)+u_p(k-1)) + (Cu-1.)*u_m(k)   ! corresponds to up-um
        cffv      = 0.5*(v_p(k)+v_p(k-1)) + (Cv-1.)*v_m(k)   ! corresponds to vp-vm
        normVel   = (cffu*cffu+cffv*cffv+cffw*cffw)          !
        epsilon   = ceps_nemo * t_p(k,ntra) * SQRT(t_p(k,ntra)) * imxld0(k)
        IF(opt==1) THEN  ! HB09 tke_p equation is simply   d tke_p / dz = E tke - D tke_p
          epsilon = 0.; normVel = 0.
        ENDIF
        CALL get_tke_p_R10( t_p(k-1,ntra),t_p(k,ntra), tke_m(k-1), tke_m(k), w_p(k-1),w_p(k),normVel,epsilon,Hz(k),cff*beta1 )
        eps(k)    = epsilon
      ENDIF
      !=====================================================================
      CALL get_vort_p_R10(vort_p(k-1),vort_p(k),fcor,a_p(k-1),a_p(k),w_p(k-1),w_p(k), &
      beta1,beta2,Hz(k),delta0,wpmin)
    ENDDO
    !=======================================================================
    ! At this point, up and vp contain up-Cu ue  and vp-Cv ve
    u_p(N) = up0; v_p(N) = vp0
    DO k=1,N-1
      u_p(k) = u_p(k) + Cu*u_m(k)
      v_p(k) = v_p(k) + Cv*v_m(k)
    ENDDO
    u_p(0) = u_p(1); v_p(0) = v_p(1)
    IF(tke_comput == 0) THEN
      DO k = 0,N
        t_p(k,ntra) = MAX( t_p(k,ntra) + tke_m(k), tkep_min )
      ENDDO
    ELSE
      DO k = 0,N
        t_p(k,ntra) = MAX( t_p(k,ntra), tkep_min )
      ENDDO
    ENDIF
  !---------------------------------------------------------------------------------------------------
  END SUBROUTINE mass_flux_R10
  !===================================================================================================

  !===================================================================================================
  SUBROUTINE mass_flux_R10_cor(u_m,v_m,t_m,tke_m,z_w,Hz,tp0,up0,vp0,wp0,mf_params,eos_params,fcor,ecor, &
    tkep_min, mxlp_min, small_ap,lin_eos,zinv,N,ntra,nparams,neos,  &
    a_p,u_p,v_p,w_p,t_p,B_p,ent,det,vort_p,eps)
    !---------------------------------------------------------------------------------------------------
    !!==========================================================================<br />
    !!                  ***  ROUTINE mass_flux_R10  ***                         <br />
    !! ** Purposes : solve mass-flux equations for Rio et al. (2010)
    !!                        entrainment/detrainement closure                  <br />
    !! \begin{align*} E &= a^{\rm p} \beta_1 \max(0, \partial_z w^{\rm p}) \\
    !!                D &= - a^{\rm p} \beta_2 \min\left( 0 , \partial_z w^{\rm p} \right) -  a^{\rm p} w^{\rm p} \delta_0  \end{align*} <br />
    !!@note we must have \( \beta_1 < 1 \) and \( 1 < \beta_2 < 2 \) otherwise unphysical values of \( a^{\rm p} \) are obtained  @endnote <br />
    !! \begin{align*}
    !! \partial_z(a^{\rm p} w^{\rm p}) &= E-D  \\
    !! \partial_z(a^{\rm p} w^{\rm p} \phi^{\rm p}) &= E \phi^{\rm e} - D \phi^{\rm p}
    !! =   E \left\{ \overline{\phi} + \underbrace{\frac{a^{\rm p}}{1-a^{\rm p}} (\overline{\phi}-\phi^{\rm p})}_{\rm mass\_flux\_small\_ap=False} \right\} - D \phi^{\rm p} \\
    !! w^{\rm p} \partial_z w^{\rm p} &= \left( 1 + \underbrace{\frac{a^{\rm p}}{1-a^{\rm p}}}_{\rm mass\_flux\_small\_ap=False} \right) \left[ b' + b \epsilon \right] (w^{\rm p})^2 + a B^{\rm p} \\
    !! \partial_z(e^{\rm p}-e) &= \left( \frac{E}{-a^{\rm p} w^{\rm p}} \right) \left[ \left( 1 - \underbrace{\frac{a^{\rm p}}{1-a^{\rm p}}}_{\rm mass\_flux\_small\_ap=False}  \right) (e^{\rm p}-e) - \frac{1}{2} \left( 1 + \underbrace{\frac{a^{\rm p}}{1-a^{\rm p}}}_{\rm mass\_flux\_small\_ap=False}  \right) \| \mathbf{v}^{\rm p} - \mathbf{v} \|^2    \right] - \partial_z e - \frac{c_\epsilon}{w^{\rm p} l_{\epsilon}} e^{\rm p} \sqrt{e^{\rm p}} \\
    !! \partial_z(a^{\rm p} w^{\rm p} [u^{\rm p}-C_u \overline{u}]) &=  E \left\{  (1-C_u) \overline{u} + \underbrace{\frac{a^{\rm p}}{1-a^{\rm p}}\left( (1-C_u) \overline{u} -[u^{\rm p}-C_u \overline{u}]   \right)}_{\rm mass\_flux\_small\_ap=False}  \right\} - D [u^{\rm p}-C_u \overline{u}]
    !! \end{align*}<br />
    !
    ! \partial_z(a^{\rm p} w^{\rm p} e^{\rm p}) &=  E \left\{ e + \frac{1}{2} \| \mathbf{v}^{\rm p} - \mathbf{v} \|^2 + \underbrace{\frac{a^{\rm p}}{1-a^{\rm p}} \left( (e-e^{\rm p}) + \frac{1}{2} \| \mathbf{v}^{\rm p} - \mathbf{v} \|^2 \right)}_{\rm mass\_flux\_small\_ap=False} \right\} - D e^{\rm p}  \\
    !!==========================================================================<br />
    USE scm_par, ONLY : grav,ceps_nemo,wpmin
    IMPLICIT NONE
    INTEGER, INTENT(IN   )                 :: N                     !! number of vertical levels
    INTEGER, INTENT(IN   )                 :: ntra                  !! number of tracers
    INTEGER, INTENT(IN   )                 :: nparams               !! number of parameters in the EDOs
    INTEGER, INTENT(IN   )                 :: neos                  !! number of parameters in the EDOs
    REAL(8), INTENT(IN   )                 :: z_w(0:N)              !! depth at cell interfaces [m]
    REAL(8), INTENT(IN   )                 :: Hz(1:N)               !! layer thickness [m]
    REAL(8), INTENT(IN   )                 :: u_m(1:N     )         !! mean zonal velocity [m/s]
    REAL(8), INTENT(IN   )                 :: v_m(1:N     )         !! mean meridional velocity [m/s]
    REAL(8), INTENT(IN   )                 :: t_m(1:N,ntra)         !! mean tracer
    REAL(8), INTENT(IN   )                 :: tke_m(0:N)            !! mean TKE [m2/s2]
    REAL(8), INTENT(IN   )                 :: up0                   !! surface value for plume zonal velocity [m/s]
    REAL(8), INTENT(IN   )                 :: vp0                   !! surface value for plume meridional velocity [m/s]
    REAL(8), INTENT(IN   )                 :: wp0                   !! surface value for plume vertical velocity [m/s]
    REAL(8), INTENT(IN   )                 :: tp0(1:ntra)           !! surface value for plume tracers
    REAL(8), INTENT(IN   )                 :: mf_params(1:nparams)  !! parameters in the ODEs
    REAL(8), INTENT(IN   )                 :: eos_params(1:neos)    !! parameters in the EOS (for lin eos only)
    REAL(8), INTENT(IN   )                 :: fcor                  !! Coriolis frequency [s-1]
    REAL(8), INTENT(IN   )                 :: ecor                  !! NT Coriolis frequency [s-1]
    REAL(8), INTENT(IN   )                 :: tkep_min
    REAL(8), INTENT(IN   )                 :: mxlp_min
    LOGICAL, INTENT(IN   )                 :: small_ap              !! (T) small area approximation (F) no approximation
    LOGICAL, INTENT(IN   )                 :: lin_eos               !!
    REAL(8), INTENT(  OUT)                 :: a_p(0:N)              !! fractional area occupied by the plume
    REAL(8), INTENT(  OUT)                 :: w_p(0:N)              !! vertical velocity in the plume [m/s]
    REAL(8), INTENT(  OUT)                 :: u_p(0:N)              !! zonal velocity in the plume [m/s]
    REAL(8), INTENT(  OUT)                 :: v_p(0:N)              !! meridional velocity in the plume [m/s]
    REAL(8), INTENT(  OUT)                 :: t_p(0:N,ntra)         !! tracer values in the plume
    REAL(8), INTENT(  OUT)                 :: B_p(0:N)              !! buoyancy forcing  [m/s2]
    REAL(8), INTENT(  OUT)                 :: ent(1:N)              !! diagnostics : entrainment [m-1]
    REAL(8), INTENT(  OUT)                 :: det(1:N)              !! diagnostics : detrainment [m-1]
    REAL(8), INTENT(  OUT)                 :: vort_p(0:N)
    REAL(8), INTENT(  OUT)                 :: eps(1:N)          !! diagnostics : TKE dissipation [m2 s-3]
    REAL(8), INTENT(INOUT)                 :: zinv                  !! depth at which w_p = wmin  [m]
    ! local variables
    REAL(8)                                :: du_m(0:N)          ! generalized up - Cu*ubar + z fperp
    REAL(8)                                :: dv_m(0:N)          ! generalized vp - Cv*vbar
    REAL(8)                                :: delta0
    INTEGER                                :: k,itrc,iter
    REAL(8)                                :: cff,cff1, rho_m(1:N), N2sfc
    REAL(8)                                :: tke_env, u_env,v_env,t_env
    REAL(8)                                :: aa,bb,bp,beta1,zbb,zbp
    REAL(8)                                :: hinv,rho_p,rho_w,idwdz0
    REAL(8)                                :: cffu,cffv,cffw,apr,normVel
    REAL(8)                                :: temp_p,salt_p, Cu, Cv, tpm, apm, beta2, cffmax
    LOGICAL                                :: found    = .false.
    REAL(8)                                :: mxld(0:N), imxld0(1:N)
    REAL(8)                                :: lup, ldwn, epsilon, dtke, rn2
    REAL(8)                                :: BpTilde,cor_u,cor_v,frc_u,frc_v
    !======================================================================
    !print*,'in mass_flux_R10_cor **',zinv
    ! unpack parameters (mf_params  = [Cent,Cdet,wp_a,wp_b,wp_bp,up_c,vp_c,bc_ap,delta0]
    !!!!! TEST modulation of plume radius by Ro (not working)
    ! REAL(8)                                :: B0, Ro
    ! B0=7.302579717857959e-08
    ! Ro = min((B0/fcor)**(1/2)/(fcor*abs(zinv)) , (B0/fcor)**(1/2)/(fcor*1000))
    ! cff= tanh(Ro**2)
    ! beta1 = mf_params(1); aa     = mf_params(3)
    ! bb    = mf_params(4); bp     = mf_params(5)/ABS(zinv*cff) !bp     = mf_params(5)/(-zinv)
    ! Cu    = mf_params(6); Cv     = mf_params(7)
    ! beta2 = mf_params(2); delta0 = mf_params(9)/ABS(zinv*cff) !delta0 = mf_params(9)/(-zinv)
    !!!!!!!!!
    beta1 = mf_params(1); aa     = mf_params(3)
    bb    = mf_params(4); bp     = mf_params(5)/ABS(zinv) !bp     = mf_params(5)/(-zinv)
    Cu    = mf_params(6); Cv     = mf_params(7)
    beta2 = mf_params(2); delta0 = mf_params(9)/ABS(zinv) !delta0 = mf_params(9)/(-zinv)
    !=======================================================================
    ! initialize plume properties with surface values
    a_p(N) = mf_params(8) ; a_p(0:N-1       ) = 0.
    w_p(N) = wp0                  ; w_p(0:N-1       ) = 0.
    !
    v_p(N) = vp0 ; v_p(0:N-1) = 0.;  u_p(N) = up0 ; u_p(0:N-1) = 0.
    !
    t_p(N,1:ntra) = tp0(1:ntra)   ; t_p(0:N-1,1:ntra) = 0.
    B_p(0:N) = 0.; ent(1:N) = 0.  ; det(1:N) = 0.  ; eps(1:N) = 0.
    !
    vort_p(N) = fcor
    !
    DO k = 2,N
      du_m(k) = u_m(k) - u_m(k-1)
      dv_m(k) = v_m(k) - v_m(k-1)
    ENDDO
    du_m(1) = du_m(2)
    dv_m(1) = dv_m(2)
    !
    t_p(N,ntra) = 0.  ! in this case t_p(ntra) contains t_p - tke
    !
    DO k = 1,N
      cff       = 0.5*(z_w(k)+z_w(k-1))
      lup       = MAX( -cff     , mxlp_min)
      ldwn      = MAX(  cff-zinv, mxlp_min)
      imxld0(k) = 1./SQRT(lup*ldwn)
    ENDDO
    imxld0(N) = 0.5*(3.*imxld0(N-1)-imxld0(N-2))
    !
    !=======================================================================
    ! Linear eos version
    !=======================================================================
    IF(lin_eos) THEN
      DO k = 1,N
        CALL eos_val_lin(t_m(k,1),t_m(k,2),eos_params,neos,rho_m(k))
      ENDDO
    ELSE
      DO k = 1,N
        CALL eos_val(t_m(k,1),t_m(k,2),0.5*(z_w(k)+z_w(k-1)),eos_params(1),rho_m(k))
      ENDDO
    ENDIF
    !=======================================================================
    N2sfc  = -(grav/eos_params(1))*(rho_m(N)-rho_m(N-1))/Hz(N)
    !=======================================================================
    DO k=N,1,-1

      ! TESTS : Modulate ent/det by vorticity
      ! fait pas grand chose...
      !
      ! cff=(1+EXP(-vort_p(k)/fcor))*0.5
      ! if (fcor/vort_p(k)>0) then
      !   cff = 1
      ! else
      !   cff=tanh((abs(fcor/vort_p(k)))**0.37)
      ! endif 
      ! cff=tanh(abs((1-abs(vort_p(k)/fcor)))**0.37+0.31)
      ! !
      ! beta1 = cff*mf_params(1)
      ! beta2 = cff*mf_params(2)
      ! ! bb    = cff*mf_params(4)
      ! ! aa    = cff*mf_params(3)
      ! delta0= (1/cff)*mf_params(9)/ABS(zinv)
      ! bp    = (1/cff)*mf_params(5)/ABS(zinv)
      !
      !
      ! Compute B_p
      temp_p = t_p(k,1); salt_p = t_p(k,2)
      !
      IF(lin_eos) THEN; CALL eos_val_lin(temp_p,salt_p,eos_params,neos,rho_p)
      ELSE; CALL eos_val(temp_p,salt_p,0.5*(z_w(k)+z_w(k-1)),eos_params(1),rho_p); ENDIF
      !! Compute \( B^{\rm p}_{k} \) \[ B^{\rm p}_{k} = - \frac{g}{\rho_0} \left( \rho^{\rm p}_{k+1/2} - \overline{\rho}_k \right) \]
      B_p(k) = - grav * ( rho_p - rho_m(k) ) / eos_params(1)
      ! Wp equation first
      !! If \( {\rm small\_ap = False} : (b',b) \rightarrow \frac{(b',b)}{1-a^{\rm p}_{k+1/2}} \) <br />
      cff    = 1.
      IF(.not.small_ap) cff = 1./(1.-a_p(k))
      zbb = cff*bb; zbp = cff*bp
      BpTilde = B_p(k) + (ecor/aa)*u_p(k)
      !! Compute \( w^{\rm p}_{k-1/2} \) :: call \(  {\rm get\_w\_p\_R10}  \)
      !  \[ (w^{\rm p})^{2}_{k+1/2} - (w^{\rm p})^{2}_{k-1/2} =
      !  h_k (b' + b \epsilon_k) \left((w^{\rm p})^{2}_{k+1/2} + (w^{\rm p})^{2}_{k-1/2})\right)
      ! + 2 a h_k B^{\rm p}_{k}\]
      CALL get_w_p_R10(w_p(k-1),w_p(k),aa,zbb,zbp,beta1,Hz(k),BpTilde,wpmin,hinv,found)
      ! diagnostics
      cff       = (w_p(k)-w_p(k-1))/(0.5*Hz(k)*(w_p(k)+w_p(k-1)))
      ent  (k)  = MAX(0., -beta1*cff)
      det  (k)  = MAX(0.,  beta2*cff) + delta0
      !
      IF(found) THEN
        zinv     = z_w(k)-hinv; found = .false.
      ENDIF
      !! Compute \( a^{\rm p}_{k-1/2} \) :: call \(  {\rm get\_a\_p\_R10}  \)
      ! \[  (a^{\rm p} w^{\rm p})_{k-\frac{1}{2}}  =   (a^{\rm p} w^{\rm p})_{k+\frac{1}{2}}
      ! - \beta_1 \left(\frac{a^{\rm p}_{k+\frac{1}{2}} + a^{\rm p}_{k-\frac{1}{2}} }{2}\right)
      ! \left(  \max(0, (\delta_z w^{\rm p})_k) + \min\left( \frac{w^{\rm p} h_k}{\beta_1} \delta_0 ,
      ! (\delta_z w^{\rm p})_k + \frac{w^{\rm p} h_k}{\beta_1} (\delta_1)_k \right) \right) \\  \]
      CALL get_a_p_R10(a_p(k-1),a_p(k),w_p(k-1),w_p(k),beta1,beta2,Hz(k),delta0,wpmin)
      cff = a_p(k)/(1.-a_p(k))
      IF(small_ap) cff = 0.
      ! Compute tracers (except TKE_p)
      DO itrc = 1,ntra-1
        t_env = t_m(k,itrc) + cff*(t_m(k,itrc)-t_p(k,itrc)) !! Compute environment \( \phi^{\rm e}_k\)
        !!\begin{align*}
        !! \phi^{\rm e}_k &= \overline{\phi}_k  \hspace{7cm} \mbox{small_ap = True} \\
        !! \phi^{\rm e}_k &= \overline{\phi}_k + \left( \frac{a^{\rm p}_{k+1/2}}{1-a^{\rm p}_{k+1/2}} \right) ( \overline{\phi}_k - \phi^{\rm p}_{k+1/2} ) \hspace{1cm} \mbox{small_ap = False}
        !! \end{align*}
        CALL get_t_p_R10(t_p(k-1,itrc),t_p(k,itrc),t_env,  &
        a_p(k-1),a_p(k),w_p(k-1),w_p(k),beta1,beta2,Hz(k),delta0,wpmin) !! Compute \( \phi^{\rm p}_{k-1/2}\) :: call \(  {\rm get\_t\_p\_R10}  \)
      ENDDO
      !=======
      u_env = u_m(k); v_env = v_m(k) ! m/s
      frc_u = Cu*du_m(k)-Hz(k)*ecor  ! m/s
      frc_v = Cv*dv_m(k)             ! m/s
      ! frc_u = Cu*du_m(k)-Hz(k)*ecor + Cu*u_p(k) ! m/s, test
      ! frc_v = Cv*dv_m(k)            + Cu*v_p(k) ! m/s, test
      !==v====
      ! Compute up
      cor_u   =   fcor*Hz(k)*v_p(k) ! m2/s2
      CALL get_u_p_R10(u_p(k-1),u_p(k),u_env,a_p(k-1),a_p(k),w_p(k-1),w_p(k),   &
                    beta1,beta2,Hz(k),delta0,frc_u,cor_u,wpmin) !! Compute \( u^{\rm p}_{k-1/2}\) :: call \(  {\rm get\_t\_p\_R10}  \)
      !if( u_p(k-1).ne.u_p(k-1) ) then
      !  print*,'NaN in get_u_p_R10 :: ',u_p(k-1),u_p(k),u_env,a_p(k-1),a_p(k)
      !  stop
      !endif
      ! Compute vp
      cor_v   = -fcor*Hz(k)*u_p(k) ! m2/s2
      CALL get_u_p_R10(v_p(k-1),v_p(k),v_env,a_p(k-1),a_p(k),w_p(k-1),w_p(k),   &
                    beta1,beta2,Hz(k),delta0,frc_v,cor_v,wpmin) !! Compute \( v^{\rm p}_{k-1/2}\) :: call \(  {\rm get\_t\_p\_R10}  \)
      !if( v_p(k-1).ne.v_p(k-1) ) then
      !  print*,'NaN in get_v_p_R10 :: ',v_p(k-1),v_p(k),v_env,a_p(k-1),a_p(k)
      !  stop
      !endif
      ! Compute TKEplume - TKEmean
      cff       = a_p(k)/(1.-a_p(k))
      IF(small_ap) cff = 0.
      cffw      = 0.5*(w_p(k)+w_p(k-1))
      cffu      = 0.5*(u_p(k)+u_p(k-1))
      cffv      = 0.5*(v_p(k)+v_p(k-1))
      normVel   = (1.+cff)*(cffu*cffu+cffv*cffv+cffw*cffw)/(1.-cff)
      dtke      = tke_m(k)-tke_m(k-1)
      epsilon   = ceps_nemo * (t_p(k,ntra)+tke_m(k))   &
                * SQRT(t_p(k,ntra)+tke_m(k)) * imxld0(k)
      CALL get_dtke_p_R10(t_p(k-1,ntra),t_p(k,ntra),dtke,normVel,w_p(k-1),w_p(k),  &
                            Hz(k),epsilon,(1.-cff)*beta1) !! Compute \( e^{\rm p}_{k-1/2}\) :: call \(  {\rm get\_dtke\_p\_R10}  \)
      eps(k)    = epsilon
      !
      CALL get_vort_p_R10(vort_p(k-1),vort_p(k),fcor,a_p(k-1),a_p(k),w_p(k-1),w_p(k), &
                                                      beta1,beta2,Hz(k),delta0,wpmin)
      !
    ENDDO
    !=======================================================================
    ! At this point, up and vp contain up-Cu ue  and vp-Cv ve
    !u_p(N) = up0; v_p(N) = vp0
    !DO k=1,N-1
    !  u_p(k) = u_p(k) + Cu*u_m(k) - z_w(k)*ecor
    !  v_p(k) = v_p(k) + Cv*v_m(k)
    !ENDDO
    !u_p(0) = u_p(1); v_p(0) = v_p(1)
    !
    DO k = 0,N
    t_p(k,ntra) = MAX( t_p(k,ntra) + tke_m(k), tkep_min )
    ENDDO
    !
    !print*,'out mass_flux_R10_cor **',zinv,u_p(0:10),v_p(0:10)
  !---------------------------------------------------------------------------------------------------
  END SUBROUTINE mass_flux_R10_cor
  !===================================================================================================



END MODULE scm_mfc

MODULE scm_par
   IMPLICIT NONE
   PUBLIC
   REAL(8), PUBLIC, PARAMETER ::   grav     = 9.81           !! Gravity of Earth
   REAL(8), PUBLIC, PARAMETER ::   vkarmn   = 0.41           !! Von Karman constant
!   REAL(8), PUBLIC, PARAMETER ::   rho0     = 1027.          !! Boussinesq reference density [kg/m3]
!   REAL(8), PUBLIC, PARAMETER ::   cp       = 3985.0d0       !! Specific heat capacity of saltwater [J/kg K]
   REAL(8), PUBLIC, PARAMETER ::   rpi      = 4.*ATAN(1.)    !! \( \pi \)
   REAL(8), PUBLIC, PARAMETER ::   wpmin    = 1.e-08         !! Minimum value of \( w^{\rm p} \) [m/s] <br />
   !
   !---------------------------------------------------------------------
   ! SCM parameters for TKE turbulent closure <br />
   !---------------------------------------------------------------------
   REAL(8), PUBLIC, PARAMETER ::   tke_min0  = 1.e-4          !! surface minimum value of tke for Dirichlet condition [m2/s2]
   ! NEMO constants
   REAL(8), PUBLIC, PARAMETER ::   ceps_nemo = 0.5 * sqrt(2.)   !! Constant \( c_\epsilon \) in NEMO
   REAL(8), PUBLIC, PARAMETER ::     cm_nemo = 0.1              !! Constant \( c_m \) in NEMO
   REAL(8), PUBLIC, PARAMETER ::     ct_nemo = 0.1              !! Constant \( c_s \) in NEMO
   REAL(8), PUBLIC, PARAMETER ::     ce_nemo = 0.1              !! Constant \( c_e \) in NEMO
   REAL(8), PUBLIC, PARAMETER ::    Ric_nemo = 2./(2.+ceps_nemo/cm_nemo) !! Critical Richardson number
   ! MNH constants
   REAL(8), PUBLIC, PARAMETER ::   ceps_mnh = 0.845   !! Constant \( c_\epsilon \) in MesoNH
   REAL(8), PUBLIC, PARAMETER ::     cm_mnh = 0.126   !! Constant \( c_m \) in MesoNH
   REAL(8), PUBLIC, PARAMETER ::     ct_mnh = 0.143   !! Constant \( c_s \) in MesoNH
   REAL(8), PUBLIC, PARAMETER ::     ce_mnh = 0.34    !! Constant \( c_e \) in MesoNH
   REAL(8), PUBLIC, PARAMETER ::    Ric_mnh = 0.143   !! Critical Richardson number
   ! RS81 constants
   REAL(8), PUBLIC, PARAMETER ::   ceps_r81 = 0.7     !! Constant \( c_\epsilon \) in Redelsperger & Sommeria 1981
   REAL(8), PUBLIC, PARAMETER ::     cm_r81 = 0.0667  !! Constant \( c_m \) in Redelsperger & Sommeria 1981
   REAL(8), PUBLIC, PARAMETER ::     ct_r81 = 0.167   !! Constant \( c_s \) in Redelsperger & Sommeria 1981
   REAL(8), PUBLIC, PARAMETER ::     ce_r81 = 0.4     !! Constant \( c_e \) in Redelsperger & Sommeria 1981
   REAL(8), PUBLIC, PARAMETER ::    Ric_r81 = 0.139   !! Critical Richardson number
   !
   REAL(8), PUBLIC, PARAMETER ::   avmolT      = 1.4D-07
   REAL(8), PUBLIC, PARAMETER ::   avmolu      = 1.3D-06
   REAL(8), PUBLIC, PARAMETER ::   mxl_min  = 10.      !(avm_bak / 0.1) / sqrt( tke_min )  !! minimum value for mixing lengths [m]
   REAL(8), PUBLIC, PARAMETER ::   mxl_min0 = 0.04    !! minimum surface value for miwing lengths [m]
   REAL(8), PUBLIC, PARAMETER ::   pdlrmin  = 0.1     !! minimum value for the inverse Prandtl number
!   REAL(8), PUBLIC, PARAMETER ::   tke_min  = 1.e-8   !( avm_bak/(cm_nemo*mxl_min) )**2
   !
   REAL(8), PUBLIC, PARAMETER ::  bshear    =   1.e-20 !! minimum shear
   REAL(8), PUBLIC, PARAMETER ::  rsmall    =   1.e-20
   !
   REAL(8), PARAMETER ::  galp    =  0.53    ! parameter for Galperin mixing length limitation
   REAL(8), PARAMETER ::  chk     =  1400./grav ! charnock coefficient
   REAL(8), PARAMETER ::  Zosmin  =  1.e-2   ! min surface roughness length
   REAL(8), PARAMETER ::  Zobmin  =  1.e-4   ! min bottom  roughness length
   !
   !=======================================================================
   ! CANUTO A stability function (Canuto & al. A 2001)
   !=======================================================================
   REAL(8),  parameter                :: ce1         =  1.44D0
   REAL(8),  parameter                :: ce2         =  1.92D0
   REAL(8),  parameter                :: ce3minus    = -0.62D0
   REAL(8),  parameter                :: ce3plus     =  1.0D0
   REAL(8),  parameter                :: sig_k       =  1.0D0
   REAL(8),  parameter                :: sig_e       =  1.3D0
   REAL(8),  parameter                :: anLimitFact =  0.5D0
   REAL(8),  parameter                :: asLimitFact =  1.0D0
   REAL(8),  parameter                :: small       =  1.0D-10
!
   REAL(8),  parameter                :: cc1       =  5.0000D0
   REAL(8),  parameter                :: cc2       =  0.8000D0
   REAL(8),  parameter                :: cc3       =  1.9680D0
   REAL(8),  parameter                :: cc4       =  1.1360D0
   REAL(8),  parameter                :: cc5       =  0.0000D0
   REAL(8),  parameter                :: cc6       =  0.4000D0
   REAL(8),  parameter                :: ct1       =  5.9500D0
   REAL(8),  parameter                :: ct2       =  0.6000D0
   REAL(8),  parameter                :: ct3       =  1.0000D0
   REAL(8),  parameter                :: ct4       =  0.0000D0
   REAL(8),  parameter                :: ct5       =  0.3333D0
   REAL(8),  parameter                :: ctt       =  0.7200D0
   REAL(8),  parameter                :: a1        =  2./3. - cc2/2.
   REAL(8),  parameter                :: a2        =  1.    - cc3/2.
   REAL(8),  parameter                :: a3        =  1.    - cc4/2.
   REAL(8),  parameter                :: a4        =          cc5/2.
   REAL(8),  parameter                :: a5        =  1./2. - cc6/2.
   REAL(8),  parameter                :: at1       =           1. - ct2
   REAL(8),  parameter                :: at2       =           1. - ct3
   REAL(8),  parameter                :: at3       =  2. *   ( 1. - ct4)
   REAL(8),  parameter                :: at4       =  2. *   ( 1. - ct5)
   REAL(8),  parameter                :: at5       =  2.*ctt*( 1. - ct5)
   REAL(8),  parameter                :: xN        =   0.5*cc1
   REAL(8),  parameter                :: xNt       =   ct1
   REAL(8),  parameter                :: d0        =   36.* xN**3. * xNt**2.
   REAL(8),  parameter                :: d1        =   84.*a5*at3 * xN**2. * xNt  + 36.*at5 * xN**3. * xNt
   REAL(8),  parameter                :: d2        =   9.*(at2**2.-at1**2.) * xN**3. - 12.*(a2**2.-3.*a3**2.) * xN * xNt**2.
   REAL(8),  parameter                :: d3        =   12.*a5*at3*(a2*at1-3.*a3*at2) * xN + 12.*a5*at3*(a3**2.-a2**2.) * xNt       &
                                                        + 12.*at5*(3.*a3**2.-a2**2.) * xN * xNt
   REAL(8),  parameter                :: d4        =   48.*a5**2.*at3**2. * xN + 36.*a5*at3*at5 * xN**2.
   REAL(8),  parameter                :: d5        =   3.*(a2**2.-3.*a3**2.)*(at1**2.-at2**2.) * xN
   REAL(8),  parameter                :: n0        =   36.*a1 * xN**2. * xNt**2.
   REAL(8),  parameter                :: n1        = - 12.*a5*at3*(at1+at2) * xN**2.              &
                                                     + 8.*a5*at3*(6.*a1-a2-3.*a3) * xN * xNt      &
                                                     + 36.*a1*at5 * xN**2. * xNt
   REAL(8),  parameter                :: n2        =   9.*a1*(at2**2.-at1**2.) * xN**2.
   REAL(8),  parameter                :: nt0       =   12.*at3 * xN**3. * xNt
   REAL(8),  parameter                :: nt1       =   12.*a5*at3**2.  * xN**2.
   REAL(8),  parameter                :: nt2       =   9.*a1*at3*(at1-at2) * xN**2. + (  6.*a1*(a2-3.*a3)                         &
                                                          - 4.*(a2**2.-3.*a3**2.) )*at3 * xN * xNt
   REAL(8),  parameter                :: cm0       =  ( (a2*a2 - 3.0*a3*a3 + 3.0*a1*xN)/(3.0*xN*xN) )**0.25
   REAL(8),  parameter                :: cm3       = cm0**3
   REAL(8),  parameter                :: cm3_inv   = 1./cm3
   REAL(8),  parameter                :: cm0inv2   = 1./cm0**2
   REAL(8),  parameter                :: anMinNum  = -(d1 + nt0) + sqrt((d1+nt0)**2. - 4.*d0*(d4+nt1))
   REAL(8),  parameter                :: anMinDen  = 2.*(d4+nt1)
   REAL(8),  parameter                :: anMin     = anMinNum / anMinDen
   !=======================================================================
   !
END MODULE scm_par

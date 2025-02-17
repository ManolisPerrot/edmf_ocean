MODULE scm_par
   IMPLICIT NONE
   PUBLIC
   REAL(8), PUBLIC, PARAMETER ::   grav     = 9.81           !! Gravity of Earth
   REAL(8), PUBLIC, PARAMETER ::   vkarmn   = 0.41           !! Von Karman constant
   REAL(8), PUBLIC, PARAMETER ::   rho0     = 1027.          !! Boussinesq reference density [kg/m3]
   REAL(8), PUBLIC, PARAMETER ::   cp       = 3985.0d0       !! Specific heat capacity of saltwater [J/kg K]
   REAL(8), PUBLIC, PARAMETER ::   rpi      = 4.*ATAN(1.)    !! \( \pi \)
   REAL(8), PUBLIC, PARAMETER ::   wpmin    = 1.e-08         !! Minimum value of \( w^{\rm p} \) [m/s] <br />
   !
   !---------------------------------------------------------------------
   ! SCM parameters for TKE turbulent closure <br />
   !---------------------------------------------------------------------
   !REAL(8), PUBLIC, PARAMETER ::   tke_min   = 1.e-8          !! minimum TKE value [m2/s2]
   REAL(8), PUBLIC, PARAMETER ::   tke_min0  = 1.e-4          !! surface minimum value of tke for Dirichlet condition [m2/s2]
   !REAL(8), PUBLIC, PARAMETER ::   avm_bak   = 1.e-4          !! background eddy-viscosity [m2/s]
   !REAL(8), PUBLIC, PARAMETER ::   avt_bak   = 1.e-5          !! background eddy-diffusivity [m2/s]
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
   REAL(8), PARAMETER ::  c1 = 5.
   REAL(8), PARAMETER ::  c2 = 0.8
   REAL(8), PARAMETER ::  c3 = 1.968
   REAL(8), PARAMETER ::  c4 = 1.136
   REAL(8), PARAMETER ::  c5 = 0.
   REAL(8), PARAMETER ::  c6 = 0.4
   REAL(8), PARAMETER ::  cb1 = 5.95
   REAL(8), PARAMETER ::  cb2 = 0.6
   REAL(8), PARAMETER ::  cb3 = 1.
   REAL(8), PARAMETER ::  cb4 = 0.
   REAL(8), PARAMETER ::  cb5 = 0.33333
   REAL(8), PARAMETER ::  cbb = 0.72
   REAL(8), PARAMETER ::  a1  = 0.66666666667 - 0.5*c2
   REAL(8), PARAMETER ::  a2  = 1.            - 0.5*c3
   REAL(8), PARAMETER ::  a3  = 1.            - 0.5*c4
   REAL(8), PARAMETER ::  a5  = 0.5           - 0.5*c6
   REAL(8), PARAMETER ::  nn  = 0.5*c1  ! N
   REAL(8), PARAMETER ::  nb  = cb1     ! Nt
   REAL(8), PARAMETER ::  ab1 = 1. - cb2
   REAL(8), PARAMETER ::  ab2 = 1. - cb3
   REAL(8), PARAMETER ::  ab3 = 2.*(1.-cb4)
   REAL(8), PARAMETER ::  ab5 = 2.*cbb*(1.-cb5)
   REAL(8), PARAMETER ::  sf_d0 = 36.0*nn*nn*nn*nb*nb  ! d0
   REAL(8), PARAMETER ::  sf_d1 = 84.0*a5*ab3*nn*nn*nb+36.0*ab5*nn*nn*nn*nb ! d1
   REAL(8), PARAMETER ::  sf_d2 = 9.0*(ab2*ab2-ab1*ab1)*nn*nn*nn           &
                                       - 12.0*(a2*a2-3.*a3*a3)*nn*nb*nb
   REAL(8), PARAMETER ::  sf_d3 = 12.0*a5*ab3*(a2*ab1-3.0*a3*ab2)* nn      &
                                      + 12.0*a5*ab3*(    a3*a3-a2*a2)* nb     &
                                      + 12.0*   ab5*(3.0*a3*a3-a2*a2)*nn*nb
   REAL(8), PARAMETER ::  sf_d4 = 48.0*a5*a5*ab3*ab3*nn + 36.0*a5*ab3*ab5*nn*nn
   REAL(8), PARAMETER ::  sf_d5 = 3.0*(a2*a2-3.0*a3*a3)*(ab1*ab1-ab2*ab2)*nn
   REAL(8), PARAMETER ::  sf_n0 = 36.0*a1*nn*nn*nb*nb
   REAL(8), PARAMETER ::  sf_n1 = - 12.0*a5*ab3*(ab1+ab2)*nn*nn                         &
             + 8.0*a5*ab3*(6.0*a1-a2-3.0*a3)*nn*nb + 36.0*a1*ab5*nn*nn*Nb
   REAL(8), PARAMETER ::  sf_n2 = 9.0*a1*(ab2*ab2-ab1*ab1)*nn*nn
   REAL(8), PARAMETER ::  sf_nb0 = 12.0*ab3*nn*nn*nn*nb   ! nt0
   REAL(8), PARAMETER ::  sf_nb1 = 12.0*a5*ab3*ab3*nn*nn  ! nt1
   REAL(8), PARAMETER ::  sf_nb2 = 9.0*a1*ab3*(ab1-ab2)*nn*nn + ( 6.0*a1*(a2-3.0*a3)     &
                                  - 4.0*(a2*a2-3.0*a3*a3) )*ab3 * nn * nb  ! nt2
   REAL(8), PARAMETER ::  lim_am0 = sf_d0*sf_n0
   REAL(8), PARAMETER ::  lim_am1 = sf_d0*sf_n1 + sf_d1*sf_n0
   REAL(8), PARAMETER ::  lim_am2 = sf_d1*sf_n1 + sf_d4*sf_n0
   REAL(8), PARAMETER ::  lim_am3 = sf_d4*sf_n1
   REAL(8), PARAMETER ::  lim_am4 = sf_d2*sf_n0
   REAL(8), PARAMETER ::  lim_am5 = sf_d2*sf_n1+sf_d3*sf_n0
   REAL(8), PARAMETER ::  lim_am6 = sf_d3*sf_n1
   !
END MODULE scm_par

MODULE scm_par
   IMPLICIT NONE
   PUBLIC
   REAL(8), PUBLIC, PARAMETER ::   grav    = 9.81           !! Gravity of Earth
   REAL(8), PUBLIC, PARAMETER ::   vkarmn  = 0.41           !! Von Karman constant
   REAL(8), PUBLIC, PARAMETER ::   rho0    = 1027.          !! Boussinesq reference density [kg/m3]
   REAL(8), PUBLIC, PARAMETER ::   cp      = 3985.0d0       !! Specific heat capacity of saltwater [J/kg K]
   REAL(8), PUBLIC, PARAMETER ::   rpi     = 4.*ATAN(1.)    !! \( \pi \)
   !REAL(8), PUBLIC, PARAMETER ::   wpmin   = 1.e-08         !! Minimum value of \( w^{\rm p} \) [m/s] <br />
   !
   !---------------------------------------------------------------------
   ! SCM parameters for TKE turbulent closure <br />
   !---------------------------------------------------------------------
   REAL(8), PUBLIC, PARAMETER ::   tke_min   = 1.e-8          !! minimum TKE value [m2/s2]
   REAL(8), PUBLIC, PARAMETER ::   tke_min0  = 1.e-4          !! surface minimum value of tke for Dirichlet condition [m2/s2]
   REAL(8), PUBLIC, PARAMETER ::   avm_bak   = 1.e-4          !! background eddy-viscosity [m2/s]
   REAL(8), PUBLIC, PARAMETER ::   avt_bak   = 1.e-5          !! background eddy-diffusivity [m2/s]
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
   REAL(8), PUBLIC, PARAMETER ::   mxl_min = (avm_bak / 0.1) / sqrt( tke_min )  !! minimum value for mixing lengths [m]
   REAL(8), PUBLIC, PARAMETER ::   mxl_min0 = 0.04    !! minimum surface value for miwing lengths [m]
   REAL(8), PUBLIC, PARAMETER ::   pdlrmin  = 0.1     !! minimum value for the inverse Prandtl number
   !
   REAL(8), PUBLIC, PARAMETER ::  bshear    =   1.e-20 !! minimum shear
   REAL(8), PUBLIC, PARAMETER ::  rsmall    =   1.e-20
   !
END MODULE scm_par

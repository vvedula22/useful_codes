!#######################################################################
!     Parameters used for Aliev-Panfilov Ventricular Myocyte Model
!#######################################################################
!     Scaling factors
!     Voltage scaling
      REAL(KIND=8) :: Vscale  = 100.0D0 ! 100.0
!     Time scaling
      REAL(KIND=8) :: Tscale  = 12.9D0  ! 12.9
!     Voltage offset parameter
      REAL(KIND=8) :: Voffset = -80.0D0 ! -80.0
!-----------------------------------------------------------------------
!     Model parameters
      REAL(KIND=8) :: alpha = 0.05D0
      REAL(KIND=8) :: a     = 0.002D0
      REAL(KIND=8) :: b     = 0.15D0
      REAL(KIND=8) :: c     = 8.0D0
      REAL(KIND=8) :: mu1   = 0.2D0
      REAL(KIND=8) :: mu2   = 0.3D0
!-----------------------------------------------------------------------
!     Activation coupling parameters
      REAL(KIND=8) :: Xrest = -80.0D0
      REAL(KIND=8) :: Xcrit = -30.0D0
      REAL(KIND=8) :: K_T   = 5D-3
      REAL(KIND=8) :: eps0  = 0.1D0
      REAL(KIND=8) :: eps1  = 1.0D0
      REAL(KIND=8) :: xi_T  = 1.0D0
!#######################################################################

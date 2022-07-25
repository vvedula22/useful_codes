!--------------------------------------------------------------------
!
!     Parameters used for Aliev-Panfilov Ventricular Myocyte Model.
!
!--------------------------------------------------------------------

!     Scaling factors
!     Voltage scaling (mV)
      REAL(KIND=RKIND) :: Vscale  = 100._RKIND
!     Time scaling (ms)
      REAL(KIND=RKIND) :: Tscale  = 12.9_RKIND
!     Voltage offset parameter (mV)
      REAL(KIND=RKIND) :: Voffset = -80._RKIND
!-----------------------------------------------------------------------
!     Model parameters
      REAL(KIND=RKIND) :: alpha = 1.E-2_RKIND
      REAL(KIND=RKIND) :: a     = 2.E-3_RKIND
      REAL(KIND=RKIND) :: b     = 0.15_RKIND
      REAL(KIND=RKIND) :: c     = 8._RKIND
      REAL(KIND=RKIND) :: mu1   = 0.2_RKIND
      REAL(KIND=RKIND) :: mu2   = 0.3_RKIND
!-----------------------------------------------------------------------
!     Electromechanics coupling parameters: active stress model
!     Resting voltage (mV)
      REAL(KIND=RKIND) :: Vrest = -80._RKIND
!     Critical voltage (mV)
      REAL(KIND=RKIND) :: Vcrit = -30._RKIND
!     Saturation potential (Pa/mV)
      REAL(KIND=RKIND) :: K_T = 5.0E3_RKIND
!     Minimum activation (ms^{-1})
      REAL(KIND=RKIND) :: eps_0 = 0.1_RKIND
!     Maximum activation (ms^{-1})
      REAL(KIND=RKIND) :: eps_i = 1._RKIND
!     Transition rate (mV^{-1})
      REAL(KIND=RKIND) :: xi_T  = 1._RKIND

!     Cm: Cell capacitance per unit surface area
      REAL(KIND=RKIND) :: Cm  = 1._RKIND
!     sV: Surface to volume ratio
      REAL(KIND=RKIND) :: sV  = 1._RKIND
!     rho: Cellular resistivity
      REAL(KIND=RKIND) :: rho = 1._RKIND
!#######################################################################

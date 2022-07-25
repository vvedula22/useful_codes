!#######################################################################
!     Parameters used for Bueno-Orovio Ventricular Myocyte Model (endo)
!#######################################################################
!     Scaling factors
!     Voltage scaling
      REAL(KIND=8) :: Vscale  = 85.7D0
!     Time scaling
      REAL(KIND=8) :: Tscale  = 1.0D0
!     Voltage offset parameter
      REAL(KIND=8) :: Voffset = -84.0D0
!-----------------------------------------------------------------------
!     Model parameters (epi, endo, myo)
      REAL(KIND=8) :: u_o(3)      = (/0.0D0,     0.0D0,    0.0D0   /)
      REAL(KIND=8) :: u_u(3)      = (/1.55D0,    1.56D0,   1.61D0  /)
      REAL(KIND=8) :: theta_v(3)  = (/0.3D0,     0.3D0,    0.3D0   /)
      REAL(KIND=8) :: theta_w(3)  = (/0.13D0,    0.13D0,   0.13D0  /)
      REAL(KIND=8) :: thetam_v(3) = (/6.0D-3,    0.2D0,    0.1D0   /)
      REAL(KIND=8) :: theta_o(3)  = (/6.0D-3,    6.0D-3,   5.0D-3  /)
      REAL(KIND=8) :: taum_v1(3)  = (/60.0D0,    75.0D0,   80.0D0  /)
      REAL(KIND=8) :: taum_v2(3)  = (/1.15D3,    10.0D0,   1.4506D0/)
      REAL(KIND=8) :: taup_v(3)   = (/1.4506D0,  1.4506D0, 1.4506D0/)
      REAL(KIND=8) :: taum_w1(3)  = (/60.0D0,    6.0D0,    70.0D0  /)
      REAL(KIND=8) :: taum_w2(3)  = (/15.0D0,    140.0D0,  8.0D0   /)
      REAL(KIND=8) :: km_w(3)     = (/65.0D0,    200.0D0,  200.0D0 /)
      REAL(KIND=8) :: um_w(3)     = (/0.03D0,    0.016D0,  0.016D0 /)
      REAL(KIND=8) :: taup_w(3)   = (/200.0D0,   280.0D0,  280.0D0 /)
      REAL(KIND=8) :: tau_fi(3)   = (/0.11D0,    0.1D0,    0.078D0 /)
      REAL(KIND=8) :: tau_o1(3)   = (/400.0D0,   470.0D0,  410.0D0 /)
      REAL(KIND=8) :: tau_o2(3)   = (/6.0D0,     6.0D0,    7.0D0   /)
      REAL(KIND=8) :: tau_so1(3)  = (/30.0181D0, 40.0D0,   91.0D0  /)
      REAL(KIND=8) :: tau_so2(3)  = (/0.9957D0,  1.2D0,    0.8D0   /)
      REAL(KIND=8) :: k_so(3)     = (/2.0458D0,  2.0D0,    2.1D0   /)
      REAL(KIND=8) :: u_so(3)     = (/0.65D0,    0.65D0,   0.6D0   /)
      REAL(KIND=8) :: tau_s1(3)   = (/2.7342D0,  2.7342D0, 2.7342D0/)
      REAL(KIND=8) :: tau_s2(3)   = (/16.0D0,    2.0D0,    2.0D0   /)
      REAL(KIND=8) :: k_s(3)      = (/2.0994D0,  2.0994D0, 2.0994D0/)
      REAL(KIND=8) :: u_s(3)      = (/0.9087D0,  0.9087D0, 0.9087D0/)
      REAL(KIND=8) :: tau_si(3)   = (/1.8875D0,  2.9013D0, 3.3849D0/)
      REAL(KIND=8) :: tau_winf(3) = (/0.07D0,    0.0273D0, 0.01D0  /)
      REAL(KIND=8) :: ws_inf(3)   = (/0.94D0,    0.78D0,   0.5D0   /)
!-----------------------------------------------------------------------
!     Activation coupling parameters
      REAL(KIND=8) :: Xrest = -84.0D0
      REAL(KIND=8) :: Xcrit = -30.0D0
      REAL(KIND=8) :: K_T   = 5D-3
      REAL(KIND=8) :: eps0  = 0.1D0
      REAL(KIND=8) :: eps1  = 1.0D0
      REAL(KIND=8) :: xi_T  = 1.0D0
!#######################################################################

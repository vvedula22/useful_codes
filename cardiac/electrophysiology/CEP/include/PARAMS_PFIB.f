!-----------------------------------------------------------------------
!
!     This module defines parameters for the Stewart's cellular
!     activation model for Purkinje fiber cells.
!
!     Reference for Stewart's electrophysiology model:
!        Stewart, P., et al., (2009). Mathematical models of the
!        electrical action potential of Purkinje fibre cells.
!        Philosophical Transactions of the Royal Society A:
!        Mathematical, Physical and Engineering Sciences, 367(1896),
!        pp.2225-2255.
!        https://doi.org/10.1098/rsta.2008.0283
!
!-----------------------------------------------------------------------

!     Default model parameters
!     R: Gas constant
      REAL(KIND=RKIND) :: Rc = 8314.472_RKIND      ! units: J/mol/K
!     T: Temperature
      REAL(KIND=RKIND) :: Tc = 310._RKIND          ! units: K
!     F: Faraday constant
!     Note: F=96.4867 C/mmol in paper, but C=96485.3415 C/mmol in CellML
      REAL(KIND=RKIND) :: Fc = 96485.3415_RKIND    ! units: C/mmol
!     Cm: Cell capacitance per unit surface area
!     Note: Cm = 2.0 uF/cm2 in paper, but C = 0.185 uF in CellML
      REAL(KIND=RKIND) :: Cm = 0.185_RKIND         ! units: uF/cm^{2}
!     V_c: Cytoplasmic volume
!     Note: V_c = 16.404 um3 in the paper
      REAL(KIND=RKIND) :: V_c = 0.016404_RKIND    ! units: um^{3}
!     V_sr: Sacroplasmic reticulum volume (value chosen from CellML)
!     Note: V_sr = 1.094 um3 in the paper
      REAL(KIND=RKIND) :: V_sr = 0.001094_RKIND    ! units: um^{3}
!     V_ss: Subspace volume (value from CellML)
!     Note: V_ss = 0.05468 um3 in the paper
      REAL(KIND=RKIND) :: V_ss = 5.468E-5_RKIND    ! units: um^{3}
!     K_o: Extracellular K concentration
      REAL(KIND=RKIND) :: K_o = 5.4_RKIND          ! units: mM
!     Na_o: Extracellular Na concentration
      REAL(KIND=RKIND) :: Na_o = 140._RKIND        ! units: mM
!     Ca_o: Extracellular Ca concentration
      REAL(KIND=RKIND) :: Ca_o = 2._RKIND          ! units: mM
!     G_Na: Maximal I_Na conductance
      REAL(KIND=RKIND) :: G_Na = 130.5744_RKIND    ! units: nS/pF
!     G_K1: Maximal I_K1 conductance
      REAL(KIND=RKIND) :: G_K1 = 0.065_RKIND       ! units: nS/pF
!     G_to: Maximal I_to conductance
      REAL(KIND=RKIND) :: G_to = 0.08184_RKIND     ! units: nS/pF
!     G_sus: Maximal I_sus conductance
      REAL(KIND=RKIND) :: G_sus = 0.0227_RKIND     ! units: nS/pF
!     G_fK: Maximal hyperpolarization current I_f conductance
      REAL(KIND=RKIND) :: G_fK = 0.0234346_RKIND     ! units: nS/pF
!     G_fNa: Maximal hyperpolarization current I_f conductance
      REAL(KIND=RKIND) :: G_fNa = 0.0145654_RKIND   ! units: nS/pF
!     G_Kr: Maximal I_Kr conductance
      REAL(KIND=RKIND) :: G_Kr = 0.0918_RKIND      ! units: nS/pF
!     G_Ks: Maximal epicardial I_Ks conductance, units: nS/pF
      REAL(KIND=RKIND) :: G_Ks = 0.2352_RKIND      ! units: nS/pF
!     p_KNa: Relative I_Ks permeability to Na
      REAL(KIND=RKIND) :: p_KNa = 0.03_RKIND       ! dimensionless
!     G_CaL: Maximal I_CaL conductance
      REAL(KIND=RKIND) :: G_CaL = 3.98E-5_RKIND    ! units: cm^{3}/uF/ms
!     K_NaCa: Maximal I_NaCa
      REAL(KIND=RKIND) :: K_NaCa = 1000._RKIND     ! units: pA/pF
!     gamma: Voltage dependent parameter of I_NaCa
      REAL(KIND=RKIND) :: gamma = 0.35_RKIND       ! dimensionless
!     K_mCa: Ca_i half-saturation constant for I_NaCa
      REAL(KIND=RKIND) :: K_mCa = 1.38_RKIND       ! units: mM
!     K_mNai: Na_i half-saturation constant for I_NaCa
      REAL(KIND=RKIND) :: K_mNai = 87.5_RKIND      ! units: mM
!     K_sat: Saturation factor for I_NaCa
      REAL(KIND=RKIND) :: K_sat = 0.1_RKIND        ! dimensionless
!     alpha: Factor enhancing outward nature of I_NaCa
      REAL(KIND=RKIND) :: alpha = 2.5_RKIND        ! dimensionless
!     p_NaK: Maximal I_NaK
      REAL(KIND=RKIND) :: p_NaK = 2.724_RKIND      ! units: pA/pF
!     K_mK: K_o half-saturation constant of I_NaK
      REAL(KIND=RKIND) :: K_mK = 1._RKIND          ! units: mM
!     K_mNa: Na_i half-saturation constant of I_NaK
      REAL(KIND=RKIND) :: K_mNa = 40._RKIND        ! units: mM
!     G_pK: Maximal I_pK conductance
      REAL(KIND=RKIND) :: G_pK = 0.0146_RKIND      ! units: nS/pF
!     G_pCa: Maximal I_pCa conductance
      REAL(KIND=RKIND) :: G_pCa = 0.1238_RKIND     ! units: pA/pF
!     K_pCa: Half-saturation constant of I_pCa
      REAL(KIND=RKIND) :: K_pCa = 5.E-4_RKIND      ! units: mM
!     G_bNa: Maximal I_bNa conductance
      REAL(KIND=RKIND) :: G_bNa = 2.9E-4_RKIND     ! units: nS/pF
!     G_bCa: Maximal I_bCa conductance
      REAL(KIND=RKIND) :: G_bCa = 5.92E-4_RKIND    ! units: nS/pF
!     Vmax_up: Maximal I_up conductance
      REAL(KIND=RKIND) :: Vmax_up = 0.006375_RKIND ! units: mM/ms
!     K_up: Half-saturation constant of I_up
      REAL(KIND=RKIND) :: K_up = 2.5E-4_RKIND      ! units: mM
!     V_rel: Maximal I_rel conductance
!     Note: V_rel=40.8 mM/ms in the paper, but V_rel=0.102 /ms in CellML
      REAL(KIND=RKIND) :: V_rel = 0.102_RKIND      ! units: mM/ms
!     k1p: R to O and RI to I, I_rel transition rate
      REAL(KIND=RKIND) :: k1p = 0.15_RKIND         ! units: mM^{-2}/ms
!     k2p: O to I and R to RI, I_rel transition rate
      REAL(KIND=RKIND) :: k2p = 0.045_RKIND        ! units: mM^{-1}/ms
!     k3: O to R and I to RI, I_rel transition rate
      REAL(KIND=RKIND) :: k3 = 0.06_RKIND          ! units: ms^{-1}
!     k4: I to O and Ri to I, I_rel transition rate
!         Note k4=1.5e-5 in paper, but k4=0.005 in CellML
      REAL(KIND=RKIND) :: k4 = 0.005_RKIND         ! units: ms^{-1}
!     EC: Ca_sr half-saturation constant of k_casr
      REAL(KIND=RKIND) :: EC = 1.5_RKIND           ! units: mM
!     max_sr: Maximum value of k_casr
      REAL(KIND=RKIND) :: max_sr = 2.5_RKIND       ! dimensionless
!     min_sr: Minimum value of k_casr
      REAL(KIND=RKIND) :: min_sr = 1._RKIND        ! dimensionless
!     V_leak: Maximal I_leak conductance
      REAL(KIND=RKIND) :: V_leak = 3.6E-4_RKIND    ! units: mM/ms
!     V_xfer: Maximal I_xfer conductance
      REAL(KIND=RKIND) :: V_xfer = 0.0038_RKIND    ! units: mM/ms
!     Buf_c: Total cytoplasmic buffer concentration
      REAL(KIND=RKIND) :: Buf_c = 0.2_RKIND        ! units: mM
!     K_bufc: Ca_i half-saturation constant for cytplasmic buffer
      REAL(KIND=RKIND) :: K_bufc = 0.001_RKIND     ! units: mM
!     Buf_sr: Total sacroplasmic buffer concentration
      REAL(KIND=RKIND) :: Buf_sr = 10._RKIND       ! units: mM
!     K_bufsr: Ca_sr half-saturation constant for subspace buffer
      REAL(KIND=RKIND) :: K_bufsr = 0.3_RKIND      ! units: mM
!     Buf_ss: Total subspace buffer concentration
      REAL(KIND=RKIND) :: Buf_ss = 0.4_RKIND       ! units: mM
!     K_bufss: Ca_ss half-saturation constant for subspace buffer
      REAL(KIND=RKIND) :: K_bufss = 2.5E-4_RKIND   ! units: mM
!     Resting potential
      REAL(KIND=RKIND) :: Vrest = -85.23_RKIND     ! units: mV
!-----------------------------------------------------------------------
!     Scaling factors
!     Voltage scaling
      REAL(KIND=RKIND) :: Vscale  = 1._RKIND
!     Time scaling
      REAL(KIND=RKIND) :: Tscale  = 1._RKIND
!     Voltage offset parameter
      REAL(KIND=RKIND) :: Voffset = 0._RKIND
!-----------------------------------------------------------------------
!     Variables
!     Reverse potentials for Na, K, Ca
      REAL(KIND=RKIND) :: E_Na
      REAL(KIND=RKIND) :: E_K
      REAL(KIND=RKIND) :: E_Ca
      REAL(KIND=RKIND) :: E_Ks

!     Cellular transmembrane currents
!     I_Na: Fast sodium current
      REAL(KIND=RKIND) :: I_Na
!     I_K1: inward rectifier outward current
      REAL(KIND=RKIND) :: I_K1
!     I_to: transient outward current
      REAL(KIND=RKIND) :: I_to
!     I_Kr: rapid delayed rectifier current
      REAL(KIND=RKIND) :: I_Kr
!     I_Ks: slow delayed rectifier current
      REAL(KIND=RKIND) :: I_Ks
!     I_CaL: L-type Ca current
      REAL(KIND=RKIND) :: I_CaL
!     I_NaCa: Na-Ca exchanger current
      REAL(KIND=RKIND) :: I_NaCa
!     I_NaK: Na-K pump current
      REAL(KIND=RKIND) :: I_NaK
!     I_pCa: plateau Ca current
      REAL(KIND=RKIND) :: I_pCa
!     I_pK: plateau K current
      REAL(KIND=RKIND) :: I_pK
!     I_bCa: background Ca current
      REAL(KIND=RKIND) :: I_bCa
!     I_lean: background Na current
      REAL(KIND=RKIND) :: I_bNa
!     I_leak: sacroplasmic reticulum Ca leak current
      REAL(KIND=RKIND) :: I_leak
!     I_up: sacroplasmic reticulum Ca pump current
      REAL(KIND=RKIND) :: I_up
!     I_rel: Ca induced Ca release current
      REAL(KIND=RKIND) :: I_rel
!     I_xfer: diffusive Ca current
      REAL(KIND=RKIND) :: I_xfer
!     I_sus: sustained current
      REAL(KIND=RKIND) :: I_sus
!     I_f: hyperpolarization-activated currents
      REAL(KIND=RKIND) :: I_f, I_fNa, I_fK
!-----------------------------------------------------------------------
!     State variables
      REAL(KIND=RKIND) :: V
      REAL(KIND=RKIND) :: K_i
      REAL(KIND=RKIND) :: Na_i
      REAL(KIND=RKIND) :: Ca_i
      REAL(KIND=RKIND) :: Ca_ss
      REAL(KIND=RKIND) :: Ca_sr
      REAL(KIND=RKIND) :: R_bar

!     Gating variables (runtime, steady state)
      REAL(KIND=RKIND) :: xr1, xr1i
      REAL(KIND=RKIND) :: xr2, xr2i
      REAL(KIND=RKIND) :: xs, xsi
      REAL(KIND=RKIND) :: m, mi
      REAL(KIND=RKIND) :: h, hi
      REAL(KIND=RKIND) :: j, ji
      REAL(KIND=RKIND) :: d, di
      REAL(KIND=RKIND) :: f, fi
      REAL(KIND=RKIND) :: f2, f2i
      REAL(KIND=RKIND) :: fCass, fCassi
      REAL(KIND=RKIND) :: s, si
      REAL(KIND=RKIND) :: r, ri
      REAL(KIND=RKIND) :: y, yi

!     Other variables
      REAL(KIND=RKIND) :: k1
      REAL(KIND=RKIND) :: k2
      REAL(KIND=RKIND) :: k_casr
      REAL(KIND=RKIND) :: O
!#######################################################################

!--------------------------------------------------------------------
!
!     Parameters for the Courtemanche-Ramirez-Nattel (CRN) model of
!     the human atrial myocyte.
!
!     Reference:
!        Courtemanche M, Ramirez RJ, Nattel S (1998). Ionic mechanisms
!        underlying human atrial action potential properties: insights
!        from a mathematical model. Am J Physiol Heart Circ Physiol,
!        275(1), H301-H321.
!        https://doi.org/10.1152/ajpheart.1998.275.1.H301
!
!  ------------------------- UNIT SYSTEM -------------------------
!
!     Units are native to the paper.  The model is kept internally
!     consistent with the original formulation, so these values are
!     not rescaled:
!
!        time           ms
!        voltage        mV
!        current        ABSOLUTE pA (the paper's pA/pF densities
!                       multiplied by Cm = 100 pF at the point of use)
!        conductance    nS/pF in this file; x Cm at use -> nS
!        concentration  mM (= mmol/L)
!        volume         um^3
!        SR fluxes      mM/ms (I_rel, I_tr, I_up, I_up_leak; these are
!                       fluxes, NOT electrical currents)
!
!     Two dimensional identities make the system close:
!        pA / pF          = mV/ms  -> dV/dt = -(SUM I + Istim)/Cm is
!                                    in mV/ms, so the clock is ms
!        pA / (C/mmol * um^3) = mM/ms -> I/(F*Vol) needs no extra
!                                    conversion factor anywhere
!
!     CAUTION: the CEP models in this solver do not share a unit
!     system.  They differ in the clock, in whether currents are
!     absolute or densities, and in the role Cm plays in dV/dt.  Do
!     not copy stimulus blocks, time steps or parameter values
!     between input decks.  The cross-model table is kept in one
!     place: example/README_CEP_UNITS.md
!
!     NOTE (Eq. 68 / Fn): the currents entering Fn must be ABSOLUTE
!     pA, not pA/pF.  Eqs. 65/66 compare Fn against 3.4175E-13, so a
!     pA/pF implementation leaves the membrane term 100x low and SR
!     release never fires.  This is why CRN is not harmonized to the
!     TTP-style densities used elsewhere in the solver.
!
!--------------------------------------------------------------------
!     Default model parameters
!     R: Universal gas constant
      REAL(KIND=RKIND) :: R = 8.3143_RKIND          ! units: J/(K*mol)
!     T: Absolute temperature
      REAL(KIND=RKIND) :: T = 310._RKIND            ! units: K
!     F: Faraday constant
      REAL(KIND=RKIND) :: F = 96.4867_RKIND         ! units: C/mmol
!     Cm: Membrane capacitance
      REAL(KIND=RKIND) :: Cm = 100._RKIND           ! units: pF
!     Vol_cell: Cell volume (not used by the equations; Table 1)
      REAL(KIND=RKIND) :: Vol_cell = 20100._RKIND   ! units: um^3
!     Vol_i: Intracellular (myoplasm) volume
      REAL(KIND=RKIND) :: Vol_i = 13668._RKIND      ! units: um^3
!     Vol_up: SR uptake compartment (NSR) volume
      REAL(KIND=RKIND) :: Vol_up = 1109.52_RKIND    ! units: um^3
!     Vol_rel: SR release compartment (JSR) volume
      REAL(KIND=RKIND) :: Vol_rel = 96.48_RKIND     ! units: um^3
!     K_o: Extracellular K concentration
      REAL(KIND=RKIND) :: K_o = 5.4_RKIND           ! units: mM
!     Na_o: Extracellular Na concentration
      REAL(KIND=RKIND) :: Na_o = 140._RKIND         ! units: mM
!     Ca_o: Extracellular Ca concentration
      REAL(KIND=RKIND) :: Ca_o = 1.8_RKIND          ! units: mM
!-----------------------------------------------------------------------
!     Maximal conductances (Table 1, nS/pF; multiplied by Cm at use)
!     g_Na: I_Na
      REAL(KIND=RKIND) :: g_Na = 7.8_RKIND          ! units: nS/pF
!     g_K1: I_K1
      REAL(KIND=RKIND) :: g_K1 = 0.09_RKIND         ! units: nS/pF
!     g_to: I_to
      REAL(KIND=RKIND) :: g_to = 0.1652_RKIND       ! units: nS/pF
!     g_Kr: I_Kr
      REAL(KIND=RKIND) :: g_Kr = 0.0294_RKIND       ! units: nS/pF
!     g_Ks: I_Ks
      REAL(KIND=RKIND) :: g_Ks = 0.129_RKIND        ! units: nS/pF
!     g_CaL: I_Ca,L
      REAL(KIND=RKIND) :: g_CaL = 0.1238_RKIND      ! units: nS/pF
!     g_bCa: I_b,Ca
      REAL(KIND=RKIND) :: g_bCa = 0.00113_RKIND     ! units: nS/pF
!     g_bNa: I_b,Na
      REAL(KIND=RKIND) :: g_bNa = 0.000674_RKIND    ! units: nS/pF
!-----------------------------------------------------------------------
!     Pumps and exchangers (Table 1, pA/pF; multiplied by Cm at use)
!     I_NaK_max: Maximal I_NaK
      REAL(KIND=RKIND) :: I_NaK_max = 0.60_RKIND    ! units: pA/pF
!     I_NaCa_max: Maximal I_NaCa
      REAL(KIND=RKIND) :: I_NaCa_max = 1600._RKIND  ! units: pA/pF
!     I_pCa_max: Maximal I_p,Ca
      REAL(KIND=RKIND) :: I_pCa_max = 0.275_RKIND   ! units: pA/pF
!     Km_Na_i: [Na]i half-saturation constant for I_NaK
      REAL(KIND=RKIND) :: Km_Na_i = 10._RKIND       ! units: mM
!     Km_K_o: [K]o half-saturation constant for I_NaK
      REAL(KIND=RKIND) :: Km_K_o = 1.5_RKIND        ! units: mM
!     Km_Na: [Na]o half-saturation constant for I_NaCa
      REAL(KIND=RKIND) :: Km_Na = 87.5_RKIND        ! units: mM
!     Km_Ca: [Ca]o half-saturation constant for I_NaCa
      REAL(KIND=RKIND) :: Km_Ca = 1.38_RKIND        ! units: mM
!     k_sat: Saturation factor for I_NaCa
      REAL(KIND=RKIND) :: k_sat = 0.1_RKIND         ! dimensionless
!     gamma: Voltage dependence parameter for I_NaCa
      REAL(KIND=RKIND) :: gamma = 0.35_RKIND        ! dimensionless
!-----------------------------------------------------------------------
!     Sarcoplasmic reticulum
!     I_up_max: Maximal I_up
      REAL(KIND=RKIND) :: I_up_max = 0.005_RKIND    ! units: mM/ms
!     K_up: [Ca]i half-saturation constant for I_up
      REAL(KIND=RKIND) :: K_up = 0.00092_RKIND      ! units: mM
!     Ca_up_max: Maximal Ca concentration in uptake compartment
      REAL(KIND=RKIND) :: Ca_up_max = 15._RKIND     ! units: mM
!     k_rel: Maximal release rate for I_rel
      REAL(KIND=RKIND) :: k_rel = 30._RKIND         ! units: 1/ms
!     tau_tr: NSR -> JSR transfer time constant
      REAL(KIND=RKIND) :: tau_tr = 180._RKIND       ! units: ms
!-----------------------------------------------------------------------
!     Ca buffers (instantaneous equilibrium)
!     Cmdn_max: Total calmodulin concentration in myoplasm
      REAL(KIND=RKIND) :: Cmdn_max = 0.05_RKIND     ! units: mM
!     Trpn_max: Total troponin concentration in myoplasm
      REAL(KIND=RKIND) :: Trpn_max = 0.07_RKIND     ! units: mM
!     Csqn_max: Total calsequestrin concentration in JSR
      REAL(KIND=RKIND) :: Csqn_max = 10._RKIND      ! units: mM
!     Km_Cmdn: [Ca]i half-saturation constant for calmodulin
      REAL(KIND=RKIND) :: Km_Cmdn = 0.00238_RKIND   ! units: mM
!     Km_Trpn: [Ca]i half-saturation constant for troponin
      REAL(KIND=RKIND) :: Km_Trpn = 0.0005_RKIND    ! units: mM
!     Km_Csqn: [Ca]rel half-saturation constant for calsequestrin
      REAL(KIND=RKIND) :: Km_Csqn = 0.8_RKIND       ! units: mM
!-----------------------------------------------------------------------
!     Kinetics
!     K_Q10: Temperature scaling factor for I_Kur and I_to kinetics
      REAL(KIND=RKIND) :: K_Q10 = 3._RKIND          ! dimensionless
!-----------------------------------------------------------------------
!     Fn_rel_scale: scale applied to the RELEASE term of Fn (Eq. 68)
!     ONLY.  This is the single deliberate deviation from the
!     published equations anywhere in this module; it is needed to
!     reproduce the paper's own Fig. 15.
!
!        1.0   Eq. 68 exactly as published.
!        0.01  DEFAULT; reproduces Fig. 15.
      REAL(KIND=RKIND) :: Fn_rel_scale = 0.01_RKIND ! dimensionless
!-----------------------------------------------------------------------
!     Resting potential, used only by the stretch-activated current
!     I_sac = Ksac*(Vrest - V)
      REAL(KIND=RKIND) :: Vrest = -81.18_RKIND      ! units: mV
!-----------------------------------------------------------------------
!     Scaling factors (declared for parity with other models; unused)
!     Voltage scaling
      REAL(KIND=RKIND) :: Vscale = 1._RKIND
!     Time scaling
      REAL(KIND=RKIND) :: Tscale = 1._RKIND
!     Voltage offset parameter
      REAL(KIND=RKIND) :: Voffset = 0._RKIND
!-----------------------------------------------------------------------
!     Variables
!     Reversal potentials for Na, K, Ca
      REAL(KIND=RKIND) :: E_Na
      REAL(KIND=RKIND) :: E_K
      REAL(KIND=RKIND) :: E_Ca

!     Cellular transmembrane currents (absolute pA)
!     I_Na: Fast Na current
      REAL(KIND=RKIND) :: I_Na
!     I_K1: Inward rectifier K current
      REAL(KIND=RKIND) :: I_K1
!     I_to: Transient outward K current
      REAL(KIND=RKIND) :: I_to
!     I_Kur: Ultrarapid delayed rectifier K current
      REAL(KIND=RKIND) :: I_Kur
!     I_Kr: Rapid delayed rectifier K current
      REAL(KIND=RKIND) :: I_Kr
!     I_Ks: Slow delayed rectifier K current
      REAL(KIND=RKIND) :: I_Ks
!     I_CaL: L-type Ca current
      REAL(KIND=RKIND) :: I_CaL
!     I_pCa: Sarcolemmal Ca pump current
      REAL(KIND=RKIND) :: I_pCa
!     I_NaK: Na-K pump current
      REAL(KIND=RKIND) :: I_NaK
!     I_NaCa: Na-Ca exchanger current
      REAL(KIND=RKIND) :: I_NaCa
!     I_bNa: Background Na current
      REAL(KIND=RKIND) :: I_bNa
!     I_bCa: Background Ca current
      REAL(KIND=RKIND) :: I_bCa
!     I_sac: Stretch-activated current
      REAL(KIND=RKIND) :: I_sac

!     Sarcoplasmic reticulum Ca fluxes (mM/ms)
!     I_rel: Ca release from the JSR (CICR)
      REAL(KIND=RKIND) :: I_rel
!     I_tr: NSR -> JSR transfer
      REAL(KIND=RKIND) :: I_tr
!     I_up: Ca uptake into the NSR
      REAL(KIND=RKIND) :: I_up
!     I_up_leak: Ca leak from the NSR
      REAL(KIND=RKIND) :: I_up_leak

!     Fn: Ca release flux trigger (Eq. 68).  Computed in CRN_GETF and
!     consumed by CRN_UPDATEG for the u and v gates; CRN_GETF must
!     always be called before CRN_UPDATEG.
      REAL(KIND=RKIND) :: Fn

!     State variables (local copies)
      REAL(KIND=RKIND) :: V
      REAL(KIND=RKIND) :: Na_i
      REAL(KIND=RKIND) :: K_i
      REAL(KIND=RKIND) :: Ca_i
      REAL(KIND=RKIND) :: Ca_up
      REAL(KIND=RKIND) :: Ca_rel

!     Gating variables (local copies).  Note the g-prefix: the gate f
!     would otherwise collide with the Faraday constant F (Fortran is
!     case-insensitive).
      REAL(KIND=RKIND) :: gm, gh, gj, gd, gf, gfca, gxr, gxs, goa,
     2   goi, gua, gui, gu, gv, gw

!#######################################################################

!--------------------------------------------------------------------
!
!     Parameters for Nygreen Model.
!
!     Reference for Nygreen model:
!        Nygren A, Fiset C, Firek L, Clark JW, Lindblad DS, Clark RB,
!        Giles WR (1998).
!        Mathematical model of an adult human atrial cell: the role
!        of K+ currents in repolarization. Circulation research, 82(1)
!        https://pubmed.ncbi.nlm.nih.gov/9440706/
!
!     Nygreen mathematical model is based on the classical formulation of
!     Hodgkin and Huxley:
!        A. L. Hodgkin and A. F. Huxley (1952).
!        A Quantitative Description of Membrane Current and its Application
!        to Conduction and Excitation in Nerve. The Journal of Physiology, 117(4)
!        https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1392413/
!
!---------------------------------------------------------------------
!     Default model parameters
!     Na_b: Na concentration in bulk (bathing) medium
      REAL(KIND=RKIND) :: Na_b = 130._RKIND           ! units: mmol/L
!     K_b: K concentration in bulk (bathing) medium
      REAL(KIND=RKIND) :: K_b = 5.4_RKIND            ! units: mmol/L
!     Ca_b: Ca concentration in bulk (bathing) medium
      REAL(KIND=RKIND) :: Ca_b = 1.8_RKIND           ! units: mmol/L
!     Mg_i: Mg concentration in the intracellular medium
      REAL(KIND=RKIND) :: Mg_i = 2.5_RKIND           ! units: mmol/L
!     E_Ca,app: Apparent reversal potential for I_Ca,L
      REAL(KIND=RKIND) :: E_Ca_app = 60._RKIND        ! units: mV
!     K_Ca: Half-maximum Ca^2+ binding concentration for f_Ca
!     Note that K_Ca=0.25mmol/L in the paper, but K_Ca=0.025mmol/L
!     in CellML repository
      REAL(KIND=RKIND) :: K_Ca = 0.025_RKIND         ! units: mmol/L
!     R: Universal Gas constant
      REAL(KIND=RKIND) :: R = 8314._RKIND            ! units: mJ/mol/K
!     T: Absolute temperature
      REAL(KIND=RKIND) :: T = 306.15_RKIND           ! units: K
!     F: Faraday's constant
      REAL(KIND=RKIND) :: F = 96487._RKIND           ! units: C/mmol
!     Cm: Membrane capacitance
      REAL(KIND=RKIND) :: Cm = 0.05_RKIND            ! units: nF
!     Vol_i: Total cytosolic volume
      REAL(KIND=RKIND) :: Vol_i = 0.005884_RKIND     ! units: nL
!     Vol_c: Volume of the extracellular cleft space
!     Vol_c = 0.136 * Vol_i
      REAL(KIND=RKIND) :: Vol_c = 0.000800224_RKIND  ! units: nL
!     Vol_d: Volume of the diffusion-restricted subsarcolemmal space
!     Vol_d = 0.02 * Vol_i
      REAL(KIND=RKIND) :: Vol_d = 0.00011768_RKIND   ! units: nL
!     Vol_rel: Volume of the sarcoplasmic reticulum release compartment
      REAL(KIND=RKIND) :: Vol_rel = 0.0000441_RKIND  ! units: nL
!     Vol_up: Volume of the sarcoplasmic reticulum uptake compartment
      REAL(KIND=RKIND) :: Vol_up = 0.0003969_RKIND   ! units: nL
!     tau_Na: Time constant of diffusion of Na from the bulk medium to
!     the extracellular cleft space
      REAL(KIND=RKIND) :: tau_Na = 14.3_RKIND        ! units: s
!     tau_Na: Time constant of diffusion of K from the bulk medium to
!     the extracellular cleft space
      REAL(KIND=RKIND) :: tau_K = 10._RKIND          ! units: s
!     tau_Ca: Time constant of diffusion of Ca from the bulk medium to
!     the extracellular cleft space
      REAL(KIND=RKIND) :: tau_Ca = 24.7_RKIND        ! units: s
!     tau_di: Time constant of diffusion of Na from the restricted
!     subsarcolemmal space to the cytosol
      REAL(KIND=RKIND) :: tau_di = 0.01_RKIND        ! units: s
!     Ibar_NaK: Maximum Na-K pump current
      REAL(KIND=RKIND) :: Ibar_NaK = 70.8253_RKIND   ! units: pA
!     K_NaK,K: Half-maximum K binding concentration for I_NaK
      REAL(KIND=RKIND) :: K_NaK_K = 1._RKIND         ! units: mmol/L
!     K_NaK,Na: Half-maximum Na binding concentration for I_NaK
      REAL(KIND=RKIND) :: K_NaK_Na = 11._RKIND       ! units: mmol/L
!     Ibar_CaP: Maximum Ca pump current
      REAL(KIND=RKIND) :: Ibar_CaP = 4._RKIND        ! units: pA
!     K_CaP: Half-maximum Ca binding concentration for I_CaP
      REAL(KIND=RKIND) :: K_CaP = 0.0002_RKIND       ! units: mmol/L
!     K_NaCa: Scaling factor for I_NaCa
      REAL(KIND=RKIND) :: K_NaCa = 0.0374842_RKIND   ! units: pA/(mmol/L)^4
!     gamma: Position of energy barrier controlling voltage dependence
!     of I_NaCa
      REAL(KIND=RKIND) :: gamma = 0.45_RKIND         ! dimensionless
!     d_NaCa: Denominator constant for I_NaCa
      REAL(KIND=RKIND) :: d_NaCa = 0.0003_RKIND      ! units: (mmol/L)^{-4}
!     phi_Na,en: Electroneutral Na influx
      REAL(KIND=RKIND) :: phi_Na_en = -1.68_RKIND    ! units: pA
!     Ibar_up: Maximum sarcoplasmic reticulum uptake current
      REAL(KIND=RKIND) :: Ibar_up = 2800._RKIND      ! units: pA
!     K_cyca: Half-maximum binding concentration for Ca_i to I_up
      REAL(KIND=RKIND) :: K_cyca = 0.0003_RKIND      ! units: mmol/L
!     K_srca: Half-maximum binding concentration for Ca_up to I_up
      REAL(KIND=RKIND) :: K_srca = 0.5_RKIND         ! units: mmol/L
!     K_xcs: Ratio of forward to back reactions for I_up
      REAL(KIND=RKIND) :: K_xcs = 0.4_RKIND          ! dimensionless
!     tau_tr: Time constant of diffusion ("translocation") of Ca from
!     sarcoplasmic reticulum uptake to release compartment
      REAL(KIND=RKIND) :: tau_tr = 0.01_RKIND        ! units: s
!     alpha_rel: Scaling factor for I_rel
      REAL(KIND=RKIND) :: alpha_rel = 200000._RKIND  ! units: pA L/mmol
!     K_rel,i: Half-activation Ca_i for I_rel
      REAL(KIND=RKIND) :: K_rel_i = 0.0003_RKIND     ! units: mmol/L
!     K_rel,d: Half-activation Ca_d for I_rel
      REAL(KIND=RKIND) :: K_rel_d = 0.003_RKIND      ! units: mmol/L
!     r_recov: Recovery rate constant for the sarcoplasmic reticulum
!     release channel
      REAL(KIND=RKIND) :: r_recov = 0.815_RKIND      ! units: s^{-1}
!-----------------------------------------------------------------------
!     Maximum Conductance Values
!     P_Na: Permeability for I_Na
      REAL(KIND=RKIND) :: P_Na = 0.0016_RKIND        ! units: nL/s
!     Gbar_Ca,L: Maximum conductance for I_Ca,L
      REAL(KIND=RKIND) :: Gbar_CaL = 6.75_RKIND      ! units: nS
!     Gbar_t: Maximum conductance for I_t
      REAL(KIND=RKIND) :: Gbar_t = 7.5_RKIND         ! units: nS
!     Gbar_sus: Maximum conductance for I_sus
      REAL(KIND=RKIND) :: Gbar_sus = 2.75_RKIND      ! units: nS
!     Gbar_K,s: Maximum conductance for I_K,s
      REAL(KIND=RKIND) :: Gbar_Ks = 1._RKIND         ! units: nS
!     Gbar_K,r: Maximum conductance for I_K,r
      REAL(KIND=RKIND) :: Gbar_Kr = 0.5_RKIND        ! units: nS
!     Gbar_K1: Maximum conductance for I_K1
      REAL(KIND=RKIND) :: Gbar_K1 = 3._RKIND         ! units: nS
!     Gbar_B,Na: Maximum conductance for I_B,Na
      REAL(KIND=RKIND) :: Gbar_BNa = 0.060599_RKIND  ! units: nS
!     Gbar_B,Ca: Maximum conductance for I_B,Ca
      REAL(KIND=RKIND) :: Gbar_BCa = 0.078681_RKIND  ! units: nS
!-----------------------------------------------------------------------
!     Resting potential
      REAL(KIND=RKIND) :: Vrest = -75.0_RKIND        ! units: mV
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

!     Cellular transmembrane currents
!     I_Na: Na current
      REAL(KIND=RKIND) :: I_Na
!     I_Ca,L: L-type Ca current
      REAL(KIND=RKIND) :: I_CaL
!     I_t: Transient outward K current
      REAL(KIND=RKIND) :: I_t
!     I_sus: Sustained outward K current
      REAL(KIND=RKIND) :: I_sus
!     I_K1: Inward rectifier K current
      REAL(KIND=RKIND) :: I_K1
!     I_B,Na: Background Inward Na current
      REAL(KIND=RKIND) :: I_BNa
!     I_B,Ca: Background Inward Ca current
      REAL(KIND=RKIND) :: I_BCa
!     I_NaK: Na-K pump current
      REAL(KIND=RKIND) :: I_NaK
!     I_CaP: Sarcolemmal Ca pump current
      REAL(KIND=RKIND) :: I_CaP
!     I_NaCa: Na-Ca exchange current
      REAL(KIND=RKIND) :: I_NaCa
!     Delayed rectifier K currents
      REAL(KIND=RKIND) :: I_Ks, I_Kr
!     Other Ca currents due to sacroplasmic reticulum
      REAL(KIND=RKIND) :: I_di, I_up, I_tr, I_rel
!-----------------------------------------------------------------------
!     State variables
      REAL(KIND=RKIND) :: V
      REAL(KIND=RKIND) :: Na_c
      REAL(KIND=RKIND) :: K_c
      REAL(KIND=RKIND) :: Ca_c
      REAL(KIND=RKIND) :: Na_i
      REAL(KIND=RKIND) :: K_i
      REAL(KIND=RKIND) :: Ca_i
      REAL(KIND=RKIND) :: Ca_d
      REAL(KIND=RKIND) :: Ca_up
      REAL(KIND=RKIND) :: Ca_rel

!     F_1: Relative amount of "inactive precursor" in I_rel
      REAL(KIND=RKIND) :: F_1
!     F_2: Relative amount of "activator" in I_rel
      REAL(KIND=RKIND) :: F_2

!     Occupancy variables
      REAL(KIND=RKIND) :: O_C
      REAL(KIND=RKIND) :: O_TC
      REAL(KIND=RKIND) :: O_TMgC
      REAL(KIND=RKIND) :: O_TMgMg
      REAL(KIND=RKIND) :: O_Calse

!     Gating variables
      REAL(KIND=RKIND) :: gm, mbar, tau_m
      REAL(KIND=RKIND) :: gh_1, hbar_1, tau_h1
      REAL(KIND=RKIND) :: gh_2, hbar_2, tau_h2
      REAL(KIND=RKIND) :: gd_L, dbar_L, tau_dL
      REAL(KIND=RKIND) :: gf_L1, fbar_L1, tau_fL1
      REAL(KIND=RKIND) :: gf_L2, fbar_L2, tau_fL2
      REAL(KIND=RKIND) :: gr, rbar, tau_r
      REAL(KIND=RKIND) :: gs, sbar, tau_s
      REAL(KIND=RKIND) :: gr_sus, rbar_sus, tau_rsus
      REAL(KIND=RKIND) :: gs_sus, sbar_sus, tau_ssus
      REAL(KIND=RKIND) :: gn, nbar, tau_n
      REAL(KIND=RKIND) :: gp_a, pbar_a, tau_pa

!     Other local variables
      REAL(KIND=RKIND) :: p_i, r_act, r_inact, f_Ca

!#######################################################################

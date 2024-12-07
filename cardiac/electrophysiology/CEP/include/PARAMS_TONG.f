!--------------------------------------------------------------------
!
!     Parameters for the Tong model for an isolated uterine smooth
!     muscle cell.
!
!        Tong, W. C., Choi, C. Y., Karche, S., Holden, A. V., Zhang, H.,
!        & Taggart, M. J. (2011). A computational model of the ionic
!        currents, Ca2+ dynamics and action potentials underlying
!        contraction of isolated uterine smooth muscle. PloS one, 6(4),
!        e18685. doi:10.1371/journal.pone.0018685
!
!     Most parameters are extracted from CellML repository:
!       https://models.cellml.org/e/263/Tong_Choi_Kharche_Holden_Zhang_Taggart_2011.cellml/view
!
!--------------------------------------------------------------------

!     Default model parameters
!     T: Temperature
      REAL(KIND=RKIND) :: T = 308.0_RKIND         ! [K]
!     F: Faraday constant
      REAL(KIND=RKIND) :: F = 96485.0_RKIND       ! [C/mol]
!     R: Gas constant
      REAL(KIND=RKIND) :: R = 8314._RKIND         ! [J/K/kmol]
!     Na+ valency
      REAL(KIND=RKIND) :: z_Na = 1.0_RKIND        ! [dimensionless]
!     K+ valency
      REAL(KIND=RKIND) :: z_K = 1.0_RKIND         ! [dimensionless]
!     Ca2+ valency
      REAL(KIND=RKIND) :: z_Ca = 2.0_RKIND        ! [dimensionless]
!     Intracellular Na+ concentration
      REAL(KIND=RKIND) :: Na_i = 4.0_RKIND        ! [mM]
!     Extracellular Na+ concentration
      REAL(KIND=RKIND) :: Na_o = 130.0_RKIND      ! [mM]
!     Intracellular K+ concentration
      REAL(KIND=RKIND) :: K_i = 140.0_RKIND       ! [mM]
!     Extracellular K+ concentration
      REAL(KIND=RKIND) :: K_o = 6.0_RKIND         ! [mM]
!     Intracellular Cl- concentration
      REAL(KIND=RKIND) :: Cl_i = 46.0_RKIND       ! [mM]
!     Extracellular Cl- concentration
      REAL(KIND=RKIND) :: Cl_o = 130.0_RKIND      ! [mM]
!     Extracellular Ca+ concentration
      REAL(KIND=RKIND) :: Ca_o = 2.5_RKIND        ! [mM]
!     Extracellular Mg2+ concentration
      REAL(KIND=RKIND) :: Mg_o = 0.5_RKIND        ! [mM]
!     Specific membrane capacitance
      REAL(KIND=RKIND) :: C_m = 1.0_RKIND         ! [uF/cm^2]
!     Cell surface area to volume ratio
      REAL(KIND=RKIND) :: AV_c = 4.0_RKIND        ! [1/m]
!     Proportion of free Ca2+ ions
!        Defined as buff in CellML & src code
      REAL(KIND=RKIND) :: bet = 0.015_RKIND       ! [dimensionless]
!--------------------------------------------------------------------
!     Plasma membrane Ca-ATPase
!     PMCA Ca2+ flux
      REAL(KIND=RKIND) :: J_PMCA = 3.5E-7_RKIND   ! [mM/ms]
!     Half-saturation concentration for PMCA
!        Defined in uM, converted to mM
      REAL(KIND=RKIND) :: K_mPMCA = 5.0E-4_RKIND  ! [mM]
!     Hill coefficient for PMCA
      REAL(KIND=RKIND) :: n_PMCA = 2.0_RKIND      ! [dimensionless]
!--------------------------------------------------------------------
!     NaCa Exchanger
!     NaCa exchanger flux
      REAL(KIND=RKIND) :: J_NaCa = 3.5E-6_RKIND   ! [mM/ms]
!     Half-saturation concentration for NaCa exchanger allosteric factor
!     Wrong in Table S4 corrected with src code and CellML model
      REAL(KIND=RKIND) :: K_mAllo = 3.0E-4_RKIND  ! [mM]
!     Hill coefficient for NaCa exchanger allosteric factor
      REAL(KIND=RKIND) :: n_Allo = 4.0_RKIND      ! [dimensionless]
!     Saturation factor for NaCa exchanger
      REAL(KIND=RKIND) :: k_sat = 0.27_RKIND      ! [dimensionless]
!     Partition parameter for NaCa exchanger flux
!        xgamma in CellML
      REAL(KIND=RKIND) :: gam = 0.35_RKIND        ! [dimensionless]
!     Na_i dissociation constant for NaCa exchanger
      REAL(KIND=RKIND) :: K_mNai = 30.0_RKIND     ! [mM]
!     Ca_i dissociation constant for NaCa exchanger
!        Defined in uM, converted to mM checked with CellML & src code
      REAL(KIND=RKIND) :: K_mCai = 7.0E-3_RKIND   ! [mM]
!     Na_o dissociation constant for NaCa exchanger
      REAL(KIND=RKIND) :: K_mNao = 87.5_RKIND     ! [mM]
!     Ca_o dissociation constant for NaCa exchanger
      REAL(KIND=RKIND) :: K_mCao = 1.3_RKIND      ! [mM]
!--------------------------------------------------------------------
!     L-type Calcium (ICaL)
!     Maximum conductance of ICaL
      REAL(KIND=RKIND) :: g_CaL = 0.6_RKIND       ! [mS/uF]
!     Reversal potential of ICaL
      REAL(KIND=RKIND) :: E_CaL = 45.0_RKIND      ! [mV]
!     Half-saturation concentration for ICaL
!     Defined in uM, converted to mM; kmca in CellML
      REAL(KIND=RKIND) :: K_mCaL = 1.0E-3_RKIND   ! [mM]
!     Inactivation Steady State Slope Factor
      REAL(KIND=RKIND) :: f_ss_sf = 7.0_RKIND     ! [mV] / src code
!--------------------------------------------------------------------
!     T-type Calcium (ICaT)
!     Maximum conductance of ICaT
      REAL(KIND=RKIND) :: g_CaT = 0.058_RKIND     ! [mS/uF]
!     Reversal potential of ICaT
      REAL(KIND=RKIND) :: E_CaT = 42.0_RKIND      ! [mV]
!--------------------------------------------------------------------
!     Sodium Current (INa)
!     Maximum conductance of INa, can range 0 to 0.12 (double check)
      REAL(KIND=RKIND) :: g_Na = 0.0_RKIND        ! [mS/uF]
!--------------------------------------------------------------------
!     Hyperpolarization-activated Current (Ih)
!     Maximum conductance of Ih
      REAL(KIND=RKIND) :: g_h = 0.0542_RKIND      ! [mS/uF]
!     Permeability of I_h to NaK
!        Defined separately as P_Na=0.35 and P_K=1 in CellML, but their
!        ratio is the same
      REAL(KIND=RKIND) :: P_NaK = 0.35_RKIND      ! [dimensionless]
!--------------------------------------------------------------------
!     Potassium Currents (IK1, IK2, IKa, IKb, IK(Ca))
!     Maximum conductance of total IK
!        Not used in CellML; Instead g_K1, etc. expressions are solved
      REAL(KIND=RKIND) :: g_K = 0.8_RKIND         ! [mS/uF]
!     Maximum conductance of background potassium current
      REAL(KIND=RKIND) :: g_b = 0.004_RKIND       ! [mS/uF]
!     Maximum conductance of IK1
      REAL(KIND=RKIND) :: g_K1 = 0.52_RKIND       ! 0.65*g_K [mS/uF]
!     Maximum conductance of IK2
      REAL(KIND=RKIND) :: g_K2 = 0.032_RKIND      ! 0.04[mS/uF]
!     Maximum conductance of IKa
      REAL(KIND=RKIND) :: g_Ka = 0.16_RKIND       ! 0.2[mS/uF]
!     ??? g_BK = g_K; g_K is never used, but defined in Table S4
      REAL(KIND=RKIND) :: g_BK = 0.8_RKIND        ! [mS/uF]
!     Maximum conductance of calcium-activated potassium currents
!     (BK, BKaB)
      REAL(KIND=RKIND) :: g_KCa = 0.8_RKIND       ! [mS/uF]
!     Proportion of I_a in IK(Ca)
!     Value from paper source code not defined in Table S4
!     gbka in CellML
      REAL(KIND=RKIND) :: p_a = 0.2_RKIND         ! [dimensionless]
!     Proportion of I_ab1 in IK(Ca)
!        Value from the paper's source code not defined in Table S4
!        gbkab in CellML
      REAL(KIND=RKIND) :: p_b = 0.1_RKIND         ! [dimensionless]
!--------------------------------------------------------------------
!     Nonselective Cation Current (I_NSCC)
!     Maximum conductance of I_NSCC leak current
!        Used here but not defined in CellML,
!        g_L = 0.009685 mS/pF by Testrow et al.
      REAL(KIND=RKIND) :: g_L = 0.0_RKIND         ! [mS/uF]
!     Maximum conductance of I_NSCC
      REAL(KIND=RKIND) :: g_NS = 0.0123_RKIND     ! [mS/uF]
!     Permeability of I_NSCC to Ca:Cs
!        Defined separately as PnsCa=0.89 and PnsCs=1 in CellML,
!        but their ratio is the same
      REAL(KIND=RKIND) :: P_CaCs = 0.89_RKIND     ! [dimensionless]
!     Permeability of I_NSCC to Na:Cs
!        Defined separately as PnsNa = 0.9_RKIND and PnsCs=1 in CellML,
!        but their ratio is the same
      REAL(KIND=RKIND) :: P_NaCs = 0.9_RKIND      ! [dimensionless]
!     Permeability of I_NSCC to K:Cs
!        Defined separately as PnsK = 1.3_RKIND and PnsCs=1 in CellML,
!        but their ratio is the same
      REAL(KIND=RKIND) :: P_KCs = 1.3_RKIND       ! [dimensionless]
!     Half-saturation concentration for Mg_o inhibition of I_NSCC
!        Not explicitly defined in CellML, but used in I_NSCC equation
!        Value in Table S4 of the paper = 0.28
      REAL(KIND=RKIND) :: K_dMg = 0.281007_RKIND  ! [mM]
!--------------------------------------------------------------------
!     Calcium-Activate Chloride Current (I_Cl(Ca))
!     Maximium conductance of I_Cl(Ca)
      REAL(KIND=RKIND) :: g_Cl = 0.1875_RKIND     ! [mS/uF]
!--------------------------------------------------------------------
!     Na-K Exchanger Current (I_NaK)
!     Maximum conductance of I_NaK
!        ginak in CellML
      REAL(KIND=RKIND) :: g_NaK = 1.7_RKIND       ! [A/F]
!     Half-saturation concentration for K_o dependency of I_NaK
!        nakKmko in CellML
      REAL(KIND=RKIND) :: K_mK = 2.0_RKIND        ! [mM]
!     Hill coefficient for K_o dependency of I_NaK
!        not explicitly defined in CellML but used in knak equation
      REAL(KIND=RKIND) :: n_K = 1.5_RKIND
!     Half-saturation concentration for Na_i dependency of I_NaK
!        nakKmnai in CellML
      REAL(KIND=RKIND) :: K_mNa = 22.0_RKIND      ! [mM]
!     Hill coefficient for Na_i dependency of I_NaK
!        not explicitly defined in CellML but used in nnak equation
      REAL(KIND=RKIND) :: n_Na = 2.0_RKIND
!--------------------------------------------------------------------
!     Calcium-Dependent Force
!     Half-saturation Ca concentration for force
      REAL(KIND=RKIND) :: K_mF = 161.301E-6_RKIND ! [mM]
!     Hill coefficient for Ca_i dependency of force
      REAL(KIND=RKIND) :: n_F = 3.60205_RKIND     ! [dimensionless]
!     Maximal force
      REAL(KIND=RKIND) :: F_max = 3.0_RKIND       ! [uN]
!--------------------------------------------------------------------
!     Cellular transmembrane currents
      REAL(KIND=RKIND) :: I_CaL
      REAL(KIND=RKIND) :: I_Na
      REAL(KIND=RKIND) :: I_CaT
      REAL(KIND=RKIND) :: I_h
      REAL(KIND=RKIND) :: I_K1
      REAL(KIND=RKIND) :: I_K2
      REAL(KIND=RKIND) :: I_Ka
      REAL(KIND=RKIND) :: I_KCa
      REAL(KIND=RKIND) :: I_ClCa

      REAL(KIND=RKIND) :: I_b
      REAL(KIND=RKIND) :: I_NSCC
      REAL(KIND=RKIND) :: I_NaK
      REAL(KIND=RKIND) :: I_NaCa
      REAL(KIND=RKIND) :: J_Camem
      REAL(KIND=RKIND) :: J_PMCA_var
      REAL(KIND=RKIND) :: J_NaCa_var
!--------------------------------------------------------------------
!     Scaling factors
!     Voltage scaling
      REAL(KIND=RKIND) :: Vscale  = 1._RKIND
!     Time scaling
      REAL(KIND=RKIND) :: Tscale  = 1._RKIND
!     Voltage offset parameter
      REAL(KIND=RKIND) :: Voffset = 0._RKIND
!--------------------------------------------------------------------
!     State variables
      REAL(KIND=RKIND) :: V
      REAL(KIND=RKIND) :: Ca_i

!     Gating variables
!     Gating variables for I_Na
      REAL(KIND=RKIND) :: m, m_i
      REAL(KIND=RKIND) :: h, h_i
!     Gating variables for I_CaL
      REAL(KIND=RKIND) :: d, d_i
      REAL(KIND=RKIND) :: f_1, f_1i
      REAL(KIND=RKIND) :: f_2, f_2i
!     Gating variables for I_CaT
      REAL(KIND=RKIND) :: b_g, b_gi
      REAL(KIND=RKIND) :: g, g_i
!     Gating variables for I_K1
      REAL(KIND=RKIND) :: q, q_i
      REAL(KIND=RKIND) :: r_1, r_1i
      REAL(KIND=RKIND) :: r_2, r_2i
!     Gating variables for I_K2
      REAL(KIND=RKIND) :: p, p_i
      REAL(KIND=RKIND) :: k_1, k_1i
      REAL(KIND=RKIND) :: k_2, k_2i
!     Gating variables for I_Ka
      REAL(KIND=RKIND) :: s, s_i
      REAL(KIND=RKIND) :: x_g, x_gi
!     Gating variables for I_KCa, I_a, I_ab1
      REAL(KIND=RKIND) :: x_a, x_ai
      REAL(KIND=RKIND) :: x_ab1, x_ab1i
!     Gating variables for I_h
      REAL(KIND=RKIND) :: y, y_i
!     Gating variables for I_ClCa
      REAL(KIND=RKIND) :: c, c_i
!     For Ca_i dependent force
      REAL(KIND=RKIND) :: w, w_i
!#######################################################################

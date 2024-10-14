      !     Default model parameters
      ! T: Temperature
      REAL(KIND=RKIND) :: T = 308_RKIND  ! [K] - CellML
      ! F: Faraday constant
      REAL(KIND=RKIND) :: F = 96485_RKIND  ! [C/mol] - CellML
      ! R: Gas constant
      REAL(KIND=RKIND) :: R = 8314_RKIND  ! [J/K/kmol] - CellML
      ! Na+ valency
      REAL(KIND=RKIND) :: z_Na = 1_RKIND  ! [dimensionless] - CellML
      ! K+ valency
      REAL(KIND=RKIND) :: z_K = 1_RKIND  ! [dimensionless] - CellML
      ! Ca2+ valency
      REAL(KIND=RKIND) :: z_Ca = 2_RKIND  ! [dimensionless] - CellML
      ! Intracellular Na+ concentration
      REAL(KIND=RKIND) :: Na_i = 4_RKIND  ! [mM] - CellML
      ! Extracellular Na+ concentration
      REAL(KIND=RKIND) :: Na_o = 130_RKIND  ! [mM] - CellML
      ! Intracellular K+ concentration
      REAL(KIND=RKIND) :: K_i = 140_RKIND  ! [mM] - CellML
      ! Extracellular K+ concentration
      REAL(KIND=RKIND) :: K_o = 6_RKIND  ! [mM] - CellML
      ! Intracellular Cl- concentration
      REAL(KIND=RKIND) :: Cl_i = 46_RKIND  ! [mM] - CellML
      ! Extracellular Cl- concentration
      REAL(KIND=RKIND) :: Cl_o = 130_RKIND  ! [mM] - CellML
      ! Extracellular Ca+ concentration
      REAL(KIND=RKIND) :: Ca_o = 2.5_RKIND  ! [mM] - CellML
      ! Extracellular Mg2+ concentration
      REAL(KIND=RKIND) :: Mg_o = 0.5_RKIND ! [mM] - CellML
      ! Specific membrane capacitance
      REAL(KIND=RKIND) :: C_m = 1_RKIND ! [uF/cm^2] - CellML
      ! Cell surface area to volume ratio
      REAL(KIND=RKIND) :: AV_c = 4_RKIND ! [1/m] - CellML
      ! Proportion of free Ca2+ ions
      ! defined as buff in CellML & src code
      REAL(KIND=RKIND) :: bet = 0.015_RKIND ! [dimensionless] - CellML

      ! Plasma membrane Ca-ATPase
      ! PMCA Ca2+ flux
      REAL(KIND=RKIND) :: J_PMCA = 3.5e-7_RKIND ! [mM/ms] - CellML
      ! Half-saturation concentration for PMCA
      ! Defined in uM, converted to mM
      REAL(KIND=RKIND) :: K_mPMCA = 0.0005_RKIND ! [mM] - CellML
      ! Hill coefficient for PMCA
      REAL(KIND=RKIND) :: n_PMCA = 2_RKIND ! [dimensionless] - CellML

      ! NaCa Exchanger
      ! NaCa exchanger flux
      REAL(KIND=RKIND) :: J_NaCa = 3.5e-6_RKIND ! [mM/ms] - CellML
      ! Half-saturation concentration for NaCa exchanger allosteric factor
      ! Wrong in Table S4 corrected with src code and CellML model
      REAL(KIND=RKIND) :: K_mAllo = 0.0003_RKIND ! [mM] - CellML
      ! Hill coefficient for NaCa exchanger allosteric factor
      REAL(KIND=RKIND) :: n_Allo = 4_RKIND ! [dimensionless] - CellML
      ! Saturation factor for NaCa exchanger
      REAL(KIND=RKIND) :: k_sat = 0.27_RKIND ! [dimensionless] - CellML
      ! Partition parameter for NaCa exchanger flux
      ! xgamma in CellML
      REAL(KIND=RKIND) :: gam = 0.35_RKIND ! [dimensionless] - CellML
      ! Na_i dissociation constant for NaCa exchanger
      REAL(KIND=RKIND) :: K_mNai = 30.0_RKIND ! [mM] - CellML
      ! Ca_i dissociation constant for NaCa exchanger
      ! Defined in uM, converted to mM checked with CellML & src code
      REAL(KIND=RKIND) :: K_mCai = 7e-3_RKIND ! [mM] - CellML
      ! Na_o dissociation constant for NaCa exchanger
      REAL(KIND=RKIND) :: K_mNao= 87.5_RKIND ! [mM] - CellML
      ! Ca_o dissociation constant for NaCa exchanger
      REAL(KIND=RKIND) :: K_mCao = 1.3_RKIND ! [mM] - CellML

      ! L-type Calcium (ICaL)
      ! Maximum conductance of ICaL
      REAL(KIND=RKIND) :: g_CaL = 0.6_RKIND ! [mS/uF] - CellML
      ! Reversal potential of ICaL
      REAL(KIND=RKIND) :: E_CaL = 45_RKIND ! [mV] - CellML
      ! Half-saturation concentration for ICaL
      ! Defined in uM, converted to mM
      ! kmca in CellML
      REAL(KIND=RKIND) :: K_mCaL = 1e-3_RKIND ! [mM] - CellML
      ! Inactivation Steady State Slope Factor
      REAL(KIND=RKIND) :: f_ss_sf = 7_RKIND ! [mV] - CellML / src code

      ! T-type Calcium (ICaT)
      ! Maximum conductance of ICaT
      REAL(KIND=RKIND) :: g_CaT = 0.058_RKIND ! [mS/uF] - CellML
      ! Reversal potential of ICaT
      REAL(KIND=RKIND) :: E_CaT = 42_RKIND ! [mV] - CellML

      ! Sodium Current (INa)
      ! maximum conductance of INa
      ! can range 0 to 0.12 (double check)
      REAL(KIND=RKIND) :: g_Na = 0_RKIND ! [mS/uF] - CellML

      ! Hyperpolarization-activated Current (Ih)
      ! maximum conductance of Ih
      REAL(KIND=RKIND) :: g_h = 0.0542_RKIND ! [mS/uF] - CellML
      ! Permeability of I_h to Na:K
      ! defined as P_Na=0.35 and P_K=1 separately in CellML ratio is the same
      REAL(KIND=RKIND) :: P_NaK = 0.35_RKIND ! [dimensionless]

      ! Potassium Currents (IK1, IK2, IKa, IKb, IK(Ca))
      ! Maximum conductance of total IK
      ! Not used in CellML instead the g_K1, etc. expressions are solved
      REAL(KIND=RKIND) :: g_K = 0.8_RKIND ! [mS/uF]
      ! Maximum conductance of background potassium current
      REAL(KIND=RKIND) :: g_b = 0.004_RKIND ! [mS/uF] - CellML
      ! Maximum conductance of IK1
      REAL(KIND=RKIND) :: g_K1 = 0.52_RKIND ! 0.65*g_K [mS/uF] - CellML
      ! Maximum conductance of IK2
      REAL(KIND=RKIND) :: g_K2 = 0.032_RKIND ! 0.04[mS/uF] - CellML
      ! Maximum conductance of IKa
      REAL(KIND=RKIND) :: g_Ka = 0.16_RKIND ! 0.2[mS/uF] - CellML
      ! ???
      REAL(KIND=RKIND) :: g_BK = 0.8_RKIND ! g_K Never used, but defined in Table S4
      ! Maximum conductance of calcium-activated potassium currents (BK, BKaB)
      REAL(KIND=RKIND) :: g_KCa = 0.8_RKIND ! [mS/uF] - CellML
      ! proportion of I_a in IK(Ca)
      ! Value from paper source code not defined in Table S4
      ! gbka in CellML
      REAL(KIND=RKIND) :: p_a = 0.2_RKIND ! [dimensionless] - CellML
      ! proportion of I_ab1 in IK(Ca)
      ! Value from paper source code not defined in Table S4
      ! gbkab in CellML
      REAL(KIND=RKIND) :: p_b = 0.1_RKIND  ! [dimensionless] - CellML

      ! Nonselective Cation Current (I_NSCC)
      ! Maximum conductance of I_NSCC leak current
      ! Used here but not defined in CellML
      REAL(KIND=RKIND) :: g_L = 0_RKIND ! [mS/uF] ! 0.009685 mS/pF by Testrow et al.
      ! Maximum conductance of I_NSCC
      REAL(KIND=RKIND) :: g_NS = 0.0123_RKIND ! [mS/uF]
      ! Permeability of I_NSCC to Ca:Cs
      ! defined as PnsCa=0.89 and PnsCs=1 separately in CellML ratio is the same
      REAL(KIND=RKIND) :: P_CaCs = 0.89_RKIND ! [dimensionless]
      ! Permeability of I_NSCC to Na:Cs
      ! defined as PnsNa = 0.9_RKIND and PnsCs=1 separately in CellML ratio is the same
      REAL(KIND=RKIND) :: P_NaCs = 0.9_RKIND ! [dimensionless]
      ! Permeability of I_NSCC to K:Cs
      ! defined as PnsK = 1.3_RKIND and PnsCs=1 separately in CellML ratio is the same
      REAL(KIND=RKIND) :: P_KCs = 1.3_RKIND ! [dimensionless]
      ! Half-saturation concentration for Mg_o inhibition of I_NSCC
      ! Not explicitly defined in CellML, but used in I_NSCC equation as 0.281007
      REAL(KIND=RKIND) :: K_dMg = 0.281007_RKIND ! [mM] ! Table S4: 0.28

      ! Calcium-Activate Chloride Current (I_Cl(Ca))
      ! Maximium conductance of I_Cl(Ca)
      REAL(KIND=RKIND) :: g_Cl = 0.1875_RKIND ! [mS/uF] - CellML

      ! Na-K Exchanger Current (I_NaK)
      ! Maximum conductance of I_NaK
      ! ginak in CellML
      REAL(KIND=RKIND) :: g_NaK = 1.7_RKIND ! [A/F]
      ! Half-saturation concentration for K_o dependency of I_NaK
      ! nakKmko in CellML
      REAL(KIND=RKIND) :: K_mK = 2.0_RKIND ! [mM]
      ! Hill coefficient for K_o dependency of I_NaK
      ! not explicitly defined in CellML but used in knak equation as 1.5
      REAL(KIND=RKIND) :: n_K = 1.5_RKIND
      ! Half-saturation concentration for Na_i dependency of I_NaK
      ! nakKmnai in CellML
      REAL(KIND=RKIND) :: K_mNa = 22.0_RKIND ! [mM]
      ! Hill coefficient for Na_i dependency of I_NaK
      ! not explicitly defined in CellML but used in nnak equation as 2
      REAL(KIND=RKIND) :: n_Na = 2_RKIND

      ! Calcium-Dependent Force
      ! Half-saturation Ca concentration for force
      REAL(KIND=RKIND) :: K_mF = 161.301_RKIND * 1E-6 ! [mM]
      ! Hill coefficient for Ca_i dependency of force
      REAL(KIND=RKIND) :: n_F = 3.60205_RKIND ! [dimensionless]
      ! Maximal force
      REAL(KIND=RKIND) :: F_max = 3.0_RKIND ! [uN]

      !     Cellular transmembrane currents
      REAL(KIND=RKIND) :: I_CaL
      REAL(KIND=RKIND) :: I_Na
      REAL(KIND=RKIND) :: I_CaT
      REAL(KIND=RKIND) :: I_h
      REAL(KIND=RKIND) :: I_K1
      REAL(KIND=RKIND) :: I_K2
      REAL(KIND=RKIND) :: I_Ka
      REAL(KIND=RKIND) :: I_KCa
      REAL(KIND=RKIND) :: I_b
      REAL(KIND=RKIND) :: I_ClCa
      REAL(KIND=RKIND) :: I_NSCC
      REAL(KIND=RKIND) :: I_NaK
      REAL(KIND=RKIND) :: I_NaCa
      REAL(KIND=RKIND) :: J_PMCA_var
      REAL(KIND=RKIND) :: J_NaCa_var
      REAL(KIND=RKIND) :: J_Camem
      REAL(KIND=RKIND) :: I_NSNa
      REAL(KIND=RKIND) :: I_NSCa
      REAL(KIND=RKIND) :: I_NSK

      !-----------------------------------------------------------------------
!     Scaling factors
!     Voltage scaling
      REAL(KIND=RKIND) :: Vscale  = 1._RKIND
!     Time scaling
      REAL(KIND=RKIND) :: Tscale  = 1._RKIND
!     Voltage offset parameter
      REAL(KIND=RKIND) :: Voffset = 0._RKIND

      !-----------------------------------------------------------------------
      !     State variables
      REAL(KIND=RKIND) :: V
      REAL(KIND=RKIND) :: Ca_i

      !     Gating variables (runtime, steady state)
      REAL(KIND=RKIND) :: m
      REAL(KIND=RKIND) :: h
      REAL(KIND=RKIND) :: b_gate ! renamed from b to avoid interferring with temp vars in GetF and GetG
      REAL(KIND=RKIND) :: g
      REAL(KIND=RKIND) :: d
      REAL(KIND=RKIND) :: f_1
      REAL(KIND=RKIND) :: f_2
      REAL(KIND=RKIND) :: q
      REAL(KIND=RKIND) :: r_1
      REAL(KIND=RKIND) :: r_2
      REAL(KIND=RKIND) :: p
      REAL(KIND=RKIND) :: k_1
      REAL(KIND=RKIND) :: k_2
      REAL(KIND=RKIND) :: x_a
      REAL(KIND=RKIND) :: x_ab1
      REAL(KIND=RKIND) :: s
      REAL(KIND=RKIND) :: x_gate  ! renamed from x to x_Ka due to lack of case-sensitivty in Fortran
      REAL(KIND=RKIND) :: y
      REAL(KIND=RKIND) :: c
      REAL(KIND=RKIND) :: w

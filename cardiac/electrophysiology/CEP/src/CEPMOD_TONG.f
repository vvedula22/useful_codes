!-----------------------------------------------------------------------
!
!     This module defines data structures for the Tong cellular
!     activation model for uterine smooth muscle cells. There is no
!     contraction model coupled to this excitation model yet.
!
!     Reference for Tong uterine myocyte excitation model:
!        Tong, W. C., Choi, C. Y., Karche, S., Holden, A. V., Zhang, H.,
!        & Taggart, M. J. (2011). A computational model of the ionic
!        currents, Ca2+ dynamics and action potentials underlying
!        contraction of isolated uterine smooth muscle. PloS one, 6(4),
!        e18685. doi:10.1371/journal.pone.0018685
!
!     Equations given in the manuscript appendix were corrected using
!     the model from CellML repository:
!       https://models.cellml.org/e/263/Tong_Choi_Kharche_Holden_Zhang_Taggart_2011.cellml/view
!
!-----------------------------------------------------------------------
      MODULE TONGMOD
      USE TYPEMOD
      USE UTILMOD, ONLY : stdL
      IMPLICIT NONE

      PRIVATE

      INCLUDE "PARAMS_TONG.f"

      PUBLIC :: TONG_INIT
      PUBLIC :: TONG_READPARFF
      PUBLIC :: TONG_INTEGFE
      PUBLIC :: TONG_INTEGRK

      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE TONG_INIT(nX, nG, X, Xg)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX, nG
      REAL(KIND=RKIND), INTENT(OUT) :: X(nX), Xg(nG)

!     Initialize state variables. Values taken from paper's src code.
      X(1)   = -53.90915441282156_RKIND      ! V     [mV]
      X(2)   = 0.0001161881607214449_RKIND   ! Ca_i  [mM]

!     Initialize gating variables. Values taken from paper's src code.
      Xg(1)  = 0.1253518889572223_RKIND      ! m     I_Na
      Xg(2)  = 0.404599170710196_RKIND       ! h     I_Na
      Xg(3)  = 0.508117603077852_RKIND       ! b_g   I_CaT
      Xg(4)  = 0.03582573962705717_RKIND     ! g     I_CaT
      Xg(5)  = 0.01036961357784695_RKIND     ! d     I_CaL
      Xg(6)  = 0.9065941499695301_RKIND      ! f_1   I_CaL
      Xg(7)  = 0.9065967263076083_RKIND      ! f_2   I_CaL
      Xg(8)  = 0.2060363247740295_RKIND      ! q     I_K1
      Xg(9)  = 0.1922244113609531_RKIND      ! r_1   I_K1
      Xg(10) = 0.1932803618375963_RKIND      ! r_2   I_K1
      Xg(11) = 0.1174074734567931_RKIND      ! p     I_K2
      Xg(12) = 0.9968385770271651_RKIND      ! k_1   I_K2
      Xg(13) = 0.9968408069904307_RKIND      ! k_2   I_K2
      Xg(14) = 0.0003569126518797985_RKIND   ! x_a   I_BKa
      Xg(15) = 0.002220456569762898_RKIND    ! x_ab1 I_BKab
      Xg(16) = 0.0307583106982354_RKIND      ! s     I_Ka
      Xg(17) = 0.08785242843398365_RKIND     ! x_g   I_Ka
      Xg(18) = 0.002604864867063448_RKIND    ! y     I_h
      Xg(19) = 0.0003764413740731269_RKIND   ! c     I_Cl
      Xg(20) = 0.2345238135343783_RKIND      ! w     Force

      RETURN
      END SUBROUTINE TONG_INIT
!-----------------------------------------------------------------------
      SUBROUTINE TONG_READPARFF(fname)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: fname

      INTEGER fid

      fid = 1528

      OPEN(fid, FILE=TRIM(fname))
      CALL GETRVAL(fid, "T", T)
      CALL GETRVAL(fid, "F", F)
      CALL GETRVAL(fid, "R", R)
      CALL GETRVAL(fid, "z_Na", z_Na)
      CALL GETRVAL(fid, "z_K", z_K)
      CALL GETRVAL(fid, "z_Ca", z_Ca)
      CALL GETRVAL(fid, "Na_i", Na_i)
      CALL GETRVAL(fid, "Na_o", Na_o)
      CALL GETRVAL(fid, "K_i", K_i)
      CALL GETRVAL(fid, "K_o", K_o)
      CALL GETRVAL(fid, "Cl_i", Cl_i)
      CALL GETRVAL(fid, "Cl_o", Cl_o)
      CALL GETRVAL(fid, "Ca_o", Ca_o)
      CALL GETRVAL(fid, "Mg_o", Mg_o)
      CALL GETRVAL(fid, "C_m", C_m)
      CALL GETRVAL(fid, "AV_c", AV_c)
      CALL GETRVAL(fid, "bet", bet)
      CALL GETRVAL(fid, "J_PMCA", J_PMCA)
      CALL GETRVAL(fid, "K_mPMCA", K_mPMCA)
      CALL GETRVAL(fid, "n_PMCA", n_PMCA)
      CALL GETRVAL(fid, "J_NaCa", J_NaCa)
      CALL GETRVAL(fid, "K_mAllo", K_mAllo)
      CALL GETRVAL(fid, "n_Allo", n_Allo)
      CALL GETRVAL(fid, "k_sat", k_sat)
      CALL GETRVAL(fid, "gam", gam)
      CALL GETRVAL(fid, "K_mNai", K_mNai)
      CALL GETRVAL(fid, "K_mCai", K_mCai)
      CALL GETRVAL(fid, "K_mNao", K_mNao)
      CALL GETRVAL(fid, "K_mCao", K_mCao)
      CALL GETRVAL(fid, "g_CaL", g_CaL)
      CALL GETRVAL(fid, "E_CaL", E_CaL)
      CALL GETRVAL(fid, "K_mCaL", K_mCaL)
      CALL GETRVAL(fid, "f_ss_sf", f_ss_sf)
      CALL GETRVAL(fid, "g_CaT", g_CaT)
      CALL GETRVAL(fid, "E_CaT", E_CaT)
      CALL GETRVAL(fid, "g_Na", g_Na)
      CALL GETRVAL(fid, "g_h", g_h)
      CALL GETRVAL(fid, "P_NaK", P_NaK)
      CALL GETRVAL(fid, "g_K", g_K)
      CALL GETRVAL(fid, "g_b", g_b)
      CALL GETRVAL(fid, "g_K1", g_K1)
      CALL GETRVAL(fid, "g_K2", g_K2)
      CALL GETRVAL(fid, "g_Ka", g_Ka)
      CALL GETRVAL(fid, "g_BK", g_BK)
      CALL GETRVAL(fid, "g_KCa", g_KCa)
      CALL GETRVAL(fid, "p_a", p_a)
      CALL GETRVAL(fid, "p_b", p_b)
      CALL GETRVAL(fid, "g_L", g_L)
      CALL GETRVAL(fid, "g_NS", g_NS)
      CALL GETRVAL(fid, "P_CaCs", P_CaCs)
      CALL GETRVAL(fid, "P_NaCs", P_NaCs)
      CALL GETRVAL(fid, "P_KCs", P_KCs)
      CALL GETRVAL(fid, "K_dMg", K_dMg)
      CALL GETRVAL(fid, "g_Cl", g_Cl)
      CALL GETRVAL(fid, "g_NaK", g_NaK)
      CALL GETRVAL(fid, "K_mK", K_mK)
      CALL GETRVAL(fid, "n_K", n_K)
      CALL GETRVAL(fid, "K_mNa", K_mNa)
      CALL GETRVAL(fid, "n_Na", n_Na)
      CALL GETRVAL(fid, "K_mF", K_mF)
      CALL GETRVAL(fid, "n_F", n_F)
      CALL GETRVAL(fid, "F_max", F_max)
      CALL GETRVAL(fid, "Vscale", Vscale)
      CALL GETRVAL(fid, "Tscale", Tscale)
      CALL GETRVAL(fid, "Voffset", Voffset)
      CLOSE(fid)

      RETURN
      END SUBROUTINE TONG_READPARFF
!-----------------------------------------------------------------------
!     Time integration performed using Forward Euler method
      SUBROUTINE TONG_INTEGFE(nX, nG, X, Xg, dt, Istim, Ksac, RPAR)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX, nG
      REAL(KIND=RKIND), INTENT(INOUT) :: X(nX), Xg(nG), RPAR(22)
      REAL(KIND=RKIND), INTENT(IN) :: dt, Istim, Ksac

      REAL(KIND=RKIND) :: f(nX)

!     Get time derivatives of X (RHS)
      CALL TONG_GETF(nX, nG, X, Xg, f, Istim, Ksac, RPAR)

!     Updating gating variables
      CALL TONG_UPDATEG(dt, nX, nG, X, Xg)

!     Update state variables
      X = X + dt*f(:)

      RETURN
      END SUBROUTINE TONG_INTEGFE
!-----------------------------------------------------------------------
!     Time integration performed using 4th order Runge-Kutta method
      SUBROUTINE TONG_INTEGRK(nX, nG, X, Xg, dt, Istim, Ksac, RPAR)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX, nG
      REAL(KIND=RKIND), INTENT(INOUT) :: X(nX), Xg(nG), RPAR(22)
      REAL(KIND=RKIND), INTENT(IN) :: dt, Istim, Ksac

      REAL(KIND=RKIND) :: dt6, Xrk(nX), Xgr(nG), frk(nX,4)

      dt6 = dt/6._RKIND

!     RK4: 1st pass
      Xrk = X
      CALL TONG_GETF(nX, nG, Xrk, Xg, frk(:,1), Istim, Ksac, RPAR)

!     Update gating variables by half-dt
      Xgr = Xg
      CALL TONG_UPDATEG(0.5_RKIND*dt, nX, nG, X, Xgr)

!     RK4: 2nd pass
      Xrk = X + 0.5_RKIND*dt*frk(:,1)
      CALL TONG_GETF(nX, nG, Xrk, Xgr, frk(:,2), Istim, Ksac, RPAR)

!     RK4: 3rd pass
      Xrk = X + 0.5_RKIND*dt*frk(:,2)
      CALL TONG_GETF(nX, nG, Xrk, Xgr, frk(:,3), Istim, Ksac, RPAR)

!     Update gating variables by full-dt
      Xgr = Xg
      CALL TONG_UPDATEG(dt, nX, nG, X, Xgr)

!     RK4: 4th pass
      Xrk = X + dt*frk(:,3)
      CALL TONG_GETF(nX, nG, Xrk, Xgr, frk(:,4), Istim, Ksac, RPAR)

      X = X + dt6*(frk(:,1) + 2._RKIND*(frk(:,2) + frk(:,3)) + frk(:,4))
      Xg = Xgr

      RETURN
      END SUBROUTINE TONG_INTEGRK
!------------------------------------------------------
!     Compute currents and time derivatives of state variables
      SUBROUTINE TONG_GETF(nX, nG, X, Xg, dX, I_stim, K_sac, RPAR)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX, nG
      REAL(KIND=RKIND), INTENT(IN) :: X(nX), Xg(nG), I_stim, K_sac
      REAL(KIND=RKIND), INTENT(OUT) :: dX(nX)
      REAL(KIND=RKIND), INTENT(INOUT) :: RPAR(22)

      REAL(KIND=RKIND) :: den1, den2, den3, num1, num2, var1, var2
      REAL(KIND=RKIND) :: RTdivF, VFdRT, I_sac, f_Ca, E_Na, E_h, E_K,
     2   I_a, I_ab1, E_Cl, PpCa, E_NS, f_Mg, g_NSCa, g_NSNa, g_NSK,
     3   I_NSCa, f_allo, f1_NaCa, f2_NaCa

!     Local copies of state variables
      V     = X(1)
      Ca_i  = X(2)

!     Local copies of gating variables
      m     = Xg(1)
      h     = Xg(2)
      b_g   = Xg(3)
      g     = Xg(4)
      d     = Xg(5)
      f_1   = Xg(6)
      f_2   = Xg(7)
      q     = Xg(8)
      r_1   = Xg(9)
      r_2   = Xg(10)
      p     = Xg(11)
      k_1   = Xg(12)
      k_2   = Xg(13)
      x_a   = Xg(14)
      x_ab1 = Xg(15)
      s     = Xg(16)
      x_g   = Xg(17)
      y     = Xg(18)
      c     = Xg(19)
      w     = Xg(20)

!     Stretch activated currents if needed
      I_sac = K_sac * 0.0_RKIND

!     V*F/RT
      RTdivF = (R*T)/F
      VFdRT  = V/RTdivF

!     I_CaL: L-type Calcium current
      den1  = 1.0_RKIND + (Ca_i/K_mCaL)**4.0_RKIND
      f_Ca  = 1.0_RKIND / den1
      I_CaL = g_CaL*d*d*f_Ca*(0.8_RKIND*f_1 + 0.2_RKIND*f_2)*(V-E_CaL)

!     I_Na: Sodium current
      E_Na = RTdivF * LOG(Na_o/Na_i)
      I_Na = g_Na * (m**3.0_RKIND) * h * (V-E_Na)

!     I_CaT: T-type Calcium current
      I_CaT = g_CaT * (b_g*b_g) * g * (V-E_CaT)

!     I_h: Hyperpolarisation-activated current
      num1 = K_o + (P_NaK*Na_o)
      den1 = K_i + (P_NaK*Na_i)
      E_h  = RTdivF * LOG(num1/den1)
      I_h  = g_h * y * (V-E_h)

!     I_K1: Voltage-dependent potassium current
      E_K  = RTdivF * LOG(K_o/K_i)
!     Appdx S1: 0.38, 0.62 CellML, src: 0.38, 0.63
      I_K1 = g_K1*(q*q)*(0.38_RKIND*r_1 + 0.63_RKIND*r_2)*(V-E_K)

!     I_K2: Voltage-dependent potassium current
      I_K2 = g_K2*(p*p)*(0.75_RKIND*k_1 + 0.25_RKIND*k_2)*(V-E_K)

!     I_Ka: Transient potassium current
      I_Ka = g_Ka * s * x_g * (V-E_K)

!     I_K,Ca: Calcium-activated potassium current
!     Calculated as I_BKa and I_BKab in CellML & src code
!     Should expand to the same
      I_a   = x_a   * (V-E_K)
      I_ab1 = x_ab1 * (V-E_K)
      I_KCa = g_KCa * (p_a*I_a + p_b*I_ab1)

!     I_b: Background current
      I_b = g_b * (V-E_K)

!     I_Cl(Ca) or I_Cl: Calcium-activated chloride current
      E_Cl = RTdivF * LOG(Cl_i/Cl_o)
      I_ClCa = g_Cl * c * (V-E_Cl)

!     I_NSCC: Non-selective cation current
      PpCa = P_CaCs / (1.0_RKIND + EXP(VFdRT))
      num1 = P_NaCs*Na_o + P_KCs*K_o + 4.0_RKIND*PpCa*Ca_o
      den1 = P_NaCs*Na_i + P_KCs*K_i + 4.0_RKIND*PpCa*Ca_i
      E_NS = RTdivF * LOG(num1/den1)

      den1 = 1.0_RKIND + (Mg_o/K_dMg)**1.29834_RKIND
      f_Mg = 0.108043_RKIND + (0.903902_RKIND/den1)

      den1 = 5.25e-4_RKIND
      den2 = 1.0_RKIND + (150.0_RKIND/(Ca_o + 1e-8_RKIND))**2.0_RKIND
      g_NSCa = g_NS*0.015_RKIND/(den1*den2)

      den1 = 0.0123_RKIND
      den2 = 1.0_RKIND + (150.0_RKIND/(Na_o + 1e-8_RKIND))**2.0_RKIND
      g_NSNa = g_NS*0.03_RKIND/(den1*den2)

      den1 = 0.0123_RKIND
      den2 = 1.0_RKIND + (150.0_RKIND/(K_o  + 1e-8_RKIND))**2.0_RKIND
      g_NSK  = g_NS*0.0357_RKIND/(den1*den2)

      I_NSCa = g_NSCa * f_Mg * (V-E_NS)
      I_NSCC = (g_NSCa + g_NSNa + g_NSK + g_L) * f_Mg * (V-E_NS)

!     I_NaK: Sodium-potassium pump current
!     CellML: 0.1245 Table S4: 0.125
      var1 = 0.1245_RKIND*EXP(-0.1_RKIND*VFdRT)
      var2 = 0.00219_RKIND*EXP((Na_o/49.71_RKIND)+(-1.9_RKIND*VFdRT))
      den1 = 1.0_RKIND + var1 + var2
      den2 = 1.0_RKIND + ((K_mK/K_o)**n_K)
      den3 = 1.0_RKIND + ((K_mNa/Na_i)**n_Na)
      I_NaK = g_NaK / (den1*den2*den3)

!---------------------
!     Calculate fluxes
!     J_PMCA: Ca flux from plasma membrane Ca-ATPase
      den1 = 1.0_RKIND + (K_mPMCA/Ca_i)**n_PMCA
      J_PMCA_var = J_PMCA/den1

!     J_NaCa: flux from Sodium-Calcium exchanger
      den1 = 1.0_RKIND + ((K_mAllo/Ca_i)**n_Allo)
      f_allo  = 1.0_RKIND / den1
      f1_NaCa = EXP((gam-1.0_RKIND)*VFdRT)
      f2_NaCa = EXP(gam*VFdRT)

      den1 = 1.0_RKIND + (k_sat*f1_NaCa)
      den2 = ((K_mCao + Ca_o)*Na_i**3._RKIND)
     2     + (Ca_i*(K_mNao**3._RKIND + Na_o**3._RKIND))
      den2 = den2+ ((K_mNai**3._RKIND)*Ca_o*(1._RKIND+(Ca_i/K_mCai)))
      den2 = den2+ ((Na_o**3._RKIND)*K_mCai*(1._RKIND+(Na_i/K_mNai)**3))

      num1 = J_NaCa*f_allo
      num2 = (Na_i**3._RKIND)*Ca_o*f2_NaCa
      num2 = num2 - ((Na_o**3._RKIND)*Ca_i*f1_NaCa)

!     Note: Equations for I_NaCa/J_NaCa vary between all 3 sources. Both
!     CellML & src include the 0.5 factor, while Appdx S1 does not
      J_NaCa_var = -(num1*num2) / (den1*den2)

!     I_NaCa: Sodium-Calcium exchanger current
      I_NaCa = -0.5*(z_Ca*F*J_NaCa_var) / (C_m*bet*AV_c)

!     J_Ca,mem
!     Equation reproduced from source code not given in Appendix S1,
!     but given as Eq. 7 in manuscript
      J_Camem = ((AV_c*C_m*bet)/(z_Ca*F)) * (I_CaL + I_CaT + I_NSCa)

!-----------------------------
!     Compute time derivatives
!     dV/dt
      dX(1) = -(I_CaL + I_Na + I_CaT + I_h + I_K1 + I_K2 + I_Ka + I_KCa
     2        + I_ClCa + I_b + I_NSCC + I_NaK + I_NaCa + I_stim)

!     dCa_i/dt
      dX(2) = -(J_Camem + J_NaCa_var + J_PMCA_var)

!-----------------------------------------------------------------------
!     Output: start at 3 as first two values are reserved for tolerances
      RPAR(3)  = I_CaL
      RPAR(4)  = I_Na
      RPAR(5)  = I_CaT
      RPAR(6)  = I_h
      RPAR(7)  = I_K1
      RPAR(8)  = I_K2
      RPAR(9)  = I_Ka
      RPAR(10) = I_KCa
      RPAR(11) = I_ClCa
      RPAR(12) = I_b
      RPAR(13) = I_NSCC
      RPAR(14) = I_NaK
      RPAR(15) = I_NaCa
      RPAR(16) = J_Camem
      RPAR(17) = J_PMCA_var
      RPAR(18) = J_NaCa_var

!     Force is calculted separately and output directly to file
      RPAR(19) = F_max * w

      RETURN
      END SUBROUTINE TONG_GETF
!-----------------------------------------------------------------------
!     Update all the gating variables
      SUBROUTINE TONG_UPDATEG(dt, nX, nG, X, Xg)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX, nG
      REAL(KIND=RKIND), INTENT(IN) :: dt, X(nX)
      REAL(KIND=RKIND), INTENT(INOUT) :: Xg(nG)

!     Define temp vars
      REAL(KIND=RKIND) :: RTdF, VFdRT, tau, rval, den1, den2, var1, var2

!     Local copies of state variables
      V     = X(1)
      Ca_i  = X(2)

!     Local copies of gating variables
      m     = Xg(1)
      h     = Xg(2)
      b_g   = Xg(3)
      g     = Xg(4)
      d     = Xg(5)
      f_1   = Xg(6)
      f_2   = Xg(7)
      q     = Xg(8)
      r_1   = Xg(9)
      r_2   = Xg(10)
      p     = Xg(11)
      k_1   = Xg(12)
      k_2   = Xg(13)
      x_a   = Xg(14)
      x_ab1 = Xg(15)
      s     = Xg(16)
      x_g   = Xg(17)
      y     = Xg(18)
      c     = Xg(19)
      w     = Xg(20)

!     V*F/RT
      RTdF  = (R*T)/F
      VFdRT = V*F/(R*T)

!     Note: values used in gating variables vary greatly between
!     Appdx S1 & src/CellML (mostly by precision/sigfigs)

!     m gate for I_Na (CellML)
      den1   = 1.0_RKIND + EXP(-(V+35.9584_RKIND)/9.24013_RKIND)
      m_i    = 1.0_RKIND/den1
      den1   = 1.0_RKIND + EXP((V+38.0_RKIND)/10.0_RKIND)
      tau    = 0.25_RKIND + (7.0_RKIND/den1)
      Xg(1)  = m_i - ((m_i - m)*EXP(-dt/tau))

!     h gate for I_Na
      den1   = 1.0_RKIND + EXP((V+57.0_RKIND)/8.0_RKIND)
      h_i    = 1.0_RKIND/den1
      rval   = (V+47.5_RKIND)/1.5_RKIND
      den1   = 1.0_RKIND + (rval*rval)
      tau    = 0.9_RKIND + (1002.85_RKIND/den1)
      Xg(2)  = h_i - ((h_i - h)*EXP(-dt/tau))

!     b gate for I_CaT
      den1   = 1.0_RKIND + EXP(-(V+54.23_RKIND)/9.88_RKIND)
      b_gi   = 1.0_RKIND/den1
      rval   = (V+66.0_RKIND)/26.0_RKIND
      den1   = 1.0_RKIND + (rval*rval)
      tau    = 0.45_RKIND + 3.9_RKIND/den1
      Xg(3)  = b_gi - ((b_gi - b_g)*EXP(-dt/tau))

!     g gate for I_CaT (CellML)
      den1   = 1.0_RKIND + EXP((V+72.978_RKIND)/4.64_RKIND)
      g_i    = 0.02_RKIND + 0.98_RKIND/den1
      den1   = 1.0_RKIND + EXP((V-417.43_RKIND)/203.18_RKIND)
      den2   = 1.0_RKIND + EXP(-(V+61.11_RKIND)/8.07_RKIND)
      tau    = 150.0_RKIND * (1.0_RKIND - (1.0_RKIND/(den1*den2)))
      Xg(4)  = g_i - ((g_i - g)*EXP(-dt/tau))

!     d gate for I_CaL
      den1   = 1.0_RKIND + EXP(-(V+22.0_RKIND)/7.0_RKIND)
      d_i    = 1.0_RKIND/den1
      rval   = (V+29.97_RKIND)/9.0_RKIND
      den1   = 1.0_RKIND + (rval*rval)
      tau    = 2.29_RKIND + 5.7_RKIND/den1
      Xg(5)  = d_i - ((d_i - d)*EXP(-dt/tau))

!     f_1 gate for I_CaL
      den1   = 1.0_RKIND + EXP((V+38.0_RKIND)/f_ss_sf)
      f_1i   = 1.0_RKIND/den1
      tau    = 12.0_RKIND
      Xg(6)  = f_1i - ((f_1i - f_1)*EXP(-dt/tau))

!     f_2 gate for I_CaL (CellML)
      f_2i   = f_1i
      den1   = 1.0_RKIND + EXP((V+13.9629_RKIND)/45.3782_RKIND)
      den2   = 1.0_RKIND + EXP(-(V+9.49866_RKIND)/3.3945_RKIND)
      tau    = 90.9699 * (1.0_RKIND- (1.0_RKIND/(den1*den2)))
      Xg(7)  = f_2i - ((f_2i - f_2)*EXP(-dt/tau))

!     q gate for I_K1 (CellML)
      den1   = 1.0_RKIND + EXP(-(V+18.6736_RKIND)/26.6606_RKIND)
      q_i    = 0.978613_RKIND/den1
      rval   = (V+60.71_RKIND)/15.79_RKIND
      den1   = 1.0_RKIND + (rval*rval)
      tau    = 500.0_RKIND/den1
      Xg(8)  = q_i - ((q_i - q)*EXP(-dt/tau))

!     r_1 gate for I_K1 (CellML)
      den1   = 1.0_RKIND + EXP((V+63.0_RKIND)/6.3_RKIND)
      r_1i   = 1.0_RKIND/den1
      rval   = (V+62.7133_RKIND)/35.8611_RKIND
      den1   = 1.0_RKIND + (rval*rval)
      tau    = 5.0e3_RKIND/den1
      Xg(9)  = r_1i - ((r_1i - r_1)*EXP(-dt/tau))

!     r_2 gate for I_K1
      r_2i   = r_1i
      den1   = 1.0_RKIND + EXP((V+22.0_RKIND)/4.0_RKIND)
      tau    = 3.0e4_RKIND + (2.2e5_RKIND/den1)
      Xg(10) = r_2i - ((r_2i - r_2)*EXP(-dt/tau))

!     p gate for I_K2 (CellML)
      ! Equation varies greatly between CellML/src and Appdx S1
      ! Appdx S1:
      ! p_i_den = 1+EXP((-(V+0.948))/17.91)
      ! p_i = 1/p_ss_den
      ! CellML/src
      den1   = 1.0_RKIND + EXP(-(V+17.91_RKIND)/18.4_RKIND)
      p_i    = 0.948_RKIND/den1
      rval   = (V+64.1_RKIND)/28.67_RKIND
      den1   = 1.0_RKIND + (rval*rval)
      tau    = 100.0_RKIND/den1
      Xg(11) = p_i - ((p_i - p)*EXP(-dt/tau))

!     k_1 gate for I_K2
      den1   = 1.0_RKIND + EXP((V+21.2_RKIND)/5.7_RKIND)
      k_1i   = 1.0_RKIND/den1
      den1   = 1.0_RKIND + EXP((V-315.0_RKIND)/50.0_RKIND)
      den2   = 1.0_RKIND + EXP(-(V+74.9_RKIND)/8.0_RKIND)
      tau    = 1.0e6_RKIND * (1.0_RKIND - (1.0_RKIND/(den1*den2)))
      Xg(12) = k_1i - ((k_1i - k_1)*EXP(-dt/tau))

!     k_2 gate for I_K2
      k_2i   = k_1i
      den1   = 1.0_RKIND + EXP((V-132.868_RKIND)/25.3992_RKIND)
      den2   = 1.0_RKIND + EXP(-(V+24.9203_RKIND)/2.67915_RKIND)
      tau    = 2.5e6_RKIND * (1.0_RKIND - (1.0_RKIND/(den1*den2)))
      Xg(13) = k_2i - ((k_2i - k_2)*EXP(-dt/tau))

!     x_a gate for I_KCa (CellML)
      rval   = ((1000.0_RKIND*Ca_i) + 1538.29_RKIND)/739.057_RKIND
      den1   = 1.0_RKIND + (rval*rval)
      rval   = ((1000.0_RKIND*Ca_i) - 0.0630535_RKIND)/0.161942_RKIND
      den2   = 1.0_RKIND + (rval*rval)
      var1   = (8.38384_RKIND/den1) - (0.749234_RKIND/den2)

      rval   = ((1000.0_RKIND*Ca_i) + 0.237503_RKIND)/2.39278e-4_RKIND
      den1   = 1.0_RKIND + (rval)**0.42291_RKIND
      var2   = (5011.47_RKIND/den1) - 37.5137_RKIND

      den1   = 1.0_RKIND + EXP((-var1*(V-var2))/RTdF)
      x_ai   = 1.0_RKIND/den1

      rval   = -(V-158.779_RKIND)/52.1497_RKIND
      den1   = 1.0_RKIND + (rval*rval)
      tau    = 2.40914_RKIND/den1
      Xg(14) = x_ai - ((x_ai - x_a)*EXP(-dt/tau))

!     x_ab1 gate for I_KCa
      rval   = ((1000.0_RKIND*Ca_i)+228.71_RKIND)/684.946_RKIND
      den1   = 1.0_RKIND + (rval*rval)
      rval   = ((1000.0_RKIND*Ca_i)-0.219_RKIND)/0.428_RKIND
      den2   = 1.0_RKIND + (rval*rval)
      var1   = (1.4_RKIND/den1) - (0.681249_RKIND/den2)

      rval   = ((1000.0_RKIND*Ca_i)+0.401189_RKIND)/3.99115e-3_RKIND
      den1   = 1.0_RKIND + (rval)**0.668054_RKIND
      var2   = (8540.23_RKIND/den1) - 109.275_RKIND

      den1   = 1.0_RKIND + EXP((-var1*(V-var2))/RTdF)
      x_ab1i = 1.0_RKIND/den1

      rval   = (V-153.019_RKIND)/66.4952_RKIND
      den1   = 1.0_RKIND + (rval*rval)
      tau    = 13.8049_RKIND/den1
      Xg(15) = x_ab1i - ((x_ab1i - x_ab1)*EXP(-dt/tau))

!     s gate for I_Ka
      den1   = 1.0_RKIND + EXP(-(V+27.79_RKIND)/7.57_RKIND)
      s_i    = 1.0_RKIND/den1
      rval   = (V+20.52_RKIND)/35.0_RKIND
      den1   = 1.0_RKIND + (rval*rval)
      tau    = 17.0_RKIND/den1
      Xg(16) = s_i - ((s_i - s)*EXP(-dt/tau))

!     x gate for I_Ka
      den1   = 1.0_RKIND + EXP((V+69.5_RKIND)/6.0_RKIND)
      x_gi   = 0.02_RKIND + (0.98_RKIND/den1)
      rval   = (V+34.18_RKIND)/120.0_RKIND
      den1   = 1.0_RKIND + (rval*rval)
      tau    = 7.5_RKIND + (10.0_RKIND/den1)
      Xg(17) = x_gi - ((x_gi - x_g)*EXP(-dt/tau))

!     y gate for I_h (CellML)
      den1   = 1.0_RKIND + EXP((V+105.39_RKIND)/8.6553_RKIND)
      y_i    = 1.0_RKIND/den1
      den1   = (3.5e-6_RKIND*EXP(-0.0497_RKIND*V))
      den1   = den1 + (0.04003_RKIND*EXP(0.05211_RKIND*V))
      tau    = 1.0_RKIND/den1
      Xg(18) = y_i - ((y_i - y)*EXP(-dt/tau))

!     c gate for I_ClCa
      var1   = 6e-4_RKIND*EXP(2.53_RKIND*VFdRT)
      var2   = 0.1_RKIND*EXP(-5.0_RKIND*VFdRT)
      rval   = var1/Ca_i
      den1   = (rval*rval) + rval + 1.0_RKIND
      c_i    = 1.0_RKIND / (1.0_RKIND + (var2*den1))

      den1   = EXP((V+4.56_RKIND)/11.62_RKIND)
      den2   = EXP(-(V+25.5_RKIND)/11.62_RKIND)
      var1   = 210.0_RKIND/(1.0_RKIND+den1)
      var2   = 170.0_RKIND/(1.0_RKIND+den2)
      tau    = var1 + var2 - 160.0_RKIND
      Xg(19) = c_i - ((c_i - c)*EXP(-dt/tau))

!     w gate for Calcium-dependent force (VV: check expression for tau)
      den1   = 1.0_RKIND + ((K_mF/Ca_i)**n_F)
      w_i    = 1.0/den1
      den1   = 1.0_RKIND + ((Ca_i/K_mF)**n_F)
      var1   = (1.0_RKIND - 0.234845_RKIND)/den1
      tau    = 4.0e3_RKIND*(0.234845_RKIND + var1)
      Xg(20) = w_i - ((w_i - w)*EXP(-dt/tau))

      RETURN
      END SUBROUTINE TONG_UPDATEG
!------------------------------------------------------
      SUBROUTINE GETRVAL(fileId, skwrd, rVal)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: fileId
      CHARACTER(LEN=*), INTENT(IN) :: skwrd
      REAL(KIND=RKIND), INTENT(INOUT) :: rVal

      INTEGER(KIND=IKIND) :: slen, i, ios
      CHARACTER(LEN=stdL) :: sline, scmd, sval

      REWIND(fileId)
      slen = LEN(TRIM(skwrd))
      DO
         READ(fileId,'(A)',END=001) sline
         sline = ADJUSTC(sline)
         slen  = LEN(TRIM(sline))
         IF (sline(1:1).EQ.'#' .OR. slen.EQ.0) CYCLE

         DO i=1, slen
            IF (sline(i:i) .EQ. ':') EXIT
         END DO

         IF (i .GE. slen) THEN
            STOP "Error: inconsistent input file format"
         END IF

         scmd = sline(1:i-1)
         sval = sline(i+1:slen)
         sval = ADJUSTC(sval)

!        Remove any trailing comments
         slen = LEN(TRIM(sval))
         DO i=1, slen
            IF (sval(i:i) .EQ. '#') EXIT
         END DO
         sval = TRIM(ADJUSTC(sval(1:i-1)))

         IF (TRIM(skwrd) .EQ. TRIM(scmd)) THEN
            READ(sval,*,IOSTAT=ios) rval
            IF (ios .NE. 0) THEN
               STOP " Error: while reading "//TRIM(skwrd)
            END IF
            EXIT
         END IF
      END DO

 001  RETURN

! 001  STOP " Error: EOF reached while finding command <"//
!     2   TRIM(skwrd)//">"

      END SUBROUTINE GETRVAL
!-----------------------------------------------------------------------
!     Removes any leading spaces or tabs
      PURE FUNCTION ADJUSTC(str)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: str
      CHARACTER(LEN=LEN(str)) ADJUSTC

      INTEGER(KIND=IKIND) i

      DO i=1, LEN(str)
         IF (str(i:i) .NE. " " .AND. str(i:i) .NE. "  ") EXIT
      END DO
      IF (i .GT. LEN(str)) THEN
         ADJUSTC = ""
      ELSE
         ADJUSTC = str(i:)
      END IF

      RETURN
      END FUNCTION ADJUSTC
!-----------------------------------------------------------------------
      END MODULE TONGMOD
!#######################################################################

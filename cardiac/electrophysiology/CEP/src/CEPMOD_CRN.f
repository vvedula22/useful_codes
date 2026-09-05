!-----------------------------------------------------------------------
!
!     This module defines data structures for the Courtemanche-Ramirez-
!     Nattel (CRN) model of the human atrial myocyte.
!
!     Reference:
!        Courtemanche M, Ramirez RJ, Nattel S (1998). Ionic mechanisms
!        underlying human atrial action potential properties: insights
!        from a mathematical model. Am J Physiol Heart Circ Physiol,
!        275(1), H301-H321.
!        https://doi.org/10.1152/ajpheart.1998.275.1.H301
!
!     This is a line-for-line transcription of the validated Python
!     reference implementation (ct_model.py); the state split, RPAR
!     layout, singularity guards and the Fn_rel_scale switch are
!     identical.  State variables X(1:6) = V, Na_i, K_i, Ca_i, Ca_up,
!     Ca_rel are integrated by FE/RK4; gating variables Xg(1:15) = m,
!     h, j, d, f, f_Ca, xr, xs, oa, oi, ua, ui, u, v, w by the
!     Rush-Larsen exponential update.
!
!     UNITS: millisecond clock, currents in ABSOLUTE pA (the paper's
!     densities multiplied by Cm = 100 pF at use), so
!     dV/dt = -(SUM I)/Cm is in pA/pF = mV/ms.  The stimulus amplitude
!     is therefore absolute pA, not the pA/pF density a TTP deck uses.
!     See the UNIT SYSTEM block in PARAMS_CRN.f, and
!     example/README_CEP_UNITS.md for the cross-model table.
!
!     CAUTION: this module should not be converted to TTP-style
!     densities.  Eq. 68's Fn thresholds require absolute pA, and a
!     density formulation silently disables SR release.
!
!-----------------------------------------------------------------------

      MODULE CRNMOD
      USE TYPEMOD
      USE UTILMOD, ONLY : stdL
      IMPLICIT NONE

      PRIVATE

      REAL(KIND=RKIND), PARAMETER :: eps = EPSILON(eps)

!     Tolerance for detecting the removable singularities of the rate
!     expressions (appendix: "Some fractional equations require
!     evaluation of a limit to determine their values at membrane
!     potentials for which their denominator is zero.")
      REAL(KIND=RKIND), PARAMETER :: stol = 1.E-10_RKIND

      INCLUDE "PARAMS_CRN.f"

      PUBLIC :: CRN_INIT
      PUBLIC :: CRN_READPARFF
      PUBLIC :: CRN_INTEGFE
      PUBLIC :: CRN_INTEGRK

      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE CRN_INIT(nX, nG, X, Xg)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX, nG
      REAL(KIND=RKIND), INTENT(OUT) :: X(nX), Xg(nG)

!     Initial conditions: the model at rest (paper Table 2)
!     State variables
      X(1)   = -81.18_RKIND       ! V       (units: mV)
      X(2)   =  1.117E1_RKIND     ! Na_i    (units: mM)
      X(3)   =  1.39E2_RKIND      ! K_i     (units: mM)
      X(4)   =  1.013E-4_RKIND    ! Ca_i    (units: mM)
      X(5)   =  1.488_RKIND       ! Ca_up   (units: mM)
      X(6)   =  1.488_RKIND       ! Ca_rel  (units: mM)

!     Gating variables
      Xg(1)  =  2.908E-3_RKIND    ! m       (dimensionless)
      Xg(2)  =  9.649E-1_RKIND    ! h       (dimensionless)
      Xg(3)  =  9.775E-1_RKIND    ! j       (dimensionless)
      Xg(4)  =  1.367E-4_RKIND    ! d       (dimensionless)
      Xg(5)  =  9.996E-1_RKIND    ! f       (dimensionless)
      Xg(6)  =  7.755E-1_RKIND    ! f_Ca    (dimensionless)
      Xg(7)  =  3.296E-5_RKIND    ! xr      (dimensionless)
      Xg(8)  =  1.869E-2_RKIND    ! xs      (dimensionless)
      Xg(9)  =  3.043E-2_RKIND    ! oa      (dimensionless)
      Xg(10) =  9.992E-1_RKIND    ! oi      (dimensionless)
      Xg(11) =  4.966E-3_RKIND    ! ua      (dimensionless)
      Xg(12) =  9.986E-1_RKIND    ! ui      (dimensionless)
      Xg(13) =  0._RKIND          ! u       (dimensionless)
      Xg(14) =  1._RKIND          ! v       (dimensionless)
      Xg(15) =  9.992E-1_RKIND    ! w       (dimensionless)

      RETURN
      END SUBROUTINE CRN_INIT
!-----------------------------------------------------------------------
      SUBROUTINE CRN_READPARFF(fname)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: fname

      INTEGER fid

      fid = 1652

      OPEN(fid, FILE=TRIM(fname))
      CALL GETRVAL(fid, "R", R)
      CALL GETRVAL(fid, "T", T)
      CALL GETRVAL(fid, "F", F)
      CALL GETRVAL(fid, "Cm", Cm)
      CALL GETRVAL(fid, "Vol_cell", Vol_cell)
      CALL GETRVAL(fid, "Vol_i", Vol_i)
      CALL GETRVAL(fid, "Vol_up", Vol_up)
      CALL GETRVAL(fid, "Vol_rel", Vol_rel)
      CALL GETRVAL(fid, "K_o", K_o)
      CALL GETRVAL(fid, "Na_o", Na_o)
      CALL GETRVAL(fid, "Ca_o", Ca_o)

      CALL GETRVAL(fid, "g_Na", g_Na)
      CALL GETRVAL(fid, "g_K1", g_K1)
      CALL GETRVAL(fid, "g_to", g_to)
      CALL GETRVAL(fid, "g_Kr", g_Kr)
      CALL GETRVAL(fid, "g_Ks", g_Ks)
      CALL GETRVAL(fid, "g_CaL", g_CaL)
      CALL GETRVAL(fid, "g_bCa", g_bCa)
      CALL GETRVAL(fid, "g_bNa", g_bNa)

      CALL GETRVAL(fid, "I_NaK_max", I_NaK_max)
      CALL GETRVAL(fid, "I_NaCa_max", I_NaCa_max)
      CALL GETRVAL(fid, "I_pCa_max", I_pCa_max)
      CALL GETRVAL(fid, "Km_Na_i", Km_Na_i)
      CALL GETRVAL(fid, "Km_K_o", Km_K_o)
      CALL GETRVAL(fid, "Km_Na", Km_Na)
      CALL GETRVAL(fid, "Km_Ca", Km_Ca)
      CALL GETRVAL(fid, "k_sat", k_sat)
      CALL GETRVAL(fid, "gamma", gamma)

      CALL GETRVAL(fid, "I_up_max", I_up_max)
      CALL GETRVAL(fid, "K_up", K_up)
      CALL GETRVAL(fid, "Ca_up_max", Ca_up_max)
      CALL GETRVAL(fid, "k_rel", k_rel)
      CALL GETRVAL(fid, "tau_tr", tau_tr)

      CALL GETRVAL(fid, "Cmdn_max", Cmdn_max)
      CALL GETRVAL(fid, "Trpn_max", Trpn_max)
      CALL GETRVAL(fid, "Csqn_max", Csqn_max)
      CALL GETRVAL(fid, "Km_Cmdn", Km_Cmdn)
      CALL GETRVAL(fid, "Km_Trpn", Km_Trpn)
      CALL GETRVAL(fid, "Km_Csqn", Km_Csqn)

      CALL GETRVAL(fid, "K_Q10", K_Q10)
      CALL GETRVAL(fid, "Fn_rel_scale", Fn_rel_scale)
      CALL GETRVAL(fid, "Vrest", Vrest)

!     Scaling factors
      CALL GETRVAL(fid, "Vscale", Vscale)
      CALL GETRVAL(fid, "Tscale", Tscale)
      CALL GETRVAL(fid, "Voffset", Voffset)

      CLOSE(fid)

      RETURN
      END SUBROUTINE CRN_READPARFF
!-----------------------------------------------------------------------
!     Time integration performed using Forward Euler method
      SUBROUTINE CRN_INTEGFE(nX, nG, X, Xg, dt, Istim, Ksac, RPAR)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX, nG
      REAL(KIND=RKIND), INTENT(INOUT) :: X(nX), Xg(nG), RPAR(21)
      REAL(KIND=RKIND), INTENT(IN) :: dt, Istim, Ksac

      REAL(KIND=RKIND) :: f(nX)

!     Get time derivatives (RHS); also computes Fn for CRN_UPDATEG
      CALL CRN_GETF(nX, nG, X, Xg, f, Istim, Ksac, RPAR)

!     Update gating variables
      CALL CRN_UPDATEG(dt, nX, nG, X, Xg)

!     Update state variables
      X = X + dt*f(:)

      RETURN
      END SUBROUTINE CRN_INTEGFE
!-----------------------------------------------------------------------
!     Time integration performed using 4th order Runge-Kutta method
      SUBROUTINE CRN_INTEGRK(nX, nG, X, Xg, dt, Istim, Ksac, RPAR)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX, nG
      REAL(KIND=RKIND), INTENT(INOUT) :: X(nX), Xg(nG), RPAR(21)
      REAL(KIND=RKIND), INTENT(IN) :: dt, Istim, Ksac

      REAL(KIND=RKIND) :: dt6, Xrk(nX), Xgr(nG), frk(nX,4)

      dt6 = dt/6._RKIND

!     RK4: 1st pass
      Xrk = X
      CALL CRN_GETF(nX, nG, Xrk, Xg, frk(:,1), Istim, Ksac, RPAR)

!     Update gating variables by half-dt
      Xgr = Xg
      CALL CRN_UPDATEG(0.5_RKIND*dt, nX, nG, X, Xgr)

!     RK4: 2nd pass
      Xrk = X + 0.5_RKIND*dt*frk(:,1)
      CALL CRN_GETF(nX, nG, Xrk, Xgr, frk(:,2), Istim, Ksac, RPAR)

!     RK4: 3rd pass
      Xrk = X + 0.5_RKIND*dt*frk(:,2)
      CALL CRN_GETF(nX, nG, Xrk, Xgr, frk(:,3), Istim, Ksac, RPAR)

!     Update gating variables by full-dt
      Xgr = Xg
      CALL CRN_UPDATEG(dt, nX, nG, X, Xgr)

!     RK4: 4th pass
      Xrk = X + dt*frk(:,3)
      CALL CRN_GETF(nX, nG, Xrk, Xgr, frk(:,4), Istim, Ksac, RPAR)

      X = X + dt6*(frk(:,1) + 2._RKIND*(frk(:,2) + frk(:,3)) + frk(:,4))
      Xg = Xgr

      RETURN
      END SUBROUTINE CRN_INTEGRK
!-----------------------------------------------------------------------
!     Compute currents and time derivatives of state variables
      SUBROUTINE CRN_GETF(nX, nG, X, Xg, dX, I_stim, K_sac, RPAR)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX, nG
      REAL(KIND=RKIND), INTENT(IN) :: X(nX), Xg(nG), I_stim, K_sac
      REAL(KIND=RKIND), INTENT(OUT) :: dX(nX)
      REAL(KIND=RKIND), INTENT(INOUT) :: RPAR(21)

      REAL(KIND=RKIND) :: RT, gKurv, sgma, fNaK, e1, e2, n1, d1, B1,
     2   B2, I_ion

!     Local copies of state variables
      V      = X(1)
      Na_i   = X(2)
      K_i    = X(3)
      Ca_i   = X(4)
      Ca_up  = X(5)
      Ca_rel = X(6)

!     Local copies of gating variables
      gm     = Xg(1)
      gh     = Xg(2)
      gj     = Xg(3)
      gd     = Xg(4)
      gf     = Xg(5)
      gfca   = Xg(6)
      gxr    = Xg(7)
      gxs    = Xg(8)
      goa    = Xg(9)
      goi    = Xg(10)
      gua    = Xg(11)
      gui    = Xg(12)
      gu     = Xg(13)
      gv     = Xg(14)
      gw     = Xg(15)

!     Stretch-activated current
      I_sac  = K_sac * (Vrest - V)

!     R*T/F (mV)
      RT     = R * T / F

!     Equilibrium potentials (Eq. 28)
      E_Na   = RT * LOG(Na_o/Na_i)
      E_K    = RT * LOG(K_o/K_i)
      E_Ca   = 0.5_RKIND * RT * LOG(Ca_o/Ca_i)

!     I_Na: fast Na current (Eq. 29)
      I_Na   = Cm * g_Na * gm*gm*gm * gh * gj * (V - E_Na)

!     I_K1: inward rectifier K current (Eq. 35)
      I_K1   = Cm * g_K1 * (V - E_K)
     2       / (1._RKIND + EXP(0.07_RKIND*(V + 80._RKIND)))

!     I_to: transient outward K current (Eq. 36)
      I_to   = Cm * g_to * goa*goa*goa * goi * (V - E_K)

!     I_Kur: ultrarapid delayed rectifier K current (Eqs. 41, 42).
!     Note g_Kur is voltage-dependent and defined inline in the paper
!     (Eq. 42), not in Table 1.
      gKurv  = 0.005_RKIND + 0.05_RKIND
     2       / (1._RKIND + EXP(-(V - 15._RKIND)/13._RKIND))
      I_Kur  = Cm * gKurv * gua*gua*gua * gui * (V - E_K)

!     I_Kr: rapid delayed rectifier K current (Eq. 47)
      I_Kr   = Cm * g_Kr * gxr * (V - E_K)
     2       / (1._RKIND + EXP((V + 15._RKIND)/22.4_RKIND))

!     I_Ks: slow delayed rectifier K current (Eq. 50)
      I_Ks   = Cm * g_Ks * gxs*gxs * (V - E_K)

!     I_CaL: L-type Ca current (Eq. 53).  Note the fixed 65 mV
!     apparent reversal potential, not E_Ca.
      I_CaL  = Cm * g_CaL * gd * gf * gfca * (V - 65._RKIND)

!     I_NaK: Na-K pump current (Eqs. 57-59)
      sgma   = (EXP(Na_o/67.3_RKIND) - 1._RKIND) / 7._RKIND
      fNaK   = 1._RKIND / (1._RKIND
     2       + 0.1245_RKIND*EXP(-0.1_RKIND*F*V/(R*T))
     3       + 0.0365_RKIND*sgma*EXP(-F*V/(R*T)))
      I_NaK  = Cm * I_NaK_max * fNaK
     2       / (1._RKIND + (Km_Na_i/Na_i)**1.5_RKIND)
     3       * K_o / (K_o + Km_K_o)

!     I_NaCa: Na-Ca exchanger current (Eq. 60)
      e1     = EXP(gamma*F*V/(R*T))
      e2     = EXP((gamma - 1._RKIND)*F*V/(R*T))
      n1     = e1*Na_i*Na_i*Na_i*Ca_o - e2*Na_o*Na_o*Na_o*Ca_i
      d1     = (Km_Na**3 + Na_o**3) * (Km_Ca + Ca_o)
     2       * (1._RKIND + k_sat*e2)
      I_NaCa = Cm * I_NaCa_max * n1 / d1

!     Background currents (Eqs. 61, 62)
      I_bCa  = Cm * g_bCa * (V - E_Ca)
      I_bNa  = Cm * g_bNa * (V - E_Na)

!     I_pCa: sarcolemmal Ca pump current (Eq. 63)
      I_pCa  = Cm * I_pCa_max * Ca_i / (0.0005_RKIND + Ca_i)

!     SR fluxes, in mM/ms (NOT pA)
!     I_rel: Ca release from the JSR (Eq. 64)
      I_rel  = k_rel * gu*gu * gv * gw * (Ca_rel - Ca_i)

!     I_tr: NSR -> JSR transfer (Eqs. 69, 70)
      I_tr   = (Ca_up - Ca_rel) / tau_tr

!     I_up: Ca uptake into the NSR (Eq. 71) and its leak (Eq. 72)
      I_up   = I_up_max / (1._RKIND + K_up/Ca_i)
      I_up_leak = Ca_up / Ca_up_max * I_up_max

!     Fn: Ca release flux trigger (Eq. 68).  Fn_rel_scale scales the
!     RELEASE term only: 1.0 = Eq. 68 as published, 0.01 (the default)
!     reproduces the paper's Fig. 15.  See PARAMS_CRN.f.  The currents
!     here are absolute pA, deliberately.
      Fn     = Fn_rel_scale * 1.E-12_RKIND * Vol_rel * I_rel
     2       - (5.E-13_RKIND/F)
     3       * (0.5_RKIND*I_CaL - 0.2_RKIND*I_NaCa)

!-----------------------------------------------------------------------
!     Time derivatives of the state variables
!     dV/dt (Eqs. 18, 19)
      I_ion  = I_Na + I_K1 + I_to + I_Kur + I_Kr + I_Ks + I_CaL
     2       + I_pCa + I_NaK + I_NaCa + I_bNa + I_bCa
      dX(1)  = -(I_ion + I_stim)/Cm + I_sac

!     dNa_i/dt (Eq. 21)
      dX(2)  = (-3._RKIND*I_NaK - 3._RKIND*I_NaCa - I_bNa - I_Na)
     2       / (F*Vol_i)

!     dK_i/dt (Eq. 22).  I_b,K is identically zero in this model.
      dX(3)  = (2._RKIND*I_NaK - I_K1 - I_to - I_Kur - I_Kr - I_Ks)
     2       / (F*Vol_i)

!     dCa_i/dt (Eqs. 23-25), with instantaneous buffering
      B1     = (2._RKIND*I_NaCa - I_pCa - I_CaL - I_bCa)
     2       / (2._RKIND*F*Vol_i)
     3       + (Vol_up*(I_up_leak - I_up) + I_rel*Vol_rel) / Vol_i
      n1     = Ca_i + Km_Trpn
      d1     = Ca_i + Km_Cmdn
      B2     = 1._RKIND + Trpn_max*Km_Trpn/(n1*n1)
     2       + Cmdn_max*Km_Cmdn/(d1*d1)
      dX(4)  = B1 / B2

!     dCa_up/dt (Eq. 26)
      dX(5)  = I_up - I_up_leak - I_tr*Vol_rel/Vol_up

!     dCa_rel/dt (Eq. 27), with instantaneous calsequestrin buffering
      n1     = Ca_rel + Km_Csqn
      dX(6)  = (I_tr - I_rel)
     2       / (1._RKIND + Csqn_max*Km_Csqn/(n1*n1))

!     Quantities to be written to file
      RPAR(3)  = I_Na
      RPAR(4)  = I_K1
      RPAR(5)  = I_to
      RPAR(6)  = I_Kur
      RPAR(7)  = I_Kr
      RPAR(8)  = I_Ks
      RPAR(9)  = I_CaL
      RPAR(10) = I_pCa
      RPAR(11) = I_NaK
      RPAR(12) = I_NaCa
      RPAR(13) = I_bNa
      RPAR(14) = I_bCa
      RPAR(15) = I_rel
      RPAR(16) = I_tr
      RPAR(17) = I_up
      RPAR(18) = I_up_leak
      RPAR(19) = Fn
      RPAR(20) = I_stim
      RPAR(21) = I_sac

      RETURN
      END SUBROUTINE CRN_GETF
!-----------------------------------------------------------------------
!     Update all the gating variables (Rush-Larsen, Eq. 77).  Reads
!     X(1) = V and X(4) = Ca_i, plus the module-scope Fn set by the
!     preceding CRN_GETF call.
      SUBROUTINE CRN_UPDATEG(dt, n, nG, X, Xg)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: n, nG
      REAL(KIND=RKIND), INTENT(IN) :: dt, X(n)
      REAL(KIND=RKIND), INTENT(INOUT) :: Xg(nG)

      REAL(KIND=RKIND) :: a, b, am, bm, tau, xinf, e1

      V     = X(1)
      Ca_i  = X(4)

!     m: I_Na activation (Eqs. 30, 34)
      a     = V + 47.13_RKIND
      IF (ABS(a) .LT. stol) THEN
         am = 3.2_RKIND
      ELSE
         am = 0.32_RKIND*a / (1._RKIND - EXP(-0.1_RKIND*a))
      END IF
      bm    = 0.08_RKIND * EXP(-V/11._RKIND)
      tau   = 1._RKIND / (am + bm)
      xinf  = am * tau
      Xg(1) = xinf + (Xg(1) - xinf)*EXP(-dt/tau)

!     h: I_Na fast inactivation (Eqs. 31, 34)
      IF (V .LT. -40._RKIND) THEN
         am = 0.135_RKIND * EXP(-(V + 80._RKIND)/6.8_RKIND)
         bm = 3.56_RKIND*EXP(0.079_RKIND*V)
     2      + 3.1E5_RKIND*EXP(0.35_RKIND*V)
      ELSE
         am = 0._RKIND
         bm = 1._RKIND / (0.13_RKIND*(1._RKIND
     2      + EXP(-(V + 10.66_RKIND)/11.1_RKIND)))
      END IF
      tau   = 1._RKIND / (am + bm)
      xinf  = am * tau
      Xg(2) = xinf + (Xg(2) - xinf)*EXP(-dt/tau)

!     j: I_Na slow inactivation (Eqs. 32, 33, 34)
      IF (V .LT. -40._RKIND) THEN
         am = (-1.2714E5_RKIND*EXP(0.2444_RKIND*V)
     2      - 3.474E-5_RKIND*EXP(-0.04391_RKIND*V))
     3      * (V + 37.78_RKIND)
     4      / (1._RKIND + EXP(0.311_RKIND*(V + 79.23_RKIND)))
         bm = 0.1212_RKIND*EXP(-0.01052_RKIND*V)
     2      / (1._RKIND + EXP(-0.1378_RKIND*(V + 40.14_RKIND)))
      ELSE
         am = 0._RKIND
         bm = 0.3_RKIND*EXP(-2.535E-7_RKIND*V)
     2      / (1._RKIND + EXP(-0.1_RKIND*(V + 32._RKIND)))
      END IF
      tau   = 1._RKIND / (am + bm)
      xinf  = am * tau
      Xg(3) = xinf + (Xg(3) - xinf)*EXP(-dt/tau)

!     d: I_Ca,L activation (Eq. 54)
      a     = V + 10._RKIND
      IF (ABS(a) .LT. stol) THEN
         tau = 1._RKIND / (0.035_RKIND*6.24_RKIND*2._RKIND)
      ELSE
         tau = (1._RKIND - EXP(-a/6.24_RKIND))
     2       / (0.035_RKIND*a*(1._RKIND + EXP(-a/6.24_RKIND)))
      END IF
      xinf  = 1._RKIND / (1._RKIND + EXP(-a/8._RKIND))
      Xg(4) = xinf + (Xg(4) - xinf)*EXP(-dt/tau)

!     f: I_Ca,L voltage inactivation (Eq. 55)
      a     = V + 10._RKIND
      tau   = 9._RKIND / (0.0197_RKIND
     2      * EXP(-(0.0337_RKIND*0.0337_RKIND)*a*a) + 0.02_RKIND)
      xinf  = 1._RKIND / (1._RKIND + EXP((V + 28._RKIND)/6.9_RKIND))
      Xg(5) = xinf + (Xg(5) - xinf)*EXP(-dt/tau)

!     f_Ca: I_Ca,L Ca-dependent inactivation (Eq. 56)
      tau   = 2._RKIND
      xinf  = 1._RKIND / (1._RKIND + Ca_i/0.00035_RKIND)
      Xg(6) = xinf + (Xg(6) - xinf)*EXP(-dt/tau)

!     xr: I_Kr activation (Eqs. 48, 49)
      a     = V + 14.1_RKIND
      IF (ABS(a) .LT. stol) THEN
         am = 0.0003_RKIND*5._RKIND
      ELSE
         am = 0.0003_RKIND*a / (1._RKIND - EXP(-a/5._RKIND))
      END IF
      b     = V - 3.3328_RKIND
      IF (ABS(b) .LT. stol) THEN
         bm = 7.3898E-5_RKIND*5.1237_RKIND
      ELSE
         bm = 7.3898E-5_RKIND*b / (EXP(b/5.1237_RKIND) - 1._RKIND)
      END IF
      tau   = 1._RKIND / (am + bm)
      xinf  = 1._RKIND / (1._RKIND + EXP(-a/6.5_RKIND))
      Xg(7) = xinf + (Xg(7) - xinf)*EXP(-dt/tau)

!     xs: I_Ks activation (Eqs. 51, 52)
      a     = V - 19.9_RKIND
      IF (ABS(a) .LT. stol) THEN
         am = 4.E-5_RKIND*17._RKIND
         bm = 3.5E-5_RKIND*9._RKIND
      ELSE
         am = 4.E-5_RKIND*a / (1._RKIND - EXP(-a/17._RKIND))
         bm = 3.5E-5_RKIND*a / (EXP(a/9._RKIND) - 1._RKIND)
      END IF
      tau   = 0.5_RKIND / (am + bm)
      xinf  = (1._RKIND + EXP(-a/12.7_RKIND))**(-0.5_RKIND)
      Xg(8) = xinf + (Xg(8) - xinf)*EXP(-dt/tau)

!     oa: I_to activation (Eqs. 37, 38)
      am    = 0.65_RKIND / (EXP(-(V + 10._RKIND)/8.5_RKIND)
     2      + EXP(-(V - 30._RKIND)/59._RKIND))
      bm    = 0.65_RKIND / (2.5_RKIND
     2      + EXP((V + 82._RKIND)/17._RKIND))
      tau   = 1._RKIND / (am + bm) / K_Q10
      xinf  = 1._RKIND / (1._RKIND
     2      + EXP(-(V + 20.47_RKIND)/17.54_RKIND))
      Xg(9) = xinf + (Xg(9) - xinf)*EXP(-dt/tau)

!     oi: I_to inactivation (Eqs. 39, 40)
      am    = 1._RKIND / (18.53_RKIND
     2      + EXP((V + 113.7_RKIND)/10.95_RKIND))
      bm    = 1._RKIND / (35.56_RKIND
     2      + EXP(-(V + 1.26_RKIND)/7.44_RKIND))
      tau   = 1._RKIND / (am + bm) / K_Q10
      xinf  = 1._RKIND / (1._RKIND + EXP((V + 43.1_RKIND)/5.3_RKIND))
      Xg(10) = xinf + (Xg(10) - xinf)*EXP(-dt/tau)

!     ua: I_Kur activation (Eqs. 43, 44)
      am    = 0.65_RKIND / (EXP(-(V + 10._RKIND)/8.5_RKIND)
     2      + EXP(-(V - 30._RKIND)/59._RKIND))
      bm    = 0.65_RKIND / (2.5_RKIND
     2      + EXP((V + 82._RKIND)/17._RKIND))
      tau   = 1._RKIND / (am + bm) / K_Q10
      xinf  = 1._RKIND / (1._RKIND
     2      + EXP(-(V + 30.3_RKIND)/9.6_RKIND))
      Xg(11) = xinf + (Xg(11) - xinf)*EXP(-dt/tau)

!     ui: I_Kur inactivation (Eqs. 45, 46)
      am    = 1._RKIND / (21._RKIND
     2      + EXP(-(V - 185._RKIND)/28._RKIND))
      bm    = EXP((V - 158._RKIND)/16._RKIND)
      tau   = 1._RKIND / (am + bm) / K_Q10
      xinf  = 1._RKIND / (1._RKIND
     2      + EXP((V - 99.45_RKIND)/27.48_RKIND))
      Xg(12) = xinf + (Xg(12) - xinf)*EXP(-dt/tau)

!     u: I_rel activation (Eq. 65).  EXPIT is the overflow-safe
!     logistic: the sigmoid slope factor is 1.367E-15, so the argument
!     reaches ~1E3 in normal operation and a naive EXP overflows.
      tau   = 8._RKIND
      xinf  = EXPIT((Fn - 3.4175E-13_RKIND)/13.67E-16_RKIND)
      Xg(13) = xinf + (Xg(13) - xinf)*EXP(-dt/tau)

!     v: I_rel Ca-flux inactivation (Eq. 66)
      tau   = 1.91_RKIND + 2.09_RKIND
     2      * EXPIT((Fn - 3.4175E-13_RKIND)/13.67E-16_RKIND)
      xinf  = 1._RKIND
     2      - EXPIT((Fn - 6.835E-14_RKIND)/13.67E-16_RKIND)
      Xg(14) = xinf + (Xg(14) - xinf)*EXP(-dt/tau)

!     w: I_rel voltage inactivation (Eq. 67)
      a     = V - 7.9_RKIND
      IF (ABS(a) .LT. stol) THEN
         tau = 6._RKIND*0.2_RKIND/1.3_RKIND
      ELSE
         e1  = EXP(-a/5._RKIND)
         tau = 6._RKIND*(1._RKIND - e1)
     2       / ((1._RKIND + 0.3_RKIND*e1)*a)
      END IF
      xinf  = 1._RKIND - 1._RKIND
     2      / (1._RKIND + EXP(-(V - 40._RKIND)/17._RKIND))
      Xg(15) = xinf + (Xg(15) - xinf)*EXP(-dt/tau)

      RETURN
      END SUBROUTINE CRN_UPDATEG
!-----------------------------------------------------------------------
!     1/(1+EXP(-z)) evaluated without overflow.  Needed only for the
!     Fn-driven gates u and v (Eqs. 65, 66); mathematically identical
!     to the naive form inside the normal range and matches the Python
!     reference implementation bit for bit.
      PURE FUNCTION EXPIT(z) RESULT(res)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: z
      REAL(KIND=RKIND) res, e

      IF (z .GE. 0._RKIND) THEN
         res = 1._RKIND / (1._RKIND + EXP(-z))
      ELSE
         e   = EXP(z)
         res = e / (1._RKIND + e)
      END IF

      RETURN
      END FUNCTION EXPIT
!-----------------------------------------------------------------------
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
               WRITE(*,'(A)') " Error: while reading "//TRIM(skwrd)
               STOP " Error in input file.. aborting"
            END IF
            EXIT
         END IF
      END DO

 001  RETURN

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
      END MODULE CRNMOD
!#######################################################################

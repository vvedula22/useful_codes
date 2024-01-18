!-----------------------------------------------------------------------
!
!     This module defines data structures for the Nygren atrial myocyte
!     activation model for electrophysiology. A phenomenological
!     active stress model is used for excitation-contraction coupling.
!
!     Reference for Nygren model:
!        Nygren A, Fiset C, Firek L, Clark JW, Lindblad DS, Clark RB,
!        Giles WR (1998). Mathematical model of an adult human atrial
!        cell: the role of K+ currents in repolarization. Circulation
!        Research, 82(1). https://pubmed.ncbi.nlm.nih.gov/9440706/
!
!-----------------------------------------------------------------------

      MODULE NYGMOD
      USE TYPEMOD
      USE UTILMOD, ONLY : stdL
      IMPLICIT NONE

      PRIVATE

      REAL(KIND=RKIND), PARAMETER :: eps = EPSILON(eps)

      INCLUDE "PARAMS_NYG.f"

      PUBLIC :: NYG_INIT
      PUBLIC :: NYG_READPARFF
      PUBLIC :: NYG_INTEGFE
      PUBLIC :: NYG_INTEGRK

      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE NYG_INIT(nX, nG, X, Xg)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX, nG
      REAL(KIND=RKIND), INTENT(OUT) :: X(nX), Xg(nG)

      X(:)   = 0._RKIND

!     Initialize state variables
      X(1)   = -74.2525_RKIND     ! V       (units: mV)
      X(2)   =  130.011_RKIND     ! Na_c    (units: mmol/L)
      X(3)   =  5.3581_RKIND      ! K_c     (units: mmol/L)
      X(4)   =  1.8147_RKIND      ! Ca_c    (units: mmol/L)
      X(5)   =  8.5547_RKIND      ! Na_i    (units: mmol/L)
      X(6)   =  129.4350_RKIND    ! K_i     (units: mmol/L)
      X(7)   =  6.729E-5_RKIND    ! Ca_i    (units: mmol/L)
      X(8)   =  7.2495E-5_RKIND   ! Ca_d    (units: mmol/L)
      X(9)   =  0.6646_RKIND      ! Ca_up   (units: mmol/L)
      X(10)  =  0.6465_RKIND      ! Ca_rel  (units: mmol/L)
      X(11)  =  0.4284_RKIND      ! F_1     (dimensionless)
      X(12)  =  0.0028_RKIND      ! F_2     (dimensionless)

!     Initialize occupancy state variables
      X(13)  =  0.0275_RKIND      ! O_C     (dimensionless)
      X(14)  =  0.0133_RKIND      ! O_TC    (dimensionless)
      X(15)  =  0.1961_RKIND      ! O_TMgC  (dimensionless)
      X(16)  =  0.7094_RKIND      ! O_TMgMg (dimensionless)
      X(17)  =  0.4369_RKIND      ! O_Calse (dimensionless)

!     Initialize gating variables
      Xg(1)  =  3.2017E-3_RKIND   ! m       (dimensionless)
      Xg(2)  =  0.8814_RKIND      ! h_1     (dimensionless)
      Xg(3)  =  0.8742_RKIND      ! h_2     (dimensionless)
      Xg(4)  =  1.3005E-5_RKIND   ! d_L     (dimensionless)
      Xg(5)  =  0.9986_RKIND      ! f_L1    (dimensionless)
      Xg(6)  =  0.9986_RKIND      ! f_L2    (dimensionless)
      Xg(7)  =  1.0678E-3_RKIND   ! r       (dimensionless)
      Xg(8)  =  0.949_RKIND       ! s       (dimensionless)
      Xg(9)  =  1.5949E-4_RKIND   ! r_sus   (dimensionless)
      Xg(10) =  0.9912_RKIND      ! s_sus   (dimensionless)
      Xg(11) =  4.8357E-3_RKIND   ! n       (dimensionless)
      Xg(12) =  0.0001_RKIND      ! p_a     (dimensionless)

      RETURN
      END SUBROUTINE NYG_INIT
!-----------------------------------------------------------------------
      SUBROUTINE NYG_READPARFF(fname)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: fname

      INTEGER fid

      fid = 1528

      OPEN(fid, FILE=TRIM(fname))
      CALL GETRVAL(fid, "Na_b", Na_b)
      CALL GETRVAL(fid, "K_b", K_b)
      CALL GETRVAL(fid, "Ca_b", Ca_b)
      CALL GETRVAL(fid, "Mg_i", Mg_i)
      CALL GETRVAL(fid, "E_Ca_app", E_Ca_app)
      CALL GETRVAL(fid, "K_Ca", K_Ca)
      CALL GETRVAL(fid, "R", R)
      CALL GETRVAL(fid, "T", T)
      CALL GETRVAL(fid, "F", F)
      CALL GETRVAL(fid, "Cm", Cm)
      CALL GETRVAL(fid, "Vol_i", Vol_i)
      CALL GETRVAL(fid, "Vol_c", Vol_c)
      CALL GETRVAL(fid, "Vol_d", Vol_d)
      CALL GETRVAL(fid, "Vol_rel", Vol_rel)
      CALL GETRVAL(fid, "Vol_up", Vol_up)
      CALL GETRVAL(fid, "tau_Na", tau_Na)
      CALL GETRVAL(fid, "tau_K", tau_K)
      CALL GETRVAL(fid, "tau_Ca", tau_Ca)
      CALL GETRVAL(fid, "tau_di", tau_di)
      CALL GETRVAL(fid, "Ibar_NaK", Ibar_NaK)
      CALL GETRVAL(fid, "K_NaK_K", K_NaK_K)
      CALL GETRVAL(fid, "K_NaK_Na", K_NaK_Na)
      CALL GETRVAL(fid, "Ibar_CaP", Ibar_CaP)
      CALL GETRVAL(fid, "K_CaP", K_CaP)
      CALL GETRVAL(fid, "K_NaCa", K_NaCa)
      CALL GETRVAL(fid, "gamma", gamma)
      CALL GETRVAL(fid, "d_NaCa", d_NaCa)
      CALL GETRVAL(fid, "phi_Na_en", phi_Na_en)
      CALL GETRVAL(fid, "Ibar_up", Ibar_up)
      CALL GETRVAL(fid, "K_cyca", K_cyca)
      CALL GETRVAL(fid, "K_srca", K_srca)
      CALL GETRVAL(fid, "K_xcs", K_xcs)
      CALL GETRVAL(fid, "tau_tr", tau_tr)
      CALL GETRVAL(fid, "alpha_rel", alpha_rel)
      CALL GETRVAL(fid, "K_rel_i", K_rel_i)
      CALL GETRVAL(fid, "K_rel_d", K_rel_d)
      CALL GETRVAL(fid, "r_recov", r_recov)

      CALL GETRVAL(fid, "P_Na", P_Na)
      CALL GETRVAL(fid, "Gbar_CaL", Gbar_CaL)
      CALL GETRVAL(fid, "Gbar_t", Gbar_t)
      CALL GETRVAL(fid, "Gbar_sus", Gbar_sus)
      CALL GETRVAL(fid, "Gbar_Ks", Gbar_Ks)
      CALL GETRVAL(fid, "Gbar_Kr", Gbar_Kr)
      CALL GETRVAL(fid, "Gbar_K1", Gbar_K1)
      CALL GETRVAL(fid, "Gbar_BNa", Gbar_BNa)
      CALL GETRVAL(fid, "Gbar_BCa", Gbar_BCa)

!     Scaling factors
      CALL GETRVAL(fid, "Vscale", Vscale)
      CALL GETRVAL(fid, "Tscale", Tscale)
      CALL GETRVAL(fid, "Voffset", Voffset)

      CLOSE(fid)

      RETURN
      END SUBROUTINE NYG_READPARFF
!-----------------------------------------------------------------------
!     Time integration performed using Forward Euler method
      SUBROUTINE NYG_INTEGFE(nX, nG, X, Xg, dt, Istim, Ksac, RPAR)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX, nG
      REAL(KIND=RKIND), INTENT(INOUT) :: X(nX), Xg(nG), RPAR(20)
      REAL(KIND=RKIND), INTENT(IN) :: dt, Istim, Ksac

      REAL(KIND=RKIND) :: f(nX)

!     Get time derivatives (RHS)
      CALL NYG_GETF(nX, nG, X, Xg, f, Istim, Ksac, RPAR)

!     Update gating variables
      CALL NYG_UPDATEG(dt, nX, nG, X, Xg)

!     Update state variables
      X = X + dt*f(:)

      RETURN
      END SUBROUTINE NYG_INTEGFE
!-----------------------------------------------------------------------
!     Time integration performed using 4th order Runge-Kutta method
      SUBROUTINE NYG_INTEGRK(nX, nG, X, Xg, dt, Istim, Ksac, RPAR)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX, nG
      REAL(KIND=RKIND), INTENT(INOUT) :: X(nX), Xg(nG), RPAR(20)
      REAL(KIND=RKIND), INTENT(IN) :: dt, Istim, Ksac

      REAL(KIND=RKIND) :: dt6, Xrk(nX), Xgr(nG), frk(nX,4)

      dt6 = dt/6._RKIND

!     RK4: 1st pass
      Xrk = X
      CALL NYG_GETF(nX, nG, Xrk, Xg, frk(:,1), Istim, Ksac, RPAR)

!     Update gating variables by half-dt
      Xgr = Xg
      CALL NYG_UPDATEG(0.5_RKIND*dt, nX, nG, X, Xgr)

!     RK4: 2nd pass
      Xrk = X + 0.5_RKIND*dt*frk(:,1)
      CALL NYG_GETF(nX, nG, Xrk, Xgr, frk(:,2), Istim, Ksac, RPAR)

!     RK4: 3rd pass
      Xrk = X + 0.5_RKIND*dt*frk(:,2)
      CALL NYG_GETF(nX, nG, Xrk, Xgr, frk(:,3), Istim, Ksac, RPAR)

!     Update gating variables by full-dt
      Xgr = Xg
      CALL NYG_UPDATEG(dt, nX, nG, X, Xgr)

!     RK4: 4th pass
      Xrk = X + dt*frk(:,3)
      CALL NYG_GETF(nX, nG, Xrk, Xgr, frk(:,4), Istim, Ksac, RPAR)

      X = X + dt6*(frk(:,1) + 2._RKIND*(frk(:,2) + frk(:,3)) + frk(:,4))
      Xg = Xgr

      RETURN
      END SUBROUTINE NYG_INTEGRK
!-----------------------------------------------------------------------
!     Compute currents and time derivatives of state variables
      SUBROUTINE NYG_GETF(nX, nG, X, Xg, dX, I_stim, K_sac, RPAR)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX, nG
      REAL(KIND=RKIND), INTENT(IN) :: X(nX), Xg(nG), I_stim, K_sac
      REAL(KIND=RKIND), INTENT(OUT) :: dX(nX)
      REAL(KIND=RKIND), INTENT(INOUT) :: RPAR(20)

      REAL(KIND=RKIND) :: RT, n1, n2, n3, d1, d2, d3, I_sac, dOdt

!     Local copies of state variables
      V       = X(1)
      Na_c    = X(2)
      K_c     = X(3)
      Ca_c    = X(4)
      Na_i    = X(5)
      K_i     = X(6)
      Ca_i    = X(7)
      Ca_d    = X(8)
      Ca_up   = X(9)
      Ca_rel  = X(10)
      F_1     = X(11)
      F_2     = X(12)

!     Local copies of occupancy variables
      O_C     = X(13)
      O_TC    = X(14)
      O_TMgC  = X(15)
      O_TMgMg = X(16)
      O_Calse = X(17)

!     Local copies of gating variables
      gm      = Xg(1)
      gh_1    = Xg(2)
      gh_2    = Xg(3)
      gd_L    = Xg(4)
      gf_L1   = Xg(5)
      gf_L2   = Xg(6)
      gr      = Xg(7)
      gs      = Xg(8)
      gr_sus  = Xg(9)
      gs_sus  = Xg(10)
      gn      = Xg(11)
      gp_a    = Xg(12)

!     Stretch-activated currents
      I_sac = K_sac * (Vrest - V)

!     Ca diffusion surrent from the diffusion-restricted subsarcolemmal
!     space to the cytosol
      I_di = (Ca_d-Ca_i)*2._RKIND*F*Vol_d/tau_di

!     Compute reverse potentials
      RT   = R * T / F
      E_K  = RT * LOG(K_c/K_i)
      E_Na = RT * LOG(Na_c/Na_i)
      E_Ca = 0.5_RKIND * RT * LOG(Ca_c/Ca_i)

!     Some additional local variables
!     f_Ca, p_i
      f_Ca = Ca_d / (Ca_d + K_Ca)
      p_i  = 1._RKIND / (1._RKIND + EXP((V+55._RKIND)/24._RKIND))

!     r_act, r_inact
      n1   = Ca_i / (Ca_i + K_rel_i)
      n2   = Ca_d / (Ca_d + K_rel_d)
      r_act   = 203.8_RKIND*(n1**4._RKIND + n2**4._RKIND)
      r_inact = 33.96_RKIND + 339.6_RKIND*(n1**4._RKIND)

!     Now, compute all the currents
!     I_Na: Na current
      n1    = P_Na*(gm**3._RKIND)
      n2    = 0.9_RKIND*gh_1 + 0.1_RKIND*gh_2
      n3    = Na_c*V*F*(EXP((V-E_Na)/RT) - 1._RKIND)
      d1    = RT*(EXP(V/RT) - 1._RKIND)
      I_Na  = n1*n2*n3/d1

!     I_Ca: L-type Ca current
      n1    = Gbar_CaL * gd_L
      n2    = f_Ca*gf_L1 + (1._RKIND-f_Ca)*gf_L2
      n3    = (V - E_Ca_app)
      I_CaL = n1*n2*n3

!     I_t: transient outward K current
      I_t   = Gbar_t * gr * gs * (V - E_K)

!     I_sus: sustained outward K current
      I_sus = Gbar_sus * gr_sus * gs_sus * (V - E_K)

!     I_Ks: slow delayed rectifier K current
      I_Ks  = Gbar_Ks * gn * (V - E_K)

!     I_Kr: rapid delayed rectifier K current
      I_Kr  = Gbar_Kr * gp_a * p_i * (V - E_K)

!     I_K1: inwardly rectifing K current
      n1    = Gbar_K1 * (K_c**0.4457_RKIND) * (V - E_K)
      d1    = 1._RKIND + EXP(1.5_RKIND*(V-E_K+3.6_RKIND)/RT)
      I_K1  = n1 / d1

!     I_B,Na: background Na current
      I_bNa = Gbar_BNa * (V - E_Na)

!     I_B,Ca: background Ca current
      I_bCa = Gbar_BCa * (V - E_Ca)

!     I_NaK: Na-K pump current
      n1    = Ibar_NaK * K_c * (Na_i**1.5_RKIND) * (V+150._RKIND)
      d1    = K_c + K_NaK_K
      d2    = (Na_i)**1.5_RKIND + (K_NaK_Na)**1.5_RKIND
      d3    = V + 200._RKIND
      I_NaK = n1 / (d1*d2*d3)

!     I_CaP: sarcolemmal Ca pump current
      I_CaP = Ibar_CaP * Ca_i / (Ca_i + K_CaP)

!     I_NaCa: Na-Ca exchange current
      n1 = (Na_i)**3._RKIND * Ca_c * EXP(gamma*V/RT)
      n1 = n1 - ((Na_c)**3._RKIND * Ca_i * EXP((gamma-1._RKIND)*V/RT))
      n1 = K_NaCa * n1
      d1 = ((Na_i)**3._RKIND * Ca_c) + ((Na_c)**3._RKIND * Ca_i)
      d1 = 1._RKIND + d_NaCa*d1
      I_NaCa = n1 / d1

!     I_up: sacroplasmic reticulum Ca uptake current
      n1   = (Ca_i/K_cyca) - (K_xcs*K_xcs*Ca_up/K_srca)
      n1   = Ibar_up * n1
      d1   = (Ca_i + K_cyca)/K_cyca
      d2   = K_xcs*(Ca_up + K_srca)/K_srca
      I_up = n1 / (d1 + d2)

!     I_tr: sacroplasmic reticulum Ca translocation current
      I_tr = (Ca_up - Ca_rel)*2._RKIND*F*Vol_rel/tau_tr

!     I_rel: Ca induced Ca current (CICR)
      n1    = alpha_rel * (Ca_rel - Ca_i)
      n2    = F_2 / (F_2 + 0.25_RKIND)
      I_rel = n1*n2*n2

!-----------------------------------------------------------------------
!     Now compute time derivatives of the state variables:
!     dV/dt: includes stimulus and stretch-activated currents
!      dX(1)  = -(I_Na + I_CaL + I_t + I_sus + I_K1 + I_BNa + I_BCa
!     2       +   I_NaK + I_CaP + I_NaCa + I_stim)/Cm +  I_sac
!     CellML modification for dV/dt
      dX(1) = -(I_Na + I_CaL + I_t + I_sus + I_K1 + I_Kr + I_Ks + I_BNa
     2      +   I_BCa + I_NaK + I_CaP + I_NaCa + I_stim)/Cm +  I_sac

!     dNa_c/dt
      n1     = Na_b - Na_c
      d1     = tau_Na
      n2     = I_Na + I_BNa + 3._RKIND*(I_NaK + I_NaCa) + phi_Na_en
      d2     = Vol_c*F
      dX(2)  = (n1/d1) + (n2/d2)

!     dK_c/dt
      n1     = K_b - K_c
      d1     = tau_K
      n2     = I_t + I_sus + I_K1 + I_Ks + I_Kr - 2._RKIND*I_NaK
      d2     = Vol_c*F
      dX(3)  = (n1/d1) + (n2/d2)

!     dCa_c/dt
      n1     = Ca_b - Ca_c
      d1     = tau_Ca
      n2     = I_CaL + I_BCa + I_CaP - 2._RKIND*I_NaCa
      d2     = 2._RKIND*Vol_c*F
      dX(4)  = (n1/d1) + (n2/d2)

!     dNa_i/dt
      n1     = -(I_Na + I_BNa + 3._RKIND*(I_NaK + I_NaCa) + phi_Na_en)
      d1     = Vol_i*F
      dX(5)  = n1/d1

!     dK_i/dt
      n1     = -(I_t + I_sus + I_K1 + I_Ks + I_Kr - 2._RKIND*I_Nak)
      d1     = Vol_i*F
      dX(6)  = n1/d1

!     dCa_i/dt: will be updated dO/dt term later
      n1     = -(-I_di + I_BCa + I_CaP - 2._RKIND*I_NaCa + I_up - I_rel)
      d1     = 2._RKIND*Vol_i*F
      dX(7)  = n1/d1

!     dCa_d/dt
!      dX(8)  = -(I_CaL-I_di) / (2._RKIND*Vol_d*F)
!     CellML modification for dCa_d/dt
      dX(8)  = -(I_CaL+I_di) / (2._RKIND*Vol_d*F)

!     dCa_up/dt
      dX(9)  = (I_up-I_tr) / (2._RKIND*Vol_up*F)

!     dCa_rel/dt: will be updated with dO_calse/dt later
      dX(10) = (I_tr-I_rel) / (2._RKIND*Vol_rel*F)

!     dF_1/dt
      dX(11)  = r_recov*(1._RKIND - F_1 - F_2) - r_act*F_1

!     dF_2/dt
      dX(12)  = r_act*F_1 - r_inact*F_2

!-----------------------------------------------------------------------
!     Now compute time derivatives of the occupancy variables:
!     dO_C/dt
      dX(13) = 2.E5_RKIND*Ca_i*(1._RKIND-O_C) - 476._RKIND*O_C

!     dO_TC/dt
      dX(14) = 7.84E4_RKIND*Ca_i*(1._RKIND-O_TC) - 392._RKIND*O_TC

!     dO_TMgC/dt
      dX(15) = 2.E5_RKIND*Ca_i*(1._RKIND - O_TMgC - O_TMgMg)
     2       - 6.6_RKIND*O_TMgC

!     dO_TMgMg/dt
      dX(16) = 2.E3_RKIND*Mg_i*(1._RKIND - O_TMgC - O_TMgMg)
     2       - 666._RKIND*O_TMgMg

!     dO_Calse/dt
      dX(17) = 480._RKIND*Ca_rel*(1._RKIND-O_Calse) - 400._RKIND*O_Calse
      dX(10) = dX(10) - 31._RKIND*dX(17)

!     dO/dt
      dOdt  = 0.045_RKIND*dX(13) + 0.08_RKIND*dX(14) + 0.16_RKIND*dX(15)
      dX(7) = dX(7) - dOdt

!     Quantities to be written to file
      RPAR(3)  = I_Na
      RPAR(4)  = I_K1
      RPAR(5)  = I_t
      RPAR(6)  = I_sus
      RPAR(7)  = I_K1
      RPAR(8)  = I_BNa
      RPAR(9)  = I_BCa
      RPAR(10) = I_NaK
      RPAR(11) = I_CaP
      RPAR(12) = I_NaCa
      RPAR(13) = I_Ks
      RPAR(14) = I_Kr
      RPAR(15) = I_di
      RPAR(16) = I_up
      RPAR(17) = I_tr
      RPAR(18) = I_rel
      RPAR(19) = I_stim
      RPAR(20) = I_sac

      RETURN
      END SUBROUTINE NYG_GETF
!-----------------------------------------------------------------------
!     Update all the gating variables
      SUBROUTINE NYG_UPDATEG(dt, n, nG, X, Xg)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: n, nG
      REAL(KIND=RKIND), INTENT(IN) :: dt, X(n)
      REAL(KIND=RKIND), INTENT(INOUT) :: Xg(nG)

      REAL(KIND=RKIND) :: a, e1, e2

      V      = X(1)

      gm     = Xg(1)
      gh_1   = Xg(2)
      gh_2   = Xg(3)
      gd_L   = Xg(4)
      gf_L1  = Xg(5)
      gf_L2  = Xg(6)
      gr     = Xg(7)
      gs     = Xg(8)
      gr_sus = Xg(9)
      gs_sus = Xg(10)
      gn     = Xg(11)
      gp_a   = Xg(12)

!     m: activation gating variable for I_Na
      e1     = EXP(-(V+27.12_RKIND)/8.21_RKIND)
      a      = (V + 25.57_RKIND) / 28.8_RKIND
      e2     = EXP(-a*a)
      mbar   = 1._RKIND / (1._RKIND + e1)
      tau_m  = 4.2E-5_RKIND*e2 + 2.4E-5_RKIND
      Xg(1)  = mbar - (mbar - gm)*EXP(-dt/tau_m)

!     h_1: fast inactivation gating variable for I_Na
      e1     = EXP((V+63.6_RKIND)/5.3_RKIND)
      e2     = EXP((V+35.1_RKIND)/3.2_RKIND)
      hbar_1 = 1._RKIND / (1._RKIND + e1)
      tau_h1 = (0.03_RKIND/(1._RKIND + e2)) + 3.E-4_RKIND
      Xg(2)  = hbar_1 - (hbar_1 - gh_1)*EXP(-dt/tau_h1)

!     h_2: slow inactivation gating variable for I_Na
      hbar_2 = hbar_1
      tau_h2 = (0.12_RKIND/(1._RKIND + e2)) + 3.E-3_RKIND
      Xg(3)  = hbar_2 - (hbar_2 - gh_2)*EXP(-dt/tau_h2)

!     d_L: activation gating variable for I_Ca,L
      e1     = EXP(-(V+9._RKIND)/5.8_RKIND)
      a      = (V+35._RKIND)/30._RKIND
      e2     = EXP(-a*a)
      dbar_L = 1._RKIND / (1._RKIND + e1)
      tau_dL = 2.7E-3_RKIND*e2 + 2.E-3_RKIND
      Xg(4)  = dbar_L - (dbar_L - gd_L)*EXP(-dt/tau_dL)

!     f_L1: fast inactivation gating variable for I_Ca,L
      e1      = EXP((V+27.4_RKIND)/7.1_RKIND)
      a       = (V+40._RKIND)/14.4_RKIND
      e2      = EXP(-a*a)
      fbar_L1 = 1._RKIND / (1._RKIND + e1)
      tau_fL1 = 0.161_RKIND*e2 + 0.01_RKIND
      Xg(5)   = fbar_L1 - (fbar_L1 - gf_L1)*EXP(-dt/tau_fL1)

!     f_L2: slow inactivation gating variable for I_Ca,L
      a       = (V+40._RKIND)/14.2_RKIND
      e2      = EXP(-a*a)
      fbar_L2 = fbar_L1
      tau_fL2 = 1.3323_RKIND*e2 + 0.0626_RKIND
      Xg(6)   = fbar_L2 - (fbar_L2 - gf_L2)*EXP(-dt/tau_fL2)

!     r: activation gating variable for I_t
      e1    = EXP(-(V-1._RKIND)/11._RKIND)
      a     = V/30._RKIND
      e2    = EXP(-a*a)
      rbar  = 1._RKIND / (1._RKIND + e1)
      tau_r = 3.5E-3_RKIND*e2 + 1.5E-3_RKIND
      Xg(7) = rbar - (rbar - gr)*EXP(-dt/tau_r)

!     s: inactivation gating variable for I_t
      e1    = EXP((V+40.5_RKIND)/11.5_RKIND)
      a     = (V+52.45)/14.97_RKIND
      e2    = EXP(-a*a)
      sbar  = 1._RKIND / (1._RKIND + e1)
      tau_s = 0.4812_RKIND*e2 + 0.01414_RKIND
      Xg(8) = sbar - (sbar - gs)*EXP(-dt/tau_s)

!     r_sus: activation gating variable for I_sus
      e1       = EXP(-(V+4.3_RKIND)/8._RKIND)
      e2       = EXP((V+5._RKIND)/12._RKIND)
      rbar_sus = 1._RKIND / (1._RKIND + e1)
      tau_rsus = (9.E-3_RKIND/(1._RKIND+e2)) + 5.E-4_RKIND
      Xg(9)    = rbar_sus - (rbar_sus - gr_sus)*EXP(-dt/tau_rsus)

!     s_sus: inactivation gating variable for I_sus
      e1       = EXP((V+20._RKIND)/10._RKIND)
      e2       = EXP((V+60._RKIND)/10._RKIND)
      sbar_sus = (0.4_RKIND/(1._RKIND + e1)) + 0.6_RKIND
      tau_ssus = (0.047_RKIND/(1._RKIND + e2)) + 0.3_RKIND
      Xg(10)   = sbar_sus - (sbar_sus - gs_sus)*EXP(-dt/tau_ssus)

!     n: activation gating variable for I_K,s
      e1     = EXP(-(V-19.9_RKIND)/12.7_RKIND)
      a      = (V-20._RKIND)/20._RKIND
      e2     = EXP(-a*a)
      nbar   = 1._RKIND / (1._RKIND + e1)
      tau_n  = 0.7_RKIND + 0.4_RKIND*e2
      Xg(11) = nbar - (nbar - gn)*EXP(-dt/tau_n)

!     p_a: activation gating variable for I_K,r
      e1     = EXP(-(V+15._RKIND)/6._RKIND)
      a      = (V+20.1376_RKIND) / 22.1996_RKIND
      e2     = EXP(-a*a)
      pbar_a = 1._RKIND / (1._RKIND + e1)
      tau_pa = 0.03118_RKIND + 0.21718_RKIND*e2
      Xg(12) = pbar_a - (pbar_a - gp_a)*EXP(-dt/tau_pa)

      RETURN
      END SUBROUTINE NYG_UPDATEG
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
               STOP " Error: while reading "//TRIM(skwrd)
            END IF
            EXIT
         END IF
      END DO

 001  RETURN

!  001  STOP " Error: EOF reached while finding command <"//
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
      END MODULE NYGMOD
!#######################################################################

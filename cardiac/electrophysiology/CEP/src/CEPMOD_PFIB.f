!-----------------------------------------------------------------------
!
!     This module defines data structures for the Stewart's
!     cellular activation model for Purkinje fiber cells.
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

      MODULE PFIBMOD
      USE TYPEMOD
      USE UTILMOD, ONLY : stdL
      IMPLICIT NONE

      PRIVATE

      REAL(KIND=RKIND), PARAMETER :: eps = EPSILON(eps)

      INCLUDE "PARAMS_PFIB.f"

      PUBLIC :: PFIB_INIT
      PUBLIC :: PFIB_READPARFF
      PUBLIC :: PFIB_INTEGFE
      PUBLIC :: PFIB_INTEGRK

      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE PFIB_INIT(nX, nG, X, Xg)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX, nG
      REAL(KIND=RKIND), INTENT(OUT) :: X(nX), Xg(nG)

!     Note that the initial values listed here are directly from the
!     original paper. However, these are quite different from CellML
c!     Initialize state variables
c      X(1)   = -74.7890522727_RKIND      ! V      (units: mV)
c      X(2)   = 136.9896086978_RKIND      ! K_i    (units: mM)
c      X(3)   =   8.5447311020_RKIND      ! Na_i   (units: mM)
c      X(4)   =   0.0001720623_RKIND      ! Ca_i   (units: mM)
c      X(5)   =   0.0006146554_RKIND      ! Ca_ss  (units: mM)
c      X(6)   =   3.2830723338_RKIND      ! Ca_sr  (units: mM)
c      X(7)   =   0.8199969443_RKIND      ! R'     (dimensionless)

c!     Initialize gating variables
c      Xg(1)  =   0.4663168269_RKIND      ! x_r1   (dimensionless)
c      Xg(2)  =   0.3657472179_RKIND      ! x_r2   (dimensionless)
c      Xg(3)  =   0.0486609588_RKIND      ! x_s    (dimensionless)
c      Xg(4)  =   0.0145766758_RKIND      ! m      (dimensionless)
c      Xg(5)  =   0.2979720207_RKIND      ! h      (dimensionless)
c      Xg(6)  =   0.0692509548_RKIND      ! j      (dimensionless)
c      Xg(7)  =   0.0001356656_RKIND      ! d      (dimensionless)
c      Xg(8)  =   0.5943228461_RKIND      ! f      (dimensionless)
c      Xg(9)  =   0.8265709174_RKIND      ! f_2    (dimensionless)
c      Xg(10) =   0.9767040566_RKIND      ! f_cass (dimensionless)
c      Xg(11) =   0.9717098312_RKIND      ! s      (dimensionless)
c      Xg(12) =   0.0006830833_RKIND      ! r      (dimensionless)
c      Xg(13) =   0.0184308075_RKIND      ! y      (dimensionless)

!     O is treated as a gate in the paper with an initial value,
!     but not in CellML. There is no ODE for O but only an algebraic eq.
!      Xg(14) =   0.0000006152_RKIND      ! O      (dimensionless)

!     Below are the initial values from CellML code:
!     Initialize state variables
      X(1)   = -69.1370441636_RKIND      ! V      (units: mV)
      X(2)   = 136.7818941602_RKIND      ! K_i    (units: mM)
      X(3)   =   8.8042028653_RKIND      ! Na_i   (units: mM)
      X(4)   =   0.0001018782_RKIND      ! Ca_i   (units: mM)
      X(5)   =   0.0004468187_RKIND      ! Ca_ss  (units: mM)
      X(6)   =   3.1083688666_RKIND      ! Ca_sr  (units: mM)
      X(7)   =   0.9915800519_RKIND      ! R'     (dimensionless)

!     Initialize gating variables
      Xg(1)  =   0.0055028200_RKIND      ! x_r1   (dimensionless)
      Xg(2)  =   0.3132132864_RKIND      ! x_r2   (dimensionless)
      Xg(3)  =   0.0095370852_RKIND      ! x_s    (dimensionless)
      Xg(4)  =   0.0417391656_RKIND      ! m      (dimensionless)
      Xg(5)  =   0.1906787337_RKIND      ! h      (dimensionless)
      Xg(6)  =   0.2382198362_RKIND      ! j      (dimensionless)
      Xg(7)  =   0.0002879063_RKIND      ! d      (dimensionless)
      Xg(8)  =   0.9893285603_RKIND      ! f      (dimensionless)
      Xg(9)  =   0.9954748904_RKIND      ! f_2    (dimensionless)
      Xg(10) =   0.9999554296_RKIND      ! f_cass (dimensionless)
      Xg(11) =   0.9638610180_RKIND      ! s      (dimensionless)
      Xg(12) =   0.0010361809_RKIND      ! r      (dimensionless)
      Xg(13) =   0.0457562668_RKIND      ! y      (dimensionless)

      RETURN
      END SUBROUTINE PFIB_INIT
!-----------------------------------------------------------------------
      SUBROUTINE PFIB_READPARFF(fname)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: fname

      INTEGER fid

      fid = 1528

      OPEN(fid, FILE=TRIM(fname))
      CALL GETRVAL(fid, "Rc", Rc)
      CALL GETRVAL(fid, "Tc", Tc)
      CALL GETRVAL(fid, "Fc", Fc)
      CALL GETRVAL(fid, "Cm", Cm)
      CALL GETRVAL(fid, "sV", sV)
      CALL GETRVAL(fid, "rho", rho)
      CALL GETRVAL(fid, "V_c", V_c)
      CALL GETRVAL(fid, "V_sr", V_sr)
      CALL GETRVAL(fid, "V_ss", V_ss)
      CALL GETRVAL(fid, "K_o", K_o)
      CALL GETRVAL(fid, "Na_o", Na_o)
      CALL GETRVAL(fid, "Ca_o", Ca_o)
      CALL GETRVAL(fid, "G_K1", G_K1)
      CALL GETRVAL(fid, "G_to", G_to)
      CALL GETRVAL(fid, "G_sus", G_sus)
      CALL GETRVAL(fid, "G_fK", G_fK)
      CALL GETRVAL(fid, "G_fNa", G_fNa)
      CALL GETRVAL(fid, "G_Kr", G_Kr)
      CALL GETRVAL(fid, "G_Ks", G_Ks)
      CALL GETRVAL(fid, "G_Na", G_Na)
      CALL GETRVAL(fid, "p_KNa", p_KNa)
      CALL GETRVAL(fid, "G_CaL", G_CaL)
      CALL GETRVAL(fid, "K_NaCa", K_NaCa)
      CALL GETRVAL(fid, "gamma", gamma)
      CALL GETRVAL(fid, "K_mCa", K_mCa)
      CALL GETRVAL(fid, "K_mNai", K_mNai)
      CALL GETRVAL(fid, "K_sat", K_sat)
      CALL GETRVAL(fid, "alpha", alpha)
      CALL GETRVAL(fid, "p_NaK", p_NaK)
      CALL GETRVAL(fid, "K_mK", K_mK)
      CALL GETRVAL(fid, "K_mNa", K_mNa)
      CALL GETRVAL(fid, "G_pK", G_pK)
      CALL GETRVAL(fid, "G_pCa", G_pCa)
      CALL GETRVAL(fid, "K_pCa", K_pCa)
      CALL GETRVAL(fid, "G_bNa", G_bNa)
      CALL GETRVAL(fid, "G_bCa", G_bCa)
      CALL GETRVAL(fid, "Vmax_up", Vmax_up)
      CALL GETRVAL(fid, "K_up", K_up)
      CALL GETRVAL(fid, "V_rel", V_rel)
      CALL GETRVAL(fid, "k1p", k1p)
      CALL GETRVAL(fid, "k2p", k2p)
      CALL GETRVAL(fid, "k3", k3)
      CALL GETRVAL(fid, "k4", k4)
      CALL GETRVAL(fid, "EC", EC)
      CALL GETRVAL(fid, "max_sr", max_sr)
      CALL GETRVAL(fid, "min_sr", min_sr)
      CALL GETRVAL(fid, "V_leak", V_leak)
      CALL GETRVAL(fid, "V_xfer", V_xfer)
      CALL GETRVAL(fid, "Buf_c", Buf_c)
      CALL GETRVAL(fid, "K_bufc", K_bufc)
      CALL GETRVAL(fid, "Buf_sr", Buf_sr)
      CALL GETRVAL(fid, "K_bufsr", K_bufsr)
      CALL GETRVAL(fid, "Buf_ss", Buf_ss)
      CALL GETRVAL(fid, "K_bufss", K_bufss)
      CALL GETRVAL(fid, "Vrest", Vrest)

!     Scaling factors
      CALL GETRVAL(fid, "Vscale", Vscale)
      CALL GETRVAL(fid, "Tscale", Tscale)
      CALL GETRVAL(fid, "Voffset", Voffset)

      CLOSE(fid)

      RETURN
      END SUBROUTINE PFIB_READPARFF
!-----------------------------------------------------------------------
!     Time integration performed using Forward Euler method
      SUBROUTINE PFIB_INTEGFE(nX, nG, X, Xg, dt, Istim, Ksac, RPAR)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX, nG
      REAL(KIND=RKIND), INTENT(INOUT) :: X(nX), Xg(nG), RPAR(20)
      REAL(KIND=RKIND), INTENT(IN) :: dt, Istim, Ksac

      REAL(KIND=RKIND) :: f(nX)

!     Get time derivatives (RHS)
      CALL PFIB_GETF(nX, nG, X, Xg, f, Istim, Ksac, RPAR)

!     Update gating variables
      CALL PFIB_UPDATEG(dt, nX, nG, X, Xg)

!     Update state variables
      X = X + dt*f(:)

      RETURN
      END SUBROUTINE PFIB_INTEGFE
!-----------------------------------------------------------------------
!     Time integration performed using 4th order Runge-Kutta method
      SUBROUTINE PFIB_INTEGRK(nX, nG, X, Xg, dt, Istim, Ksac, RPAR)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX, nG
      REAL(KIND=RKIND), INTENT(INOUT) :: X(nX), Xg(nG), RPAR(20)
      REAL(KIND=RKIND), INTENT(IN) :: dt, Istim, Ksac

      REAL(KIND=RKIND) :: dt6, Xrk(nX), Xgr(nG), frk(nX,4)

      dt6 = dt/6._RKIND

!     RK4: 1st pass
      Xrk = X
      CALL PFIB_GETF(nX, nG, Xrk, Xg, frk(:,1), Istim, Ksac, RPAR)

!     Update gating variables by half-dt
      Xgr = Xg
      CALL PFIB_UPDATEG(0.5_RKIND*dt, nX, nG, X, Xgr)

!     RK4: 2nd pass
      Xrk = X + 0.5_RKIND*dt*frk(:,1)
      CALL PFIB_GETF(nX, nG, Xrk, Xgr, frk(:,2), Istim, Ksac, RPAR)

!     RK4: 3rd pass
      Xrk = X + 0.5_RKIND*dt*frk(:,2)
      CALL PFIB_GETF(nX, nG, Xrk, Xgr, frk(:,3), Istim, Ksac, RPAR)

!     Update gating variables by full-dt
      Xgr = Xg
      CALL PFIB_UPDATEG(dt, nX, nG, X, Xgr)

!     RK4: 4th pass
      Xrk = X + dt*frk(:,3)
      CALL PFIB_GETF(nX, nG, Xrk, Xgr, frk(:,4), Istim, Ksac, RPAR)

      X = X + dt6*(frk(:,1) + 2._RKIND*(frk(:,2) + frk(:,3)) + frk(:,4))
      Xg = Xgr

      RETURN
      END SUBROUTINE PFIB_INTEGRK
!-----------------------------------------------------------------------
!     Compute currents and time derivatives of state variables
      SUBROUTINE PFIB_GETF(nX, nG, X, Xg, dX, I_stim, K_sac, RPAR)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX, nG
      REAL(KIND=RKIND), INTENT(IN) :: X(nX), Xg(nG), I_stim, K_sac
      REAL(KIND=RKIND), INTENT(OUT) :: dX(nX)
      REAL(KIND=RKIND), INTENT(INOUT) :: RPAR(20)

      REAL(KIND=RKIND) :: RT, a, b, tau, sq5, e1, e2, e3, e4, n1, n2,
     2   d1, d2, d3, I_sac

!     Local copies of state variables
      V     = X(1)
      K_i   = X(2)
      Na_i  = X(3)
      Ca_i  = X(4)
      Ca_ss = X(5)
      Ca_sr = X(6)
      R_bar = X(7)

!     Local copies of gating variables
      xr1   = Xg(1)
      xr2   = Xg(2)
      xs    = Xg(3)
      m     = Xg(4)
      h     = Xg(5)
      j     = Xg(6)
      d     = Xg(7)
      f     = Xg(8)
      f2    = Xg(9)
      fcass = Xg(10)
      s     = Xg(11)
      r     = Xg(12)
      y     = Xg(13)

!     Stretch-activated currents
      I_sac = K_sac * (Vrest - V)

!      Diff = 1._RKIND / (1.0D1 * rho * Cm * sV)
      RT   = Rc * Tc / Fc
      E_K  = RT * LOG(K_o/K_i)
      E_Na = RT * LOG(Na_o/Na_i)
      E_Ca = 0.5_RKIND * RT * LOG(Ca_o/Ca_i)
      E_Ks = RT * LOG( (K_o + p_KNa*Na_o)/(K_i + p_KNa*Na_i) )

!     I_Na: Fast sodium current
      I_Na = G_Na * (m**3._RKIND) * h * j * (V - E_Na)

!     I_to: transient outward current
      I_to = G_to * r * s * (V - E_K)

!     I_K1: inward rectifier outward current
      e1   = EXP(0.06_RKIND*(V - E_K - 200._RKIND))
      e2   = EXP(2.E-4_RKIND*(V - E_K + 100._RKIND))
      e3   = EXP(0.1_RKIND*(V - E_K - 10._RKIND))
      e4   = EXP(-0.5_RKIND*(V - E_K))
      a    = 0.1_RKIND/(1._RKIND + e1)
      b    = (3._RKIND*e2 + e3) / (1._RKIND + e4)
      tau  = a / (a + b)
      sq5  = SQRT(K_o/5.4_RKIND)
      I_K1 = G_K1 * sq5 * tau * (V - E_K)

!     I_Kr: rapid delayed rectifier current
      I_Kr = G_Kr * sq5 * xr1 * xr2 * (V - E_K)

!     I_Ks: slow delayed rectifier current
      I_Ks = G_Ks * (xs**2._RKIND) * (V - E_Ks)

!     I_CaL: L-type Ca current
      a     = 2._RKIND*(V-15._RKIND)/RT
      b     = 2._RKIND*a*Fc * (0.25_RKIND*Ca_ss*EXP(a) - Ca_o) /
     2   (EXP(a)-1._RKIND)
      I_CaL = G_CaL * d * f * f2 * fcass * b

!     I_NaCa: Na-Ca exchanger current
      e1     = EXP(gamma*V/RT)
      e2     = EXP((gamma-1._RKIND)*V/RT)
      n1     = e1*(Na_i**3._RKIND)*Ca_o - e2*(Na_o**3._RKIND)*Ca_i*alpha
      d1     = K_mNai**3._RKIND + Na_o**3._RKIND
      d2     = K_mCa + Ca_o
      d3     = 1._RKIND + K_sat*e2
      I_NaCa = K_NaCa * n1 / (d1*d2*d3)

!     I_NaK: Na-K pump current
      e1    = EXP(-0.1_RKIND*V/RT)
      e2    = EXP(-V/RT)
      n1    = P_NaK * K_o * Na_i
      d1    = K_o + K_mK
      d2    = Na_i + K_mNa
      d3    = 1._RKIND + 0.1245_RKIND*e1 + 0.0353_RKIND*e2
      I_NaK = n1 / (d1*d2*d3)

!     I_pCa: plateau Ca current
      I_pCa = G_pCa * Ca_i / (K_pCa + Ca_i)

!     I_pK: plateau K current
      I_pK  = G_pK * (V-E_K) /
     2   (1._RKIND + EXP((25._RKIND-V)/5.98_RKIND))

!     I_bCa: background Ca current
      I_bCa = G_bCa * (V - E_Ca)

!     I_bNa: background Na current
      I_bNa = G_bNa * (V - E_Na)

!     I_leak: Sacroplasmic Reticulum Ca leak current
      I_leak = V_leak * (Ca_sr - Ca_i)

!     I_up: Sacroplasmic Reticulum Ca pump current
      I_up  = Vmax_up / (1._RKIND + (K_up/Ca_i)**2._RKIND)

!     I_rel: Ca induced Ca current (CICR)
      k_casr = max_sr - ((max_sr-min_sr)/
     2   (1._RKIND + (EC/Ca_sr)**2._RKIND) )
      k1     = k1p / k_casr
      O      = k1 * R_bar * (Ca_ss**2._RKIND) /
     2   (k3 + k1*(Ca_ss**2._RKIND))
      I_rel  = V_rel * O * (Ca_sr - Ca_ss)

!     I_xfer: diffusive Ca current between Ca subspae and cytoplasm
      I_xfer = V_xfer * (Ca_ss - Ca_i)

!-----------------------------------------------------------------------
!     Now compute time derivatives
!     dV/dt: rate of change of transmembrane voltage
      dX(1)  = -(I_Na + I_to + I_K1 + I_Kr + I_Ks + I_CaL + I_NaCa +
     2   I_NaK + I_pCa + I_pK + I_bCa + I_bNa + I_stim) + I_sac

!     dK_i/dt
      dX(2)  = -(Cm/(V_c*Fc)) * (I_K1 + I_to + I_Kr + I_Ks + I_pK -
     2   2._RKIND*I_NaK + I_stim)

!     dNa_i/dt
      dX(3)  = -(Cm/(V_c*Fc)) * (I_Na + I_bNa +
     2   3._RKIND*(I_NaK + I_NaCa))

!     dCa_i/dt
      n1     = (I_leak - I_up)*V_sr/V_c + I_xfer
      n2     = -(Cm/(V_c*Fc)) * (I_bCa + I_pCa - 2._RKIND*I_Naca)
     2   / 2._RKIND
      d1     = 1._RKIND + K_bufc*Buf_c/(Ca_i + K_bufc)**2._RKIND
      dX(4)  = (n1 + n2)/d1

!     dCa_ss: rate of change of Ca_ss
      n1     = (-I_CaL*Cm/(2._RKIND*Fc) + I_rel*V_sr - V_c*I_xfer)/V_ss
      d1     = 1._RKIND + K_bufss*Buf_ss/(Ca_ss + K_bufss)**2._RKIND
      dX(5)  = n1 / d1

!     dCa_sr: rate of change of Ca_sr
      n1     = I_up - I_leak - I_rel
      d1     = 1._RKIND + K_bufsr*Buf_sr/(Ca_sr + K_bufsr)**2._RKIND
      dX(6)  = n1 / d1

!     Rbar: ryanodine receptor
      k2     = k2p * k_casr
      dX(7)  = -k2*Ca_ss*R_bar + k4*(1._RKIND - R_bar)

!     Quantities to be written to file
      RPAR(3)  = I_Na
      RPAR(4)  = I_K1
      RPAR(5)  = I_to
      RPAR(6)  = I_Kr
      RPAR(7)  = I_Ks
      RPAR(8)  = I_CaL
      RPAR(9)  = I_NaCa
      RPAR(10) = I_NaK
      RPAR(11) = I_pCa
      RPAR(12) = I_pK
      RPAR(13) = I_bCa
      RPAR(14) = I_bNa
      RPAR(15) = I_leak
      RPAR(16) = I_up
      RPAR(17) = I_rel
      RPAR(18) = I_xfer
      RPAR(19) = I_stim
      RPAR(20) = I_sac

      RETURN
      END SUBROUTINE PFIB_GETF
!-----------------------------------------------------------------------
!     Update all the gating variables
      SUBROUTINE PFIB_UPDATEG(dt, n, nG, X, Xg)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: n, nG
      REAL(KIND=RKIND), INTENT(IN) :: dt, X(n)
      REAL(KIND=RKIND), INTENT(INOUT) :: Xg(nG)

      REAL(KIND=RKIND) :: a, b, c, tau

      V     = X(1)
      Ca_ss = X(5)

      xr1   = Xg(1)
      xr2   = Xg(2)
      xs    = Xg(3)
      m     = Xg(4)
      h     = Xg(5)
      j     = Xg(6)
      d     = Xg(7)
      f     = Xg(8)
      f2    = Xg(9)
      fcass = Xg(10)
      s     = Xg(11)
      r     = Xg(12)
      y     = Xg(13)

!     xr1: activation gate for I_Kr
      xr1i   = 1._RKIND/(1._RKIND + EXP(-(26._RKIND+V)/7._RKIND))
      a      = 450._RKIND/(1._RKIND + EXP(-(45._RKIND+V)/10._RKIND))
      b      = 6._RKIND/(1._RKIND + EXP((30._RKIND+V)/11.5_RKIND))
      tau    = a*b
      Xg(1)  = xr1i - (xr1i - xr1)*EXP(-dt/tau)

!     xr2: inactivation gate for I_Kr
      xr2i   = 1._RKIND /(1._RKIND + EXP((88._RKIND+V)/24._RKIND))
      a      = 3._RKIND /(1._RKIND + EXP(-(60._RKIND+V)/20._RKIND))
      b      = 1.12_RKIND/(1._RKIND + EXP(-(60._RKIND-V)/20._RKIND))
      tau    = a*b
      Xg(2)  = xr2i - (xr2i - xr2)*EXP(-dt/tau)

!     xs: activation gate for I_Ks
      xsi    = 1._RKIND/(1._RKIND + EXP(-(5._RKIND+V)/14._RKIND))
      a      = 1400._RKIND/SQRT(1._RKIND + EXP((5._RKIND-V)/6._RKIND))
      b      = 1._RKIND/(1._RKIND + EXP((V-35._RKIND)/15._RKIND))
      tau    = a*b + 80._RKIND
      Xg(3)  = xsi - (xsi - xs)*EXP(-dt/tau)

!     m: activation gate for I_Na
      mi     = 1._RKIND/( (1._RKIND +
     2   EXP(-(56.86_RKIND+V)/9.03_RKIND))**2._RKIND )
      a      = 1._RKIND/(1._RKIND + EXP(-(60._RKIND+V)/5._RKIND))
      b      = 0.1_RKIND/(1._RKIND + EXP((35._RKIND+V)/5._RKIND))
     2       + 0.1_RKIND/(1._RKIND + EXP((V-50._RKIND)/200._RKIND))
      tau    = a*b
      Xg(4)  = mi - (mi - m)*EXP(-dt/tau)

!     h: fast inactivation gate for I_Na
      hi     = 1._RKIND/( (1._RKIND
     2       + EXP((71.55_RKIND+V)/7.43_RKIND))**2._RKIND )
      IF (V .GE. -40._RKIND) THEN
         a   = 0._RKIND
         b   = 0.77_RKIND/(0.13_RKIND*(1._RKIND
     2       + EXP(-(10.66_RKIND+V)/11.1_RKIND)))
      ELSE
         a   = 5.7E-2_RKIND*EXP(-(80._RKIND+V)/6.8_RKIND)
         b   = 2.7_RKIND*EXP(0.079_RKIND*V)
     2       + 310000._RKIND*EXP(0.3485_RKIND*V)
      END IF
      tau    = 1._RKIND / (a + b)
      Xg(5)  = hi - (hi - h)*EXP(-dt/tau)

!     j: slow inactivation gate for I_Na
      ji     = 1._RKIND/( (1._RKIND
     2       + EXP((71.55_RKIND+V)/7.43_RKIND))**2._RKIND )
      IF (V .GE. -40._RKIND) THEN
         a   = 0._RKIND
         b   = 0.6_RKIND*EXP(5.7E-2_RKIND*V)
     2       / (1._RKIND + EXP(-0.1_RKIND*(V+32._RKIND)))
      ELSE
         a   = -(25428._RKIND*EXP(0.2444_RKIND*V)
     2       + 6.948E-6_RKIND*EXP(-0.04391_RKIND*V)) * (V+37.78_RKIND)
     3       / (1._RKIND + EXP(0.311_RKIND*(79.23_RKIND+V)))
         b   = 0.02424_RKIND*EXP(-0.01052_RKIND*V) /
     2         (1._RKIND + EXP(-0.1378_RKIND*(40.14_RKIND+V)))
      END IF
      tau    = 1._RKIND / (a + b)
      Xg(6)  = ji - (ji - j)*EXP(-dt/tau)

!     d: activation gate for I_CaL
      di     = 1._RKIND/(1._RKIND + EXP(-(8._RKIND+V)/7.5_RKIND))
      a      = 1.4_RKIND/(1._RKIND + EXP(-(35._RKIND+V)/13._RKIND))
     2       + 0.25_RKIND
      b      = 1.4_RKIND/(1._RKIND + EXP((5._RKIND+V)/5._RKIND))
      c      = 1._RKIND/(1._RKIND + EXP((50._RKIND-V)/20._RKIND))
      tau    = a*b + c
      Xg(7)  = di - (di - d)*EXP(-dt/tau)

!     f: slow inactivation gate for I_CaL
      fi     = 1._RKIND/(1._RKIND + EXP((20._RKIND+V)/7._RKIND))
      a      = 1102.5_RKIND*EXP(-((V+27._RKIND)**2._RKIND)/225._RKIND)
      b      = 200._RKIND/(1._RKIND + EXP((13._RKIND-V)/10._RKIND))
      c      = 180._RKIND/(1._RKIND + EXP((30._RKIND+V)/10._RKIND))
     2       + 20._RKIND
      tau    = a + b + c
c!     for spiral wave breakup
c      IF (V .GT. 0._RKIND) tau = tau*2._RKIND
      Xg(8)  = fi - (fi - f)*EXP(-dt/tau)

!     f2: fast inactivation gate for I_CaL
      f2i    = 0.67_RKIND/(1._RKIND + EXP((35._RKIND+V)/7._RKIND))
     2       + 0.33_RKIND
      a      = 562._RKIND*EXP(-((27._RKIND+V)**2._RKIND) /240._RKIND)
      b      = 31._RKIND/(1._RKIND + EXP((25._RKIND-V)/10._RKIND))
      c      = 80._RKIND/(1._RKIND + EXP((30._RKIND+V)/10._RKIND))
      tau    = a + b + c
      Xg(9)  = f2i - (f2i - f2)*EXP(-dt/tau)

!     fCass: inactivation gate for I_CaL into subspace
      c      = 1._RKIND/(1._RKIND + (Ca_ss/0.05_RKIND)**2._RKIND)
      fcassi = 0.6_RKIND*c  + 0.4_RKIND
      tau    = 80._RKIND*c + 2._RKIND
      Xg(10) = fcassi - (fcassi - fcass)*EXP(-dt/tau)

!     s: inactivation gate for I_to
      IF (i.EQ.1 .OR. i.EQ.3) THEN
         si  = 1._RKIND/(1._RKIND + EXP((20._RKIND+V)/5._RKIND))
         tau = 85._RKIND*EXP(-((V+45._RKIND)**2._RKIND) /320._RKIND)
     2       + 5._RKIND/(1._RKIND+EXP((V-20._RKIND)/5._RKIND))
     3       + 3._RKIND
      ELSE IF (i .EQ. 2) THEN
         si  = 1._RKIND/(1._RKIND + EXP((28._RKIND+V)/5._RKIND))
         tau = 1000._RKIND*EXP(-((V+67._RKIND)**2._RKIND) /1000._RKIND)
     2       + 8._RKIND
      END IF
      Xg(11) = si - (si - s)*EXP(-dt/tau)

!     r: activation gate for I_to
      ri     = 1._RKIND/(1._RKIND + EXP((20._RKIND-V)/6._RKIND))
      tau    = 9.5_RKIND*EXP(-((V+40._RKIND)**2._RKIND) /1800._RKIND)
     2       + 0.8_RKIND
      Xg(12) = ri - (ri - r)*EXP(-dt/tau)

      RETURN
      END SUBROUTINE PFIB_UPDATEG
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
      SUBROUTINE PARSESTR(strng, toks, ntoks)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: strng
      CHARACTER(LEN=*), DIMENSION(1024), INTENT(OUT) :: toks
      INTEGER(KIND=IKIND), INTENT(OUT) :: ntoks

      INTEGER(KIND=IKIND) :: i, j, slen
      CHARACTER(LEN=stdL) :: dlmtr, token

      dlmtr = ''
      token = ''

      dlmtr = '< (,=")>'
      ntoks = 1
      slen  = LEN(TRIM(strng))

      ntoks = 0
      i = 0
      DO WHILE (i .LT. slen)
         i = i + 1
         IF (INDEX(dlmtr,strng(i:i)) .NE. 0) CYCLE
         DO j=i+1, slen
            IF (INDEX(dlmtr,strng(j:j)) .NE. 0) EXIT
         END DO
         IF (j .LE. slen) THEN
            ntoks = ntoks + 1
            toks(ntoks) = strng(i:j-1)
            i = j-1
         ELSE
            EXIT
         END IF
      END DO

      RETURN
      END SUBROUTINE PARSESTR
!-----------------------------------------------------------------------
      END MODULE TTPMOD
!#######################################################################

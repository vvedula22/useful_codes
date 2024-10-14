! Created referencing CEPMOD_TTP.f
      MODULE TONGMOD
      USE TYPEMOD
      USE UTILMOD, ONLY : stdL
      IMPLICIT NONE

      PRIVATE

      INCLUDE "PARAMS_TONG.f"

      PUBLIC :: TONG_INIT
      PUBLIC :: TONG_READPARFF
      PUBLIC :: TONG_INTEGFE

      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE TONG_INIT(nX, nG, X, Xg)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX, nG
      REAL(KIND=RKIND), INTENT(OUT) :: X(nX), Xg(nG)

!       Initialize state variables
!       values from paper src code
      X(1) = -53.90915441282156 ! V [mV]
      X(2) = 0.0001161881607214449 ! Ca_i [mM]
      X(3) = 0 ! F [uN]

!       Initialize gating variables
!       values from paper src code
      Xg(1) = 0.1253518889572223 ! m [unitless] - I_Na
      Xg(2) = 0.404599170710196 ! h [unitless] - I_Na
      Xg(3) = 0.508117603077852 ! b [unitless] - I_CaT
      Xg(4) = 0.03582573962705717 ! g [unitless] - I_CaT
      Xg(5) = 0.01036961357784695 ! d [unitless] - I_CaL
      Xg(6) = 0.9065941499695301 ! f_1 [unitless] - I_CaL
      Xg(7) = 0.9065967263076083 ! f_2 [unitless] - I_CaL
      Xg(8) = 0.2060363247740295 ! q [unitless] - I_K1
      Xg(9) = 0.1922244113609531 ! r_1 [unitless] - I_K1
      Xg(10) = 0.1932803618375963 ! r_2 [unitless] - I_K1
      Xg(11) = 0.1174074734567931 ! p [unitless] - I_K2
      Xg(12) = 0.9968385770271651 ! k_1 [unitless] - I_K2
      Xg(13) = 0.9968408069904307 ! k_2 [unitless] - I_K2
      Xg(14) = 0.0003569126518797985 ! x_a [unitless] - I_BKa
      Xg(15) = 0.002220456569762898 ! x_ab1 [unitless] - I_BKab
      Xg(16) = 0.0307583106982354 ! s [unitless] - I_Ka
      Xg(17) = 0.08785242843398365 ! x [unitless] - I_Ka
      Xg(18) = 0.002604864867063448 ! y [unitless] - I_h
      Xg(19) = 0.0003764413740731269 ! c [unitless] - I_Cl
      Xg(20) = 0.2345238135343783 ! w [unitless] - Force

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
!   Time integration performed using Forward Euler method
      SUBROUTINE TONG_INTEGFE(nX, nG, X, Xg, dt, I_stim, RPAR)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX, nG
      REAL(KIND=RKIND), INTENT(INOUT) :: X(nX), Xg(nG), RPAR(21)
      REAL(KIND=RKIND), INTENT(IN) :: dt, I_stim

      REAL(KIND=RKIND) :: dX(nX), dXg(nG)

!   Get time derivatives of X (RHS)
      CALL TONG_GETF(nX, nG, X, Xg, dX, I_stim, RPAR)

!   Get time derivatives of Xg
      CALL TONG_GETG(nX, nG, X, Xg, dXg)

!   Update state & gating variables
      X = X + dt*dX(:)
      Xg = Xg + dt*dXg(:)
      X(3) = F_max * Xg(20) ! Calculate force separately

      RETURN
      END SUBROUTINE TONG_INTEGFE
!------------------------------------------------------
      SUBROUTINE TONG_GETF(nX, nG, X, Xg, dX, I_stim, RPAR)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX, nG
      REAL(KIND=RKIND), INTENT(IN) :: X(nX), Xg(nG), I_stim
      REAL(KIND=RKIND), INTENT(OUT) :: dX(nX)
      REAL(KIND=RKIND), INTENT(INOUT) :: RPAR(21)
     
! Define temp vars
      REAL(KIND=RKIND) :: RTdivF, den1, den2, den3, num1, num2
      REAL(KIND=RKIND) :: var1, var2, var3, E_Na, E_h, E_K, E_Cl, E_Ns
      REAL(KIND=RKIND) :: g_NSNa, g_NSCa, g_NSK, I_NSL, I_a, I_ab1
   
! local copies of state variables
      V = X(1)
      Ca_i = X(2)
   
! local copies of gating variables
      m = Xg(1)
      h = Xg(2)
      b_gate = Xg(3)
      g = Xg(4)
      d = Xg(5)
      f_1 = Xg(6)
      f_2 = Xg(7)
      q = Xg(8)
      r_1 = Xg(9)
      r_2 = Xg(10)
      p = Xg(11)
      k_1 = Xg(12)
      k_2 = Xg(13)
      x_a = Xg(14)
      x_ab1 = Xg(15)
      s = Xg(16)
      x_gate = Xg(17)
      y = Xg(18)
      c = Xg(19)
      w = Xg(20)
   
      RTdivF = (R*T)/F
      
! I_CaL
      den1 = 1.0 + (Ca_i/K_mCaL)**4.0
      var1 = 1 / den1
      I_CaL = g_CaL * d**2.0 * var1 * (0.8*f_1 + 0.2*f_2)*(V-E_CaL)
      
! I_Na
      E_Na = RTdivF*log(Na_o/Na_i)
      I_Na = g_Na * m**3 * h * (V-E_Na)
   
! I_CaT
      I_CaT = g_CaT * b_gate**2 * g * (V-E_CaT)
      
! I_h
      num1 = K_o + (P_NaK*Na_o)
      den1 = K_i + (P_NaK*Na_i)
      E_h = RTdivF*log(num1/den1)
      I_h = g_h * y * (V-E_h)
   
!I_K1
      E_K = RTdivF * log(K_o/K_i)
      I_K1 = g_K1*q**2*(0.38*r_1 + 0.63*r_2)*(V-E_K) ! Appdx S1: 0.38, 0.62 CellML, src: 0.38, 0.63
   
!I_K2
      I_K2 = g_K2*p**2*(0.75*k_1 + 0.25*k_2)*(V-E_K)
   
!I_Ka
      I_Ka = g_Ka*s*x_gate*(V-E_K)
      
!I_K,Ca
      ! Calculated as I_BKa and I_BKab in CellML & src code
      ! Should expand to the same
      I_a = x_a * (V-E_K)
      I_ab1 = x_ab1 * (V-E_K)
      I_KCa = g_KCa * (p_a*I_a + p_b*I_ab1)
   
!I_b
      I_b = g_b*(V-E_K)
   
!I_Cl(Ca) or I_Cl
      E_Cl = RTdivF*log(Cl_i/Cl_o)
      I_ClCa = g_Cl*c*(V-E_Cl)
   
!I_NSCC
      var1 = P_CaCs/(1+exp((V*F)/(R*T)))
      num1 = P_NaCs*Na_o + P_KCs*K_o + 4*var1*Ca_o
      den1 = P_NaCs*Na_i + P_KCs*K_i + 4*var1*Ca_i
      E_NS = RTdivF * log(num1/den1)
   
      den1 = 1 + (Mg_o/K_dMg)**1.29834 ! CellML: 1.29834 Appdx S1: 1.3
      var1 = 0.108043 + (0.903902/den1) ! CellML: 0.108043, 0.903902 Appdx S1: 0.1, 0.9
   
      den1 = 5.25e-4
      den2 = 1+(150/(Ca_o+1e-8))**2
      var2 = (1/den1) *(0.03/den2) ! gs_cao in CellML
      g_NSCa = g_NS*0.5*var2 ! in CellML variable gnsCa=0.5 is used instead of 0.5
      
      den1 = 0.0123
      den2 = 1+(150/(Na_o+1e-8))**2
      var2 = (1/den1) *(0.03/den2)
      g_NSNa = g_NS*var2 ! in CellML variable gnsNa=1 is included
   
      den1 = 0.0123
      den2 = 1+(150/(K_o+1e-8))**2
      var2 = (1/den1) *(0.03/den2)
      g_NSK = g_NS*1.19*var2 ! in CellML variable gnsK=1.19 is used instead of 1.19
   
      I_NSCa = g_NSCa*var1*(V-E_NS)
      I_NSNa = g_NSNa*var1*(V-E_NS)
      I_NSK = g_NSK*var1*(V-E_NS)
      I_NSL = g_L*var1*(V-E_NS)
      I_NSCC = I_NSCa+I_NSNa+I_NSK+I_NSL    
   
!I_NaK
      ! CellML: 0.1245 Table S4: 0.125
      den1 = 1+0.1245*exp(-0.1*V*F/(R*T))+0.00219*exp(Na_o/49.71) * 
     1 exp(-1.9*V*F/(R*T))
      var1 = 1/den1
      den1 = 1+((K_mK/K_o)**n_K)
      var2 = 1/den1
      den1 = 1+((K_mNa/Na_i)**n_Na)
      var3 = 1/den1
      I_NaK = g_NaK*var1*var2*var3
   
! Calculate fluxes
   
!J_PMCA
      den1 = 1 + (K_mPMCA/Ca_i)**n_PMCA
      J_PMCA_var = J_PMCA/den1
      
!J_NaCa, I_NaCa
      den1 = 1+((K_mAllo/Ca_i)**n_Allo)
      var1 = 1/den1
      var2 = exp((gam-1)*((V*F)/(R*T)))
      var3 = exp(gam*((V*F)/(R*T)))
      den1 = 1+k_sat*var2 ! 1st coefficient
      den2 = (K_mCao*Na_i**3) + ((K_mNao**3)*Ca_i) + (Ca_o*Na_i**3) + 
     1 (Na_o**3*Ca_i) ! 2nd coefficient 1st half of sum
      den3 = ((K_mNai**3)*Ca_o)*(1+(Ca_i/K_mCai)) + (Na_o**3*K_mCai) * 
     1 (1+(Na_i/K_mNai)**3) ! 2nd coefficient 2nd half of sum
      num1 = J_NaCa*var1
      num2 = (Na_i**3*Ca_o*var3)-(Na_o**3*Ca_i*var2)
      ! Equations for I_NaCa/J_NaCa vary between all 3 sources both CellML & src
      ! include the 0.5 factor, while Appdx S1 does not
      J_NaCa_var = -(num1*num2) / (den1*(den2+den3))
      I_NaCa = -0.5 * ((z_Ca*F*J_NaCa_var)/(C_m*bet)) * (1/AV_c)
   
!J_Ca,mem
      ! Equation reproduced from source code not given in Appendix S1 
      J_Camem = ((AV_c*C_m*bet)/(z_Ca*F))*(I_CaL+I_CaT+I_NSCa)
      
      ! Compute time derivatives
      
!dV/dt
      dX(1) = -(I_CaL+I_Na+I_CaT+I_h+I_K1+I_K2+I_Ka+I_KCa+I_b+I_ClCa+
     1 I_NSCC+I_NaK+I_NaCa+I_stim)
      
!dCa_i/dt
      dX(2) = -(J_Camem+J_NaCa_var+J_PMCA_var)
   
   ! Start at 3 to let first 2 outputs be V and Ca_i
      RPAR(3) = I_CaL
      RPAR(4) = I_Na
      RPAR(5) = I_CaT
      RPAR(6) = I_h
      RPAR(7) = I_K1
      RPAR(8) = I_K2
      RPAR(9) = I_Ka
      RPAR(10) = I_KCa
      RPAR(11) = I_b
      RPAR(12) = I_ClCa
      RPAR(13) = I_NSCC
      RPAR(14) = I_NaK
      RPAR(15) = I_NaCa
      RPAR(16) = J_PMCA_var
      RPAR(17) = J_NaCa_var
      RPAR(18) = J_Camem
      RPAR(19) = I_NSNa
      RPAR(20) = I_NSCa
      RPAR(21) = I_NSK
      
      RETURN
      END SUBROUTINE TONG_GETF
      !------------------------------------------------------
      SUBROUTINE TONG_GETG(nX, nG, X, Xg, dXg)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX, nG
      REAL(KIND=RKIND), INTENT(IN) :: X(nX), Xg(nG)
      REAL(KIND=RKIND), INTENT(INOUT) :: dXg(nG)
   
   !     Define temp vars
      REAL(KIND=RKIND) :: den1, den2, ss, tau, var1, var2
      
   !     Local copies of state variables
      V = X(1)
      Ca_i = X(2)
   
   !     Local copies of gating variables
      m = Xg(1)
      h = Xg(2)
      b_gate = Xg(3)
      g = Xg(4)
      d = Xg(5)
      f_1 = Xg(6)
      f_2 = Xg(7)
      q = Xg(8)
      r_1 = Xg(9)
      r_2 = Xg(10)
      p = Xg(11)
      k_1 = Xg(12)
      k_2 = Xg(13)
      x_a = Xg(14)
      x_ab1 = Xg(15)
      s = Xg(16)
      x_gate = Xg(17)
      y = Xg(18)
      c = Xg(19)
      w = Xg(20)
   
   ! Note: values used in gating variables vary greatly between Appdx S1 &
   ! src/CellML (mostly by precision/sigfigs)
   
   ! I_CaL
      ! d
      den1 = 1 + exp(-(V+22)/7)
      ss = 1/den1
      
      den1 = 1+((V+29.97)/9)**2
      tau = 2.29 + 5.7/den1
      
      dXg(5) = (ss-d)/tau
      
      ! f_1
      den1 = 1 + exp((V+38)/f_ss_sf) ! inactivation steady state, slope factor: normally 7
      ss = 1/den1 ! also used for f_2
      
      tau = 12 ! [ms]
      
      dXg(6) = (ss-f_1)/tau
      
      ! f_2
      den1 = (1+exp((V+13.9629)/45.3782)) ! CellML: 13.9629, 45.3782 Appdx S1: 13.96, 45.38
      den2 = (1+exp(-(V+9.49866)/3.3945)) ! CellML: 9.49866, 3.3945 Appdx S1: 9.5, 3.39
      tau = 90.9699 * (1-(1/(den1*den2))) ! CellML: 90.9699 Appdx S1: 90.97
      
      dXg(7) = (ss-f_2)/tau
      
      ! I_Na
      ! m
      den1 = 1+exp((-(V+35.9584))/9.24013) ! CellML: 35.9584, 9.24013 Appdx S1: 35.96, 9.24
      ss = 1/den1
      
      den1 = 1+exp((V+38)/10)
      tau = 0.25+(7/den1)
      
      dXg(1)=(ss-m)/tau
      
      ! h
      den1 = 1+exp((V+57)/8)
      ss = 1/den1
      
      den1 = 1+((V+47.5)/1.5)**2
      tau = 0.9+(1002.85/den1)
      
      dXg(2) = (ss-h)/tau
   
      
      !I_CaT
      ! b
      den1 = 1+exp((-(V+54.23))/9.88)
      ss = 1/den1
      
      den1 = 1+((V+66)/26)**2
      tau = 0.45 + 3.9/den1
      
      dXg(3) = (ss-b_gate)/tau
      
      ! g
      den1 = 1+exp((V+72.978)/4.64) ! CellML: 72.978 Appdx S1: 72.98
      ss = 0.02 + 0.98/den1
      
      den1 = 1+exp((V-417.43)/203.18)
      den2 = 1+exp((-(V+61.11))/8.07)
      tau = 150 - 150/(den1*den2)
      
      dXg(4) = (ss-g)/tau
   
      
      !I_h
      ! y
      den1 = 1+exp((V+105.39)/8.6553) ! CellML: 8.6553 Appdx S1: 8.66
      ss = 1/den1
      
      den1 = (3.5e-6*exp(-0.0497*V)) + (0.04003*exp(0.05211*V))  ! CellML: 0.04003, 0.05211 Appdx S1: 0.04, 0.0521
      tau = 1/den1
      
      dXg(18) = (ss-y)/tau
   
      
      !I_K1
      den1 = 1+exp((-(V+18.6736))/26.6606) ! CellML: 18.6736, 26.6606 Appdx S1: 18.67, 26.66
      ss = 0.978613/den1 ! CellML: 0.978613, Appdx S1: 1
      
      den1 = 1+((V+60.71)/15.79)**2
      tau = 500.0/den1
      
      dXg(8) = (ss-q)/tau
      
      ! r_1
      den1 = 1+exp((V+63)/6.3)
      ss = 1/den1 ! also used for r_2
      
      den1 = 1+((V+62.7133)/35.8611)**2 ! CellML: 62.7133, 35.8611 Appdx S1: 62.71, 35.86
      tau = 5e3/den1 ! CellML/src: 5e3, Appdx S1: 5e4
      
      dXg(9) = (ss-r_1)/tau
      
      ! r_2
      den1 = 1+exp((V+22)/4)
      tau = 3.0e4 + (2.2e5/den1)
      
      dXg(10) = (ss-r_2)/tau
   
      
      ! I_K2
      ! p_ss equation varies greatly between CellML/src and Appdx S1
      ! Appdx S1:
      ! p_ss_den = 1+exp((-(V+0.948))/17.91)
      ! p_ss = 1/p_ss_den
      ! CellML/src
      
      ! p    
      den1 = 1+exp(-(V+17.91)/18.4)
      ss = 0.948/den1
      
      den1 = 1+((V+64.1)/28.67)**2
      tau = 100/den1
      
      dXg(11) = (ss-p)/tau
      
      ! k_1
      den1 = 1+exp((V+21.2)/5.7)
      ss = 1/den1 ! also used for k_2
      
      den1 = 1+exp((V-315)/50)
      den2 = 1+exp((-(V+74.9))/8)
      tau = 1e6*(1 - (1/(den1 * den2)))
      
      dXg(12) = (ss-k_1)/tau
      
      ! k_2
      den1 = 1+exp((V-132.868)/25.3992) ! CellML: 132.868, 25.3992 Appdx S1: 132.87, 25.4
      den2 = 1+exp((-(V+24.9203))/2.67915) ! CellML: 24.9203, 2.67915 Appdx S1: 24.92, 2.68
      tau = 2.5e6*(1-(1/(den1*den2)))
      
      dXg(13) = (ss-k_2)/tau
      
      
      ! IKa
      ! s
      den1 = 1+exp((-(V+27.79))/7.57)
      ss = 1/den1
      
      den1 = 1+((V+20.52)/35)**2
      tau = 17/den1
      
      dXg(16) = (ss-s)/tau
      
      ! x
      den1 = 1+exp((V+69.5)/6)
      ss = 0.02+(0.98/den1)
      
      den1 = 1+((V+34.18)/120)**2
      tau = 7.5 + (10/den1)
      
      dXg(17) = (ss-x_gate)/tau
   
      
      ! I_KCa
      !BKa
      ! x_a    
      den1 = (((1000*Ca_i)+1538.29)/739.057)**2 ! CellML: 739.057 Appdx S1: 739.06
      den2 = (((1000*Ca_i)-0.0630535)/0.161942)**2 ! CellML: 0.0630535, 0.161942 Appdx S1: 0.063, 0.162
      var1 = (8.38/(1+den1))-(0.749/(1+den2)) ! CellML: 8.38384, 0.749234 Appdx S1: 8.38, 0.749
      
      den1 = (((1000*Ca_i) + 0.237503)/2.39278e-4)**0.42291 ! CellML: 0.237503, ~2.39278e-4, 0.42291 Appdx S1: 0.238, 2.39e-4, 0.423
      var2 = (5011.47/(1+den1))-37.5137 ! CellML: 37.5137 Appdx S1: 37.51
      
      den1 = exp((-var1*F*(V-var2))/(R*T))
      ss = 1/(1+den1)
      
      den1 = ((V-158.779)/(-52.1497))**2 ! CellML: 158.779, 52.1497 Appdx S1: 158.78, 52.15
      tau = 2.40914/(1+den1) ! CellML: 2.40914 Appdx S1: 2.41
      
      dXg(14) = (ss-x_a)/tau
      
      !BKab
      ! x_ab1
      den1 = (((1000*Ca_i)+228.71)/684.946)**2 ! CellML: 684.946 Appdx S1: 684.95
      den2 = (((1000*Ca_i)-0.219)/0.428)**2  ! CellML: 0.218988, 0.428335
      var1 = (1.4/(1+den1))-(0.681249/(1+den2)) ! CellML: 0.681249 Appdx S1: 0.681
      
      den1 = (((1000*Ca_i)+0.401189)/3.99115e-3)**0.668054 ! CelLML: 0.401189, ~3.99115e-3, 0.668054 Appdx S1: 0.401, 3.99e-3, 0.668
      var2 = (8540.23/(1+den1))-109.275 ! CellML: 109.275 Appdx S1: 109.28
      
      den1 = exp((-var1*F*(V-var2))/(R*T))
      ss = 1/(1+den1)
      
      den1 = ((V-153.019)/66.4952)**2 ! CellML: 153.019, 66.4952 Appdx S1: 153.02, 66.5
      tau = 13.8049/(1+den1) ! CellML: 13.8049 Appdx S1: 13.8
      
      dXg(15) = (ss-x_ab1)/tau
   
      ! I_ClCa (I_Cl in CellML/src)
      ! c
      var1 = 0.0006*exp((2.53*F*V)/(R*T))
      var2 = 0.1*exp((-5*F*V)/(R*T))
      den1 = (var1/Ca_i)**2 + (var1/Ca_i) + 1
      ss = 1/(1+(var2*den1))
      
      den1 = exp((V+4.56)/11.62)
      den2 = exp((-(V+25.5))/11.62)
      tau = (210.0/(1+den1)) + (170/(1+den2)) - 160.0
      
      dXg(19) = (ss-c)/tau
   
      
      ! Calcium-Dependent Force
      ! w
      den1 = 1+((K_mF/Ca_i)**n_F) ! also used for tau denominator
      ss = 1.0/den1
      
      tau = 4000 * (0.234845 + (1.0-0.234845)/den1) ! Src Code: 0.234845 Appdx S1: 0.235
      
      dXg(20) = (ss-w)/tau
   
      RETURN
      END SUBROUTINE TONG_GETG
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
      SUBROUTINE GETRVEC(fileId, skwrd, nd, rVec)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: fileId, nd
      CHARACTER(LEN=*), INTENT(IN) :: skwrd
      REAL(KIND=RKIND), INTENT(INOUT) :: rVec(nd)
  
      INTEGER(KIND=IKIND) :: slen, i, ios, nt
      CHARACTER(LEN=stdL) :: sline, scmd, sval
      CHARACTER(LEN=stdL), DIMENSION(250) :: tokList
  
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
            CALL PARSESTR(sval, tokList, nt)
            IF (nt .NE. nd) THEN
               DO i=1, nt
                  WRITE(*,'(I2,2X,A)') i, TRIM(tokList(i))
               END DO
               STOP " Error: Unexpected token length "//TRIM(skwrd)
            END IF

            DO i=1, nt
               READ(tokList(i),*,IOSTAT=ios) rvec(i)
               IF (ios .NE. 0) THEN
                  STOP " Error: while reading "//TRIM(skwrd)
               END IF
            END DO
            EXIT
         END IF
      END DO
  
 001  RETURN
  
! 001  STOP " Error: EOF reached while finding command <"//
!     2   TRIM(skwrd)//">"
  
      END SUBROUTINE GETRVEC
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
      END MODULE TONGMOD
!#######################################################################
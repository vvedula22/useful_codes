!
! Copyright (c) Stanford University, The Regents of the University of
!               California, and others.
!
! All Rights Reserved.
!
! See Copyright-SimVascular.txt for additional details.
!
! Permission is hereby granted, free of charge, to any person obtaining
! a copy of this software and associated documentation files (the
! "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish,
! distribute, sublicense, and/or sell copies of the Software, and to
! permit persons to whom the Software is furnished to do so, subject
! to the following conditions:
!
! The above copyright notice and this permission notice shall be included
! in all copies or substantial portions of the Software.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!-----------------------------------------------------------------------
!
!     This routine embodies formulation for solving electrophysiology
!     cellular activation model at the nodal level.
!
!-----------------------------------------------------------------------

!     Initialize Cardiac Electrophysiology Model
      SUBROUTINE CEPINIT()
      USE CEPMOD
      IMPLICIT NONE

      SELECT CASE (cep%cepType)
      CASE (cepModel_AP)
!        Clear log file
         WRITE(oFile,'(A)') "log_AP.txt"
         CALL CTXT(oFile)

         WRITE(*,'(4X,A)') "Initializing Aliev-Panfilov model"
         CALL AP_INIT(nX, X)

      CASE (cepModel_BO)
!        Clear log file
         IF (cep%imyo .EQ. 1) THEN
            WRITE(oFile,'(A)') "log_BO_epi.txt"
         ELSE IF (cep%imyo .EQ. 2) THEN
            WRITE(oFile,'(A)') "log_BO_endo.txt"
         ELSE IF (cep%imyo .EQ. 3) THEN
            WRITE(oFile,'(A)') "log_BO_myo.txt"
         END IF
         CALL CTXT(oFile)

!        Initialize state variables
         WRITE(*,'(4X,A)') "Initializing Bueno-Orovio model"
         CALL BO_INIT(nX, X)

      CASE (cepModel_FN)
!        Clear log file
         WRITE(oFile,'(A)') "log_FN.txt"
         CALL CTXT(oFile)

!        Initialize state variables
         WRITE(*,'(4X,A)') "Initializing Fitzhugh-Nagumo model"
         CALL FN_INIT(nX, X)

      CASE (cepModel_TTP)
!        Clear log file
         IF (cep%imyo .EQ. 1) THEN
            WRITE(oFile,'(A)') "log_TTP_epi.txt"
         ELSE IF (cep%imyo .EQ. 2) THEN
            WRITE(oFile,'(A)') "log_TTP_endo.txt"
         ELSE IF (cep%imyo .EQ. 3) THEN
            WRITE(oFile,'(A)') "log_TTP_myo.txt"
         END IF
         CALL CTXT(oFile)

!        Initialize state variables
         WRITE(*,'(4X,A)') "Initializing tenTusscher-Panfilov model"
         CALL TTP_INIT(cep%imyo, nX, nG, X, Xg)

      END SELECT

      RETURN
      END SUBROUTINE CEPINIT
!#######################################################################
!     State variable integration
      SUBROUTINE CEPINTEG()
      USE CEPMOD
      IMPLICIT NONE

      INTEGER, PARAMETER :: iwrite = 20
      INTEGER i, j, k, no, cntr
      REAL(KIND=8) :: cTS, Istim, Ksac, eX
      LOGICAL flag

      INTEGER, ALLOCATABLE :: IPAR(:)
      REAL(KIND=8), ALLOCATABLE :: Xo(:), RPAR(:)

!     Output to files
      no = nX
      IF (cep%cepType .EQ. cepModel_BO) THEN
         no = no + 3
         ALLOCATE(RPAR(5))
      ELSE IF (cep%cepType .EQ. cepModel_TTP) THEN
         no = no + 16
         ALLOCATE(RPAR(18))
      ELSE
         ALLOCATE(RPAR(2))
      END IF
      RPAR(:) = 0D0

      IF (cep%emCpld) no = no + 2
      ALLOCATE(Xo(no))
      Xo = 0.0D0

      IF (cep%odeS%tIntType .EQ. tIntType_CN2) THEN
         ALLOCATE(IPAR(2))
         IPAR(1) = cep%odeS%maxItr
         IPAR(2) = 0
         RPAR(1) = cep%odeS%absTol
         RPAR(2) = cep%odeS%relTol
      END IF

!     Feedback coefficient for stretch-activated-currents
      Ksac = 0D0

      IF (iProg) THEN
         WRITE(*,'(4X,A)',ADVANCE='NO') "Time integration progress: "
      END IF

      SELECT CASE (cep%cepType)
      CASE (cepModel_AP)
!        Time integrator
         SELECT CASE (cep%odes%tIntType)
         CASE (tIntType_FE)
!           Time loop
            j   = 1
            k   = 1
            cTS = 0D0
            DO i=1, nTS
!              Get stimulus current
               CALL GET_STIM_AMP(cTS, cep, Istim, flag)
               IF (flag) RETURN

!              Integrate local state variables
               CALL AP_INTEGFE(nX, X, cTS, dt, Istim, Ksac)
               Xo(1:nX) = X(:)

!              Perform NaN check
               IF (ISNAN(X(1))) STOP "NaN occurence (X). Aborted!"

!              Electromechanics excitation-activation
               IF (cep%emCpld) THEN
                  CALL AP_ACTVSTRS(X(1), dt, cep%Tact, eX)
                  IF (ISNAN(cep%Tact)) THEN
                     STOP "NaN occurence (Tact). Aborted!"
                  END IF
                  Xo(no-1) = eX
                  Xo(no)   = cep%Tact
               END IF

!              Time step advance
               cTS = cTS + dt

!              Write state variables to a log file
               IF (MOD(i,iwrite) .EQ. 0)
     2            CALL WTXT(oFile, no, Xo, cTS)

!           Display progress
               IF (i.EQ.j .AND. iProg) THEN
                  WRITE(*,'(A)',ADVANCE='NO') (k-1)*20//"%  "
                  k = k + 1
                  j = NINT(REAL((k-1)*nTS,KIND=8)/5.0D0)
               END IF
            END DO

         CASE (tIntType_RK4)
!           Time loop
            j   = 1
            k   = 1
            cTS = 0D0
            DO i=1, nTS
!              Get stimulus current
               CALL GET_STIM_AMP(cTS, cep, Istim, flag)
               IF (flag) RETURN

!              Integrate local state variables
               CALL AP_INTEGRK(nX, X, cTS, dt, Istim, Ksac)
               Xo(1:nX) = X(:)

!              Perform NaN check
               IF (ISNAN(X(1))) STOP "NaN occurence (X). Aborted!"

!              Electromechanics excitation-activation
               IF (cep%emCpld) THEN
                  CALL AP_ACTVSTRS(X(1), dt, cep%Tact, eX)
                  IF (ISNAN(cep%Tact)) THEN
                     STOP "NaN occurence (Tact). Aborted!"
                  END IF
                  Xo(no-1) = eX
                  Xo(no)   = cep%Tact
               END IF

!              Time step advance
               cTS = cTS + dt

!              Write state variables to a log file
               IF (MOD(i,iwrite) .EQ. 0)
     2            CALL WTXT(oFile, no, Xo, cTS)

!           Display progress
               IF (i.EQ.j .AND. iProg) THEN
                  WRITE(*,'(A)',ADVANCE='NO') (k-1)*20//"%  "
                  k = k + 1
                  j = NINT(REAL((k-1)*nTS,KIND=8)/5.0D0)
               END IF
            END DO

         CASE (tIntType_CN2)
!           Time loop
            j   = 1
            k   = 1
            cTS = 0D0
            DO i=1, nTS
!              Get stimulus current
               CALL GET_STIM_AMP(cTS, cep, Istim, flag)
               IF (flag) RETURN

!              Integrate local state variables
               CALL AP_INTEGCN2(nX, X, cTS, dt, Istim, Ksac, IPAR, RPAR)
               Xo(1:nX) = X(:)

!              Perform NaN check
               IF (ISNAN(X(1))) STOP "NaN occurence (X). Aborted!"

!              Electromechanics excitation-activation
               IF (cep%emCpld) THEN
                  CALL AP_ACTVSTRS(X(1), dt, cep%Tact, eX)
                  IF (ISNAN(cep%Tact)) THEN
                     STOP "NaN occurence (Tact). Aborted!"
                  END IF
                  Xo(no-1) = eX
                  Xo(no)   = cep%Tact
               END IF

!              Time step advance
               cTS = cTS + dt

!              Write state variables to a log file
               IF (MOD(i,iwrite) .EQ. 0)
     2            CALL WTXT(oFile, no, Xo, cTS)

!           Display progress
               IF (i.EQ.j .AND. iProg) THEN
                  WRITE(*,'(A)',ADVANCE='NO') (k-1)*20//"%  "
                  k = k + 1
                  j = NINT(REAL((k-1)*nTS,KIND=8)/5.0D0)
               END IF
            END DO
         END SELECT

      CASE (cepModel_BO)
!        Time integrator
         SELECT CASE (cep%odes%tIntType)
         CASE (tIntType_FE)
!           Time loop
            j   = 1
            k   = 1
            cTS = 0D0
            DO i=1, nTS
!              Get stimulus current
               CALL GET_STIM_AMP(cTS, cep, Istim, flag)
               IF (flag) RETURN

!              Integrate local state variables
               CALL BO_INTEGFE(cep%imyo, nX, X, cTS, dt, Istim, Ksac,
     2            RPAR)
               Xo(1:nX) = X(:)
               Xo(nX+1:nX+3) = RPAR(3:5)

!              Perform NaN check
               IF (ISNAN(X(1))) STOP "NaN occurence (X). Aborted!"

!              Electromechanics excitation-activation
               IF (cep%emCpld) THEN
                  CALL BO_ACTVSTRS(X(4), dt, cep%Tact, eX)
                  IF (ISNAN(cep%Tact)) THEN
                     STOP "NaN occurence (Tact). Aborted!"
                  END IF
                  Xo(no-1) = eX
                  Xo(no)   = cep%Tact
               END IF

!              Time step advance
               cTS = cTS + dt

!              Write state variables to a log file
               IF (MOD(i,iwrite) .EQ. 0)
     2            CALL WTXT(oFile, no, Xo, cTS)

!           Display progress
               IF (i.EQ.j .AND. iProg) THEN
                  WRITE(*,'(A)',ADVANCE='NO') (k-1)*20//"%  "
                  k = k + 1
                  j = NINT(REAL((k-1)*nTS,KIND=8)/5.0D0)
               END IF
            END DO

         CASE (tIntType_RK4)
!           Time loop
            j   = 1
            k   = 1
            cTS = 0D0
            DO i=1, nTS
!              Get stimulus current
               CALL GET_STIM_AMP(cTS, cep, Istim, flag)
               IF (flag) RETURN

!              Integrate local state variables
               CALL BO_INTEGRK(cep%imyo, nX, X, cTS, dt, Istim, Ksac,
     2            RPAR)
               Xo(1:nX) = X(:)
               Xo(nX+1:nX+3) = RPAR(3:5)

!              Perform NaN check
               IF (ISNAN(X(1))) STOP "NaN occurence (X). Aborted!"

!              Electromechanics excitation-activation
               IF (cep%emCpld) THEN
                  CALL BO_ACTVSTRS(X(4), dt, cep%Tact, eX)
                  IF (ISNAN(cep%Tact)) THEN
                     STOP "NaN occurence (Tact). Aborted!"
                  END IF
                  Xo(no-1) = eX
                  Xo(no)   = cep%Tact
               END IF

!              Time step advance
               cTS = cTS + dt

!              Write state variables to a log file
               IF (MOD(i,iwrite) .EQ. 0)
     2            CALL WTXT(oFile, no, Xo, cTS)

!           Display progress
               IF (i.EQ.j .AND. iProg) THEN
                  WRITE(*,'(A)',ADVANCE='NO') (k-1)*20//"%  "
                  k = k + 1
                  j = NINT(REAL((k-1)*nTS,KIND=8)/5.0D0)
               END IF
            END DO

         CASE (tIntType_CN2)
!           Time loop
            j   = 1
            k   = 1
            cTS = 0D0
            DO i=1, nTS
!              Get stimulus current
               CALL GET_STIM_AMP(cTS, cep, Istim, flag)
               IF (flag) RETURN

!              Integrate local state variables
               CALL BO_INTEGCN2(cep%imyo, nX, X, cTS, dt, Istim, Ksac,
     2            IPAR, RPAR)
               Xo(1:nX) = X(:)
               Xo(nX+1:nX+3) = RPAR(3:5)

!              Perform NaN check
               IF (ISNAN(X(1))) STOP "NaN occurence (X). Aborted!"

!              Electromechanics excitation-activation
               IF (cep%emCpld) THEN
                  CALL BO_ACTVSTRS(X(4), dt, cep%Tact, eX)
                  IF (ISNAN(cep%Tact)) THEN
                     STOP "NaN occurence (Tact). Aborted!"
                  END IF
                  Xo(no-1) = eX
                  Xo(no)   = cep%Tact
               END IF

!              Time step advance
               cTS = cTS + dt

!              Write state variables to a log file
               IF (MOD(i,iwrite) .EQ. 0)
     2            CALL WTXT(oFile, no, Xo, cTS)

!           Display progress
               IF (i.EQ.j .AND. iProg) THEN
                  WRITE(*,'(A)',ADVANCE='NO') (k-1)*20//"%  "
                  k = k + 1
                  j = NINT(REAL((k-1)*nTS,KIND=8)/5.0D0)
               END IF
            END DO
         END SELECT

      CASE (cepModel_FN)
!        Time integrator
         SELECT CASE (cep%odes%tIntType)
         CASE (tIntType_FE)
!           Time loop
            j   = 1
            k   = 1
            cTS = 0D0
            DO i=1, nTS
!              Get stimulus current
               CALL GET_STIM_AMP(cTS, cep, Istim, flag)
               IF (flag) RETURN

!              Integrate local state variables
               CALL FN_INTEGFE(nX, X, cTS, dt, Istim)
               Xo(1:nX) = X(:)

!              Perform NaN check
               IF (ISNAN(X(1))) STOP "NaN occurence (X). Aborted!"

!              Time step advance
               cTS = cTS + dt

!              Write state variables to a log file
               IF (MOD(i,iwrite) .EQ. 0)
     2            CALL WTXT(oFile, no, Xo, cTS)

!           Display progress
               IF (i.EQ.j .AND. iProg) THEN
                  WRITE(*,'(A)',ADVANCE='NO') (k-1)*20//"%  "
                  k = k + 1
                  j = NINT(REAL((k-1)*nTS,KIND=8)/5.0D0)
               END IF
            END DO

         CASE (tIntType_RK4)
!           Time loop
            j   = 1
            k   = 1
            cTS = 0D0
            DO i=1, nTS
!              Get stimulus current
               CALL GET_STIM_AMP(cTS, cep, Istim, flag)
               IF (flag) RETURN

!              Integrate local state variables
               CALL FN_INTEGRK(nX, X, cTS, dt, Istim)
               Xo(1:nX) = X(:)

!              Perform NaN check
               IF (ISNAN(X(1))) STOP "NaN occurence (X). Aborted!"

!              Time step advance
               cTS = cTS + dt

!              Write state variables to a log file
               IF (MOD(i,iwrite) .EQ. 0)
     2            CALL WTXT(oFile, no, Xo, cTS)

!           Display progress
               IF (i.EQ.j .AND. iProg) THEN
                  WRITE(*,'(A)',ADVANCE='NO') (k-1)*20//"%  "
                  k = k + 1
                  j = NINT(REAL((k-1)*nTS,KIND=8)/5.0D0)
               END IF
            END DO

         CASE (tIntType_CN2)
!           Time loop
            j   = 1
            k   = 1
            cTS = 0D0
            DO i=1, nTS
!              Get stimulus current
               CALL GET_STIM_AMP(cTS, cep, Istim, flag)
               IF (flag) RETURN

!              Integrate local state variables
               CALL FN_INTEGCN2(nX, X, cTS, dt, Istim, IPAR, RPAR)
               Xo(1:nX) = X(:)

!              Perform NaN check
               IF (ISNAN(X(1))) STOP "NaN occurence (X). Aborted!"

!              Time step advance
               cTS = cTS + dt

!              Write state variables to a log file
               IF (MOD(i,iwrite) .EQ. 0)
     2            CALL WTXT(oFile, no, Xo, cTS)

!           Display progress
               IF (i.EQ.j .AND. iProg) THEN
                  WRITE(*,'(A)',ADVANCE='NO') (k-1)*20//"%  "
                  k = k + 1
                  j = NINT(REAL((k-1)*nTS,KIND=8)/5.0D0)
               END IF
            END DO
         END SELECT

      CASE (cepModel_TTP)
!        Time integrator
         SELECT CASE (cep%odes%tIntType)
         CASE (tIntType_FE)
!           Time loop
            j   = 1
            k   = 1
            cTS = 0D0
            DO i=1, nTS
!              Get stimulus current
               CALL GET_STIM_AMP(cTS, cep, Istim, flag)
               IF (flag) RETURN

!              Integrate local state variables
               CALL TTP_INTEGFE(cep%imyo, nX, nG, X, Xg, cTS, dt, Istim,
     2            Ksac, RPAR)
               Xo(1:nX) = X(1:nX)
               Xo(nX+1:nX+16) = RPAR(3:18)

!              Perform NaN check
               IF (ISNAN(X(1))) STOP "NaN occurence (X). Aborted!"

!              Electromechanics excitation-activation
               IF (cep%emCpld) THEN
                  CALL TTP_ACTVSTRS(X(4), dt, cep%Tact, eX)
                  IF (ISNAN(cep%Tact)) THEN
                     STOP "NaN occurence (Tact). Aborted!"
                  END IF
                  Xo(no-1) = eX
                  Xo(no)   = cep%Tact
               END IF

!              Time step advance
               cTS = cTS + dt

!              Write state variables to a log file
               IF (MOD(i,iwrite) .EQ. 0)
     2            CALL WTXT(oFile, no, Xo, cTS)

!           Display progress
               IF (i.EQ.j .AND. iProg) THEN
                  WRITE(*,'(A)',ADVANCE='NO') (k-1)*20//"%  "
                  k = k + 1
                  j = NINT(REAL((k-1)*nTS,KIND=8)/5.0D0)
               END IF
            END DO

         CASE (tIntType_RK4)
!           Time loop
            j   = 1
            k   = 1
            cTS = 0D0
            DO i=1, nTS
!              Get stimulus current
               CALL GET_STIM_AMP(cTS, cep, Istim, flag)
               IF (flag) RETURN

!              Integrate local state variables
               CALL TTP_INTEGRK(cep%imyo, nX, nG, X, Xg, cTS, dt, Istim,
     2            Ksac, RPAR)
               Xo(1:nX) = X(1:nX)
               Xo(nX+1:nX+16) = RPAR(3:18)

!              Perform NaN check
               IF (ISNAN(X(1))) STOP "NaN occurence (X). Aborted!"

!              Electromechanics excitation-activation
               IF (cep%emCpld) THEN
                  CALL TTP_ACTVSTRS(X(4), dt, cep%Tact, eX)
                  IF (ISNAN(cep%Tact)) THEN
                     STOP "NaN occurence (Tact). Aborted!"
                  END IF
                  Xo(no-1) = eX
                  Xo(no)   = cep%Tact
               END IF

!              Time step advance
               cTS = cTS + dt

!              Write state variables to a log file
               IF (MOD(i,iwrite) .EQ. 0)
     2            CALL WTXT(oFile, no, Xo, cTS)

!           Display progress
               IF (i.EQ.j .AND. iProg) THEN
                  WRITE(*,'(A)',ADVANCE='NO') (k-1)*20//"%  "
                  k = k + 1
                  j = NINT(REAL((k-1)*nTS,KIND=8)/5.0D0)
               END IF
            END DO

         CASE (tIntType_CN2)
!           Time loop
            j   = 1
            k   = 1
            cTS = 0D0
            DO i=1, nTS
!              Get stimulus current
               CALL GET_STIM_AMP(cTS, cep, Istim, flag)
               IF (flag) RETURN

!              Integrate local state variables
               CALL TTP_INTEGCN2(cep%imyo, nX, nG, X, Xg, cTS, dt,
     2            Istim, Ksac, IPAR, RPAR)
               Xo(1:nX) = X(1:nX)
               Xo(nX+1:nX+16) = RPAR(3:18)

!              Perform NaN check
               IF (ISNAN(X(1))) STOP "NaN occurence (X). Aborted!"

!              Electromechanics excitation-activation
               IF (cep%emCpld) THEN
                  CALL TTP_ACTVSTRS(X(4), dt, cep%Tact, eX)
                  IF (ISNAN(cep%Tact)) THEN
                     STOP "NaN occurence (Tact). Aborted!"
                  END IF
                  Xo(no-1) = eX
                  Xo(no)   = cep%Tact
               END IF

!              Time step advance
               cTS = cTS + dt

!              Write state variables to a log file
               IF (MOD(i,iwrite) .EQ. 0)
     2            CALL WTXT(oFile, no, Xo, cTS)

!           Display progress
               IF (i.EQ.j .AND. iProg) THEN
                  WRITE(*,'(A)',ADVANCE='NO') (k-1)*20//"%  "
                  k = k + 1
                  j = NINT(REAL((k-1)*nTS,KIND=8)/5.0D0)
               END IF
            END DO
         END SELECT
      END SELECT

      WRITE(*,'(A)')
      IF (cep%odeS%tIntType .EQ. tIntType_CN2) THEN
         WRITE(*,'(4X,A)') "-------------------------------------------"
         WRITE(*,'(4X,A)') " No of Newton-Raphson fails: "//STR(IPAR(2))
         WRITE(*,'(4X,A)') "-------------------------------------------"
      END IF

      IF (ALLOCATED(IPAR)) DEALLOCATE(IPAR)
      IF (ALLOCATED(RPAR)) DEALLOCATE(RPAR)
      IF (ALLOCATED(Xo)) DEALLOCATE(Xo)

      RETURN
      END SUBROUTINE CEPINTEG
!-----------------------------------------------------------------------
      SUBROUTINE GET_STIM_AMP(time, cepL, stim_amp, flag)
      USE CEPMOD
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: time
      TYPE(cepModelType), INTENT(INOUT) :: cepL
      REAL(KIND=8), INTENT(OUT) :: stim_amp
      LOGICAL, INTENT(OUT) :: flag

      flag = .FALSE.
#ifdef S1S2REST
      IF (cepL%S1S2%cntr .LT. cepL%S1S2%nrep) THEN
#endif
         stim_amp = 0D0
         IF (time.GE.cepL%Istim%Ts .AND. time.LE.cepL%Istim%Te) THEN
            stim_amp = cepL%Istim%A
         END IF

         IF (time .GE. cepL%Istim%Te) THEN
            stim_amp = 0D0
            cepL%Istim%Ts = cepL%Istim%Ts + cepL%Istim%CL
            cepL%Istim%Te = cepL%Istim%Ts + cepL%Istim%Td
#ifdef S1S2REST
            cepL%S1S2%cntr = cepL%S1S2%cntr + 1
#endif
         END IF

#ifdef S1S2REST
      ELSE IF (cepL%S1S2%cntr .EQ. cepL%S1S2%nrep) THEN

         IF (time.GE.cepL%Istim%Ts .AND. time.LE.cepL%Istim%Te) THEN
            stim_amp = cepL%Istim%A
         END IF

         IF (time .GE. cepL%Istim%Te) THEN
            stim_amp = 0D0
            cepL%Istim%Ts = cepL%Istim%Ts + cepL%S1S2%APD + cepL%S1S2%DI
            cepL%Istim%Te = cepL%Istim%Ts + cepL%Istim%Td
            cepL%S1S2%cntr = cepL%S1S2%cntr + 1
            WRITE(*,'(4X,A,F12.2,4X,A,F8.2)')
     2         "Ts_S2: ", cepL%Istim%Ts, " DI: ", cepL%S1S2%DI
         END IF

      ELSE IF (cep%S1S2%cntr .EQ. cep%S1S2%nrep+1) THEN

         IF (time.GE.cepL%Istim%Ts .AND. time.LE.cepL%Istim%Te) THEN
            stim_amp = cepL%S1S2%Istim_A
         END IF

         IF (time .GE. cepL%Istim%Te) THEN
            stim_amp = 0D0
            cepL%Istim%Ts = cepL%Istim%Ts + cepL%Istim%CL
            cepL%Istim%Te = cepL%Istim%Ts + cepL%Istim%Td
            cepL%S1S2%cntr = 0

!           Adjust DI
            IF (cepL%S1S2%DI .GT. 1000.0D0) THEN
               cepL%S1S2%DI = cepL%S1S2%DI - 1000.0D0
            ELSE IF (cepL%S1S2%DI .GT. 300.0D0) THEN
               cepL%S1S2%DI = cepL%S1S2%DI - 100.0D0
            ELSE IF (cepL%S1S2%DI .GT. 150.0D0) THEN
               cepL%S1S2%DI = cepL%S1S2%DI - 20.0D0
            ELSE IF (cepL%S1S2%DI .GT. 50.0D0) THEN
               cepL%S1S2%DI = cepL%S1S2%DI - 10.0D0
            ELSE IF (cepL%S1S2%DI .GT. 5.0D0) THEN
               cepL%S1S2%DI = cepL%S1S2%DI - 5.0D0
            ELSE
               flag = .TRUE.
               WRITE(*,'(A)')
               WRITE(*,'(4X,A)') "S1S2 restitution protocol complete!"
               RETURN
            END IF
         END IF
      END IF
#endif

      RETURN
      END SUBROUTINE GET_STIM_AMP
!####################################################################

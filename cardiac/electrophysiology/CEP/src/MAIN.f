!#######################################################################
      PROGRAM CEPINTEG
      USE CEPMOD
      IMPLICIT NONE

      LOGICAL :: flag = .FALSE.
      INTEGER fid, i, j, k, ic, nWTF
      INTEGER c1, c2, cmax, crate
      REAL(KIND=8) cTS, Ts, Te, wtime, Istim, epsX
      CHARACTER(LEN=stdL) fName

      INTEGER, ALLOCATABLE :: IPAR(:)
      REAL(KIND=8), ALLOCATABLE :: WTF(:), RPAR(:)

      i = IARGC()
      IF (i .EQ. 0) THEN
         STOP "Input file required specifying options"
      ELSE IF (i .GT. 1) THEN
         STOP "Too many arguments"
      END IF
      CALL GETARG(1, fName)

      CALL READINPUTS(fName)

      nWTF = nX
      IF (cepType .EQ. cepType_TTP) THEN
         nWTF = nWTF + 16
         ALLOCATE(RPAR(18))
      ELSE IF (cepType .EQ. cepType_BO) THEN
         nWTF = nWTF + 3
         ALLOCATE(RPAR(5))
      ELSE
         ALLOCATE(RPAR(2))
      END IF
      RPAR(:) = 0D0

      IF (emCpl) nWTF = nWTF + 2
      ALLOCATE(WTF(nWTF))
      WTF(:) = 0D0

      IF (tIntType .EQ. tIntType_CN2) THEN
         ALLOCATE(IPAR(2))
         flag    = .TRUE.
         IPAR(1) = MAXITR
         IPAR(2) = 0

         RPAR(1) = ATOL
         RPAR(2) = RTOL
      END IF

      cTS  = 0D0
      j    = 1
      k    = 1
      CALL SYSTEM_CLOCK(count_max=cmax)
      CALL SYSTEM_CLOCK(count_rate=crate)
      CALL SYSTEM_CLOCK(c1)

      SELECT CASE (cepType)
      CASE (cepType_AP)
!        Clear log file
         WRITE(fName,'(A)') "log_AP.txt"
         CALL CTXT(fName)

!        Initialize state variables
         WRITE(*,'(4X,A)') "Initializing Aliev-Panfilov model"
         CALL AP_INIT(nX, X)

!        Loop over time
         WRITE(*,'(4X,A)',ADVANCE='NO') "Time integration progress: "
         DO i=1, nTS
!           Apply external stimulus
            ic = FLOOR(cTS/BCL)
            Ts = Tstim_start + REAL(ic,KIND=8)*BCL
            Te = Ts + Tstim_dur
            IF (cTS.GE.Ts-eps .AND. cTS.LE.Te+eps) THEN
               X(1) = X(1) + stim_amp
            END IF

!           Time integration
            SELECT CASE (tIntType)
            CASE (tIntType_FE)
               CALL AP_INTEGFE(nX, X, cTS, dt)

            CASE (tIntType_RK4)
               CALL AP_INTEGRK(nX, X, cTS, dt)

            CASE (tIntType_CN2)
               CALL AP_INTEGCN2(nX, X, cTS, dt, IPAR, RPAR)
            END SELECT
            WTF(1:nX) = X(:)

!           Perform NaN check
            IF (ISNAN(X(1))) STOP "NaN occurence (X). Aborted!"

!           Activation Coupling
            IF (emCpl) THEN
               CALL AP_ACTCPL(nX, X, dt, epsX, Tact)
               IF (ISNAN(Tact)) STOP "NaN occurence (Tact). Aborted!"
               WTF(nWTF-1) = epsX
               WTF(nWTF)   = Tact
            END IF

!           Time step advance
            cTS = cTS + dt

!           Write state variables to a log file
            CALL WTXT(fName, nWTF, WTF, cTS)

!           Display progress
            IF (i .EQ. j) THEN
               WRITE(*,'(A)',ADVANCE='NO') (k-1)*20//"%  "
               k = k + 1
               j = NINT(REAL((k-1)*nTS,KIND=8)/5.0D0)
            END IF
         END DO

      CASE (cepType_FN)
!        Clear log file
         WRITE(fName,'(A)') "log_FN.txt"
         CALL CTXT(fName)

!        Initialize state variables
         WRITE(*,'(4X,A)') "Initializing Fitzhugh-Nagumo model"
         CALL FN_INIT(nX, X)

!        Loop over time
         WRITE(*,'(4X,A)',ADVANCE='NO') "Time integration progress: "
         DO i=1, nTS
!           Apply external stimulus
            ic = FLOOR(cTS/BCL)
            Ts = Tstim_start + REAL(ic,KIND=8)*BCL
            Te = Ts + Tstim_dur
            IF (cTS.GE.Ts-eps .AND. cTS.LE.Te+eps) THEN
               X(1) = X(1) + stim_amp
            END IF

!           Time integration
            SELECT CASE (tIntType)
            CASE (tIntType_FE)
               CALL FN_INTEGFE(nX, X, cTS, dt)

            CASE (tIntType_RK4)
               CALL FN_INTEGRK(nX, X, cTS, dt)

            CASE (tIntType_CN2)
               CALL FN_INTEGCN2(nX, X, cTS, dt, IPAR, RPAR)
            END SELECT
            WTF(1:nX) = X(:)

!           Perform NaN check
            IF (ISNAN(X(1))) STOP "NaN occurence. Aborted!"

!           Time step advance
            cTS = cTS + dt

!           Write state variables to a log file
            CALL WTXT(fName, nWTF, WTF, cTS)

!           Display progress
            IF (i .EQ. j) THEN
               WRITE(*,'(A)',ADVANCE='NO') (k-1)*20//"%  "
               k = k + 1
               j = NINT(REAL((k-1)*nTS,KIND=8)/5.0D0)
            END IF
         END DO

      CASE (cepType_TTP)
!        Clear log file
         IF (mzone .EQ. 1) THEN
            WRITE(fName,'(A)') "log_TTP_epi.txt"
         ELSE IF (mzone .EQ. 2) THEN
            WRITE(fName,'(A)') "log_TTP_endo.txt"
         ELSE IF (mzone .EQ. 3) THEN
            WRITE(fName,'(A)') "log_TTP_myo.txt"
         END IF
         CALL CTXT(fName)

!        Initialize state variables
         WRITE(*,'(4X,A)') "Initializing tenTusscher-Panfilov model"
         CALL TTP_INIT(mzone, nX, X)

!        Loop over time
         WRITE(*,'(4X,A)',ADVANCE='NO') "Time integration progress: "
         DO i=1, nTS
!           Apply external stimulus
            ic = FLOOR(cTS/BCL)
            Ts = Tstim_start + REAL(ic,KIND=8)*BCL
            Te = Ts + Tstim_dur
            IF (cTS.GE.Ts-eps .AND. cTS.LE.Te+eps) THEN
               Istim = -stim_amp
            ELSE
               Istim = 0.0D0
            END IF

!           Time integration
            SELECT CASE (tIntType)
            CASE (tIntType_FE)
               CALL TTP_INTEGFE(mzone, nX, X, cTS, dt, Istim, RPAR)

            CASE (tIntType_RK4)
               CALL TTP_INTEGRK(mzone, nX, X, cTS, dt, Istim, RPAR)

            CASE (tIntType_CN2)
               CALL TTP_INTEGCN2(mzone, nX, X, cTS, dt, Istim, IPAR,
     2            RPAR)

            END SELECT
            WTF(1:nX) = X(:)
            WTF(nX+1:nX+16) = RPAR(3:18)

!           Perform NaN check
            IF (ISNAN(X(1))) STOP "NaN occurence (X). Aborted!"

!           Activation Coupling
            IF (emCpl) THEN
               CALL TTP_ACTCPL(nX, X, dt, epsX, Tact)
               IF (ISNAN(Tact)) STOP "NaN occurence (Tact). Aborted!"
               WTF(nWTF-1) = epsX
               WTF(nWTF)   = Tact
            END IF

!           Time step advance
            cTS = cTS + dt

!           Write state variables to a log file
            CALL WTXT(fName, nWTF, WTF, cTS)

!           Display progress
            IF (i .EQ. j) THEN
               WRITE(*,'(A)',ADVANCE='NO') (k-1)*20//"%  "
               k = k + 1
               j = NINT(REAL((k-1)*nTS,KIND=8)/5.0D0)
            END IF
         END DO

      CASE (cepType_BO)
!        Clear log file
         IF (mzone .EQ. 1) THEN
            WRITE(fName,'(A)') "log_BO_epi.txt"
         ELSE IF (mzone .EQ. 2) THEN
            WRITE(fName,'(A)') "log_BO_endo.txt"
         ELSE IF (mzone .EQ. 3) THEN
            WRITE(fName,'(A)') "log_BO_myo.txt"
         END IF
         CALL CTXT(fName)

!        Initialize state variables
         WRITE(*,'(4X,A)') "Initializing Bueno-Orovio model"
         CALL BO_INIT(nX, X)

!        Loop over time
         WRITE(*,'(4X,A)',ADVANCE='NO') "Time integration progress: "
         DO i=1, nTS
!           Apply external stimulus
            ic = FLOOR(cTS/BCL)
            Ts = Tstim_start + REAL(ic,KIND=8)*BCL
            Te = Ts + Tstim_dur
            IF (cTS.GE.Ts-eps .AND. cTS.LE.Te+eps) THEN
               Istim = -stim_amp
            ELSE
               Istim = 0.0D0
            END IF

!           Time integration
            SELECT CASE (tIntType)
            CASE (tIntType_FE)
               CALL BO_INTEGFE(mzone, nX, X, cTS, dt, Istim, RPAR)

            CASE (tIntType_RK4)
               CALL BO_INTEGRK(mzone, nX, X, cTS, dt, Istim, RPAR)

            CASE (tIntType_CN2)
               CALL BO_INTEGCN2(mzone, nX, X, cTS, dt, Istim, IPAR,RPAR)
            END SELECT
            WTF(1:nX) = X(:)
            WTF(nX+1:nX+3) = RPAR(3:5)

!           Perform NaN check
            IF (ISNAN(X(1))) STOP "NaN occurence (X). Aborted!"

!           Activation Coupling
            IF (emCpl) THEN
               CALL BO_ACTCPL(nX, X, dt, epsX, Tact)
               IF (ISNAN(Tact)) STOP "NaN occurence (Tact). Aborted!"
               WTF(nWTF-1) = epsX
               WTF(nWTF)   = Tact
            END IF

!           Time step advance
            cTS = cTS + dt

!           Write state variables to a log file
            CALL WTXT(fName, nWTF, WTF, cTS)

!           Display progress
            IF (i .EQ. j) THEN
               WRITE(*,'(A)',ADVANCE='NO') (k-1)*20//"%  "
               k = k + 1
               j = NINT(REAL((k-1)*nTS,KIND=8)/5.0D0)
            END IF
         END DO

      END SELECT
      WRITE(*,'(A)')

      CALL SYSTEM_CLOCK(c2)
      wtime = REAL(c2-c1,KIND=8)/REAL(crate,KIND=8)
      WRITE(*,'(4X,A,F8.2,A)') "Time elapsed: ", wtime, "s"

      IF (flag) THEN
         WRITE(*,'(A)') "----------------------------------------------"
         WRITE(*,'(A)') " No of Newton-Raphson fails: "//STR(ipar(2))
         WRITE(*,'(A)') "----------------------------------------------"
      END IF

      RETURN
      END PROGRAM CEPINTEG
!#######################################################################
      SUBROUTINE READINPUTS(fName)
      USE CEPMOD
      IMPLICIT NONE
      CHARACTER(LEN=stdL), INTENT(IN) :: fName

      INTEGER fid, itmp
      CHARACTER(LEN=stdL) :: sTmp

      fid = 1011
      WRITE(*,'(4X,A)') "Reading inputs"
      OPEN(fid, FILE=TRIM(fName))
      READ(fid,*,END=101) ! CEP Model Type !
      READ(fid,*,END=101) sTmp
      SELECT CASE (TRIM(sTmp))
      CASE ("ap", "AP")
         cepType = cepType_AP
         nX = 2
      CASE ("fn", "FN")
         cepType = cepType_FN
         nX = 2
      CASE ("ttp", "TTP")
         cepType = cepType_TTP
         nX = 19
      CASE ("bo", "BO")
         cepType = cepType_BO
         nX = 4
      CASE DEFAULT
         STOP "ERROR: Unknown electrophysiology model"
      END SELECT

      READ(fid,*,END=101) ! Myocardium zone
      READ(fid,*,END=101) mzone

      READ(fid,*,END=101) ! Time integrator !
      READ(fid,*,END=101) sTmp
      SELECT CASE (TRIM(sTmp))
      CASE ("fe", "FE", "Euler")
         tIntType = tIntType_FE
      CASE ("rk4", "RK4", "Runge")
         tIntType = tIntType_RK4
      CASE ("cn", "cn2", "CN", "CN2")
         tIntType = tIntType_CN2
      CASE DEFAULT
         STOP "ERROR: Unknown time integration scheme"
      END SELECT

      READ(fid,*,END=101) ! Basic cycle length !
      READ(fid,*,END=101) BCL
      READ(fid,*,END=101) ! No. of time steps !
      READ(fid,*,END=101) nTS
      READ(fid,*,END=101) ! time increment !
      READ(fid,*,END=101) dt
      READ(fid,*,END=101) ! External stimulus start time !
      READ(fid,*,END=101) Tstim_start
      READ(fid,*,END=101) ! External stimulus duration !
      READ(fid,*,END=101) Tstim_dur
      READ(fid,*,END=101) ! External stimulus amplitude!
      READ(fid,*,END=101) stim_amp
      READ(fid,*,END=101) ! Electro-Mechanics coupling !
      READ(fid,*,END=101) itmp
      IF (itmp .NE. 0) emCpl = .TRUE.
 101  CLOSE(fid)

      IF (emCpl .AND. cepType .EQ. cepType_FN) THEN
         STOP "ERROR: EM coupling is not allowed for Fitzhugh-Nagumo"//
     2      " model"
      END IF

      IF (mzone.LT.1 .OR. mzone.GT.3) THEN
         STOP "ERROR: invalid myocardial zone specified"
      END IF

      IF (mzone .GT. 1) THEN
         IF (cepType.NE.cepType_TTP .AND. cepType.NE.cepType_BO) THEN
            STOP "ERROR: mid-myocardium and endocardium zones are "//
     2         " allowed for tenTuscher-Panfilov and Bueno-Orovio "//
     3         " models only"
         END IF
      END IF

      ALLOCATE(X(nX))
      X = 0.0D0

      RETURN
      END SUBROUTINE READINPUTS
!#######################################################################
      SUBROUTINE CTXT(fName)
      USE CEPMOD, ONLY : stdL
      IMPLICIT NONE
      CHARACTER(LEN=stdL), INTENT(IN) :: fName

      INTEGER fid

      fid = 1011
      OPEN(fid, FILE=TRIM(fname), STATUS='UNKNOWN')
      CLOSE(fid,STATUS='DELETE')

      OPEN(fid, FILE=TRIM(fname), STATUS='NEW')
      CLOSE(fid)

      RETURN
      END SUBROUTINE CTXT
!-----------------------------------------------------------------------
      SUBROUTINE WTXT(fname, n, X, t)
      USE CEPMOD, ONLY : stdL
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      REAL(KIND=8), INTENT(IN) :: t, X(n)
      CHARACTER(LEN=stdL), INTENT(IN) :: fName

      INTEGER i, fid

      fid = 1011
      OPEN(fid, FILE=TRIM(fName), STATUS='OLD', POSITION='APPEND')
      WRITE(fid,'(2X,F10.4)',ADVANCE='NO') t
      DO i=1, n
         WRITE(fid,'(1X,1pE18.9)',ADVANCE='NO') X(i)
      END DO
      WRITE(fid,'(A)')
      CLOSE(fid)

      RETURN
      END SUBROUTINE WTXT
!#######################################################################


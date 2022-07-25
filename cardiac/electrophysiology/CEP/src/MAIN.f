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
!--------------------------------------------------------------------
!
!
!--------------------------------------------------------------------

      PROGRAM CEPMAIN
      USE CEPMOD
      IMPLICIT NONE

      LOGICAL flag
      INTEGER i, c1, c2, cmax, crate
      REAL(KIND=8) t1, t2, wtime
      CHARACTER(LEN=stdL) fName

      i = IARGC()
      IF (i .EQ. 0) THEN
         STOP "Input file required specifying options"
      ELSE IF (i .GT. 1) THEN
         STOP "Too many arguments"
      END IF
      CALL GETARG(1, fName)

      INQUIRE(FILE=TRIM(fName), EXIST=flag)
      IF (.NOT.flag) THEN
         WRITE(*,'(A)') "ERROR: Input file <"//TRIM(fName)//
     2   "> does not exist"
         STOP
      END IF

      CALL READINPUTS(fName)

c      CALL SYSTEM_CLOCK(count_max=cmax)
c      CALL SYSTEM_CLOCK(count_rate=crate)
c      CALL SYSTEM_CLOCK(c1)
      CALL CPU_TIME(t1)

      CALL CEPINIT()

      CALL CEPINTEG()

      DEALLOCATE(X, Xg)

c      CALL SYSTEM_CLOCK(c2)
c      wtime = REAL(c2-c1,KIND=8)/REAL(crate,KIND=8)
      CALL CPU_TIME(t2)
      wtime = t2 - t1

      WRITE(*,'(A)')
      WRITE(*,'(4X,A,F8.2,A)') "Time elapsed: ", wtime, "s"
      WRITE(*,'(A)')

      RETURN
      END PROGRAM CEPMAIN
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
      CALL TO_LOWER(sTmp)
      SELECT CASE (TRIM(sTmp))
      CASE ("ap")
         cep%cepType = cepModel_AP
         nX = 2
         nG = 0

      CASE ("bo")
         cep%cepType = cepModel_BO
         nX = 4

      CASE ("fn")
         cep%cepType = cepModel_FN
         nX = 2
         nG = 0

      CASE ("ttp")
         cep%cepType = cepModel_TTP
         nX = 7
         nG = 12

      CASE DEFAULT
         STOP "ERROR: Unknown electrophysiology model"
      END SELECT

      READ(fid,*,END=101) ! Myocardium zone: 1-epi; 2-endo; 3-myo
      READ(fid,*,END=101) cep%imyo

      READ(fid,*,END=101) ! Time integrator !
      READ(fid,*,END=101) sTmp
      CALL TO_LOWER(sTmp)
      SELECT CASE (TRIM(sTmp))
      CASE ("fe", "euler")
         cep%odeS%tIntType = tIntType_FE

      CASE ("rk", "rk4", "runge")
         cep%odeS%tIntType = tIntType_RK4

      CASE ("cn", "cn2")
         cep%odeS%tIntType = tIntType_CN2

      CASE DEFAULT
         STOP "ERROR: Unknown time integration scheme"
      END SELECT

!     General time integration parameters
      READ(fid,*,END=101) ! No. of time steps !
      READ(fid,*,END=101) nTS
      READ(fid,*,END=101) ! time increment !
      READ(fid,*,END=101) dt
!     Stimulus parameters
      READ(fid,*,END=101) ! External stimulus amplitude!
      READ(fid,*,END=101) cep%Istim%A
      READ(fid,*,END=101) ! External stimulus start time !
      READ(fid,*,END=101) cep%Istim%Ts
      READ(fid,*,END=101) ! External stimulus duration !
      READ(fid,*,END=101) cep%Istim%Td
      READ(fid,*,END=101) ! Basic cycle length !
      READ(fid,*,END=101) cep%Istim%CL
      cep%Istim%Te = cep%Istim%Ts + cep%Istim%Td

      READ(fid,*,END=101) ! Electro-Mechanics coupling !
      READ(fid,*,END=101) itmp
      IF (itmp .NE. 0) cep%emCpld = .TRUE.

      READ(fid,*,END=101) ! Display progress !
      READ(fid,*,END=101) itmp
      IF (itmp .NE. 0) iProg = .TRUE.

#ifdef S1S2REST
      IF (cep%cepType .NE. cepModel_TTP)
     2   STOP "ERROR: S1S2 protocol applicable for TTP model only"
!     Read S1-S2 protocol inputs
      READ(fid,*,END=101) ! S1-S2 inputs
      READ(fid,*,END=101) ! S1 repeats before S2 stimulus
      READ(fid,*,END=101) cep%S1S2%nrep
      READ(fid,*,END=101) ! Initial Diastolic interval
      READ(fid,*,END=101) cep%S1S2%DI
      READ(fid,*,END=101) ! Basic APD
      READ(fid,*,END=101) cep%S1S2%APD
      READ(fid,*,END=101) ! S2 amplitude
      READ(fid,*,END=101) cep%S1S2%Istim_A
      iProg = .FALSE.
#endif

 101  CLOSE(fid)

      IF (cep%emCpld .AND. cep%cepType .EQ. cepModel_FN) THEN
         STOP "ERROR: EM coupling is not allowed for Fitzhugh-Nagumo"//
     2      " model"
      END IF

      IF (cep%imyo.LT.1 .OR. cep%imyo.GT.3) THEN
         STOP "ERROR: invalid myocardial zone specified"
      END IF

      IF (cep%imyo .GT. 1) THEN
         IF (cep%cepType.NE.cepModel_TTP .AND.
     2       cep%cepType.NE.cepModel_BO) THEN
            STOP "ERROR: mid-myocardium and endocardium zones are "//
     2         " allowed for tenTuscher-Panfilov and Bueno-Orovio "//
     3         " models only"
         END IF
      END IF

      ALLOCATE(X(nX), Xg(nG))
      X  = 00D0
      Xg = 00D0

      RETURN
      END SUBROUTINE READINPUTS
!#######################################################################


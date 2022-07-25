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
!     This module defines data structures for cardiac electrophysiology
!     model equation. It also interfaces with individual modules for
!     the cellular activation model.
!
!-----------------------------------------------------------------------

      MODULE CEPMOD
      USE UTILMOD
      USE APMOD
      USE FNMOD
      USE TTPMOD
      USE BOMOD
      IMPLICIT NONE

!     Type of cardiac electrophysiology models: Aliev-Panfilov model,
!     Bueno-Orovio-Cherry-Fenton model, Fitzhugh-Nagumo model,
!     tenTusscher-Panfilov 2006 model
      INTEGER, PARAMETER :: cepModel_NA = 100, cepModel_AP = 101,
     2   cepModel_BO = 102, cepModel_FN = 103, cepModel_TTP = 104

!     Time integration scheme: Forward-Euler, Runge-Kutta 4th order,
!     Crank-Nicholson
      INTEGER, PARAMETER :: tIntType_NA  = 200, tIntType_FE = 201,
     2   tIntType_RK4 = 202, tIntType_CN2 = 203

!     Time integration scheme and related parameters
      TYPE odeType
!        Time integration method type
         INTEGER :: tIntType = tIntType_NA
!        Max. iterations for Newton-Raphson method
         INTEGER :: maxItr = 5
!        Absolute tolerance
         REAL(KIND=8) :: absTol = 1D-8
!        Relative tolerance
         REAL(KIND=8) :: relTol = 1D-4
      END TYPE odeType

!     External stimulus type
      TYPE stimType
!        start time
         REAL(KIND=8) :: Ts = 0D0
!        duration of stimulus
         REAL(KIND=8) :: Td = 0D0
!        end time
         REAL(KIND=8) :: Te = 0D0
!        cycle length
         REAL(KIND=8) :: CL = 0D0
!        stimulus amplitude
         REAL(KIND=8) :: A = 0D0
      END TYPE stimType

!     S1S2 protocol data type
      TYPE S1S2type
!        Num S1 repeats before S2 stimulus
         INTEGER :: nrep = 0
!        Counter to track S1 repeats
         INTEGER :: cntr = 1
!        Diastolic interval
         REAL(KIND=8) :: DI = 0D0
!        Action potential duration
         REAL(KIND=8) :: APD = 0D0
!        S2 stimulus
         REAL(KIND=8) :: Istim_A
      END TYPE S1S2type

!     Cardiac electrophysiology model type
      TYPE cepModelType
!        EM coupling
         LOGICAL :: emCpld = .FALSE.
!        Type of cardiac electrophysiology model
         INTEGER :: cepType = cepModel_NA
!        Myocardium zone id: 1-epi; 2-endo; 3-myo
         INTEGER :: imyo = 1
!        Constant for stretch-activated-currents
         REAL(KIND=8) :: Ksac = 0D0
!        Activation force
         REAL(KIND=8) :: Tact = 0D0
!        External stimulus
         TYPE(stimType) :: Istim
!        Time integration options
         TYPE(odeType) :: odeS
!        S1S2 type
         TYPE(S1S2type) :: S1S2
      END TYPE cepModelType

!     Display progress
      LOGICAL iProg
!     No. of time steps
      INTEGER nTS
!     Number of state variables
      INTEGER :: nX = 0
!     Number of gating variables
      INTEGER :: nG = 0
!     Time step increment
      REAL(KIND=8) :: dt = 0D0
!     Output log file
      CHARACTER(LEN=stdL) :: oFile
!     State variables
      REAL(KIND=8), ALLOCATABLE :: X(:)
!     Gating variables
      REAL(KIND=8), ALLOCATABLE :: Xg(:)
!     Cardiac electrophysiology type
      TYPE(cepModelType) :: cep

      END MODULE CEPMOD
!#######################################################################

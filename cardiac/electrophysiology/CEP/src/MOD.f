!#######################################################################
      MODULE CEPMOD
      USE UTILMOD
      USE APMOD
      USE FNMOD
      USE TTPMOD
      USE BOMOD
      IMPLICIT NONE

!-----------------------------------------------------------------------
!     Parameters
!     Electrophysiology model
!     Aliev-Panfilov model, Fitzhugh-Nagumo model,
!     tenTusscher-Panfilov model (2006), Bueno-Orovio model
      INTEGER, PARAMETER :: cepType_NA  = 100
      INTEGER, PARAMETER :: cepType_AP  = 101
      INTEGER, PARAMETER :: cepType_FN  = 102
      INTEGER, PARAMETER :: cepType_TTP = 103
      INTEGER, PARAMETER :: cepType_BO  = 104

!     Time integration scheme
!     Forward-Euler, Runge-Kutta 4th order, Crank-Nicholson,
      INTEGER, PARAMETER :: tIntType_NA  = 200
      INTEGER, PARAMETER :: tIntType_FE  = 201
      INTEGER, PARAMETER :: tIntType_RK4 = 202
      INTEGER, PARAMETER :: tIntType_CN2 = 203

!     Newton-Raphson parameters
!     Max. iterations
      INTEGER, PARAMETER :: MAXITR = 5
!     Absolute tolerance
      REAL(KIND=8), PARAMETER :: ATOL = 1D-8
!     Relative tolerance
      REAL(KIND=8), PARAMETER :: RTOL = 1D-4

!-----------------------------------------------------------------------
!     User-defined inputs

!     Type of electrophysiology model
      INTEGER :: cepType = cepType_NA
!     Type of time integrator
      INTEGER :: tIntType = tIntType_NA
!     Myocardium zone index: 1-epi; 2-myo; 3-endo
      INTEGER :: mzone = 1
!     No. of time steps
      INTEGER :: nTS = 0
!     Basic cycle length
      REAL(KIND=8) :: BCL = 0D0
!     Time step increment
      REAL(KIND=8) :: dt = 0D0
!     Time step to apply external stimulus
      REAL(KIND=8) :: Tstim_start = 0D0
!     Stimulus duration
      REAL(KIND=8) :: Tstim_dur = 0D0
!     Stimulus amplitude
      REAL(KIND=8) :: stim_amp = 0D0
!     Electromechanics - activation
!     Activation force
      REAL(KIND=8) :: Tact = 0D0
!     EM coupling
      LOGICAL :: emCpl = .FALSE.

!-----------------------------------------------------------------------
!     No. of state variables
      INTEGER :: nX = 0
!     State variables
      REAL(KIND=8), ALLOCATABLE :: X(:)
!-----------------------------------------------------------------------

      END MODULE CEPMOD
!#######################################################################

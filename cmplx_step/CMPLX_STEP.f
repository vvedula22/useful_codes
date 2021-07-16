!#######################################################################
      MODULE COM

      INTEGER, PARAMETER :: IK1  = SELECTED_INT_KIND(2)
      INTEGER, PARAMETER :: IK2  = SELECTED_INT_KIND(4)
      INTEGER, PARAMETER :: IK4  = SELECTED_INT_KIND(9)
      INTEGER, PARAMETER :: IK8  = SELECTED_INT_KIND(18)

      INTEGER, PARAMETER :: RK4  = SELECTED_REAL_KIND(6, 37)
      INTEGER, PARAMETER :: RK8  = SELECTED_REAL_KIND(15, 307)
      INTEGER, PARAMETER :: RK16 = SELECTED_REAL_KIND(33, 4931)

      INTEGER, PARAMETER :: IK = IK4
      INTEGER, PARAMETER :: RK = RK16

      CONTAINS
         FUNCTION GETF(x) RESULT(f)
         IMPLICIT NONE
         REAL(KIND=RK), INTENT(IN) :: x
         REAL(KIND=RK) :: f

         f = EXP(x)/SQRT(SIN(x)**3._RK + COS(x)**3._RK)

         RETURN
         END FUNCTION GETF
!-----------------------------------------------------------------------
         FUNCTION GETFDER(f, x) RESULT(df)
         IMPLICIT NONE
         REAL(KIND=RK), INTENT(IN) :: x, f
         REAL(KIND=RK) :: df

         df = f - 1.5_RK*f*SIN(x)*COS(x)*(SIN(x)-COS(x))/
     2      (SIN(x)**3._RK + COS(x)**3._RK)

         RETURN
         END FUNCTION GETFDER
!-----------------------------------------------------------------------
         FUNCTION GETERR(fn, fe) RESULT(eps)
         IMPLICIT NONE
         REAL(KIND=RK), INTENT(IN) :: fn, fe
         REAL(KIND=RK) :: eps

         eps = ABS(fn - fe)/ABS(fe)

         RETURN
         END FUNCTION GETERR
!-----------------------------------------------------------------------
         FUNCTION GETFCMPLX(z) RESULT(f)
         IMPLICIT NONE
         COMPLEX(KIND=RK), INTENT(IN) :: z
         COMPLEX(KIND=RK) :: f

         f = EXP(z)/SQRT(SIN(z)**3._RK + COS(z)**3._RK)

         RETURN
         END FUNCTION GETFCMPLX
!-----------------------------------------------------------------------
      END MODULE COM

!#######################################################################

      PROGRAM CMPLX_STEP
      USE COM
      IMPLICIT NONE
      INTEGER(KIND=IK), PARAMETER :: fid = 1056
      REAL(KIND=RK), PARAMETER :: x = 1.5

      INTEGER(KIND=IK) :: i
      CHARACTER(LEN=256) :: fname
      REAL(KIND=RK) :: fr, dfe, dx, dfn, dfc, en, ec
      COMPLEX(KIND=RK) :: z

      fr  = GETF(x)
      dfe = GETFDER(fr, x)

      WRITE(*,'(A)') REPEAT('=',72)
      WRITE(*,'(4X,A,3X,1pE24.16)') "f @ x=1.5  = ", fr
      WRITE(*,'(4X,A,3X,1pE24.16)') "df (exact) = ", dfe
      WRITE(*,'(A)') REPEAT('=',72)

      WRITE(fname,'(A,I3.3,A)') "error_fp", RK, ".dat"
      OPEN(fid, FILE=TRIM(fname))
      WRITE(fid,'(A)') "VARIABLES=h, FDIFF, CMPLX"
      DO i=1, 48
         dx  = 10**(-REAL(i,KIND=RK))
         dfn = (GETF(x+dx) - fr)/dx
         en  = GETERR(dfn, dfe)

         z   = CMPLX(x, dx, KIND=RK)
         dfc = IMAG(GETFCMPLX(z)) / dx
         ec  = GETERR(dfc, dfe)

         WRITE(*,101) dx, en, ec
         WRITE(fid,101) dx, en, ec

 101     FORMAT(8X,1pE12.3,2(3X,1pE24.16))
      END DO
      CLOSE(fid)

      RETURN
      END PROGRAM
!-----------------------------------------------------------------------


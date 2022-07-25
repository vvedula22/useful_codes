!####################################################################
!     Matrix operations
      MODULE MATFUN
      IMPLICIT NONE

      CONTAINS
!--------------------------------------------------------------------
!     Create a second order identity matrix of rank nd
      FUNCTION MAT_ID(nd) RESULT(A)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nd
      REAL(KIND=8) :: A(nd,nd)

      INTEGER :: i

      A = 0D0
      DO i=1, nd
         A(i,i) = 1D0
      END DO

      RETURN
      END FUNCTION MAT_ID
!--------------------------------------------------------------------
!     Trace of second order matrix of rank nd
      FUNCTION MAT_TRACE(A, nd) RESULT(trA)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nd
      REAL(KIND=8), INTENT(IN) :: A(nd,nd)
      REAL(KIND=8) :: trA

      INTEGER :: i

      trA = 0D0
      DO i=1, nd
         trA = trA + A(i,i)
      END DO

      RETURN
      END FUNCTION MAT_TRACE
!--------------------------------------------------------------------
!     Create a matrix from outer product of two vectors
      FUNCTION MAT_DYADPROD(u, v, nd) RESULT(A)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nd
      REAL(KIND=8), INTENT(IN) :: u(nd), v(nd)
      REAL(KIND=8) :: A(nd,nd)

      INTEGER :: i, j

      DO j=1, nd
         DO i=1, nd
            A(i,j) = u(i)*v(j)
         END DO
      END DO

      RETURN
      END FUNCTION MAT_DYADPROD
!--------------------------------------------------------------------
!     Computes the determinant of a square matrix
      RECURSIVE FUNCTION MAT_DET(A, nd) RESULT(D)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nd
      REAL(KIND=8), INTENT(IN) :: A(nd,nd)

      INTEGER :: i, j, n
      REAL(KIND=8) :: D, Am(nd-1,nd-1)

      D = 0D0
      IF (nd .EQ. 2) THEN
         D = A(1,1)*A(2,2) - A(1,2)*A(2,1)
      ELSE
         DO i=1, nd
            n = 0
            DO j=1, nd
               IF (i .EQ. j) THEN
                  CYCLE
               ELSE
                  n = n + 1
                  Am(:,n) = A(2:nd,j)
               END IF
            END DO ! j
            D = D + ( (-1D0)**REAL(1+i,KIND=8) * A(1,i) *
     2                 MAT_DET(Am,nd-1) )
         END DO ! i
      END IF ! nd.EQ.2

      RETURN
      END FUNCTION MAT_DET
!--------------------------------------------------------------------
!     This function uses LAPACK library to compute eigen values of a
!     square matrix, A of dimension nd
      FUNCTION MAT_EIG(A, nd)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nd
      REAL(KIND=8), INTENT(IN) :: A(nd,nd)
      COMPLEX*16 :: Amat(nd,nd), MAT_EIG(nd), b(nd), DUMMY(1,1),
     2   WORK(2*nd)

      INTEGER :: i, j, iok

      Amat = (0D0, 0D0)
      DO j=1, nd
         DO i=1, nd
            Amat(i,j) = CMPLX(A(i,j))
         END DO
      END DO

      CALL ZGEEV('N', 'N', nd, Amat, nd, b, DUMMY, 1, DUMMY, 1, WORK,
     2   2*nd, WORK, iok)
      IF (iok .NE. 0) THEN
         STOP "Failed to compute eigen values"
      ELSE
         MAT_EIG(:) = b(:)
      END IF

      RETURN
      END FUNCTION MAT_EIG
!--------------------------------------------------------------------
!     This function computes inverse of a square matrix
      FUNCTION MAT_INV(A, nd) RESULT(Ainv)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nd
      REAL(KIND=8), INTENT(IN) :: A(nd,nd)
      REAL(KIND=8) :: Ainv(nd,nd)

      REAL(KIND=8), PARAMETER :: epsil = EPSILON(epsil)

      REAL(KIND=8) :: d

      IF (nd .EQ. 2) THEN
         d = MAT_DET(A, nd)
         IF (ABS(d) .LT. 1D2*epsil) THEN
            STOP "Singular matrix detected to compute inverse"
         END IF

         Ainv(1,1) =  A(2,2)/d
         Ainv(1,2) = -A(1,2)/d

         Ainv(2,1) = -A(2,1)/d
         Ainv(2,2) =  A(1,1)/d

         RETURN
      ELSE IF (nd .EQ. 3) THEN
         d = MAT_DET(A, nd)
         IF (ABS(d) .LT. 1D2*epsil) THEN
            STOP "Singular matrix detected to compute inverse"
         END IF

         Ainv(1,1) = (A(2,2)*A(3,3)-A(2,3)*A(3,2)) / d
         Ainv(1,2) = (A(1,3)*A(3,2)-A(1,2)*A(3,3)) / d
         Ainv(1,3) = (A(1,2)*A(2,3)-A(1,3)*A(2,2)) / d

         Ainv(2,1) = (A(2,3)*A(3,1)-A(2,1)*A(3,3)) / d
         Ainv(2,2) = (A(1,1)*A(3,3)-A(1,3)*A(3,1)) / d
         Ainv(2,3) = (A(1,3)*A(2,1)-A(1,1)*A(2,3)) / d

         Ainv(3,1) = (A(2,1)*A(3,2)-A(2,2)*A(3,1)) / d
         Ainv(3,2) = (A(1,2)*A(3,1)-A(1,1)*A(3,2)) / d
         Ainv(3,3) = (A(1,1)*A(2,2)-A(1,2)*A(2,1)) / d

         RETURN
      ELSE IF (nd.GT.3 .AND. nd.LT.10) THEN
         d = MAT_DET(A, nd)
         IF (ABS(d) .LT. 1D2*epsil) THEN
            STOP "Singular matrix detected to compute inverse"
         END IF

         Ainv = MAT_INV_GE(A, nd)

         RETURN
      ELSE
         Ainv = MAT_INV_LP(A, nd)
      END IF

      RETURN
      END FUNCTION MAT_INV
!--------------------------------------------------------------------
!     This function computes inverse of a square matrix using
!     Gauss Elimination method
      FUNCTION MAT_INV_GE(A, nd) RESULT(Ainv)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nd
      REAL(KIND=8), INTENT(IN) :: A(nd,nd)
      REAL(KIND=8) :: Ainv(nd,nd)

      INTEGER :: i, j, k
      REAL(KIND=8) :: d, B(nd,2*nd)

!     Auxillary matrix
      B = 0D0
      DO i=1, nd
         DO j=1, nd
            B(i,j) = A(i,j)
         END DO
         B(i,nd+i) = 1D0
      END DO

!     Pivoting
      DO i=nd, 2, -1
         IF (B(i,1) .GT. B(i-1,1)) THEN
            DO j=1, 2*nd
               d = B(i,j)
               B(i,j) = B(i-1,j)
               B(i-1,j) = d
            END DO
         END IF
      END DO

!     Do row-column operations and reduce to diagonal
      DO i=1, nd
         DO j=1, nd
            IF (j .NE. i) THEN
               d = B(j,i)/B(i,i)
               DO k=1, 2*nd
                  B(j,k) = B(j,k) - d*B(i,k)
               END DO
            END IF
         END DO
      END DO

!     Unit matrix
      DO i=1, nd
         d = B(i,i)
         DO j=1, 2*nd
            B(i,j) = B(i,j)/d
         END DO
      END DO

!     Inverse
      DO i=1, nd
         DO j=1, nd
            Ainv(i,j) = B(i,j+nd)
         END DO
      END DO

      RETURN
      END FUNCTION MAT_INV_GE
!--------------------------------------------------------------------
!     This function computes inverse of a square matrix using
!     Lapack functions (DGETRF + DGETRI)
      FUNCTION MAT_INV_LP(A, nd) RESULT(Ainv)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nd
      REAL(KIND=8), INTENT(IN) :: A(nd, nd)
      REAL(KIND=8) :: Ainv(nd,nd)

      INTEGER iok, IPIV(nd)
      REAL(KIND=8) :: Ad(nd,nd), WORK(2*nd)

      Ad = A
      CALL DGETRF(nd, nd, Ad, nd, IPIV, iok)
      IF (iok .NE. 0) THEN
         STOP "ERROR: Singular matrix detected to compute inverse"
      END IF

      CALL DGETRI(nd, Ad, nd, IPIV, WORK, 2*nd, iok)
      IF (iok .NE. 0) THEN
         STOP "ERROR: Matrix inversion failed"
      END IF

      Ainv(:,:) = Ad(:,:)

      RETURN
      END FUNCTION MAT_INV_LP
!--------------------------------------------------------------------
      SUBROUTINE MAT_DISP(A, m, n)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: m
      INTEGER, INTENT(IN), OPTIONAL :: n
      REAL(KIND=8), INTENT(IN) :: A(:,:)

      INTEGER :: i, j, nr, nc

      nr = m
      nc = nr
      IF (PRESENT(n)) nc = n

      DO i=1, nr
         WRITE(*,'(3X,A)',ADVANCE='NO')
         DO j=1, nc
            WRITE(*,'(1X,1pE15.6,1X,A)',ADVANCE='NO') A(i,j), ' | '
         END DO
         WRITE(*,'(A)')
      END DO

      RETURN
      END SUBROUTINE MAT_DISP
!--------------------------------------------------------------------
      END MODULE MATFUN
!####################################################################

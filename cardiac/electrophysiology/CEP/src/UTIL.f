!#######################################################################
      MODULE UTILMOD
      IMPLICIT NONE

!     This is the standard length for all the strings used in this code
      INTEGER, PARAMETER :: stdL = 400
!     minimum value for a double number that is not considered zero
      REAL(KIND=8), PARAMETER :: eps = EPSILON(eps)

      INTERFACE NORM
         MODULE PROCEDURE NORMS, NORMV
      END INTERFACE NORM

      INTERFACE STR
         MODULE PROCEDURE :: DTSTR, NDTSTR, VDTSTR, RTSTR, NRTSTR,
     2      ITSTR, NITSTR
      END INTERFACE STR

      INTERFACE OPERATOR(//)
         MODULE PROCEDURE :: CNCSL, CNCSI, CNCIS, CNCSR
      END INTERFACE OPERATOR(//)

      CONTAINS
!####################################################################
!     This function will compute second NORM of a vector
      PURE FUNCTION NORMS(U, V)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: U(:)
      REAL(KIND=8), INTENT(IN), OPTIONAL :: V(:)
      REAL(KIND=8) NORMS

      INTEGER i, n

      n = SIZE(U)
      NORMS = 0D0
      IF (PRESENT(V)) THEN
         DO i=1, n
            NORMS = NORMS + U(i)*V(i)
         END DO
      ELSE
         DO i=1, n
            NORMS = NORMS + U(i)*U(i)
         END DO
      END IF

      RETURN
      END FUNCTION NORMS
!--------------------------------------------------------------------
      PURE FUNCTION NORMV(U, V)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: U(:,:)
      REAL(KIND=8), INTENT(IN), OPTIONAL :: V(:,:)
      REAL(KIND=8) NORMV

      INTEGER m, n, i, j

      m = SIZE(U,1)
      n = SIZE(U,2)
      NORMV = 0D0
      IF (PRESENT(V)) THEN
         SELECT CASE(m)
         CASE(1)
            DO i=1, n
               NORMV = NORMV + U(1,i)*V(1,i)
            END DO
         CASE(2)
            DO i=1, n
               NORMV = NORMV + U(1,i)*V(1,i) + U(2,i)*V(2,i)
            END DO
         CASE(3)
            DO i=1, n
               NORMV = NORMV + U(1,i)*V(1,i) + U(2,i)*V(2,i)
     2                       + U(3,i)*V(3,i)
            END DO
         CASE(4)
            DO i=1, n
               NORMV = NORMV + U(1,i)*V(1,i) + U(2,i)*V(2,i)
     2                       + U(3,i)*V(3,i) + U(4,i)*V(4,i)
            END DO
         CASE DEFAULT
            DO i=1, n
               NORMV = NORMV + SUM(U(:,i)*V(:,i))
            END DO
         END SELECT
      ELSE
         SELECT CASE(m)
         CASE(1)
            DO i=1, n
               NORMV = NORMV + U(1,i)*U(1,i)
            END DO
         CASE(2)
            DO i=1, n
               NORMV = NORMV + U(1,i)*U(1,i) + U(2,i)*U(2,i)
            END DO
         CASE(3)
            DO i=1, n
               NORMV = NORMV + U(1,i)*U(1,i) + U(2,i)*U(2,i)
     2                       + U(3,i)*U(3,i)
            END DO
         CASE(4)
            DO i=1, n
               NORMV = NORMV + U(1,i)*U(1,i) + U(2,i)*U(2,i)
     2                       + U(3,i)*U(3,i) + U(4,i)*U(4,i)
            END DO
         CASE DEFAULT
            DO i=1, n
               DO j=1, m
                  NORMV = NORMV + U(j,i)*U(j,i)
               END DO
            END DO
         END SELECT
      END IF

      RETURN
      END FUNCTION NORMV
!####################################################################
!     This routine does the cross product for a two given vector of
!     V1 and V2.
      PURE FUNCTION CROSS(V) RESULT(U)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: V(:,:)
      REAL(KIND=8) U(SIZE(V,1))

      IF (SIZE(V,1) .EQ. 2) THEN
         U(1) =  V(2,1)
         U(2) = -V(1,1)
      ELSE
         U(1) = V(2,1)*V(3,2) - V(3,1)*V(2,2)
         U(2) = V(3,1)*V(1,2) - V(1,1)*V(3,2)
         U(3) = V(1,1)*V(2,2) - V(2,1)*V(1,2)
      END IF

      RETURN
      END FUNCTION CROSS
!####################################################################
!     These functions produce string from numbers. Following is for
!     doubles
      PURE FUNCTION DTSTR(dVal) RESULT(string)

      IMPLICIT NONE

      INTEGER, PARAMETER :: l = 8

      REAL(KIND=8), INTENT(IN) :: dVal
      CHARACTER(LEN=l) string

      string = NDTSTR(dVal,l)

      RETURN
      END FUNCTION DTSTR
!--------------------------------------------------------------------
!     Similar to last one, but with a specified length
      PURE FUNCTION NDTSTR(dVal,l) RESULT(string)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: l
      REAL(KIND=8), INTENT(IN) :: dVal
      CHARACTER(LEN=l) string

      INTEGER ex, pos, i, j, abex, k, exex
      REAL(KIND=8) absn

      IF (l .EQ. 0) RETURN

      IF (dVal .NE. dVal) THEN
         IF (l .GE. 3) THEN
            string = ""
            string(l-2:l) = "NaN"
         ELSE
            string = "NaN"(1:l)
         END IF
         RETURN
      END IF

      absn = ABS(dVal)
      IF (absn .GT. HUGE(absn)) THEN
         IF (l .GE. 8) THEN
            string = ""
            string(l-7:l) = "Infinity"
         ELSE
            string = "Infinity"(1:l)
         END IF
         RETURN
      END IF

      IF (absn .LT. TINY(absn)) THEN
         pos = 1
         IF (l .GE. 3) THEN
            string(1:3) = "0.0"
            pos = 4
         END IF
         DO i=pos, l
            string(i:i) = "0"
         END DO
         RETURN
      END IF

!     This is the exponent of the number
      IF (absn .NE. 0D0) THEN
         ex = FLOOR(LOG10(absn))
      ELSE
         ex = 0
      END IF
      abex = ABS(ex)
!     How many digits exponent has
      IF (ex .NE. 0) THEN
         exex = FLOOR(LOG10(REAL(abex,8))) + 1
      ELSE
         exex = 0
      END IF

!     Number of digits of exponents and at least one for number
      i = exex + 1
!     Negative sign of the number
      IF (dVal .LT. 0D0) i = i + 1
!     For negative sign on exponent
      IF (ex .LT. 0) i = i + 1
!     That is for 'E'
      IF (ex .NE. 0) i = i + 1
!     If the sum of above is greater than the availale slots, it can
!     not be represnted with in the length. So we replace it with stars
      IF (i .GT. l) THEN
         DO i=1, l
            string(i:i) = "*"
         END DO
         RETURN
      END IF

!     Constructing the exponent first
!     pos is the position of the charecter to be written into string
      IF (ex .NE. 0) THEN
!     Writing the digits first
         DO pos=l, l-exex+1,-1
            k = MODULO(abex,10) + 1
            string(pos:pos) = "0123456789"(k:k)
            abex = abex/10
         END DO
!     Then negative sign if necessary
         IF (ex .LT. 0) THEN
            string(pos:pos) = "-"
            pos = pos - 1
         END IF
         string(pos:pos) = "E"
         pos = pos - 1
      ELSE
         pos = l
      END IF

!     l-i is the useful length remaining, beside the first number
      IF (l-i .GE. 1) THEN
         absn = absn*(1D1**(-ex/2))
         absn = absn*(1D1**(l-i-1-ex+ex/2))
         j = pos
         DO pos=j,j-l+i+2,-1
            k = FLOOR(MODULO(absn,1D1)) + 1
            string(pos:pos) = "0123456789"(k:k)
            absn = absn/1D1
         END DO
         string(pos:pos) = "."
         pos = pos - 1
         k = FLOOR(MODULO(absn,1D1)) + 1
         string(pos:pos) = "0123456789"(k:k)
      ELSE ! l-i .EQ. 0
         absn = absn*(1D1**(l-i-ex))
         k = FLOOR(MODULO(absn,1D1)) + 1
         string(pos:pos) = "0123456789"(k:k)
      END IF
      IF (dVal .LT. 0D0) string(1:1) = "-"

      RETURN
      END FUNCTION NDTSTR
!--------------------------------------------------------------------
!     Similar to last one, but with a specified length
      PURE FUNCTION VDTSTR(dVal) RESULT(string)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: dVal(:)
      CHARACTER(LEN=9*SIZE(dVal)) string

      INTEGER i, n

      n = SIZE(dVal)
      string = ""
      DO i=1, n
         string = TRIM(string)//" "//NDTSTR(dVal(i),8)
      END DO

      RETURN
      END FUNCTION VDTSTR
!--------------------------------------------------------------------
!     This is for real numbers
      PURE FUNCTION RTSTR(rVal) RESULT(string)

      IMPLICIT NONE

      INTEGER, PARAMETER :: l = 7

      REAL, INTENT(IN) :: rVal
      CHARACTER(LEN=l) string

      string = NDTSTR(REAL(rVal,8),l)

      RETURN
      END FUNCTION RTSTR
!--------------------------------------------------------------------
!     This is for real numbers with specified length
      PURE FUNCTION NRTSTR(rVal,l) RESULT(string)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: l
      REAL, INTENT(IN) :: rVal
      CHARACTER(LEN=l) string

      string = NDTSTR(REAL(rVal,8), l)

      RETURN
      END FUNCTION NRTSTR
!--------------------------------------------------------------------
!     This is for integers
      PURE FUNCTION ITSTR(iVal) RESULT(string)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: iVal
      INTEGER n
      CHARACTER(LEN=2 - MAX(0,SIGN(1,ival)) +
     2   MAXVAL((/(MIN(ABS(iVal)/10**n,1)*n,n=1,9)/))) string

      INTEGER absn, j, k, is

      absn = ABS(iVal)
      IF (absn .EQ. iVal) THEN
         is = 1
      ELSE
         is = 2
         string(1:1) = "-"
      END IF
      DO j=LEN(string),is,-1
         k = MODULO(absn,10) + 1
         string(j:j) = "0123456789"(k:k)
         absn = absn/10
      END DO

      RETURN
      END FUNCTION ITSTR
!--------------------------------------------------------------------
!     Similar to last one, but with minimum length of l
      PURE FUNCTION NITSTR(iVal, l) RESULT(string)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: iVal, l
      INTEGER n
      CHARACTER(LEN=MAX(2 - MAX(0,SIGN(1,ival)) +
     2   MAXVAL((/(MIN(ABS(iVal)/10**n,1)*n,n=1,9)/)),l)) string

      string = ITSTR(iVal)
      string = ADJUSTR(string)

      RETURN
      END FUNCTION NITSTR
!####################################################################
!     Produces a color
      PURE FUNCTION CLR(iStr,clId) RESULT(oStr)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: iStr
      INTEGER, INTENT(IN), OPTIONAL :: clId
      CHARACTER(LEN=LEN(TRIM(iStr))+9) oStr
!     Colors are: 1: White, 2: Red, 3: Green, 4: Yellow, 5: Blue,
!     6: Magenta, 7: Cyan
      CHARACTER(LEN=2),PARAMETER :: clCdL(7) = (/"29","31","32","33",
     2   "34","35","36"/)

      CHARACTER(LEN=2) clCd

      IF (PRESENT(clId)) THEN
         IF (clId .LE. 0) THEN
            clCd = clCdL(1)
         ELSE
            clCd = clCdL(1+MOD(clId-1,SIZE(clCdL)))
         END IF
      ELSE
         clCd = clCdL(2)
      END IF

      oStr = CHAR(27)//'['//clCd//'m'//TRIM(iStr)//CHAR(27)//'[0m'

      RETURN
      END FUNCTION CLR
!--------------------------------------------------------------------
!     This is for removing color
      PURE FUNCTION RMCLR(iStr) RESULT(oStr)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: iStr
      CHARACTER(LEN=LEN(iStr)) oStr

      INTEGER i, j, l

!     Searching and removing all the color codes
      l    = LEN(TRIM(iStr))
      i    = 0
      j    = 0
      oStr = ""
      DO
         i = i + 1
         IF (i .GT. l-3) EXIT
         IF (iStr(i:i+1) .EQ. CHAR(27)//"[") THEN
            IF (iStr(i+3:i+3) .EQ. "m") THEN
               i = i + 3
               CYCLE
            ELSE IF (i .GT. l-4) THEN
               EXIT
            ELSE IF (iStr(i+4:i+4) .EQ. "m") THEN
               i = i + 4
               CYCLE
            END IF
         END IF
         j = j + 1
         oStr(j:j) = iStr(i:i)
      END DO
      j = j + 1
      oStr(j:j+l-i) = iStr(i:l)

      RETURN
      END FUNCTION RMCLR
!####################################################################
      FUNCTION CNCSL(sVal,lVal)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: sVal
      LOGICAL, INTENT(IN) :: lVal

      CHARACTER(LEN=LEN(sVal) + 10) CNCSL

      IF (lVal) THEN
         CNCSL = sVal//CLR("T",3)
      ELSE
         CNCSL = sVal//CLR("F",3)
      END IF

      RETURN
      END FUNCTION CNCSL
!--------------------------------------------------------------------
!     Attaches strings and integer
      FUNCTION CNCSI(sVal,iVal)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: sVal
      INTEGER, INTENT(IN) :: iVal

      INTEGER n
      CHARACTER(LEN=LEN(sVal) + 11 - MAX(0,SIGN(1,ival)) +
     2   MAXVAL((/(MIN(ABS(iVal)/10**n,1)*n,n=1,9)/))) CNCSI

      CNCSI = sVal//CLR(STR(iVal),3)

      RETURN
      END FUNCTION CNCSI
!--------------------------------------------------------------------
      FUNCTION CNCIS(iVal,sVal)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: iVal
      CHARACTER(LEN=*), INTENT(IN) :: sVal

      INTEGER n
      CHARACTER(LEN=LEN(sVal) + 2 - MAX(0,SIGN(1,ival)) +
     2   MAXVAL((/(MIN(ABS(iVal)/10**n,1)*n,n=1,9)/))) CNCIS

      CNCIS = STR(iVal)//sVal

      RETURN
      END FUNCTION CNCIS
!--------------------------------------------------------------------
      FUNCTION CNCSR(sVal,rVal)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: sVal
      REAL(KIND=8), INTENT(IN) :: rVal

      CHARACTER(LEN=LEN(sVal) + 17) CNCSR

      CNCSR = sVal//CLR(STR(rVal),3)

      RETURN
      END FUNCTION CNCSR
!####################################################################
      END MODULE UTILMOD

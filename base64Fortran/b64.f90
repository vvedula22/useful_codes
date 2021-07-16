!**************************************************

		module params
		logical, parameter :: debug=.true.
		character, parameter :: eol=achar(10)
		character(len=64) :: b64Table = &
		"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789/"

		integer, parameter :: stdout=6
		integer, parameter :: stdl=256

		integer, parameter :: IK1 = selected_int_kind(2)  ! 8 bits
		integer, parameter :: IK2 = selected_int_kind(4)  ! 16 bits
		integer, parameter :: IK4 = selected_int_kind(9)  ! 32 bits
		integer, parameter :: IK8 = selected_int_kind(18) ! 64 bits
		integer, parameter :: IK = IK4  ! default

		integer(IK8), parameter :: MAXIK8 = huge(1_IK8)
		integer(IK4), parameter :: MAXIK4 = huge(1_IK4)
		integer(IK2), parameter :: MAXIK2 = huge(1_IK2)
		integer(IK1), parameter :: MAXIK1 = huge(1_IK1)
		integer(IK ), parameter :: MAXIK  = huge(1_IK )
		
		integer(IK ), parameter :: BITIK8 = bit_size(MAXIK8)
		integer(IK ), parameter :: BITIK4 = bit_size(MAXIK4)
		integer(IK ), parameter :: BITIK2 = bit_size(MAXIK2)
		integer(IK ), parameter :: BITIK1 = bit_size(MAXIK1)
		integer(IK ), parameter :: BITIK  = bit_size(MAXIK )

		integer(IK ), parameter :: BYTIK8 = bit_size(MAXIK8)/8_IK
		integer(IK ), parameter :: BYTIK4 = bit_size(MAXIK4)/8_IK
		integer(IK ), parameter :: BYTIK2 = bit_size(MAXIK2)/8_IK
		integer(IK ), parameter :: BYTIK1 = bit_size(MAXIK1)/8_IK
		integer(IK ), parameter :: BYTIK  = bit_size(MAXIK )/8_IK 

		integer, parameter :: RK4 = selected_real_kind(6,37)    ! 32 bits
		integer, parameter :: RK8 = selected_real_kind(15,307)  ! 64 bits
		integer, parameter :: RK16 = selected_real_kind(33,4931) ! 128 bits
		integer, parameter :: RK = RK4  ! default

		real(RK8), parameter :: MAXRK8 = huge(1._RK8)
		real(RK4), parameter :: MAXRK4 = huge(1._RK4)
		real(RK ), parameter :: MAXRK  = MAXRK4
		
		integer(IK ) :: BITRK8
		integer(IK ) :: BITRK4
		integer(IK ) :: BITRK 
		integer(IK ) :: BITCHR
		integer(IK ) :: BYTRK8
		integer(IK ) :: BYTRK4
		integer(IK ) :: BYTRK 
		integer(IK ) :: BYTCHR

		interface bit_size
			module procedure :: bit_size_r16, bit_size_r8, &
								bit_size_r4, bit_size_chr
		end interface
		
		interface byte_size
			module procedure :: byte_size_r16, byte_size_r8, &
								byte_size_r4, byte_size_chr
		end interface
		
		contains

			!==========================================

			function bit_size_r16(r) result(bits)
			implicit none
			real(RK16), intent(in) :: r
			integer(IK) :: bits
			integer(IK1) :: mold(1)
			
			bits = size(transfer(r,mold),dim=1,kind=IK)*8_IK
			return
			end function bit_size_r16

			!==========================================

			function bit_size_r8(r) result(bits)
			implicit none
			real(RK8), intent(in) :: r
			integer(IK) :: bits
			integer(IK1) :: mold(1)
			
			bits = size(transfer(r,mold),dim=1,kind=IK)*8_IK
			return
			end function bit_size_r8

			!==========================================

			function bit_size_r4(r) result(bits)
			implicit none
			real(RK4), intent(in) :: r
			integer(IK) :: bits
			integer(IK1) :: mold(1)
			
			bits = size(transfer(r,mold),dim=1,kind=IK)*8_IK
			return
			end function bit_size_r4

			!==========================================

			function bit_size_chr(r) result(bits)
			implicit none
			character(*), intent(in) :: r
			integer(IK) :: bits
			integer(IK1) :: mold(1)
			
			bits = size(transfer(r,mold),dim=1,kind=IK)*8_IK
			return
			end function bit_size_chr

			!==========================================

			function byte_size_r16(r) result(bits)
			implicit none
			real(RK16), intent(in) :: r
			integer(IK) :: bits
			
			bits = bit_size_r16(r)/8_IK
			return
			end function byte_size_r16

			!==========================================

			function byte_size_r8(r) result(bits)
			implicit none
			real(RK8), intent(in) :: r
			integer(IK) :: bits
			
			bits = bit_size_r8(r)/8_IK
			return
			end function byte_size_r8

			!==========================================

			function byte_size_r4(r) result(bits)
			implicit none
			real(RK4), intent(in) :: r
			integer(IK) :: bits
			
			bits = bit_size_r4(r)/8_IK
			return
			end function byte_size_r4

			!==========================================

			function byte_size_chr(r) result(bits)
			implicit none
			character(*), intent(in) :: r
			integer(IK) :: bits
			
			bits = bit_size_chr(r)/8_IK
			return
			end function byte_size_chr

			!==========================================

			pure subroutine encode_bits(bits,npadd,code)
			implicit none
			integer(IK1), intent(in) :: bits(:)
			integer(IK4), intent(in) :: npadd
			character(len=*), intent(out) :: code
		
			integer(IK1) :: sixb(1:4)
			integer(IK ) :: c
			integer(IK ) :: e
			integer(IK ) :: Nb
		
			Nb = size(bits,dim=1,kind=IK )
			c = 1_IK
			do e=1_IK, Nb, 3_IK
				sixb(:) = 0_IK1
				call mvbits(bits(e),2,6,sixb(1),0)
				call mvbits(bits(e),0,2,sixb(2),4)
				if ( e+1_IK .le. Nb ) then
					call mvbits(bits(e+1),4,4,sixb(2),0)
					call mvbits(bits(e+1),0,4,sixb(3),2)
				end if
				if ( e+2_IK .le. Nb ) then
					call mvbits(bits(e+2),6,2,sixb(3),0)
					call mvbits(bits(e+2),0,6,sixb(4),0)
				end if
				sixb(:) = sixb(:) + 1_IK1
				code(c  :c  ) = b64Table(sixb(1):sixb(1))
				code(c+1:c+1) = b64Table(sixb(2):sixb(2))
				code(c+2:c+2) = b64Table(sixb(3):sixb(3))
				code(c+3:c+3) = b64Table(sixb(4):sixb(4))
				c = c+4_IK
			end do
		
			if ( npadd.gt.0 ) &
				code(len(code)-npadd+1:) = repeat('=',npadd)

			end subroutine encode_bits
		
			!==========================================

			subroutine decode_bits(code,bits)
			use allFun
			implicit none
			character(len=*), intent(in) :: code
			integer(IK1), intent(out) :: bits(:)
		
			integer(IK1) :: sixb(1:4)
			integer(IK ) :: c
			integer(IK ) :: e
			integer(IK ) :: Nb
			integer(IK ) :: i
			
			Nb = size(bits,dim=1,kind=IK )
			e = 1_IK
			do c=1_IK, len(code), 4_IK
				sixb(:) = 0_IK1
				sixb(1) = index(b64Table,code(c  :c  )) - 1
				sixb(2) = index(b64Table,code(c+1:c+1)) - 1
				sixb(3) = index(b64Table,code(c+2:c+2)) - 1
				sixb(4) = index(b64Table,code(c+3:c+3)) - 1
				
				call mvbits(sixb(1),0,6,bits(e),2)
				call mvbits(sixb(2),4,2,bits(e),0)
				
				if ( e+1_IK .le. Nb ) then
					call mvbits(sixb(2),0,4,bits(e+1),4)
					call mvbits(sixb(3),2,4,bits(e+1),0)
				end if
				
				if ( e+2_IK .le. Nb ) then
					call mvbits(sixb(3),0,2,bits(e+2),6)
					call mvbits(sixb(4),0,6,bits(e+2),0)
				end if
				e = e+3_IK
			end do
			
			end subroutine decode_bits

			!==========================================

		end module params

!**************************************************

		module allFun
		use params, only: stdL,eol
		implicit none
		
		interface STR
			module procedure :: ITSTR, RTSTR, DTSTR, NDTSTR
		end interface STR

		contains
		
			!==========================================

			pure function ITSTR(iVal) result(str)
			implicit none
			integer, intent(in) :: iVal
			integer :: n,ist,j,k,slen
			character(len=256) str
			
			str = ''
			if ( iVal.lt.0 ) then
				str(1:1) = '-'
				ist = 2
			else
				ist = 1
			end if
			
			slen = 2 - max(0,sign(1,iVal)) + &
				maxval( (/ (min(abs(iVal)/10**n,1)*n,n=1, 9) /) )
			
			n = iVal
			do j=slen, ist, -1
				k = modulo(abs(n),10) + 1
				str(j:j) = '0123456789'(k:k)
				n = n/10
			end do
			
			return
			end function ITSTR
			
			!==========================================
			
			pure function RTSTR(rVal) result(str)
			implicit none
			integer, parameter :: l=8
			real*4, intent(in) :: rVal
			character(len=l) :: str
			
			str = NDTSTR(dble(rVal),l)
			
			return
			end function RTSTR
			
			!==========================================
			
			pure function DTSTR(dVal) result(str)
			implicit none
			integer, parameter :: l=8
			double precision, intent(in) :: dVal
			character(len=l) :: str
			
			str = NDTSTR(dVal,l)
			
			return
			end function DTSTR
			
			!==========================================
			
			pure function NDTSTR(dVal,l) result(str)
			implicit none
			integer, intent(in) :: l
			double precision, intent(in) :: dVal
			character(len=l) :: str
			
			integer :: i,j,k,ipos,cnt,ex,abex,nex
			double precision :: absd
			
			! check NaN !
			if ( dVal.ne.dVal ) then
				if ( l.ge.3) then
					str = ''
					str(l-2:l) = 'NaN'
				else
					str = 'NaN'(1:l)
				end if
				return
			end if
			
			! check infinity !
			absd = abs(dVal)
			if ( absd.gt.huge(absd) ) then
				if ( l.ge.8 ) then
					str = ''
					str(l-7:l) = 'Infinity'
				else
					str = 'Infinity'(1:l)
				end if
				return
			end if
			
			! check zero !
			if ( absd.lt.tiny(absd) ) then
				ipos = 1
				if ( l.ge.3 ) then
					str(1:3) = '0.0'
					ipos = 4
				end if
				do i=ipos, l
					str(i:i) = '0'
				end do
				return
			end if
			
			! count exp and number !
			ex = 0; nex = 0
			if ( absd.ne.0d0 ) ex = floor(log10(absd)) ! exponent !
			abex = abs(ex)
			
			! no of digits in exponent !
			if ( ex.ne.0 ) nex = floor(log10(real(abex,8))) + 1
			
			cnt = nex+1  ! atleast 1 number before exponent !
			if ( dVal.lt.0d0 ) cnt = cnt+1 ! for - sign of number !
			if ( ex.lt.0 ) cnt = cnt+1 ! for - sign of exponent !
			if ( ex.ne.0 ) cnt = cnt+1 ! for letter E !
			
			if ( l.lt.cnt ) then   ! check total sum !
				do i=1, l
					str(i:i) = '*'
				end do
				return
			end if
			
			if ( ex.ne.0 ) then
				do ipos=l, l-nex+1, -1
					k = modulo(abex,10) + 1
					str(ipos:ipos) = '0123456789'(k:k)
					abex = abex/10
				end do
				if ( ex.lt.0 ) then
					str(ipos:ipos) = '-'
					ipos = ipos-1
				end if
				str(ipos:ipos) = 'E'
				ipos = ipos-1
			else
				ipos = l
			end if
			
			if ( l.gt.cnt ) then
				absd = absd*(1D1**(l-cnt-1-ex))
				do i=ipos, ipos-(l-cnt)+2, -1
					k = modulo(absd,1D1) + 1
					str(i:i) = '0123456789'(k:k)
					absd = absd/1D1
				end do
				ipos = i
				str(ipos:ipos) = '.'
				ipos = ipos-1
				k = modulo(absd,1D1) + 1
				str(ipos:ipos) = '0123456789'(k:k)
			else ! l.eq.cnt
				absd = absd*(1D1**(l-cnt-ex))
				k = floor(modulo(absd,1D1)) + 1
				str(ipos:ipos) = '0123456789'(k:k)
			end if
			if ( dVal.lt.0.0d0 ) str(1:1) = '-'

			return
			end function NDTSTR

			!==========================================
			
		end module allFun

!*******************************************************

		program base64
		use params
		use allFun
		implicit none
		character(:), allocatable :: code,s1,s2
		integer(IK4) :: i,n,m
		real(RK4) :: r,d
		integer(IK4), dimension(:), allocatable :: i4arr1,i4arr2
		integer(IK8), dimension(:), allocatable :: i8arr1,i8arr2
		real(RK4), dimension(:), allocatable :: r4arr1, r4arr2
		interface
			subroutine b64_encode_I4(n,code)
			use params
			use allFun
			implicit none
			integer(IK4), intent(in) :: n
			character(len=:), allocatable, intent(out) :: code
			end subroutine b64_encode_I4
			
			subroutine b64_decode_I4(code,n)
			use params
			use allFun
			implicit none
			character(len=*), intent(in) :: code
			integer(IK4), intent(out) :: n
			end subroutine b64_decode_I4

			subroutine b64_encode_I4_a(n,code)
			use params
			use allFun
			implicit none
			integer(IK4), intent(in) :: n(:)
			character(len=:), allocatable, intent(out) :: code
			end subroutine b64_encode_I4_a
			
			subroutine b64_decode_I4_a(code,n)
			use params
			use allFun
			implicit none
			character(len=*), intent(in) :: code
			integer(IK4), intent(out) :: n(:)
			end subroutine b64_decode_I4_a

			subroutine b64_encode_R4(n,code)
			use params
			use allFun
			implicit none
			real(RK4), intent(in) :: n
			character(len=:), allocatable, intent(out) :: code
			end subroutine b64_encode_R4
			
			subroutine b64_decode_R4(code,n)
			use params
			use allFun
			implicit none
			character(len=*), intent(in) :: code
			real(RK4), intent(out) :: n
			end subroutine b64_decode_R4

			subroutine b64_encode_R4_a(n,code)
			use params
			use allFun
			implicit none
			real(RK4), intent(in) :: n(:)
			character(len=:), allocatable, intent(out) :: code
			end subroutine b64_encode_R4_a
			
			subroutine b64_decode_R4_a(code,n)
			use params
			use allFun
			implicit none
			character(len=*), intent(in) :: code
			real(RK4), intent(out) :: n(:)
			end subroutine b64_decode_R4_a

			subroutine b64_encode_str(n,code)
			use params
			use allFun
			implicit none
			character(len=*), intent(in) :: n
			character(len=:), allocatable, intent(out) :: code
			end subroutine b64_encode_str
			
			subroutine b64_decode_str(code,n)
			use params
			use allFun
			implicit none
			character(len=*), intent(in) :: code
			character(len=*), intent(out) :: n
			end subroutine b64_decode_str

			subroutine b64_encode_I8_R4_a(i,r,code)
			use params
			use allFun
			implicit none
			integer(IK8), intent(in) :: i(:)
			real(RK4), intent(in) :: r(:)
			character(len=:), allocatable, intent(out) :: code
			end subroutine b64_encode_I8_R4_a
			
			subroutine b64_decode_I8_R4_a(code,ni,nr,i,r)
			use params
			use allFun
			implicit none
			character(len=*), intent(in) :: code
			integer(IK4), intent(in) :: ni,nr
			integer(IK8), intent(out) :: i(:)
			real(RK4), intent(out) :: r(:)
			end subroutine b64_decode_I8_R4_a

		end interface
		
		BYTRK4 = byte_size_r4(MAXRK4)
		BYTRK8 = byte_size_r8(MAXRK8)
		
		allocate(i4arr1(22)); i4arr1 = 0
		allocate(r4arr1(27)); r4arr1 = 0.0
		
		! Int*4 encoding !
		write(stdout,'(A)')
		write(stdout,'(4X,A)') "Int32 encoding:"
		i = 1000000
		call b64_encode_I4(i,code)
		write(stdout,'(8X,A)') &
			trim(STR(i))//" (base64 encoding) ---> "//code
			
		call b64_decode_I4(code,n)
		write(stdout,'(8X,A)') &
			code//" (base64 decoding) ---> "//trim(STR(n))
		write(stdout,'(A)')
		
		! Real*4 encoding !
		write(stdout,'(4X,A)') "Float32 encoding.."
		r = 126354684231.78594
		call b64_encode_R4(r,code)
		write(stdout,'(8X,A)') &
			STR(r)//" (base64 encoding) ---> "//code
			
		call b64_decode_R4(code,d)
		write(stdout,'(8X,A)') &
			code//" (base64 decoding) ---> "// STR(d)
		write(stdout,'(A)')
		
		! Char* encoding !
		write(stdout,'(4X,A)') "String encoding.."
		s1 = "any carnal pleasur"
		call b64_encode_str(s1,code)
		write(stdout,'(8X,A)') &
			s1//" (base64 encoding) ---> "//code
		
		s2 = repeat(" ",len(code)*3_IK/4_IK)
		call b64_decode_str(code,s2)
		write(stdout,'(8X,A)') &
			code//" (base64 decoding) ---> "//s2
		
		code = "c3VyZS4="
		s2 = repeat(" ",len(code)*3/4)
		call b64_decode_str(code,s2)
		write(stdout,'(8X,A)') &
			code//" (base64 decoding) ---> "//s2
		write(stdout,'(A)')
		
		! Int*4 arrays !
		n = size(i4arr1,dim=1,kind=IK4)
		write(stdout,'(4X,A)') "Int32[] encoding: "
		write(stdout,'(8X,A)',advance="no") "Raw Data:       "
		do i=1, n
			i4arr1(i) = i**2
			write(stdout,'(A)',advance="no") " "//trim(STR(i4arr1(i)))
		end do
		write(stdout,'(A)')
		call b64_encode_I4_a(i4arr1,code)
		write(stdout,'(8X,A)') "Base64 Encoding: "//code
		
		allocate(i4arr2(n)); i4arr2=0
		call b64_decode_I4_a(code,i4arr2)
		write(stdout,'(8X,A)',advance="no") "Decoded Data:   "
		do i=1, n
			write(stdout,'(A)',advance="no") " "//trim(STR(i4arr2(i)))
		end do
		write(stdout,'(A)')
		write(stdout,'(A)')
		
		! Real*4 encoding !
		n = size(r4arr1,dim=1,kind=IK4)
		write(stdout,'(4X,A)') "Float32[] encoding: "
		write(stdout,'(8X,A)',advance="no") "Raw Data:       "
		do i=1, n
			r4arr1(i) = -14.0 + real(i,kind=RK4)
			write(stdout,'(A)',advance="no") " "//STR(r4arr1(i))
		end do
		write(stdout,'(A)')
		call b64_encode_R4_a(r4arr1,code)
		write(stdout,'(8X,A)') "Base64 Encoding: "//code
		
		allocate(r4arr2(n)); r4arr2=0.0
		call b64_decode_R4_a(code,r4arr2)
		write(stdout,'(8X,A)',advance="no") "Decoded Data:   "
		do i=1, n
			write(stdout,'(A)',advance="no") " "//STR(r4arr2(i))
		end do
		write(stdout,'(A)')
		write(stdout,'(A)')

		! Int8[]Float32[] compound encoding !
		n = size(r4arr1,dim=1,kind=IK4)
		allocate(i8arr1(1)); i8arr1(1) = int(n*BYTRK4,kind=IK8)
		m = size(i8arr1,dim=1,kind=IK4)
		write(stdout,'(4X,A)') "Int8[]Float32[] encoding: "
		write(stdout,'(8X,A)',advance="no") "Raw Data:       "
		do i=1, m
			write(stdout,'(A)',advance="no") " "// &
			trim(STR(int(i8arr1(i),kind=IK4)))
		end do
		do i=1, n
			r4arr1(i) = -14.0 + real(i,kind=RK4)
			write(stdout,'(A)',advance="no") " "//STR(r4arr1(i))
		end do
		write(stdout,'(A)')
		call b64_encode_I8_R4_a(i8arr1,r4arr1,code)
		write(stdout,'(8X,A)') "Base64 Encoding: "//code
		
		allocate(i8arr2(m)); i8arr2=0
		if ( allocated(r4arr2) ) deallocate(r4arr2)
		allocate(r4arr2(n)); r4arr2=0.0
		call b64_decode_I8_R4_a(code,m,n,i8arr2,r4arr2)
		write(stdout,'(8X,A)',advance="no") "Decoded Data:   "
		do i=1, m
			write(stdout,'(A)',advance="no") " "//&
				trim(STR(int(i8arr2(i),kind=IK4)))
		end do
		do i=1, n
			write(stdout,'(A)',advance="no") " "//STR(r4arr2(i))
		end do
		write(stdout,'(A)')
		write(stdout,'(A)')

		end program base64
		
!*******************************************************

		subroutine b64_encode_I4(n,code)
		use params
		use allFun
		implicit none
		integer(IK4), intent(in) :: n
		character(len=:), allocatable, intent(out) :: code
		integer(IK1), allocatable :: nIK1(:)
		integer(IK4) :: npadd,i
		
		allocate(nIK1(1: ((BYTIK4+2)/3)*3)); nIK1 = 0_IK1
		code = repeat(" ", ((BYTIK4+2)/3)*4)
		nIK1 = transfer(n,nIK1)
		npadd = mod(BYTIK4,3_IK4)
		if ( npadd>0_IK4 ) npadd = 3_IK4 - npadd
		
		call encode_bits(nIK1,npadd,code)
		
		end subroutine b64_encode_I4
		
!*******************************************************

		subroutine b64_decode_I4(code,n)
		use params
		use allFun
		implicit none
		character(len=*), intent(in) :: code
		integer(IK4), intent(out) :: n
		integer(IK1), allocatable :: nIK1(:)
		
		allocate(nIK1(1:BYTIK4)); nIK1 = 0_IK1
		
		call decode_bits(code,nIK1)
		n = transfer(nIK1,n)
		
		end subroutine b64_decode_I4
		
!*******************************************************

		subroutine b64_encode_I4_a(n,code)
		use params
		use allFun
		implicit none
		integer(IK4), intent(in) :: n(:)
		character(len=:), allocatable, intent(out) :: code
		integer(IK1), allocatable :: nIK1(:)
		integer(IK4) :: nsize,npadd,i
		
		nsize = size(n,dim=1,kind=IK4)*BYTIK4
		allocate(nIK1(1: ((nsize+2)/3)*3)); nIK1 = 0_IK1
		code = repeat(" ", ((nsize+2)/3)*4)
		nIK1 = transfer(n,nIK1)
		npadd = mod(nsize,3_IK4)
		if ( npadd>0_IK4 ) npadd = 3_IK4 - npadd
		
		call encode_bits(nIK1,npadd,code)
		
		end subroutine b64_encode_I4_a
		
!*******************************************************

		subroutine b64_decode_I4_a(code,n)
		use params
		use allFun
		implicit none
		character(len=*), intent(in) :: code
		integer(IK4), intent(out) :: n(:)
		integer(IK1), allocatable :: nIK1(:)
		integer(IK4) :: nsize
		
		nsize = len(code)*3_IK4/4_IK4
		allocate(nIK1(1:nsize)); nIK1 = 0_IK1
		
		call decode_bits(code,nIK1)
		n = transfer(nIK1,n)
		
		end subroutine b64_decode_I4_a
		
!*******************************************************

		subroutine b64_encode_R4(n,code)
		use params
		use allFun
		implicit none
		real(RK4), intent(in) :: n
		character(len=:), allocatable, intent(out) :: code
		integer(IK1), allocatable :: nIK1(:)
		integer(IK ) :: npadd,i
		
		allocate(nIK1(1: ((BYTRK4+2)/3)*3)); nIK1 = 0_IK1
		code = repeat(" ", ((BYTRK4+2)/3)*4)
		nIK1 = transfer(n,nIK1)
		npadd = mod(BYTRK4,3_IK)
		if ( npadd>0_IK ) npadd = 3_IK - npadd
		
		call encode_bits(nIK1,npadd,code)
		
		end subroutine b64_encode_R4
		
!*******************************************************

		subroutine b64_decode_R4(code,n)
		use params
		use allFun
		implicit none
		character(len=*), intent(in) :: code
		real(RK4), intent(out) :: n
		integer(IK1), allocatable :: nIK1(:)
		
		allocate(nIK1(1:BYTRK4)); nIK1 = 0_IK1
		
		call decode_bits(code,nIK1)
		n = transfer(nIK1,n)
		
		end subroutine b64_decode_R4
		
!*******************************************************

		subroutine b64_encode_R4_a(n,code)
		use params
		use allFun
		implicit none
		real(RK4), intent(in) :: n(:)
		character(len=:), allocatable, intent(out) :: code
		integer(IK1), allocatable :: nIK1(:)
		integer(IK4) :: nsize,npadd,i
		
		nsize = size(n,dim=1,kind=IK4)*BYTRK4
		allocate(nIK1(1: ((nsize+2)/3)*3)); nIK1 = 0_IK1
		code = repeat(" ", ((nsize+2)/3)*4)
		nIK1 = transfer(n,nIK1)
		npadd = mod(nsize,3_IK4)
		if ( npadd>0_IK4 ) npadd = 3_IK4 - npadd
		
		call encode_bits(nIK1,npadd,code)
		
		end subroutine b64_encode_R4_a
		
!*******************************************************

		subroutine b64_decode_R4_a(code,n)
		use params
		use allFun
		implicit none
		character(len=*), intent(in) :: code
		real(RK4), intent(out) :: n(:)
		integer(IK1), allocatable :: nIK1(:)
		integer(IK4) :: nsize
		
		nsize = len(code)*3_IK4/4_IK4
		allocate(nIK1(1:nsize)); nIK1 = 0_IK1
		
		call decode_bits(code,nIK1)
		n = transfer(nIK1,n)
		
		end subroutine b64_decode_R4_a
		
!*******************************************************

		subroutine b64_encode_str(n,code)
		use params
		use allFun
		implicit none
		character(len=*), intent(in) :: n
		character(len=:), allocatable, intent(out) :: code
		integer(IK1), allocatable :: nIK1(:)
		integer(IK4) :: npadd,i
		
		BYTCHR = byte_size(n)
		allocate(nIK1(1: ((BYTCHR+2)/3)*3)); nIK1 = 0_IK1
		code = repeat(" ", ((BYTCHR+2)/3)*4)
		nIK1 = transfer(n,nIK1)
		npadd = mod(BYTCHR,3_IK)
		if ( npadd>0_IK ) npadd = 3_IK - npadd
		
		call encode_bits(nIK1,npadd,code)
		
		end subroutine b64_encode_str
		
!*******************************************************

		subroutine b64_decode_str(code,n)
		use params
		use allFun
		implicit none
		character(len=*), intent(in) :: code
		character(len=*), intent(out) :: n
		integer(IK1), allocatable :: nIK1(:)
		
		BYTCHR = byte_size(code)*3_IK/4_IK
		allocate(nIK1(1:BYTCHR)); nIK1 = 0_IK1
		
		call decode_bits(code,nIK1)
		n = transfer(nIK1,n)
		
		end subroutine b64_decode_str
		
!*******************************************************

		subroutine b64_encode_I8_R4_a(n,d,code)
		use params
		use allFun
		implicit none
		integer(IK8), intent(in) :: n(:)
		real(RK4), intent(in) :: d(:)
		character(len=:), allocatable, intent(out) :: code
		integer(IK1), allocatable :: p1(:),p2(:),nIK1(:)
		integer(IK4) :: np1,np2,nsize,npadd,i
		
		np1 = size(n,dim=1,kind=IK4)*BYTIK8
		allocate(p1(1:np1)); p1 = 0_IK1
		p1 = transfer(n,p1)
		
		np2 = size(d,dim=1,kind=IK4)*BYTRK4
		allocate(p2(1:np2)); p2 = 0_IK1
		p2 = transfer(d,p2)
		
		nsize = np1+np2
		allocate(nIK1(1: ((nsize+2)/3)*3)); nIK1 = [p1, p2];
		code = repeat(" ", ((nsize+2)/3)*4)
		npadd = mod(nsize,3_IK4)
		if ( npadd>0_IK4 ) npadd = 3_IK4 - npadd
		
		call encode_bits(nIK1,npadd,code)
		
		end subroutine b64_encode_I8_R4_a
		
!*******************************************************

		subroutine b64_decode_I8_R4_a(code,m,n,iarr,rarr)
		use params
		use allFun
		implicit none
		character(len=*), intent(in) :: code
		integer(IK4), intent(in) :: m,n
		integer(IK8), intent(out) :: iarr(:)
		real(IK4), intent(out) :: rarr(:)
		integer(IK1), allocatable :: nIK1(:),p1(:),p2(:)
		integer(IK4) :: nsize,np1,np2
		
		nsize = len(code)*3_IK4/4_IK4
		allocate(nIK1(1:nsize)); nIK1 = 0_IK1
		
		call decode_bits(code,nIK1)
		
		np1 = m*BYTIK8
		np2 = n*BYTRK4
		
		allocate(p1(np1),p2(np2))
		p1(:)=0_IK1; p2(:)=0_IK1
		
		p1(1:np1) = nIK1(1:np1)
		p2(1:np2) = nIK1(np1+1:nsize)
		
		iarr = transfer(p1,iarr)
		rarr = transfer(p2,rarr)
		
		end subroutine b64_decode_I8_R4_a
		
!*******************************************************










	program ge
	implicit none
	integer :: fid, i, j, m, n, iok
	integer, allocatable :: piv(:)
	real(kind=8), allocatable :: A(:,:), Ai(:,:), Al(:,:), work(:)
	character(len=256) :: fname
	
	
	interface
	   function matinvs(A,n)
	   integer, intent(in) :: n
	   real(kind=8), intent(in) :: A(n,n)
	   real(kind=8) :: matinvs(n,n)
	   end function matinvs
	end interface
	
	fid = 10
	i = iargc()
	if (i .ne. 1) then
	   write(*,*) "Error: at least one input argument needed"
	   STOP
	end if
	call getarg(1,fname)
	
	open(fid,file=trim(fname))
	read(fid,*) m, n
	if (m .ne. n) then
	   write(*,*) "Error: inverse is possible only for square matrices"
	   stop
	end if
	
	allocate(A(m,n), Ai(m,n), Al(m,n))
	do i=1, m
	   read(fid,*) (A(i,j), j=1, n)
	end do
	close(fid)
	call dispMat(A, m, n, "Loaded matrix:")
	
	Ai = matinvs(A, m)
	call dispMat(Ai, m, n, "Inverse (Gauss Elimination): ")
	
	allocate(piv(n))
	call DGETRF(m, n, A, m, piv, iok)
	if (iok .ne. 0) then
	   write(*,'(4X,A)') "Error: LAPACK computing pivots"
	   STOP
	end if
	
	Al = A
	allocate(work(2*n))
	call DGETRI(n, Al, m, piv, work, 2*n, iok)
	if (iok .ne. 0) then
	   write(*,'(4X,A)') "Error: LAPACK computing inverse"
	   STOP
	end if
	call dispMat(Al, m, n, "Inverse (Lapack): ")
	
	end program ge
	
!-----------------------------------------------------------------------

	function matinvs(A, n)
	implicit none
	integer, intent(in) :: n
	real(kind=8), intent(in) :: A(n,n)
	real(kind=8) :: matinvs(n,n)
	
	integer :: m, ipv, i, j
	real(kind=8) :: pivot, tmp, B(n,2*n), detMat
	
	if (abs(detMat(A, n)) .lt. 1d3*epsilon(pivot)) then
	   write(*,'(4X,A)') "WARNING: possible singular matrix"
	end if
	
	B = 0.0d0
	do i=1, n
	   do j=1, n
	      B(i,j) = A(i,j)
	   end do
	   B(i,n+i) = 1.0d0
	end do
	
	do m=1, n
	   ipv = m
	   pivot = abs(B(m,m))
	   do i=m+1, n
	      if (abs(B(i,m)) .gt. pivot) then
	         ipv = i
	         pivot = abs(B(i,m))
	      end if
	   end do
	   if (abs(pivot) .lt. epsilon(pivot)) then
	      write(*,*) "Error: zero pivot. possible singular matrix"
	      stop
	   end if
	   
	   if (ipv .gt. m) then
	      do j=m, 2*n
	         tmp = B(m,j)
	         B(m,j) = B(ipv,j)
	         B(ipv,j) = tmp
	      end do
	   end if
	   
	   do i=1, n
	      if (i .ne. m) then
	         tmp = B(i,m)/B(m,m)
             B(i,m) = 0d0
	         do j=m+1, 2*n
	            B(i,j) = B(i,j) - tmp*B(m,j)
	         end do
	      end if
	   end do
	end do
	
	do m=1, n
	   do j=n+1, 2*n
	      B(m,j) = B(m,j)/B(m,m)
	   end do
	end do
	
	do i=1, n
	   do j=1, n
	      matinvs(i,j) = B(i,j+n)
	   end do
	end do
	
	end function matinvs
	
!-----------------------------------------------------------------------

	recursive function detMat(A, n) result(D)
	implicit none
	integer, intent(in) :: n
	real(kind=8), intent(in) :: A(n,n)
	integer :: i, j, m
	real(kind=8) :: D, Am(n-1, n-1)
	
	D = 0.0d0
	if (n .eq. 2) then
	   D = A(1,1)*A(2,2) - A(1,2)*A(2,1)
	else
	   do i=1, n
	      m = 0
	      do j=1, n
	         if (i .eq. j) then
	            cycle
	         else
	            m = m + 1
	            Am(:,m) = A(2:n,j)
	         end if
	      end do
	      D = D + ( (-1D0)**real(1+i,kind=8) * A(1,i) * detMat(Am,n-1) )
	   end do
	end if
	
	return
    end function detMat
    
!-----------------------------------------------------------------------

	subroutine dispMat(A, m, n, msg)
	implicit none
	integer, intent(in) :: m, n
	real(kind=8), intent(in) :: A(m,n)
	character(len=*), intent(in) :: msg
	integer :: i,j
	character(len=256) :: stmp
	
	stmp = ""
	if (len(msg) .gt. 0) stmp(1:len(msg)) = msg
	write(*,*)
	if (len(msg) .gt. 0) write(*,'(4X,A)') trim(stmp)
	do i=1, m
	   write(*,'(8X,A)',advance='no')
	   do j=1, n
	      write(*,'(A,F15.6)',advance='no') " ", A(i,j)
	      if (m.ne.n .and. j.eq.n/2) write(*,'(4X,A)',advance='no') "|"
	   end do
	   write(*,'(A)')
	end do
	
	end subroutine dispMat
	
!-----------------------------------------------------------------------


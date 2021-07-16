      
      module variables
         integer, parameter :: nTS = 50000, stdL = 400
         real(kind=8), parameter :: dt = 1D-4
         real(kind=8), parameter :: pi = 3.1415926535897932384626D0
         type fcType
            integer :: imthd = 0
            integer :: n = 0
            real(kind=8) qi
            real(kind=8) qs
            real(kind=8) T
            real(kind=8) ti
            real(kind=8), allocatable :: i(:)
            real(kind=8), allocatable :: r(:)
         end type fcType
      end module variables
!--------------------------------------------------------------------
      program fft_interp
      use variables
      implicit none
      
      integer :: fid, i, j, np
      real(kind=8) :: q, dq
      type(fcType) :: gt
      character(len=stdL) :: fname
      
      real(kind=8), allocatable :: time(:)
      
      i = IARGC()
      if (i .eq. 0) then
         write(*,'(A)') "ERROR: Input file name not specified"
         STOP
      else if (i .gt. 1) then
         write(*,'(A)') "ERROR: Too many arguments"
         STOP
      end if
      
      call getarg(1,fname)
      
      fid = 100
      open(fid, file=trim(fname))
      read(fid,*) i, j
      gt%n = j
      allocate(gt%i(j))
      allocate(gt%r(j))
      call fft(fid, i, gt)
      close(fid)
      
      allocate(time(nTS))
      do i=1, nTS
         time(i) = real(i-1,kind=8) * dt
      end do
      
      open(fid, file='interp_out.dat')
      do i=1, nTS
         call ifft(gt, time(i), q, dq)
         write(fid,'(F9.4, 2(X,1pE15.6))') &
            time(i), q, dq
      end do
      close(fid)
      
      return
      end program fft_interp
!--------------------------------------------------------------------
      subroutine fft(fid, np, gt)
      use variables
      implicit none
      
      integer, intent(in) :: fid, np
      type (fcType), intent(inout) :: gt
      
      integer :: i, n
      real(kind=8) :: tmp, kn, ko, s
      
      real(kind=8), allocatable :: t(:), q(:)
      
      allocate (t(np), q(np))
      do i=1, np
         read (fid,*) t(i), q(i)
      end do
      
      gt%ti = t(1)
      gt%T  = t(np) - t(1)
      gt%qi = q(1)
      gt%qs = (q(np) - q(1))/gt%T
      
      do i=1, np
         t(i) = t(i) - gt%ti
         q(i) = q(i) - gt%qi - gt%qs*t(i)
      end do
      
      do n=1, gt%n
         tmp = real(n-1,8)
         gt%r(n) = 0D0
         gt%i(n) = 0D0
         do i=1, np-1
            ko = 2D0*pi*tmp*t(i)/gt%T
            kn = 2D0*pi*tmp*t(i+1)/gt%T
            s  = (q(i+1) - q(i))/(t(i+1) - t(i))
            
            if (n .EQ. 1) then
               gt%r(n) = gt%r(n) + 5D-1*(t(i+1)-t(i))*(q(i+1)+q(i))
            else
               gt%r(n) = gt%r(n) + s*(COS(kn) - COS(ko))
               gt%i(n) = gt%i(n) - s*(SIN(kn) - SIN(ko))
            end if
         end do
         
         if (n .EQ. 1) then
            gt%r(n) = gt%r(n)/gt%T
         else
            gt%r(n) = 5D-1*gt%r(n)*gt%T/(pi*pi*tmp*tmp)
            gt%i(n) = 5D-1*gt%i(n)*gt%T/(pi*pi*tmp*tmp)
         end if
      end do
      
      return
      end subroutine fft
!--------------------------------------------------------------------
!     This is to calculate flow rate and flow acceleration (IFFT)
      subroutine ifft(gt, time, Y, dY)
      use variables
      implicit none
      
      type(fcType), intent(in) :: gt
      real(kind=8), intent(in) :: time
      real(kind=8), intent(out) :: Y, dY
      
      integer i
      real(kind=8) t, tmp, K, kd
      
      t    = DMOD(time - gt%ti, gt%T)
      tmp  = 2D0*pi/gt%T
      Y    = gt%qi + t*gt%qs
      dY   = gt%qs
      do i=1, gt%n
         kd = tmp*real(i-1,8)
         K  = t*kd
         Y  =  Y +  gt%r(i)*COS(K) - gt%i(i)*SIN(K)
         dY = dY - (gt%r(i)*SIN(K) + gt%i(i)*COS(K))*kd
      end do
      
      return 
      end subroutine ifft
!-------------------------------------------------------------------- 

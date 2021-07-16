!**********************************************************************

      module variables

      integer, parameter :: stdL = 256
      type bsType
         integer :: n, p, nc, ne, ns, eNoN
         integer, allocatable :: inn(:), ien(:,:)
         real(kind=8), allocatable :: u(:)
      end type bsType

      end module variables

!**********************************************************************

      program bsIGA
      use variables
      implicit none

      type(bsType) :: bs
      integer :: i, a, e
      character(len=stdL) :: fname

      i = iargc()
      if (i .eq. 0) then
         write(*,'(2X,A)') "Error: missing input file argument"
         stop
      else if (i .gt. 1) then
         write(*,'(2X,A)') "Error: too many arguments"
         stop
      end if
      call getarg(1,fname)

      write(*,'(A)') "===================================="
      write(*,'(2X,A)') "Reading knot vector from file: "//trim(fname)
      call readKnotV(bs, fname)

      call constBSplines(bs)

      write(*,'(2X,A)') "Knot vector:"
      write(*,'(4X,A,I2)') "size: ", bs%n
      write(*,'(4X,A,I2)') "order: ", bs%p
      write(*,'(4X,A,I2)') "num control points: ", bs%nc
      write(*,'(4X,A,I2)') "span: ", bs%ne

      write(*,'(4X,A)',advance='no') "Knot vector:  {"
      do i=1, bs%n
         write(*,'(X,F4.1)',advance='no') bs%u(i)
         if (i .lt. bs%n) then
            write(*,'(A)',advance='no') ","
         end if
      end do
      write(*,'(2X,A)') "}"

      write(*,'(4X,A)',advance='no') "INN: {"
      do e=1, bs%ne
         write(*,'(X,I2)',advance='no') bs%inn(e)
         if (e .lt. bs%ne) then
            write(*,'(A)',advance='no') ","
         end if
      end do
      write(*,'(2X,A)') "}"

      write(*,'(4X,A)') "IEN: {"
      do e=1, bs%ne
         write(*,'(6X,A)',advance='no')
         do a=1, bs%eNoN
            write(*,'(X,I2)',advance='no') bs%ien(a,e)
         end do
         write(*,'(A)')
      end do
      write(*,'(4X,A)') "}"

      write(*,'(2X,A)') "Constructing shape functions.."
      do e=1, bs%ne
         call constBSShpFn(bs, e)
      end do

      return
      end program bsIGA

!**********************************************************************

      subroutine readKnotV(bs, ctmp)
      use variables
      implicit none

      type(bsType), intent(inout) :: bs
      character(len=stdL), intent(in) :: ctmp

      integer :: fid, i

      fid = 101
      open(fid,file=trim(ctmp))
      read(fid,*) ! knot vector size
      read(fid,*) bs%n

      allocate(bs%u(bs%n))
      bs%u = 0D0
      read(fid,*) ! knot vector
      read(fid,*) (bs%u(i), i=1, bs%n)

      read(fid,*) ! num sampling points
      read(fid,*) bs%ns
      close(fid)

      return
      end subroutine readKnotV

!**********************************************************************

      subroutine constBSplines(bs)
      use variables
      implicit none

      type(bsType), intent(inout) :: bs

      integer :: i, a, e
      logical :: flag

!     check knot vector and count knot spans (elements)
      bs%p  = 0
      bs%ne = 0
      flag  = .true.
      do i=1, bs%n-1
         if (bs%u(i) .gt. bs%u(i+1)) then
            write(*,'(2X,A)') &
               "Error: knot vector must be non-decreasing"
            stop
         end if
         if (flag) then
            if (bs%u(i+1) .eq. bs%u(i)) then
               bs%p = bs%p + 1
            else
               if (i .eq. 1) then
                  write(*,'(2X,A)') "Error: violating p > 0"
                  stop
               end if
               flag = .false.
            end if
         end if
         if (bs%u(i+1) .ne. bs%u(i)) bs%ne = bs%ne + 1
         if (i .ge. bs%n-bs%p) then
            if (bs%u(i) .ne. bs%u(i+1)) then
               write(*,'(2X,A)') "Error: knot vector must be open. "// &
                  "Expecting p+1 similar knots at the end"
               stop
            end if
         end if
      end do
      bs%nc = bs%n - bs%p - 1
      bs%eNoN = bs%p + 1

      allocate(bs%inn(bs%ne), bs%ien(bs%eNoN, bs%ne))
      e = 0
      do i=bs%p+1, bs%n-bs%p-1
         if (bs%u(i) .ne. bs%u(i+1)) then
            e = e + 1
            bs%inn(e) = i
            do a=1, bs%eNoN
               bs%ien(a,e) = i + a - bs%p - 1
            end do
         end if
      end do

      return
      end subroutine constBSplines

!**********************************************************************

      subroutine constBSShpFn(bs, e)
      use variables
      implicit none

      type(bsType), intent(in) :: bs
      integer, intent(in) :: e

      integer :: p, i, s, ni, fid
      integer :: j, k, ik, pk, s1, s2, j1, j2
      real(kind=8) :: ui, saved, tmp, c, d, du
      character(len=stdL) :: fname

      real(kind=8), allocatable :: l(:), r(:), NdU(:,:,:), ders(:,:,:),&
         a(:,:)

      allocate(l(bs%p), r(bs%p), NdU(0:bs%p,0:bs%p,bs%ns), &
         ders(0:bs%p,0:bs%p,bs%ns), a(0:bs%p,2))

      ni = bs%inn(e)
      du = (bs%u(ni+1) - bs%u(ni))/2D0
      c  = 2D0*du / real(bs%ns-1, kind=8)
      do s=1, bs%ns
         ui = bs%u(ni) + c*real(s-1, kind=8)

         NdU(0,0,s) = 1D0
         do p=1, bs%p
            l(p)  = ui - bs%u(ni+1-p)
            r(p)  = bs%u(ni+p) - ui
            saved = 0D0
            do i=0, p-1
               NdU(i,p,s) = (r(i+1) + l(p-i))
               tmp    = NdU(p-1,i,s)/NdU(i,p,s)
               NdU(p,i,s) = saved + tmp*r(i+1)
               saved  = tmp*l(p-i)
            end do
            NdU(p,p,s)  = saved
         end do
         ders(:,0,s) = NdU(bs%p,:,s)

         do i=0, bs%p
            s1 = 1; s2 = 2
            a(0,1) = 1D0
            do k=1, bs%p
               d = 0D0
               ik = i - k
               pk = bs%p - k
               if (i .ge. k) then
                  a(0,s2) = a(0,s1)/NdU(ik,pk+1,s)
                  d = a(0,s2) * NdU(pk,ik,s)
               end if
               if (ik .ge. -1) then
                  j1 = 1
               else
                  j1 = -ik
               end if
               if (i-1 .le. pk) then
                  j2 = k-1
               else
                  j2 = bs%p - i
               end if
               do j=j1, j2
                  a(j,s2) = (a(j,s1)-a(j-1,s1))/NdU(ik+j,pk+1,s)
                  d = d + a(j,s2)*NdU(pk,ik+j,s)
               end do
               if (i .le. pk) then
                  a(k,s2) = -a(k-1,s1)/NdU(i,pk+1,s)
                  d = d + a(k,s2) * NdU(pk,i,s)
               end if
               ders(i,k,s) = d
               j  = s1
               s1 = s2
               s2 = j
            end do
         end do

         i = bs%p
         do k=1, bs%p
            do j=0, bs%p
               ders(j,k,s) = ders(j,k,s)*real(i,kind=8)*(du**k)
            end do
            i = i * (bs%p-k)
         end do
      end do

!======================================================================
!     Write to file
      fid = 101
      write(fname,'(A,I2.2,A)') "bsN_e", e, ".dat"
      open(fid,file=trim(fname))
      write(fid,'(A)',advance='no') "Variables=x"
      do p=1, bs%p+1
         write(fid,'(A,I2.2)',advance='no') ", N", p
      end do
      write(fid,'(A)')

      do s=1, bs%ns
         ui = bs%u(ni) + c*real(s-1, kind=8)
         write(fid,'(F12.6)',advance='no') ui
         do p=0, bs%p
            write(101,'(2X,F12.6)',advance='no') ders(p,0,s)
         end do
         write(fid,'(A)')
      end do
      close(fid)

      write(fname,'(A,I2.2,A)') "bsNx_e", e, ".dat"
      open(fid,file=trim(fname))
      write(fid,'(A)',advance='no') "Variables=x"
      do p=1, bs%p+1
         write(fid,'(A,I2.2)',advance='no') ", N", p
      end do
      write(fid,'(A)')

      do s=1, bs%ns
         ui = bs%u(ni) + c*real(s-1, kind=8)
         write(fid,'(F12.6)',advance='no') ui
         do p=0, bs%p
            write(101,'(2X,F12.6)',advance='no') ders(p,1,s)
         end do
         write(fid,'(A)')
      end do
      close(fid)

      write(fname,'(A,I2.2,A)') "bsNxx_e", e, ".dat"
      open(fid,file=trim(fname))
      write(fid,'(A)',advance='no') "Variables=x"
      do p=1, bs%p+1
         write(fid,'(A,I2.2)',advance='no') ", N", p
      end do
      write(fid,'(A)')

      do s=1, bs%ns
         ui = bs%u(ni) + c*real(s-1, kind=8)
         write(fid,'(F12.6)',advance='no') ui
         do p=0, bs%p
            write(101,'(2X,F12.6)',advance='no') ders(p,2,s)
         end do
         write(fid,'(A)')
      end do
      close(fid)
!======================================================================

      deallocate(l, r, NdU, ders)

      return
      end subroutine constBSShpFn

!**********************************************************************

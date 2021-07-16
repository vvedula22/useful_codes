!**************************************************

      module variables
      use stdParams
      use genUtils
      use vtkXMLMod
      implicit none

      type mshType
         integer :: nNo
         integer :: nEl
         integer :: eNoN
         integer :: vtkType
         integer, allocatable :: IEN(:,:)
         real(kind=8), allocatable :: x(:,:)
      end type mshType

      type, extends(mshType) :: faceType
         integer, allocatable :: gE(:)
         integer, allocatable :: gN(:)
         real(kind=8), allocatable :: xC(:,:)
         real(kind=8), allocatable :: nV(:,:)
      end type faceType

      type imageType
         integer :: iLims(2,3)
         integer :: pLims(2,3)
         real(kind=8) :: x0(3)
         real(kind=8) :: dx(3)
         integer, dimension(:), allocatable :: ival
      end type imageType

      interface destroy
         module procedure destroyMesh, destroyFace, destroyImage
      end interface destroy

      integer, parameter :: ICOORD = 1, JCOORD = 2
      integer, parameter :: nsd  = 3
      integer, parameter :: ngl  = 1
      integer, parameter :: fidP = 65

      logical :: ilog
      integer :: nxG, nyG, nz
      integer :: nx, ny
      integer :: imin, imax, jmin, jmax, kmin, kmax
      integer :: ill, iul, jll, jul, kll, kul
      integer :: istat
      real(kind=8) :: xout, yout, zout
      character(len=strL) :: fLog

      integer, allocatable :: iblank(:,:,:)
      real(kind=8), allocatable :: xc(:)
      real(kind=8), allocatable :: yc(:)
      real(kind=8), allocatable :: zc(:)

      type(faceType)  :: fa
      type(imageType) :: im

      contains
!-------------------------------------------------
         pure function cross(u, v) result(n)
         implicit none
         real(kind=8), intent(in) :: u(:), v(:)
         real(kind=8) n(nsd)

         n(1) = u(2)*v(3) - u(3)*v(2)
         n(2) = u(3)*v(1) - u(1)*v(3)
         n(3) = u(1)*v(2) - u(2)*v(1)

         return
         end function cross
!-------------------------------------------------
         subroutine destroyMesh(lM)
         implicit none
         type(mshType), intent(inout) :: lM

         if (allocated(lM%IEN))  deallocate(lM%IEN)
         if (allocated(lM%x))    deallocate(lM%x)

         return
         end subroutine destroyMesh
!-------------------------------------------------
         subroutine destroyFace(lFa)
         implicit none
         type(faceType), intent(inout) :: lFa

         if (allocated(lFa%IEN))  deallocate(lFa%IEN)
         if (allocated(lFa%x))    deallocate(lFa%x)
         if (allocated(lFa%gN))   deallocate(lFa%gN)
         if (allocated(lFa%gE))   deallocate(lFa%gE)
         if (allocated(lFa%xC))   deallocate(lFa%xC)
         if (allocated(lFa%nV))   deallocate(lFa%nV)

         return
         end subroutine destroyFace
!-------------------------------------------------
         subroutine destroyImage(lIm)
         implicit none
         type(imageType), intent(inout) :: lIm

         if (allocated(lIm%ival))  deallocate(lIm%ival)

         return
         end subroutine destroyImage
!-------------------------------------------------
      end module variables

!**************************************************

      module mpiVariables
      use mpi
      logical, parameter :: periods(2)=.false., reorder=.true.
      integer, parameter :: Xdir=0, Ydir=1
      integer, parameter :: master = 0
      integer, parameter :: mplog  = MPI_LOGICAL
      integer, parameter :: mpint  = MPI_INTEGER
      integer, parameter :: mpreal = MPI_DOUBLE_PRECISION
      integer, parameter :: mpchar = MPI_CHARACTER
      logical :: isMas, isSeq
      integer :: nProcs, np(2), myRank, myCoords(2), myCommCart, myCartRank
      integer :: myNorth, mySouth, myEast, myWest
      integer :: myVecWE, myVecSN, root=0
      integer :: myIs, myIe, myJs, myJe, ierr
      integer :: status(mpi_status_size)
      end module mpiVariables

!**************************************************

      program iblank_MPI
      use variables
      use mpiVariables
      implicit none
      integer :: i, fid
      real(kind=8) :: time
      character(len=strL) :: fname, fIm, fFa

      if (nsd .ne. 3) write(*,ftab4) &
         "ERROR: iblank is computed for 3D problems only"
      call initMPI

      i = iargc()
      if (i .ne. 1) then
         write(*,ftab4) "ERROR: one input argument needed"
         STOP
      end if
      call getarg(1, fname)

      ilog = .false.
      open(fid,file=trim(fname))
      read(fid,*)
      read(fid,'(A)') fFa
      read(fid,*)
      read(fid,'(A)') fIm
      read(fid,*)
      read(fid,*) i
      if (i .ne. 0) ilog = .true.
      close(fid)

      call readFace(fa, fFa)

      call readImage(im, fIm)

      call initParDomain

      call initGrid

      time = mpi_wtime()
      call setIblankBody
      time = mpi_wtime() - time

      if (ilog) then
         write(fidP,ftab1, advance='no') 'CPU time for iblank: '
         write(fidP,'(F7.3,A)') time, ' s'
      end if
      write(*,ftab1, advance='no') 'CPU time for iblank: '
      write(*,'(F7.3,A)') time, ' s'

      call writeIblankBody

      if (ilog) close(fidP)
      call destroy(im)
      call destroy(fa)
      call mpi_finalize(ierr)

      end program iblank_MPI

!**************************************************

      subroutine initMPI
      use variables
      use mpiVariables
      implicit none

      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,myRank,ierr)
      call mpi_comm_size(mpi_comm_world,nProcs,ierr)

      isMas = .false.
      if (myRank .eq. master) isMas = .true.

      isSeq = .false.
      if (nProcs .eq. 1) isSeq = .true.

      if (ilog) then
         write(fLog,'("log.",I5.5)') myRank
         open(fidP,file=trim(fLog),status='unknown')
      end if

      np(2) = floor(sqrt(real(nProcs,kind=8)))
      np(1) = nProcs / np(2)
      if(np(1)*np(2) .ne. nProcs) then
         write(*,ftab4) 'ERROR in number of input processors...'
         if (ilog) write(fidP,ftab4) &
            'ERROR in number of input processors...'
         call stopSim()
      end if

      return
      end subroutine initMPI

!**************************************************

      subroutine stopSim
      use variables
      use mpiVariables
      implicit none

      if (ilog) then
         write(fidP,ftab4) 'Terminating...'
         close(fidP)
      end if
      write(*,ftab4) 'Terminating...'
      call mpi_finalize(ierr)
      STOP

      end subroutine stopSim

!**************************************************

      subroutine readFace(lFa, fname)
      use variables
      use mpiVariables
      implicit none
      type(faceType), intent(inout) :: lFa
      character(len=strL), intent(in) :: fname
      type(vtkXMLType) :: vtp

      if (ilog) write(fidP,ftab1) "Reading face data.."
      istat = 0
      do
         call loadVTK(vtp, trim(fname), istat)
         if (istat .lt.0) exit

         call getVTK_numPoints(vtp, lFa%nNo,istat)
         if (istat .lt.0) exit

         call getVTK_numElems(vtp, lFa%nEl, istat)
         if (istat .lt.0) exit

         call getVTK_nodesPerElem(vtp, lFa%eNoN, istat)
         if (istat .lt.0) exit

         allocate(lFa%x(nsd,lFa%nNo), lFa%IEN(lFa%eNoN,lFa%nEl))
         allocate(lFa%xC(nsd,lFa%nEl), lFa%nV(nsd,lFa%nEl))
         call getVTK_pointCoords(vtp,lFa%x,istat)
         if (istat .lt.0) exit

         call getVTK_elemIEN(vtp,lFa%IEN,istat)
         lFa%IEN = lFa%IEN + 1
         exit
      end do

      call flushVTK(vtp)

      if (istat .lt. 0) then
         write(*,ftab4) "ERROR: vtk file read error"
         if (ilog) write(fidP,ftab4) "ERROR: vtk file read error"
         call stopSim()
      end if

      call faceInit(fa)

      return
      end subroutine readFace

!**************************************************

      subroutine readImage(lIm, fname)
      use variables
      use mpiVariables
      implicit none
      type(imageType), intent(inout) :: lIm
      character(len=strL), intent(in) :: fname
      type(vtkXMLType) :: vti

      if (ilog) write(fidP,ftab1) "Reading image data.."
      istat = 0
      do
         call loadVTK(vti, trim(fname), istat)
         if (istat .lt. 0) exit

         call getVTK_imageExtent(vti, lIm%iLims, istat)
         if (istat .lt. 0) exit

         call getVTK_imageOrigin(vti, lIm%x0, istat)
         if (istat .lt. 0) exit

         call getVTK_imageSpacing(vti, lIm%dx, istat)
         if (istat .lt. 0) exit

         call getVTK_pieceExtent(vti, lIm%pLims, istat)
         if (istat .lt. 0) exit

         nxG = lIm%pLims(2,1) - lIm%pLims(1,1) + 1
         nyG = lIm%pLims(2,2) - lIm%pLims(1,2) + 1
         nz  = lIm%pLims(2,3) - lIm%pLims(1,3) + 1
         allocate(lIm%ival(nxG*nyG*nz))
         call getVTK_pointData(vti, 'ImageScalars', lIm%ival, istat)

         exit
      end do

      call flushVTK(vti)

      if (istat .lt. 0) then
         if (ilog) write(fidP,ftab4) "ERROR: vtk file read error"
         write(*,ftab4) "ERROR: vtk file read error"
         call stopSim()
      end if

      return
      end subroutine readImage

!**************************************************

      subroutine writeImage(lIm, fname)
      use variables
      use mpiVariables
      implicit none
      type(imageType), intent(in) :: lIm
      character(len=strL) :: fname
      type(vtkXMLType) :: vti

      istat = 0
      do
         call vtkInitWriter(vti, trim(fName), istat)
         if ( istat.lt.0 ) exit

         call putVTK_imageExtent(vti, lIm%iLims, istat)
         if (istat .lt. 0) exit

         call putVTK_imageOrigin(vti, lIm%x0, istat)
         if (istat .lt. 0) exit

         call putVTK_imageSpacing(vti, lIm%dx, istat)
         if (istat .lt. 0) exit

         call putVTK_pieceExtent(vti, lIm%pLims, istat)
         if (istat .lt. 0) exit

         call putVTK_pointData(vti, 'ImageScalars', lIm%ival, istat)
         if (istat .lt. 0) exit

         call vtkWriteToFile(vti, istat)
         exit
      end do

      call flushVTK(vti)

      if (istat .lt. 0) then
         if (ilog) write(fidP,ftab4) "ERROR: vtk file write error"
         write(*,ftab4) "ERROR: vtk file write error"
         call stopSim()
      end if

      return
      end subroutine writeImage

!**************************************************

      subroutine faceInit(lFa)
      use variables
      implicit none
      type(faceType), intent(inout) :: lFa

      integer :: e, a, Ac
      real(kind=8) :: rtmp, u(nsd), v(nsd), nV(nsd)
      real(kind=8), allocatable :: xl(:,:)

      if (ilog) write(fidP,ftab1) &
         'Calculating body normals and centroids..'
      call checkNormalsDir(lFa)

      allocate(xl(nsd,lFa%eNoN))
      do e=1, lFa%nEl
         do a=1, lFa%eNoN
            Ac = lFa%IEN(a,e)
            xl(:,a) = lFa%x(:,Ac)
         end do

         u(:)  = xl(:,2) - xl(:,1)
         v(:)  = xl(:,3) - xl(:,1)
         nV(:) = CROSS(u, v)
         rtmp = sqrt(sum( nV(:)**2 ))

         lFa%nV(:,e) = nV(:)/rtmp
         lFa%xC(:,e) = sum(xl,dim=2)/real(lFa%eNoN,kind=8)
      end do
      deallocate(xl)

      return
      end subroutine faceInit

!**************************************************

      subroutine checkNormalsDir(lFa)
      use variables
      implicit none
      type(faceType), intent(inout) :: lFa
      integer :: a, e, Ac, Ec
      real(kind=8) :: u(nsd), v(nsd), r0(nsd), nV(nsd), s, smin, rtmp, &
         xCent(nsd)
      real(kind=8), allocatable :: xl(:,:)

      if (ilog) write(fidP,ftab2) 'Checking normals direction...'
      r0(1) = -1D1
      r0(2) = -1D1
      r0(3) = -1D1
      smin  = HUGE(u)
      do a=1, lFa%nNo
         s = sqrt( sum( (lFa%x(:,a)-r0(:))**2 ))
         if (s .lt. smin) then
            smin = s
            Ac   = a
         end if
      end do

      E_LOOP : do e=1, lFa%nEl
         do a=1, lFa%eNoN
            if (lFa%IEN(a,e) .eq. Ac) then
               Ec = e
               EXIT E_LOOP
            end if
         end do
      end do E_LOOP

      if (ilog) then
         write(fidP,ftab2) 'Closest Node: '//trim(STR(Ac))
         write(fidP,ftab2) 'Closest Element: '//trim(STR(Ec))
      end if

      allocate(xl(nsd,lFa%eNoN))
      do a=1, lFa%eNoN
         Ac = lFa%IEN(a,Ec)
         xl(:,a) = lFa%x(:,Ac)
      end do

      u(:)  = xl(:,2) - xl(:,1)
      v(:)  = xl(:,3) - xl(:,1)
      nV(:) = CROSS(u, v)
      xCent(:) = sum(xl,dim=2)/real(lFa%eNoN,kind=8)

      r0(:) = xCent(:) - r0(:)
      rtmp  = nV(1)*xCent(1) + nV(2)*xCent(2) + nV(3)*xCent(3)

      if (ilog) write(fidP,ftab2) 'Vector product='//trim(STR(rtmp))
      if(rtmp .lt. 0D0) then
         if (ilog) then
            write(fidP,ftab3) &
               "Closest element's normal direction is inwards"
            write(fidP,ftab4) "Reordering elements..."
         end if
         do e=1, lFa%nEl
            Ac = lFa%IEN(2,e)
            lFa%IEN(2,e) = lFa%IEN(3,e)
            lFa%IEN(3,e) = Ac
         end do
      end if

      return
      end subroutine checkNormalsDir

!**************************************************

      subroutine initParDomain
      use variables
      use mpiVariables
      implicit none

      call mpi_cart_create(mpi_comm_world, 2, np, periods, reorder, &
         myCommCart, ierr)
      call mpi_comm_rank(myCommcart, myCartRank, ierr)
      call mpi_cart_coords(myCommCart, myCartRank, 2, myCoords, ierr)
      call mpi_cart_shift(myCommCart, Xdir, 1, myWest, myEast, ierr)
      call mpi_cart_shift(myCommCart, Ydir, 1, mySouth, myNorth, ierr)

      call dom_decomp(nxG, ICOORD, nx, myIs, myIe)
      call dom_decomp(nyG, JCOORD, ny, myJs, myJe)

      if (ilog) then
         write(fidP,ftab1) 'Domain decomposition..'
         write(fidP,ftab1) '-------------------------------------------'
         write(fidP,ftab1) 'Rank  Coords(1) Coords(2) myIs myJs  nx  ny'
         write(fidP,'(4X,3(2X,I3),4(2X,I5))') &
           myCartRank, myCoords(1), myCoords(2), myIs, myJs, nx, ny
         write(fidP,ftab1)
      end if

      imin = 1
      imax = nx
      jmin = 1
      jmax = ny
      kmin = 1-ngl
      kmax = nz+ngl
      if(myCoords(1) .eq. 0      ) imin = 1  - ngl
      if(myCoords(1) .eq. np(1)-1) imax = nx + ngl
      if(myCoords(2) .eq. 0      ) jmin = 1  - ngl
      if(myCoords(2) .eq. np(2)-1) jmax = ny + ngl

      ill = 1  - ngl
      iul = nx + ngl
      jll = 1  - ngl
      jul = ny + ngl
      kll = 1
      kul = nz
      if(myCoords(1) .eq. 0      ) ill = 1
      if(myCoords(1) .eq. np(1)-1) iul = nx
      if(myCoords(2) .eq. 0      ) jll = 1
      if(myCoords(2) .eq. np(2)-1) jul = ny

      allocate(xc(1-ngl:nxG+ngl))
      allocate(yc(1-ngl:nyG+ngl))
      allocate(zc(1-ngl:nz+ngl))
      allocate(iblank(1-ngl:nx+ngl,1-ngl:ny+ngl,1-ngl:nz+ngl))

      xc = 0D0
      yc = 0D0
      zc = 0D0
      iblank = 0

      return
      end subroutine initParDomain

!**************************************************

      subroutine dom_decomp(nG, dir, n, ll, ul)
      use mpiVariables
      implicit none
      integer, intent(in) :: nG, dir
      integer, intent(out) :: n, ll, ul
      integer :: res

      n = nG/np(dir)
      res = mod(nG, np(dir))

      if(myCoords(dir) .lt. res) then
         n  = n+1
         ll = myCoords(dir)*n + 1
      else
         ll = myCoords(dir)*n + res + 1
      end if
      ul = ll+n-1

      return
      end subroutine dom_decomp

!**************************************************

      subroutine initGrid
      use variables
      use mpiVariables
      implicit none
      integer :: i, j, k

      if (ilog) write(fidP,ftab1) "Initializing image voxel grid.."
      do i=1, nxG
         xc(i) = im%x0(1) + real(i-1,kind=8)*im%dx(1)
      end do

      do j=1, nyG
         yc(j) = im%x0(2) + real(j-1,kind=8)*im%dx(2)
      end do

      do k=1, nz
         zc(k) = im%x0(3) + real(k-1,kind=8)*im%dx(3)
      end do

      do i=1, ngl
        xc(1-i) = 2.0d0*xc(1) - xc(1+i)
        xc(nxG+i) = 2.0d0*xc(nxG) - xc(nxG-i)
      end do

      do j=1, ngl
        yc(1-j) = 2.0d0*yc(1) - yc(1+j)
        yc(nyG+j) = 2.0d0*yc(nyG) - yc(nyG-j)
      end do

      do k=1, ngl
        zc(1-k) = 2.0d0*zc(1) - zc(1+k)
        zc(nz+k) = 2.0d0*zc(nz) - zc(nz-k)
      end do

      return
      end subroutine initGrid

!**************************************************

      subroutine setIblankBody
      use variables
      use mpiVariables
      implicit none
      integer, parameter :: UNDECIDED=100000, DECIDED=-UNDECIDED
      integer :: i, j, k, e, a, Ac
      integer :: iCellMin, iCellMax, iCellMinG, iCellMaxG
      integer :: jCellMin, jCellMax, jCellMinG, jCellMaxG
      integer :: kCellMin, kCellMax
      integer :: sumBlank, coronaSize, cElement
      logical :: decidedYet
      real(kind=8) :: dotNorm, xMin(3), xMax(3)

      integer, allocatable :: iblankUndecided(:,:,:), iblankTemp(:,:,:)
      real(kind=8), allocatable :: xl(:,:)

      if (ilog) write(fidP,ftab1) 'Setting up iblank'

      allocate(iblankUndecided(1-ngl:nx+ngl,1-ngl:ny+ngl,1-ngl:nz+ngl))
      allocate(iblankTemp(1-ngl:nx+ngl,1-ngl:ny+ngl,1-ngl:nz+ngl))
      allocate(xl(nsd,fa%eNoN))

      iblankUndecided(ill:iul,jll:jul,kll:kul) = UNDECIDED
      iblankTemp = 0
      coronaSize = 1

      do e=1, fa%nEl
         do a=1, fa%eNoN
            Ac = fa%IEN(a,e)
            xl(:,a) = fa%x(:,Ac)
         end do

         xMin(1) = minval(xl(1,:))
         xMax(1) = maxval(xl(1,:))
         xMin(2) = minval(xl(2,:))
         xMax(2) = maxval(xl(2,:))
         xMin(3) = minval(xl(3,:))
         xMax(3) = maxval(xl(3,:))

         if ( (xMax(1) .le. xc(1  )) .or. &
              (xMin(1) .ge. xc(nxG)) .or. &
              (xMax(2) .le. yc(1  )) .or. &
              (xMin(2) .ge. yc(nyG)) .or. &
              (xMax(3) .le. zc(1  )) .or. &
              (xMin(3) .ge. zc(nz )) ) cycle

        iCellMinG = 0; iCellMaxG = 0
        jCellMinG = 0; jCellMaxG = 0
        kCellMin  = 0; kCellMax  = 0

!       i-direction
!       -----------
        do i=0, nxG
           if ( (xc(i)-xMin(1))*(xc(i+1)-xMin(1)) .lt. 0D0 ) &
              iCellMinG = (i  )-(coronaSize-1)
           if ( (xc(i)-xMax(1))*(xc(i+1)-xMax(1)) .lt. 0D0 ) &
              iCellMaxG = (i+1)+(coronaSize-1)
        end do
        if (iCellMinG.le.0)                         iCellMinG = 1
        if (iCellMaxG.eq.0 .or. iCellMaxG.ge.nxG+1) iCellMaxG = nxG
        iCellMin = iCellMinG - (myIs-1)
        iCellMax = iCellMaxG - (myIs-1)
        if (iCellMax.lt.ill .or. iCellMin.gt.iul) cycle
        iCellMin = max(iCellMin, ill)
        iCellMax = min(iCellMax, iul)

!       j-direction
!       -----------
        do j=0, nyG
           if ( (yc(j)-xMin(2))*(yc(j+1)-xMin(2)) .lt. 0D0 ) &
              jCellMinG = (j  )-(coronaSize-1)
           if ( (yc(j)-xMax(2))*(yc(j+1)-xMax(2)) .lt. 0D0 ) &
              jCellMaxG = (j+1)+(coronaSize-1)
        end do
        if (jCellMinG.le.0)                         jCellMinG = 1
        if (jCellMaxG.eq.0 .or. jCellMaxG.ge.nyG+1) jCellMaxG = nyG
        jCellMin = jCellMinG - (myJs-1)
        jCellMax = jCellMaxG - (myJs-1)
        if (jCellMax.lt.jll .or. jCellMin.gt.jul) cycle
        jCellMin = max(jCellMin, jll)
        jCellMax = min(jCellMax, jul)

!       k-direction
!       ------------
        do k=0, nz
           if ( (zc(k)-xMin(3))*(zc(k+1)-xMin(3)) .lt. 0D0 ) &
              kCellMin = (k  )-(coronaSize-1)
           if ( (zc(k)-xMax(3))*(zc(k+1)-xMax(3)) .lt. 0D0 ) &
              kCellMax = (k+1)+(coronaSize-1)
        end do
        if (kCellMin.le.0)                       kCellMin = 1
        if (kCellMax.eq.0 .or. kCellMax.ge.nz+1) kCellMax = nz
        if (kCellMax.lt.1 .or. kCellMin.gt.nz) cycle

        iblankUndecided(iCellMin:iCellMax, &
                        jCellMin:jCellMax, &
                        kCellMin:kCellMax) = DECIDED
      end do ! e

      do k=kll, kul
        do j=jll, jul
          decidedYet=.false.
          do i=ill, iul
            dotNorm = -100.0d0
            if ( iblankUndecided(i,j,k) .eq. DECIDED ) then
              decidedYet=.true.
              call search_vertex_dotNorm(i, j, k, cElement, dotNorm)
              if (dotNorm .ge. 0.0d0) iblankTemp(i,j,k) = 1
            end if
          end do ! i

          if (.not.decidedYet) then
            i=ill
            iblankUndecided(i,j,k) = DECIDED
            call search_vertex_dotNorm(i, j, k, cElement, dotNorm)
            if (dotNorm .ge. 0.0d0) then
              iblankTemp(i,j,k) = 1
            end if
          end if
        end do ! j
      end do ! k

      do k=kll, kul
        do j=jll, jul
          do i=ill, iul
            if ( iblankUndecided(i,j,k).eq.DECIDED ) then
              iblankTemp(ill,j,k) = iblankTemp(i,j,k)
              exit
            end if
          end do
        end do
      end do

      do k=kll, kul
        do j=jll, jul
          do i=ill+1, iul
            if ( iblankUndecided(i,j,k) .eq. UNDECIDED ) THEN
              iblankTemp(i,j,k) = iblankTemp(i-1,j,k)
            end if
          end do
        end do
      end do

      do k=kll, kul
        do j=jll, jul
          do i=ill, iul
            iblank(i,j,k) = iblankTemp(i,j,k)
          end do
        end do
      end do

      call find_iblankHoles_outerboundary

      call find_isolatedIblank

      deallocate(iblankUndecided, iblankTemp, xl)

      return
      end subroutine setIblankBody

!**************************************************

      subroutine search_vertex_dotNorm(iCell, jCell, kCell, closestElement, dotNorm)
      use variables
      use mpiVariables
      implicit none
      integer, intent(in) :: iCell,jCell,kCell
      integer, intent(out) :: closestElement
      real(kind=8), intent(out) :: dotNorm
      integer, parameter :: NSIZE=200,MSIZE=20
      integer, parameter :: nCheck=3

      integer :: m,n,nc
      integer, dimension(:), allocatable :: NeighElemInd
      real(kind=8), dimension(:), allocatable :: distMarkerS
      integer :: nBodySelect,numNeighElement
      integer :: elemInd,node1,node2,node3
      integer :: shortestProbe
      integer,dimension(1)     :: iDummy(1)
      integer,dimension(MSIZE) :: cElement,closestVert
      real(kind=8) :: xCell,yCell,zCell
      real(kind=8) :: dsIntercept,xM,yM,zM
      real(kind=8) :: areaDiffMin
      real(kind=8) :: distBIElem, distBIElemMin
      real(kind=8) :: planeConst,distanceToPlane,distPointToPlane
      real(kind=8) :: side12,side23,side31,side14,side24,side34
      real(kind=8) :: area123,area124,area234,area314
      real(kind=8) :: semiPerimeter123,semiPerimeter124,semiPerimeter234,semiPerimeter314
      real(kind=8) :: epsiArea,areaDiff
      real(kind=8) :: xBI, yBI, zBI
      real(kind=8) :: xBITemp, yBITemp, zBITemp
      real(kind=8),dimension(3)     :: xVert, yVert, zVert
      real(kind=8),dimension(MSIZE) :: xBIT,yBIT,zBIT,dist

      xCell = xc(myIs+iCell-1)
      yCell = yc(myJs+jCell-1)
      zCell = zc(kCell)

      dotNorm = -1.0d0

      allocate(distMarkerS(fa%nNo),neighElemInd(NSIZE))
      do m=1, fa%nNo
         distMarkerS(m) = (fa%x(1,m)-xCell)**2 + &
                          (fa%x(2,m)-yCell)**2 + &
                          (fa%x(3,m)-zCell)**2
      end do

      do nc=1, nCheck
         iDummy = minloc(distMarkerS(1:fa%nNo))
         closestVert(nc) = iDummy(1)
         distMarkerS(closestVert(nc)) = 1.0E+20
      end do

      do nc=1, nCheck
         numNeighElement = 0
         do m=1, fa%nEl
            if ( fa%IEN(1,m) .eq. closestVert(nc) .or. &
                 fa%IEN(2,m) .eq. closestVert(nc) .or. &
                 fa%IEN(3,m) .eq. closestVert(nc) ) then
               numNeighElement = numNeighElement + 1
               NeighElemInd(numNeighElement) = m
            end if
         end do ! m

         distBIElemMin = 1.0E+16
         areaDiffMin   = 1.0E+16
         epsiArea      = 1.0E-4
         closestElement = 0

         do n=1, numNeighElement
            elemInd = NeighElemInd(n)
            node1   = fa%IEN(1,elemInd)
            node2   = fa%IEN(2,elemInd)
            node3   = fa%IEN(3,elemInd)
            call check_BIInsideTriangle(elemInd, node1, node2, node3, &
               xCell, yCell, zCell, xBITemp, yBITemp, zBITemp, &
               area123, areaDiff)
            if (dabs(areaDiff) .lt. epsiArea*area123) then
               xBI = xBITemp
               yBI = yBITemp
               zBI = zBITemp
               closestElement = elemInd
               dist(nc) = 0.0d0
               exit
            else
               call calc_BIOutsideTriangle(elemInd, node1, node2, &
                  node3, xBITemp, yBITemp, zBITemp, distBIElem, &
                  area123, areaDiff)
               if (distBIElem .le. distBIElemMin) then
                  distBIElemMin = distBIElem
                  closestElement = elemInd
                  xBI = xBITemp
                  yBI = yBITemp
                  zBI = zBITemp
               end if
               dist(nc) = distBIElemMin
            end if
         end do ! n

         xBIT(nc) = xBI
         yBIT(nc) = yBI
         zBIT(nc) = zBI
         cElement(nc) = closestElement
      end do ! nc

      iDummy = minloc(dist(1:nCheck))
      shortestProbe = iDummy(1)
      closestElement = cElement(shortestProbe)

      dotNorm = (xCell - fa%xc(1,closestElement)) &
                        *fa%nV(1,closestElement)  &
              + (yCell - fa%xc(2,closestElement)) &
                        *fa%nV(2,closestElement)  &
              + (zCell - fa%xc(3,closestElement)) &
                        *fa%nV(3,closestElement)

      deallocate(neighElemInd,distMarkerS)

      return
      end subroutine search_vertex_dotNorm

!**************************************************

      subroutine check_BIInsideTriangle(elemInd, node1, node2, node3, &
         xCell, yCell, zCell, xBITemp, yBITemp, zBITemp, area123, &
         areaDiff)
      use variables
      use mpiVariables
      implicit none
      integer, intent(in) :: elemInd, node1, node2, node3
      real(kind=8), intent(in)  :: xCell, yCell, zCell
      real(kind=8), intent(out) :: xBITemp, yBITemp, zBITemp, area123, &
         areaDiff

      integer :: iside
      real(kind=8) :: planeConst, distanceToPlane, distPointToPlane
      real(kind=8) :: side12, side23, side31, side14, side24, side34
      real(kind=8) :: area124, area234, area314
      real(kind=8) :: semiPerimeter123, semiPerimeter124
      real(kind=8) :: semiPerimeter234, semiPerimeter314

      real(kind=8) :: denom, num1, num2
      real(kind=8) :: VX, VY, VZ, PX, PY, PZ, RR, Rx, Ry, Rz
      real(kind=8) :: UX, UY, UZ, WX, WY, WZ, SS, TT, PIx, PIy, PIz, tol

      planeConst = - fa%nV(1,elemInd)*fa%x(1,node1) &
                   - fa%nV(2,elemInd)*fa%x(2,node1) &
                   - fa%nV(3,elemInd)*fa%x(3,node1)

      distanceToPlane = - ( fa%nV(1,elemInd)*xCell  &
                          + fa%nV(2,elemInd)*yCell  &
                          + fa%nV(3,elemInd)*zCell  &
                          + planeConst )

      xBITemp = xCell + fa%nV(1,elemInd)*distanceToPlane
      yBITemp = yCell + fa%nV(2,elemInd)*distanceToPlane
      zBITemp = zCell + fa%nV(3,elemInd)*distanceToPlane

      tol = 1.e-10

      RX = xBITemp
      RY = yBITemp
      RZ = zBITemp

      UX = fa%x(1,node2) - fa%x(1,node1)
      UY = fa%x(2,node2) - fa%x(2,node1)
      UZ = fa%x(3,node2) - fa%x(3,node1)

      VX = fa%x(1,node3) - fa%x(1,node1)
      VY = fa%x(2,node3) - fa%x(2,node1)
      VZ = fa%x(3,node3) - fa%x(3,node1)

      WX = RX - fa%x(1,node1)
      WY = RY - fa%x(2,node1)
      Wz = RZ - fa%x(3,node1)

      denom = (UX*VX+UY*VY+UZ*VZ)**2 - (UX*UX+UY*UY+UZ*UZ)*(VX*VX+VY*VY+VZ*VZ)
      num1 = (UX*VX+UY*VY+UZ*VZ)*(WX*VX+WY*VY+WZ*VZ) - (VX*VX+VY*VY+VZ*VZ)*(WX*UX+WY*UY+WZ*UZ)
      num2 = (UX*VX+UY*VY+UZ*VZ)*(WX*UX+WY*UY+WZ*UZ) - (UX*UX+UY*UY+UZ*UZ)*(WX*VX+WY*VY+WZ*VZ)
      SS = num1/denom
      TT = num2/denom

      if( (SS.ge.(0.0-tol)) .and. &
          (TT.ge.(0.0-tol)) .and. &
          (SS+TT.le.(1.0+tol)) ) then
         area123  = 1.0
         areadiff = 0.0
      else
         area123  = 1.0
         areadiff = 1.0
      end if

      return
      end subroutine check_BIInsideTriangle

!**************************************************

      subroutine calc_BIOutsideTriangle(elemInd, node1, node2, node3, &
         xBITemp, yBITemp, zBITemp, distBIElem, area1, aread)
      use variables
      use mpiVariables
      implicit none
      integer, intent(in) :: elemInd, node1, node2, node3
      real(kind=8), intent(in) :: xBITemp, yBITemp, zBITemp
      real(kind=8), intent(out) :: distBIElem

      real(kind=8) :: area1,  aread
      integer :: iside,  i
      integer :: isideSelect, mside
      integer :: nodeSelect1, nodeSelect2
      real(kind=8) :: aCrossbVectMagn, dotVal
      real(kind=8) :: distIntBI, distIntCG, distIntMin
      real(kind=8) :: distNorm, distNode1BINorm, distNode2BINorm
      real(kind=8) :: distVert1BI, distVert2BI
      real(kind=8) :: magnitude12, magnitudeBICG, magnitude, projectedLength
      real(kind=8) :: node12x, node12y, node12z
      real(kind=8) :: vec01x, vec01y, vec01z
      real(kind=8) :: vec12x, vec12y, vec12z
      real(kind=8) :: xBINorm, yBINorm, zBINorm
      real(kind=8) :: xCG, yCG, zCG
      real(kind=8), dimension(3) :: aVect, bVect, cVect, aCrossbVect, cCrossbVect
      real(kind=8), dimension(3) :: xInt, yInt, zInt
      real(kind=8), dimension(3) :: xVert, yVert, zVert
      real(kind=8), dimension(3) :: vect1, vect2, vect3, vect4
      real(kind=8), dimension(3,3) :: vectInt

      xVert(1)  = fa%x(1,node1)
      yVert(1)  = fa%x(2,node1)
      zVert(1)  = fa%x(3,node1)

      xVert(2)  = fa%x(1,node2)
      yVert(2)  = fa%x(2,node2)
      zVert(2)  = fa%x(3,node2)

      xVert(3)  = fa%x(1,node3)
      yVert(3)  = fa%x(2,node3)
      zVert(3)  = fa%x(3,node3)

      xCG       = fa%xc(1,elemInd)
      yCG       = fa%xc(2,elemInd)
      zCG       = fa%xc(3,elemInd)

      vect1(1:3) = (/xBITemp, yBITemp, zBITemp/)
      vect2(1:3) = (/xCG, yCG, zCG/)

      aVect(1:3) = vect2(1:3)-vect1(1:3)

      do iside=1, 3
        mside = iside +1
        if(iside.eq.3) mside = 1

        vect3(1:3) =(/xVert(iside), yVert(iside), zVert(iside)/)
        vect4(1:3) =(/xVert(mside), yVert(mside), zVert(mside)/)

        bVect(1:3)  = vect4(1:3) - vect3(1:3)
        cVect(1:3)  = vect3(1:3) - vect1(1:3)

        aCrossbVect = CROSS(aVect, bVect)
        cCrossbVect = CROSS(cVect, bVect)

        aCrossbVectMagn = aCrossbVect(1)**2 + aCrossbVect(2)**2 +aCrossbVect(3)**2
        dotVal = dot_product(cCrossbVect,aCrossbVect)
        vectInt(1:3,iside) = vect1(1:3) + aVect(1:3) *dotVal/aCrossbVectMagn
      end do

      magnitudeBICG = dsqrt( (vect1(1)-vect2(1))**2 &
                           + (vect1(2)-vect2(2))**2 &
                           + (vect1(3)-vect2(3))**2 )

      distIntMin = 1.0E+16
      isideSelect = -1000

      do iside = 1,3
         distIntBI = dsqrt( (vect1(1)-vectInt(1,iside))**2 &
                          + (vect1(2)-vectInt(2,iside))**2 &
                          + (vect1(3)-vectInt(3,iside))**2 )/magnitudeBICG

         distIntCG = dsqrt( (vect2(1)-vectInt(1,iside))**2 &
                          + (vect2(2)-vectInt(2,iside))**2 &
                          + (vect2(3)-vectInt(3,iside))**2 )/magnitudeBICG

         if ( distIntBI.le.1.0d0 .and. distIntCG.le.1.0d0) then
            distIntMin  = DMIN1(distIntBI,distIntCG)
            isideSelect = iside
         end if
      end do

      select case(isideSelect)
      case(1)
         nodeSelect1 = 1
         nodeSelect2 = 2
      case(2)
         nodeSelect1 = 2
         nodeSelect2 = 3
      case(3)
         nodeSelect1 = 3
         nodeSelect2 = 1
      end select

      vec12x = xVert(nodeSelect2) - xVert(nodeSelect1)
      vec12y = yVert(nodeSelect2) - yVert(nodeSelect1)
      vec12z = zVert(nodeSelect2) - zVert(nodeSelect1)

      magnitude12 = dsqrt(vec12x**2 + vec12y**2 + vec12z**2)

      vec01x = xVert(nodeSelect1) - xBITemp
      vec01y = yVert(nodeSelect1) - yBITemp
      vec01z = zVert(nodeSelect1) - zBITemp

      distNorm = dsqrt(  (vec12y*vec01z - vec12z*vec01y)**2  &
                       + (vec12z*vec01x - vec01z*vec12x)**2  &
                       + (vec12x*vec01y - vec01x*vec12y)**2  )/magnitude12

      node12x = xVert(nodeSelect2) - xVert(nodeSelect1)
      node12y = yVert(nodeSelect2) - yVert(nodeSelect1)
      node12z = zVert(nodeSelect2) - zVert(nodeSelect1)

      magnitude = dsqrt(node12x**2 + node12y**2 + node12z**2)

      node12x = node12x/magnitude
      node12y = node12y/magnitude
      node12z = node12z/magnitude

      projectedLength = (xBITemp - xVert(nodeSelect1))*node12x  &
                       +(yBITemp - yVert(nodeSelect1))*node12y  &
                       +(zBITemp - zVert(nodeSelect1))*node12z

      xBINorm = xVert(nodeSelect1) + projectedLength*node12x
      yBINorm = yVert(nodeSelect1) + projectedLength*node12y
      zBINorm = zVert(nodeSelect1) + projectedLength*node12z

      distNode1BINorm = dsqrt( (xVert(nodeSelect1)-xBINorm)**2 &
                             + (yVert(nodeSelect1)-yBINorm)**2 &
                             + (zVert(nodeSelect1)-zBINorm)**2 )/magnitude

      distNode2BINorm = dsqrt( (xVert(nodeSelect2)-xBINorm)**2 &
                             + (yVert(nodeSelect2)-yBINorm)**2 &
                             + (zVert(nodeSelect2)-zBINorm)**2 )/magnitude

      if ( distNode1BINorm.le.1.0d0 .and. distNode2BINorm.le.1.0d0) then
         distBIElem  =  distNorm
      else
         distVert1BI = dsqrt( (xVert(nodeSelect1)-xBITemp)**2 &
                            + (yVert(nodeSelect1)-yBITemp)**2 &
                            + (zVert(nodeSelect1)-zBITemp)**2 )

         distVert2BI = dsqrt( (xVert(nodeSelect2)-xBITemp)**2 &
                            + (yVert(nodeSelect2)-yBITemp)**2 &
                            + (zVert(nodeSelect2)-zBITemp)**2 )

         if (distVert1BI <= distVert2BI) then
            distBIElem  = distVert1BI
         else
            distBIElem  = distVert2BI
         end if
      end if

      return
      end subroutine calc_BIOutsideTriangle

!**************************************************

      subroutine find_iblankHoles_outerboundary
      use variables
      use mpiVariables
      implicit none
      integer :: i, j, k, iG, jG
      logical :: flag

      do k=kll, kul
         do j=jll, jul
            jG = myJs+j-1
            do i=ill, iul
               iG = myIs+i-1
               if(iG.lt.1 .or. iG.gt.nxG-1 .or. &
                  jG.lt.1 .or. jG.gt.nyG-1) cycle

               flag = .false.
               if (iblank(i,j,k).eq.0) then
                  if (iG .eq. nxG-1) then
                     if (iblank(i-1,j,k) .eq. 1) flag = .true.
                  end if
                  if (iG .eq. 1) then
                     if (iblank(i+1,j,k) .eq. 1) flag = .true.
                  end if
                  if (jG .eq. nyG-1) then
                     if (iblank(i,j-1,k) .eq. 1) flag = .true.
                  end if
                  if (jG .eq. 1) then
                     if (iblank(i,j+1,k) .eq. 1) flag = .true.
                  end if
                  if (k .eq. nz-1) then
                     if (iblank(i,j,k-1) .eq. 1) flag = .true.
                  end if
                  if (k .eq. 1) then
                     if (iblank(i,j,k+1) .eq. 1) flag = .true.
                  end if
               end if

               if (flag) then
                  if (ilog) then
                     write(fidP,ftab3) 'Found iblank hole at '// &
                        trim(STR(myIs+i-1))//' '//trim(STR(myJs+j-1))//&
                        ' '//trim(STR(k))
                     write(fidP,ftab3) 'This cell is blanked..'
                  end if
                  iblank(i,j,k) = 1
               end if
            end do
         end do
      end do

      return
      end subroutine find_iblankHoles_outerboundary

!**************************************************

      subroutine find_isolatedIblank
      use variables
      use mpiVariables
      implicit none
      integer :: i1, i2, j1, j2, k1, k2
      integer :: i, j, k, scanRange, sumIblank

      do k=kll, kul
         do j=jll, jul
            do i=ill, iul
               if (iBlank(i,j,k).eq.1) then
                  scanRange = 1
                  i1 = max(ill,i-scanRange)
                  i2 = min(iul,i+scanRange)
                  j1 = max(jll,j-scanRange)
                  j2 = min(jul,j+scanRange)
                  k1 = max(kll,k-scanRange)
                  k2 = min(kul,k+scanRange)

                  sumIblank = sum(iBlank(i1:i2, &
                                         j1:j2, &
                                         k1:k2) )

                  if(sumIblank .le. 1) then
                     if (ilog) then
                        write(fidP,ftab3) &
                           'Found Isolated blank Cell at '// &
                           trim(STR(myIs+i-1))//' '//&
                           trim(STR(myJs+j-1))//' '//trim(STR(k))
                        write(fidP,ftab3) 'This cell is unblanked'
                     end if
                     iBlank(i,j,k) = 0
                  end if
              end if
            end do
         end do
      end do

      return
      end subroutine find_isolatedIblank

!**************************************************

      subroutine writeIblankBody
      use variables
      use mpiVariables
      implicit none
      integer :: i, j, k, a, ip
      integer :: nxl, nyl, nzl, isl, iel, jsl, jel
      integer :: i1l, i2l, j1l, j2l, k1l, k2l
      character(len=strL) :: fname, cmd
      integer, allocatable :: tmpI(:,:,:)

      if (.not.isSeq) then
         call parCommVar(iblank, nx, ny, nz, ngl)
         write(fname,'(A)') "temp."//trim(STR(myRank))
         open(10,file=trim(fname),form="unformatted")
         write(10) nx, ny, nz, myIs, myIe, myJs, myJe
         write(10) imin, imax, jmin, jmax, kmin, kmax
         write(10) (((iblank(i,j,k),i=imin,imax),j=jmin,jmax),k=kmin,kmax)
         close(10)
         call mpi_barrier(mpi_comm_world, ierr)

         if (myRank .eq. master) then
            allocate(tmpI(1-ngl:nxG+ngl,1-ngl:nyG+ngl,1-ngl:nz+ngl))
            do ip=1, nProcs
               write(fname,'(A)') "temp."//trim(STR(ip-1))
               open(10,file=trim(fName),form="unformatted")
               read(10) nxl, nyl, nzl, isl, iel, jsl, jel
               read(10) i1l, i2l, j1l, j2l, k1l, k2l
               read(10) (((tmpI(isl+i-1,jsl+j-1,k),i=i1l,i2l),j=j1l,j2l),k=k1l,k2l)
               close(10)
            end do

            a = 0
            do k=1, nz
               do j=1, nyG
                  do i=1, nxG
                     a = a + 1
                     im%iVal(a) = tmpI(i,j,k)
                  end do
               end do
            end do
            deallocate(tmpI)
         end if

         i = 0
         call mpi_allreduce(i, j, 1, mpint, mpi_sum, mpi_comm_world, ierr)

         write(fname,'(A)') "temp."//trim(STR(myRank))
         write(cmd,'(A)') "rm  "//trim(fName)
         call system(trim(cmd))
         if (myRank .ne. master) then
            return
         end if
      else
         a = 0
         do k=1, nz
            do j=1, nyG
               do i=1, nxG
                  a = a + 1
                  im%iVal(a) = iblank(i,j,k)
               end do
            end do
         end do
      end if

      if (ilog) write(fidP,ftab1) 'Writing iblank to file...'
      write(fname,'(A)') "mask.vti"
      call writeImage(im, fname)

      return
      end subroutine writeIblankBody

!**************************************************

      subroutine parCommVar(u, IL, JL, KL, glN)
      use mpiVariables
      implicit none
      integer, intent(in) :: IL, JL, KL, glN
      integer, intent(inout) :: u(1-glN:IL+glN,1-glN:JL+glN,1-glN:KL+glN)

      call mpi_type_vector((KL+2*glN)*(JL+2*glN), glN, IL+2*glN, mpint,&
         myVecWE, ierr)
      call mpi_type_commit(myVecWE, ierr)

      call mpi_type_vector(KL+2*glN, (IL+2*glN)*glN, (IL+2*glN)*(JL+2*glN), &
         mpint, myVecSN, ierr)
      call mpi_type_commit(myVecSN, ierr)

      ! west to east !
      call mpi_sendrecv(u(IL+1-glN, 1-glN, 1-glN), 1, myVecWE, myEast, 0,  &
                        u(   1-glN, 1-glN, 1-glN), 1, myVecWE, myWest, 0,  &
                        myCommCart, status, ierr)

      ! east to west !
      call mpi_sendrecv(u(1   , 1-glN, 1-glN), 1, myVecWE, myWest, 1,  &
                        u(IL+1, 1-glN, 1-glN), 1, myVecWE, myEast, 1,  &
                        myCommCart, status, ierr)

      ! south to north !
      call mpi_sendrecv(u(1-glN, JL+1-glN, 1-glN), 1, myVecSN, myNorth, 0,  &
                        u(1-glN,    1-glN, 1-glN), 1, myVecSN, mySouth, 0,  &
                        myCommCart, status, ierr)

      ! north to south !
      call mpi_sendrecv(u(1-glN, 1   , 1-glN), 1, myVecSN, mySouth, 1,  &
                        u(1-glN, JL+1, 1-glN), 1, myVecSN, myNorth, 1,  &
                        myCommCart, status, ierr)

      call mpi_type_free(myVecWE, ierr)
      call mpi_type_free(myVecSN, ierr)

      return
      end subroutine parCommVar

!******************************************************************************

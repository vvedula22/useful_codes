!--------------------------------------------------------------------
!
!     This code computes oscillatory shear index (OSI) from surface
!     wall shear stress (WSS) data in vtp format.
!
!--------------------------------------------------------------------

!     Common variables are defined in this module
      module variables
      use genUtils
      use stdParams
      use vtkXMLMod
      implicit none

      type faceType
!        Number of face nodes
         integer :: nNo
!        Number of face elements
         integer :: nEl
!        Element number of nodes
         integer :: eNoN
!        Number of Gauss points
         integer :: nG
!        VTK type of the element
         integer :: vtkType
!        Element connectivity array
         integer, allocatable :: IEN(:,:)
!        iblank field defined in any region of interest
         integer, allocatable :: iblank(:)

!        Total surface area
         real(kind=8) :: area
!        Gauss point weights
         real(kind=8), allocatable :: w(:)
!        Gauss point coordinates in parent configuration
         real(kind=8), allocatable :: xi(:,:)
!        Element shape functions
         real(kind=8), allocatable :: N(:,:)
!        Element shape function derivatives
         real(kind=8), allocatable :: Nx(:,:,:)
!        Position coordinates
         real(kind=8), allocatable :: x(:,:)
!        Element normal
         real(kind=8), allocatable :: nV(:,:)
!        Face displacement field (if any)
         real(kind=8), allocatable :: dx(:,:)
!        Surface wall shear stress field
         real(kind=8), allocatable :: wss(:,:)
      end type faceType

!     Number of spatial dimensions (=2 or 3)
      integer :: nsd

!     File sequence controls: nstart - starting time step; nend -
!     ending time step; nfreq - frequency
      integer :: nstart, nend, nfreq

!     String to store ROI file path
      character(len=strL) :: froi

!     Vector-averaged WSS
      real(kind=8), allocatable :: awss(:,:)

!     Time-averaged WSS
      real(kind=8), allocatable :: tawss(:)

!     OSI
      real(kind=8), allocatable :: osi(:)

      type(faceType) :: fa0, fa

      contains

         pure function cross(U) result(V)
         implicit none
         real(kind=8), intent(in) :: U(:,:)
         real(kind=8) V(size(U,1))

         if (size(U,1) .eq. 2) then
            V(1) =  U(2,1)
            V(2) = -U(1,1)
         else
            V(1) = U(2,1)*U(3,2) - U(3,1)*U(2,2)
            V(2) = U(3,1)*U(1,2) - U(1,1)*U(3,2)
            V(3) = U(1,1)*U(2,2) - U(2,1)*U(1,2)
         end if

         return
         end function cross

      end module variables

!**************************************************

      program osi_calc
      use variables
      implicit none
      integer :: i, a, e, cnt, ntot, fid, ntime, iaawss
      real(kind=8) :: rtmp, aawss, dt, time
      character(len=strL) :: fname
      character(len=strL) :: prefix

      interface
         subroutine readVTP(lFa, fname, ipass)
         use variables
         implicit none
         type(faceType), intent(inout) :: lFa
         character(len=strL), intent(in) :: fname
         integer, intent(in), optional :: ipass
         end subroutine readVTP
      end interface

      fid = 10
      i = iargc()
      if (i .ne. 1) then
         write(stdout,ftab4) "Error: at least one input argument needed"
         STOP
      end if
      call getarg(1,fname)

!     Read input data from file
      iaawss = 0
      open(fid,file=trim(fname))
      read(fid,*)
      read(fid,*) nsd
      read(fid,*)
      read(fid,*) nstart, nend, nfreq
      read(fid,*)
      read(fid,*) dt
      read(fid,*)
      read(fid,'(A)') prefix
      read(fid,*)
      read(fid,*) iaawss
      read(fid,*)
      read(fid,'(A)') froi
      close(fid)

      if (nfreq .ne. 0) then
         ntot = (nend - nstart)/nfreq
      else
         ntot = 1
      end if

!     If iaawss > 0, area-averaged WSS is output to file as a function
!     of time.
      if (iaawss .gt. 0) then
         write(fname,'(A)') "aawss_"//trim(STR(nstart))//"_"// &
         trim(STR(nend))//".dat"
         open(101,file=trim(fname))
         write(101,'(A)') 'Variables=t, AWSS'
      end if

      write(stdout,'(A)') "############################################"
      write(stdout,ftab1)
      write(stdout,ftab1) "Computing OSI for time range: "// &
         trim(STR(nstart))//" to "//trim(STR(nend))
      write(stdout,ftab1)

      cnt = 0
      ntime = nstart
      time  = DBLE(ntime) * dt

      write(fname,'(A)') trim(prefix)//trim(STR(cnt))//".vtp"
      write(stdout,ftab2) "Reading file "//trim(fname)//", nTime: "// &
         trim(STR(ntime))
      call readVTP(fa0, fname, 1)

      open(fid,file=trim(froi))
      do e=1, fa0%nEl
         read(fid,*,end=101) a
         fa0%iblank(a+1) = 0
 101     exit
      end do
      close(fid)

      fa%nNo  = fa0%nNo
      fa%nEl  = fa0%nEl
      fa%eNoN = fa0%eNoN
      allocate(fa%x(nsd,fa%nNo), fa%IEN(fa%eNoN,fa%nEl), &
         fa%dx(nsd,fa%nNo), fa%wss(nsd,fa%nNo), fa%iblank(fa%nEl))
      call selecteleb(fa)
      fa%x   = fa0%x
      fa%IEN = fa0%IEN
      fa%dx  = fa0%dx
      fa%wss = fa0%wss
      fa%iblank = fa0%iblank
      call getNormals(fa)

      allocate(awss(nsd,fa%nNo), tawss(fa%nNo))
      awss  = fa%wss
      tawss = 0D0
      do a=1, fa%nNo
         tawss(a) = tawss(a) + sqrt(sum(fa%wss(:,a)**2))
      end do
      if (iaawss .gt. 0) then
         call calcAAWSS(fa, aawss)
         write(101,'(2X,F9.4,2X,1pE18.6)') time, aawss
      end if

!     Loop to compute vector-averaged and time-averaged WSS
      do cnt=1, ntot
         ntime = nstart + cnt*nfreq
         time  = DBLE(ntime)*dt
         write(fname,'(A)') trim(prefix)//trim(STR(cnt))//".vtp"
         write(stdout,ftab2) "Reading file "//trim(fname)// &
            ", nTime: "//trim(STR(ntime))
         call readVTP(fa, fname)
         call getNormals(fa)
         awss(:,:) = awss(:,:) + fa%wss(:,:)
         do a=1, fa%nNo
            rtmp = sqrt(sum(fa%wss(:,a)**2))
            tawss(a) = tawss(a) + rtmp
         end do
         if (iaawss .gt. 0) then
            call calcAAWSS(fa, aawss)
            write(101,'(2X,F9.4,2X,1pE18.6)') time, aawss
         end if
      end do
      awss(:,:) = awss(:,:)/real(ntot,kind=8)
      tawss(:)  = tawss(:) /real(ntot,kind=8)
      if (iaawss .gt. 0) close(101)

!     Calculate OSI
      allocate(osi(fa%nNo))
      osi = 0D0
      do a=1, fa%nNo
         rtmp = sqrt(sum(awss(:,a)**2))
         osi(a) = 5D-1*(1D0 - rtmp/tawss(a))
      end do

!     OSI output to a vtp file
      write(fname,'(A)') "osi_"//trim(STR(nstart))//"_"// &
         trim(STR(nend))//".vtp"
      write(stdout,ftab1)
      write(stdout,ftab1) "Writing OSI to file "//trim(fname)
      fa%x   = fa0%x
      fa%IEN = fa0%IEN
      fa%dx  = fa0%dx
      fa%wss = fa0%wss
      call writeVTP(fa, fname)

      write(stdout,ftab1)
      write(stdout,'(A)') "############################################"
      write(stdout,ftab1)

      end program osi_calc

!**************************************************
!     Read data from a vtp file
      subroutine readVTP(lFa, fname, ipass)
      use variables
      implicit none
      type(faceType), intent(inout) :: lFa
      character(len=strL), intent(in) :: fname
      integer, intent(in), optional :: ipass

      type(vtkXMLType) :: vtp
      integer :: vtkType, istat, iopt
      integer :: e, nNo, eNoN, nEl
      real(kind=8), allocatable :: tmpX(:)

      iopt = 0
      if (present(ipass)) iopt = ipass

      istat = 0
      do while (istat .eq. 0)

         call loadVTK(vtp, trim(fname), istat)
         if (istat .lt. 0) exit

         call getVTK_numPoints(vtp, nNo, istat)
         if (istat .lt. 0) exit

         call getVTK_numElems(vtp, nEl, istat)
         if (istat .lt. 0) exit

         call getVTK_nodesPerElem(vtp, eNoN, istat)
         if (istat .lt. 0) exit

         if (iopt .eq. 1) then
            lFa%nNo = nNo
            lFa%nEl = nEl
            lFa%eNon = eNoN
            allocate(lFa%ien(eNoN,nEl), lFa%x(nsd,nNo), lFa%dx(nsd,nNo), &
               lFa%wss(nsd,nNo), lFa%iblank(nEl))
            lFa%iblank = 1
            call selecteleb(lFa)
         else
            if (nNo.ne.lFa%nNo .or. nEl.ne.lFa%nEl .or. &
               eNoN.ne.lFa%eNoN) then
               write(stdout,ftab4) "Error: inconsistent degrees of "// &
                  "freedom"
               istat=-1; exit
            end if
         end if

         call getVTK_pointCoords(vtp, lFa%x, istat)
         if (istat .lt. 0) exit

         call getVTK_elemIEN(vtp, lFa%ien, istat)
         if (istat .lt. 0) exit

         call getVTK_pointData(vtp, 'FS_WSS', lFa%wss, istat)
         if (istat .lt. 0) exit

         call getVTK_pointData(vtp, 'MS_Displacement', lFa%dx, istat)
         if (istat .lt. 0) then
            lFa%dx = 0D0
            istat = 0
         end if

         call flushVTK(vtp)
         exit
      end do

      if (istat .lt. 0) then
         write(stdout,ftab4) "Error: vtk file read error"
         STOP
      end if



      return
      end subroutine readVTP

!**************************************************
!     Write data to a vtp file
      subroutine writeVTP(lFa, fname)
      use variables
      implicit none
      type(faceType), intent(in) :: lFa
      character(len=strL), intent(in) :: fname

      type(vtkXMLType) :: vtp
      integer :: istat, vtkType

      if (nsd .EQ. 3) then
         select case (lFa%eNoN)
         case(3)
            vtkType = 5   ! Tri !
         case(4)
            vtkType = 9   ! Quad !
         end select
      else
         select case (lFa%eNoN)
         case(2)
            vtkType = 3   ! Line !
         case(3)
            vtkType = 21   ! Quadr-Edge !
         end select
      end if

      istat = 0
      do while (istat .eq. 0)
         call vtkInitWriter(vtp, fname, istat)
         if (istat .lt. 0) exit

         call putVTK_pointCoords(vtp, lFa%x, istat)
         if (istat .lt. 0) exit

         call putVTK_elemIEN(vtp, lFa%IEN, vtkType, istat)
         if (istat .lt. 0) exit

         call putVTK_pointData(vtp,"FS_OSI", osi, istat)
         if (istat .lt. 0) exit

         call putVTK_pointData(vtp,"FS_AWSS", awss, istat)
         if (istat .lt. 0) exit

         call putVTK_pointData(vtp,"FS_TAWSS", tawss, istat)
         if (istat .lt. 0) exit

         call putVTK_pointData(vtp,"MS_Displacement", lFa%dx, istat)
         if (istat .lt. 0) exit

         call vtkWriteToFile(vtp, istat)
         if (istat .lt. 0) exit

         call flushVTK(vtp)
         exit
      end do

      if (istat .lt. 0) then
         write(stdout,ftab4) "Error: vtk write error"
         STOP
      end if

      return
      end subroutine writeVTP

!**************************************************
!     Compute element normals for a given face
      subroutine getNormals(lFa)
      use variables
      implicit none
      type(faceType), intent(inout) :: lFa

      integer :: a, e, g, Ac
      real(kind=8) :: xl(nsd,lFa%eNoN), n(nsd), Jac, sA

      if (allocated(lfa%nV)) then
         deallocate(lFa%nV)
         allocate(lFa%nV(nsd,lFa%nNo))
      else
         allocate(lFa%nV(nsd,lFa%nNo))
      end if

      lFa%nV = 0D0
      lFa%area = 0D0
      do e=1, lFa%nEl
         do a=1, lFa%eNoN
            Ac = lFa%IEN(a,e) + 1
            xl(:,a) = lFa%x(:,Ac) + lFa%dx(:,Ac)
         end do

         do g=1, lFa%nG
            call GNNB(lFa, e, g, n)
            Jac = SQRT(SUM( n(:)**2 ))
            sA = 0D0
            do a=1, lFa%eNoN
               Ac = lFa%IEN(a,e) + 1
               lFa%nV(:,Ac) = lFa%nV(:,Ac) + n*lFa%N(a,g)*lFa%w(g)
               sA = sA + lFa%N(a,g)
            end do
            if (lFa%iblank(e) .EQ. 1) &
               lFa%area = lFa%area + Jac*lFa%w(g)*sA
         end do
      end do
      write(stdout,ftab3) "Face area: "//STR(lFa%area)

      do a=1, lFa%nNo
         Jac = SQRT(SUM( lFa%nV(:,a)**2 ))
         if (Jac .lt. 1D1*eps) cycle
         lFa%nV(:,a) = lFa%nV(:,a) / Jac
      end do

      return
      end subroutine getNormals

!**************************************************
!     Compute area-averaged WSS
      subroutine calcAAWSS(lFa, r)
      use variables
      implicit none
      type(faceType), intent(inout) :: lFa
      real(kind=8), intent(out) :: r

      integer :: a, e, g, Ac
      real(kind=8) :: s, n(nsd), Jac, sHat

      r = 0D0
      do e=1, lFa%nEl
         if (lFa%iblank(e) .ne. 1) cycle
         do g=1, lFa%nG
            call GNNB(lFa, e, g, n)
            Jac = SQRT(SUM( n(:)**2 ))
            sHat = 0D0
            do a=1, lFa%eNoN
               Ac = lFa%IEN(a,e) + 1
               s = SQRT(SUM( lFa%wss(:,Ac)**2 ))
               sHat = sHat + s*lFa%N(a,g)
            end do
            r = r + Jac*lFa%w(g)*sHat
         end do
      end do
      r = r/lFa%area

      return
      end subroutine calcAAWSS

!**************************************************
!     Choose element type and initialize shape function & Gauss points
      subroutine selecteleb(lFa)
      use variables
      implicit none
      type(faceType), intent(inout) :: lFa

      integer :: g

      lFa%nG = lFa%eNoN
      allocate(lFa%w(lFa%nG), lFa%xi(nsd-1,lFa%nG), &
         lFa%N(lFa%eNoN,lFa%nG), lFa%Nx(nsd-1,lFa%eNoN,lFa%nG))

      call getGP(nsd-1, lFa%nG, lFa%w, lFa%xi)

      do g=1, lFa%nG
         call getShpF(nsd-1, lFa%eNoN, lFa%xi(:,g), lFa%N(:,g), &
            lFa%Nx(:,:,g))
      end do

      return
      end subroutine selecteleb

!**************************************************
!     Get Gauss point
      subroutine getGP(insd, nG, w, xi)
      implicit none
      integer, intent(in) :: insd, nG
      real(kind=8), intent(out) :: w(nG), xi(insd,nG)

      real(kind=8) :: s, t

      w = 1D0/6D0
      s = 2D0/3D0
      t = 1D0/6D0
      xi(1,1) = s; xi(2,1) = t
      xi(1,2) = t; xi(2,2) = s
      xi(1,3) = t; xi(2,3) = t

      return
      end subroutine getGP

!**************************************************
!     Get shape functions for the face element
      subroutine getShpF(insd, eNoN, xi, N, Nxi)
      implicit none
      integer, intent(in) :: insd, eNoN
      real(kind=8), intent(out) :: xi(insd), N(eNoN), Nxi(insd,eNoN)

      N(1) = xi(1)
      N(2) = xi(2)
      N(3) = 1D0 - xi(1) - xi(2)

      Nxi(1,1) =  1D0
      Nxi(2,1) =  0D0
      Nxi(1,2) =  0D0
      Nxi(2,2) =  1D0
      Nxi(1,3) = -1D0
      Nxi(2,3) = -1D0

      return
      end subroutine getShpF

!**************************************************
!     Get element normal
      subroutine GNNB(lFa, e, g, n)
      use variables
      implicit none
      type(faceType), intent(in) :: lFa
      integer, intent(in) :: e, g
      real(kind=8), intent(out) :: n(nsd)

      integer :: a, i, Ac
      real(kind=8) :: xXi(nsd,nsd-1), xl(nsd,lFa%eNoN)

      do a=1, lFa%eNoN
         Ac = lFa%IEN(a,e) + 1
         xl(:,a) = lFa%x(:,Ac) + lFa%dx(:,Ac)
      end do

      xXi = 0D0
      do a=1, lFa%eNoN
         do i=1, nsd-1
            xXi(:,i) = xXi(:,i) + lFa%Nx(i,a,g)*xl(:,a)
         end do
      end do
      n = CROSS(xXi)

      return
      end subroutine GNNB

 !**************************************************

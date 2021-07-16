!**************************************************

      module params

      integer, parameter :: eType_NA = 100, eType_LIN = 101, &
         eType_TRI = 102, eType_TET = 103, eType_BIL = 104, &
         eType_QUD = 105, eType_BIQ = 106, eType_BRK = 107, &
         eType_NRB = 108, eType_WDG = 109

      end module params

!**************************************************

      module variables
      use genUtils
      use stdParams
      use vtkXMLMod
      use params
      implicit none

      interface destroy
         module procedure destroyMesh, destroyFace
      end interface destroy

      real(kind=8), parameter :: leps = 1D-15

      type mshType
         integer :: nNo
         integer :: nEl
         integer :: eNoN
         integer :: nG
         integer :: eType = eType_NA
         integer :: vtkType
         integer, allocatable :: IEN(:,:)
         real(kind=8), allocatable :: w(:)
         real(kind=8), allocatable :: xi(:,:)
         real(kind=8), allocatable :: N(:,:)
         real(kind=8), allocatable :: Nx(:,:,:)
         real(kind=8), allocatable :: x(:,:)
         real(kind=8), allocatable :: u(:,:)
         real(kind=8), allocatable :: dx(:,:)
         real(kind=8), allocatable :: wss(:,:)
         real(kind=8), allocatable :: wssg(:)
      end type mshType

      type, extends(mshType) :: faceType
         integer, allocatable :: gE(:)
         integer, allocatable :: gN(:)
         real(kind=8), allocatable :: nV(:,:)
      end type faceType

      integer :: nsd
      integer :: nstart, nend, nfreq
      character(len=strL) :: froi, febc

      integer, allocatable :: iblank(:), incNds(:)
      real(kind=8), allocatable :: awss(:,:), tawss(:), osi(:), &
         mwss(:,:), tawssg(:)

      type(mshType)  :: msh
      type(faceType) :: fa

      contains
!-------------------------------------------------
         subroutine destroyMesh(lM)
         implicit none
         type(mshType), intent(inout) :: lM

         if (allocated(lM%IEN))  deallocate(lM%IEN)
         if (allocated(lM%w))    deallocate(lM%w)
         if (allocated(lM%xi))   deallocate(lM%xi)
         if (allocated(lM%N))    deallocate(lM%N)
         if (allocated(lM%Nx))   deallocate(lM%Nx)
         if (allocated(lM%x))    deallocate(lM%x)
         if (allocated(lM%u))    deallocate(lM%u)
         if (allocated(lM%dx))   deallocate(lM%dx)
         if (allocated(lM%wss))  deallocate(lM%wss)
         if (allocated(lM%wssg)) deallocate(lM%wssg)

         return
         end subroutine destroyMesh
!-------------------------------------------------
         subroutine destroyFace(lFa)
         implicit none
         type(faceType), intent(inout) :: lFa

         if (allocated(lFa%IEN))  deallocate(lFa%IEN)
         if (allocated(lFa%w))    deallocate(lfa%w)
         if (allocated(lFa%xi))   deallocate(lFa%xi)
         if (allocated(lFa%N))    deallocate(lFa%N)
         if (allocated(lFa%Nx))   deallocate(lFa%Nx)
         if (allocated(lFa%x))    deallocate(lFa%x)
         if (allocated(lFa%u))    deallocate(lFa%u)
         if (allocated(lFa%dx))   deallocate(lFa%dx)
         if (allocated(lFa%wss))  deallocate(lFa%wss)
         if (allocated(lFa%wssg)) deallocate(lFa%wssg)

         return
         end subroutine destroyFace
!-------------------------------------------------
      end module variables

!**************************************************

      program wss_post_2d
      use variables
      implicit none

      integer :: i, a, Ac, fid, ntime, iaawss, nNo
      integer :: iprobe
      real(kind=8) :: rtmp, aawss, aawssg, dt, time
      character(len=strL) :: fname
      character(len=strL) :: prefix
      character(len=strL) :: fhdr

      fid = 10
      i = iargc()
      if (i .ne. 1) then
         write(stdout,ftab4) "Error: at least one input argument needed"
         STOP
      end if
      call getarg(1,fname)

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
      read(fid,*) iprobe
      read(fid,*)
      read(fid,'(A)') froi
      read(fid,*)
      read(fid,'(A)') febc
      read(fid,*)
      read(fid,'(A)') fhdr
      close(fid)

      ntime = nstart
      if (ntime .lt. 100) then
         write(fname,'(A,I3.3,A)') trim(prefix), ntime, ".vtu"
      else
         write(fname,'(A)') trim(prefix)//trim(STR(ntime))//".vtu"
      end if

      call readVTU(msh, fname)
      call selectele(msh)
      call loadFace(msh, fa)
      nNo = fa%nNo

      if (iaawss .gt. 0) then
         write(fname,'(A)') "aawss_"//trim(STR(nstart))//"_"// &
         trim(STR(nend))//".dat"
         open(101,file=trim(fname))
         write(101,'(A)') 'Variables=t, AWSS, AAWSSG'
      end if

      if (iprobe.gt.0 .and. iprobe.lt.msh%nNo) then
         write(fname,'(A)') "probe_N"//trim(STR(iprobe))//"_"// &
            trim(STR(nstart))//"_"//trim(STR(nend))//".dat"
         open(102,file=trim(fname))
         write(102,'(A)') 'Variables=t, u, v'
      end if

      write(stdout,'(A)') "############################################"
      write(stdout,ftab1)
      write(stdout,ftab1) "Computing WSS metrics for time range: "// &
         trim(STR(nstart))//" to "//trim(STR(nend))
      write(stdout,ftab1)

      allocate(awss(nsd,nNo), tawss(nNo), mwss(2,nNo), tawssg(nNo))
      awss   = 0D0
      tawss  = 0D0
      tawssg = 0D0
      mwss(1,:) = HUGE(rtmp)
      mwss(2,:) = TINY(rtmp)

      do ntime = nstart, nend, nfreq
         time  = DBLE(ntime)*dt
         if (ntime .lt. 100) then
            write(fname,'(A,I3.3,A)') trim(prefix), ntime, ".vtu"
         else
            write(fname,'(A)') trim(prefix)//trim(STR(ntime))//".vtu"
         end if

         write(stdout,ftab2) "Reading file "//trim(fname)// &
            ", nTime: "//trim(STR(ntime))
         call readVTU(msh, fname)

         do a=1, nNo
            Ac = fa%gN(a)
            fa%wss(:,a) = msh%wss(:,Ac)
         end do
         call calcWSSG(msh, fa)

         aawss  = 0D0
         aawssg = 0D0
         do a=1, nNo
            awss(:,a) = awss(:,a) + fa%wss(:,a)
            rtmp  = sqrt(sum(fa%wss(:,a)**2))
            rtmp  = rtmp*real(1-iblank(a),kind=8) + &
                     leps*real(iblank(a),kind=8)
            aawss = aawss + rtmp
            tawss(a) = tawss(a) + rtmp
            if (rtmp .lt. mwss(1,a)) mwss(1,a) = rtmp
            if (rtmp .gt. mwss(2,a)) mwss(2,a) = rtmp

            rtmp  = fa%wssg(a)*real(1-iblank(a),kind=8) + &
               leps*real(iblank(a),kind=8)
            aawssg    = aawssg + rtmp
            tawssg(a) = tawssg(a) + rtmp
         end do
         rtmp   = real(nNo-SUM(iblank), kind=8)
         aawss  = aawss / rtmp
         aawssg = aawssg / rtmp

         if (iaawss .gt. 0) then
            write(101,'(2X,F9.4,2(2X,1pE13.6))') time, aawss, aawssg
         end if

         if (iprobe.gt.0 .and. iprobe.lt.msh%nNo) then
            write(102,'(2X,F9.4,2(2X,1pE13.6))') time, msh%u(1,iprobe), &
               msh%u(2,iprobe)
         end if
      end do
      if (iaawss .gt. 0) close(101)
      if (iprobe.gt.0 .and. iprobe.lt.msh%nNo) close(102)

      write(stdout,'(A)')
      i = (nend - nstart)/nfreq + 1
      awss(:,:) = awss(:,:) /real(i,kind=8)
      tawss(:)  = tawss(:)  /real(i,kind=8)
      tawssg(:) = tawssg(:) /real(i,kind=8)

      allocate(osi(nNo))
      osi = 0D0
      do a=1, nNo
         if (iblank(a) .ne. 0) cycle
         rtmp   = sqrt(sum(awss(:,a)**2))
         osi(a) = (1D0 - rtmp/tawss(a))/2D0
      end do

      write(fname,'(A)') trim(fhdr)//"_wss_"//trim(STR(nstart))//"_"// &
         trim(STR(nend))//".vtp"
      write(stdout,ftab1) "Writing WSS metrics to file "//trim(fname)
      call writeVTP(fa, fname)
      call destroy(fa)

      write(fname,'(A)') trim(fhdr)//"_wss_"//trim(STR(nstart))//"_"// &
         trim(STR(nend))//".dat"
      open(fid,file=trim(fname))
      write(fid,'(A)') "Variables=a, maxWSS, FWSS, TAWSS, OSI, TAWSSG"
      do Ac=1, msh%nNo
         a = incNds(Ac)
         if (a .eq. 0) cycle
         if (iblank(a) .eq. 0) then
            write(fid,'(2X,I4,3(2X,F12.6),2X,F6.4,2X,1pE13.6)') Ac, &
               mwss(2,a), mwss(2,a)-mwss(1,a), tawss(a), osi(a), tawssg(a)
         end if
      end do
      close(fid)

      write(stdout,ftab1)
      write(stdout,'(A)') "############################################"
      write(stdout,ftab1)

      end program wss_post_2d

!**************************************************

      subroutine readVTU(lM, fname)
      use variables
      implicit none
      type(mshType), intent(inout) :: lM
      character(len=strL), intent(in) :: fname

      type(vtkXMLType) :: vtu
      integer :: nNo, eNoN, nEl, istat

      real(kind=8), allocatable :: tmpX(:,:)

      istat = 0
      do while (istat .eq. 0)

         call loadVTK(vtu, trim(fname), istat)
         if (istat .lt. 0) exit

         call getVTK_numPoints(vtu, nNo, istat)
         if (istat .lt. 0) exit

         call getVTK_numElems(vtu, nEl, istat)
         if (istat .lt. 0) exit

         call getVTK_nodesPerElem(vtu, eNoN, istat)
         if (istat .lt. 0) exit

         if (.not.allocated(lM%x)) then
            lM%nNo = nNo
            lM%nEl = nEl
            lM%eNon = eNoN
            allocate(lM%ien(eNoN,nEl), lM%x(nsd,nNo), lM%u(nsd,nNo), &
               lM%dx(nsd,nNo), lM%wss(nsd,nNo), lM%wssg(nNo))
         end if

         allocate(tmpX(maxnsd,nNo))
         call getVTK_pointCoords(vtu, tmpX, istat)
         if (istat .lt. 0) exit
         lM%x(:,:) = tmpX(1:nsd,:)

         call getVTK_elemIEN(vtu, lM%IEN, istat)
         if (istat .lt. 0) exit

         call getVTK_pointData(vtu, 'FS_Velocity', tmpX, istat)
         if (istat .lt. 0) exit
         lM%u(:,:) = tmpX(1:nsd,:)

         call getVTK_pointData(vtu, 'FS_WSS', tmpX, istat)
         if (istat .lt. 0) exit
         lM%wss(:,:) = tmpX(1:nsd,:)
         lM%wssg(:)  = 0D0

         call getVTK_pointData(vtu, 'MS_Displacement', tmpX, istat)
         if (istat .lt. 0) then
            lM%dx = 0D0
            istat = 0
         else
            lM%dx(:,:) = tmpX(1:nsd,:)
         end if
         deallocate(tmpX)

         call flushVTK(vtu)
         exit
      end do

      if (istat .lt. 0) then
         write(stdout,ftab4) "Error: vtk file read error"
         STOP
      end if

      return
      end subroutine readVTU

!**************************************************

      subroutine loadFace(lM, lFa)
      use variables
      implicit none
      type(mshType), intent(in) :: lM
      type(faceType), intent(inout) :: lFa

      integer :: a, e, Ac, nNo

      lFa%nEl  = 0
      lFa%eNoN = lM%eNoN-1
      open(55, file=trim(febc))
      do
         read(55,*,end=110)
         lFa%nEl = lFa%nEl + 1
      end do
 110  allocate(lFa%IEN(lFa%eNoN,lFa%nEl), lFa%gE(lFa%nEl))
      rewind(55)
      do e=1, lFa%nEl
         read(55,*) lFa%gE(e), a, lFa%IEN(:,e)
      end do
      close(55)

      allocate(incNds(lM%nNo))
      incNds = 0
      nNo = 0
      do e=1, lFa%nEl
         do a=1, lFa%eNoN
            Ac = lFa%IEN(a,e)
            if (incNds(Ac) .eq. 0) then
               nNo = nNo + 1
               incNds(Ac) = nNo
            end if
        end do
      end do
      lFa%nNo = nNo

      allocate(lFa%gN(nNo))
      do Ac=1, lM%nNo
         if (incNds(Ac) .ne. 0) then
            a = incNds(Ac)
            lFa%gN(a) = Ac
         end if
      end do

!     Reorder IEN array to point local face nodes. IEN: 0 - (nEl-1)
      do e=1, lFa%nEl
         do a=1, lFa%eNoN
            Ac = lFa%IEN(a,e)
            Ac = incNds(Ac)
            lFa%IEN(a,e) = Ac - 1
         end do
      end do

      allocate(lFa%x(nsd,nNo), lFa%u(nsd,nNo), lFa%wss(nsd,nNo), &
         lFa%wssg(nNo), lFa%nV(nsd,nNo), iblank(nNo))
      do a=1, nNo
         Ac = lFa%gN(a)
         lFa%x(:,a)   = lM%x(:,Ac)
         lFa%u(:,a)   = lM%u(:,Ac)
         lFa%wss(:,a) = lM%wss(:,Ac)
         lFa%wssg(a)  = lM%wssg(Ac)
      end do

      iblank = 1
      open(55,file=trim(froi))
      do a=1, nNo
         read(55,*,end=111) Ac
         Ac = incNds(Ac)
         iblank(Ac) = 0
      end do
 111  close(55)

      call selecteleb(lFa)

      call calcFaceNormals(lM, lFa)

      return
      end subroutine loadFace

!**************************************************

      subroutine selectele(lM)
      use variables
      implicit none
      type(mshType), intent(inout) :: lM

      integer :: g

      if (nsd .eq. 3) then
         select case (lM%eNoN)
         case (8)
            lM%eType   = eType_BRK
            lM%nG      = 8
            lM%vtkType = 12
         case (6)
            lM%eType   = eType_WDG
            lM%nG      = 6
            lM%vtkType = 13
         case (4)
            lM%eType   = eType_TET
            lM%nG      = 4
            lM%vtkType = 10
         case default
            write(stdout,ftab4) &
               "Error: unknown combination of nsd & eNoN"
            STOP
         end select
      else
         select case (lM%eNoN)
         case (3)
            lM%eType   = eType_TRI
            lM%nG      = 3
            lM%vtkType = 5
         case (4)
            lM%eType   = eType_BIL
            lM%nG      = 4
            lM%vtkType = 9
         case (9)
            lM%eType   = eType_BIQ
            lM%nG      = 9
            lM%vtkType = 28
         case default
            write(stdout,ftab4) &
               "Error: unknown combination of nsd & eNoN"
            STOP
         end select
      end if

      allocate(lM%w(lM%nG), lM%xi(nsd,lM%nG), lM%N(lM%eNoN,lM%nG), &
         lM%Nx(nsd,lM%eNoN,lM%nG))

      call getGP(nsd, lM%eType, lM%nG, lM%w, lM%xi)

      do g=1, lM%nG
         call getShpF(nsd, lM%eType, lM%eNoN, lM%xi(:,g), lM%N(:,g), &
            lM%Nx(:,:,g))
      end do

      return
      end subroutine selectele

!**************************************************

      subroutine selecteleb(lFa)
      use variables
      implicit none
      type(faceType), intent(inout) :: lFa

      integer :: g, insd

      insd = nsd-1
      if (insd .eq. 2) then
         select case (lFa%eNoN)
         case (4)
            lFa%eType   = eType_BIL
            lFa%nG      = 4
         case (3)
            lFa%eType   = eType_TRI
            lFa%nG      = 3
         case default
            write(stdout,ftab4) &
               "Error: unknown combination of nsd & eNoN"
            STOP
         end select
      else if (insd .eq. 1) then
         select case (lFa%eNoN)
         case (2)
            lFa%eType   = eType_LIN
            lFa%nG      = 2
         case (3)
            lFa%eType   = eType_QUD
            lFa%nG      = 3
         case default
            write(stdout,ftab4) &
               "Error: unknown combination of nsd & eNoN"
            STOP
         end select
      end if

      allocate(lFa%w(lFa%nG), lFa%xi(insd,lFa%nG), &
         lFa%N(lFa%eNoN,lFa%nG), lFa%Nx(insd,lFa%eNoN,lFa%nG))

      call getGP(insd, lFa%eType, lFa%nG, lFa%w, lFa%xi)

      do g=1, lFa%nG
         call getShpF(insd, lFa%eType, lFa%eNoN, lFa%xi(:,g), &
            lFa%N(:,g), lFa%Nx(:,:,g))
      end do

      return
      end subroutine selecteleb

!**************************************************

      subroutine getGP(insd, eType, nG, w, xi)
      use params
      implicit none
      integer, intent(in) :: insd, eType, nG
      real(kind=8), intent(out) :: w(nG), xi(insd,nG)

      real(kind=8) s, t, lz, uz

!     3D elements
      select case (eType)
      case (eType_BRK)
         w = 1D0
         s =  1D0/SQRT(3D0)
         t = -1D0/SQRT(3D0)
         xi(1,1) = s; xi(2,1) = s; xi(3,1) = s
         xi(1,2) = t; xi(2,2) = s; xi(3,2) = s
         xi(1,3) = t; xi(2,3) = s; xi(3,3) = t
         xi(1,4) = s; xi(2,4) = s; xi(3,4) = t
         xi(1,5) = s; xi(2,5) = t; xi(3,5) = s
         xi(1,6) = t; xi(2,6) = t; xi(3,6) = s
         xi(1,7) = t; xi(2,7) = t; xi(3,7) = t
         xi(1,8) = s; xi(2,8) = t; xi(3,8) = t
      case (eType_TET)
         w = 1D0/24D0
         s = (5D0 + 3D0*SQRT(5D0))/2D1
         t = (5D0 -     SQRT(5D0))/2D1
         xi(1,1) = s; xi(2,1) = t; xi(3,1) = t
         xi(1,2) = t; xi(2,2) = s; xi(3,2) = t
         xi(1,3) = t; xi(2,3) = t; xi(3,3) = s
         xi(1,4) = t; xi(2,4) = t; xi(3,4) = t
      case (eType_WDG)
         w  =  1D0/6D0
         s  =  2D0/3D0
         t  =  1D0/6D0
         uz =  1D0/SQRT(3D0)
         lz = -1D0/SQRT(3D0)
         xi(1,1) = s; xi(2,1) = t; xi(3,1) = lz
         xi(1,2) = t; xi(2,2) = s; xi(3,2) = lz
         xi(1,3) = t; xi(2,3) = t; xi(3,3) = lz
         xi(1,4) = s; xi(2,4) = t; xi(3,4) = uz
         xi(1,5) = t; xi(2,5) = s; xi(3,5) = uz
         xi(1,6) = t; xi(2,6) = t; xi(3,6) = uz

!     2D elements
      case (eType_TRI)
         w = 1D0/6D0
         s = 2D0/3D0
         t = 1D0/6D0
         xi(1,1) = s; xi(2,1) = t
         xi(1,2) = t; xi(2,2) = s
         xi(1,3) = t; xi(2,3) = t
      case (eType_BIL)
         w = 1D0
         s =  1D0/SQRT(3D0)
         t = -1D0/SQRT(3D0)
         xi(1,1) = s; xi(2,1) = s
         xi(1,2) = t; xi(2,2) = s
         xi(1,3) = t; xi(2,3) = t
         xi(1,4) = s; xi(2,4) = t
      case (eType_BIQ)
         w(1) = 25D0/81D0; w(2) = 25D0/81D0; w(3) = 25D0/81D0
         w(4) = 25D0/81D0; w(5) = 40D0/81D0; w(6) = 40D0/81D0
         w(7) = 40D0/81D0; w(8) = 40D0/81D0; w(9) = 64D0/81D0
         s    = SQRT(6D-1)
         xi(1,1) =  -s; xi(2,1) =  -s
         xi(1,2) =   s; xi(2,2) =  -s
         xi(1,3) =   s; xi(2,3) =   s
         xi(1,4) =  -s; xi(2,4) =   s
         xi(1,5) = 0D0; xi(2,5) =  -s
         xi(1,6) =   s; xi(2,6) = 0D0
         xi(1,7) = 0D0; xi(2,7) =   s
         xi(1,8) =  -s; xi(2,8) = 0D0
         xi(1,9) = 0D0; xi(2,9) = 0D0

!     1D elements
      case (eType_LIN)
         w = 1D0
         s = 1D0/SQRT(3D0)
         xi(1,1) = -s
         xi(1,2) =  s
      case (eType_QUD)
         w(1) = 5D0/9D0; w(2) = 5D0/9D0; w(3) = 8D0/9D0
         s = SQRT(6D-1)
         xi(1,1) = -s
         xi(1,2) =  s
         xi(1,3) = 0D0
      end select

      return
      end subroutine getGP

!**************************************************

      subroutine getShpF(insd, eType, eNoN, xi, N, Nxi)
      use params
      implicit none
      integer, intent(in) :: insd, eType, eNoN
      real(kind=8), intent(out) :: xi(insd), N(eNoN), Nxi(insd,eNoN)

      real(kind=8) :: s, t, mx, my, ux, uy, uz, lx, ly, lz

!     3D elements
      select case (eType)
      case (eType_BRK)
         ux = 1D0 + xi(1); lx = 1D0 - xi(1)
         uy = 1D0 + xi(2); ly = 1D0 - xi(2)
         uz = 1D0 + xi(3); lz = 1D0 - xi(3)
         N(1) = ux*uy*uz/8D0
         N(2) = lx*uy*uz/8D0
         N(3) = lx*uy*lz/8D0
         N(4) = ux*uy*lz/8D0
         N(5) = ux*ly*uz/8D0
         N(6) = lx*ly*uz/8D0
         N(7) = lx*ly*lz/8D0
         N(8) = ux*ly*lz/8D0

         Nxi(1,1) =  uy*uz/8D0
         Nxi(2,1) =  ux*uz/8D0
         Nxi(3,1) =  ux*uy/8D0
         Nxi(1,2) = -uy*uz/8D0
         Nxi(2,2) =  lx*uz/8D0
         Nxi(3,2) =  lx*uy/8D0
         Nxi(1,3) = -uy*lz/8D0
         Nxi(2,3) =  lx*lz/8D0
         Nxi(3,3) = -lx*uy/8D0
         Nxi(1,4) =  uy*lz/8D0
         Nxi(2,4) =  ux*lz/8D0
         Nxi(3,4) = -ux*uy/8D0
         Nxi(1,5) =  ly*uz/8D0
         Nxi(2,5) = -ux*uz/8D0
         Nxi(3,5) =  ux*ly/8D0
         Nxi(1,6) = -ly*uz/8D0
         Nxi(2,6) = -lx*uz/8D0
         Nxi(3,6) =  lx*ly/8D0
         Nxi(1,7) = -ly*lz/8D0
         Nxi(2,7) = -lx*lz/8D0
         Nxi(3,7) = -lx*ly/8D0
         Nxi(1,8) =  ly*lz/8D0
         Nxi(2,8) = -ux*lz/8D0
         Nxi(3,8) = -ux*ly/8D0

      case (eType_TET)
         N(1) = xi(1)
         N(2) = xi(2)
         N(3) = xi(3)
         N(4) = 1D0 - xi(1) - xi(2) - xi(3)

         Nxi(1,1) =  1D0
         Nxi(2,1) =  0D0
         Nxi(3,1) =  0D0
         Nxi(1,2) =  0D0
         Nxi(2,2) =  1D0
         Nxi(3,2) =  0D0
         Nxi(1,3) =  0D0
         Nxi(2,3) =  0D0
         Nxi(3,3) =  1D0
         Nxi(1,4) = -1D0
         Nxi(2,4) = -1D0
         Nxi(3,4) = -1D0

      case (eType_WDG)
         ux = xi(1) ; uy = xi(2) ; uz = 1D0 - ux - uy
         s = (1D0 + xi(3))/2D0; t = (1D0 - xi(3))/2D0
         N(1) = ux*t
         N(2) = uy*t
         N(3) = uz*t
         N(4) = ux*s
         N(5) = uy*s
         N(6) = uz*s

         Nxi(1,1) =  t
         Nxi(2,1) =  0D0
         Nxi(3,1) = -ux/2D0
         Nxi(1,2) =  0D0
         Nxi(2,2) =  t
         Nxi(3,2) = -uy/2D0
         Nxi(1,3) = -t
         Nxi(2,3) = -t
         Nxi(3,3) = -uz/2D0
         Nxi(1,4) =  s
         Nxi(2,4) =  0D0
         Nxi(3,4) =  ux/2D0
         Nxi(1,5) =  0D0
         Nxi(2,5) =  s
         Nxi(3,5) =  uy/2D0
         Nxi(1,6) = -s
         Nxi(2,6) = -s
         Nxi(3,6) =  uz/2D0

!     2D elements
      case (eType_TRI)
         N(1) = xi(1)
         N(2) = xi(2)
         N(3) = 1D0 - xi(1) - xi(2)

         Nxi(1,1) =  1D0
         Nxi(2,1) =  0D0
         Nxi(1,2) =  0D0
         Nxi(2,2) =  1D0
         Nxi(1,3) = -1D0
         Nxi(2,3) = -1D0

      case (eType_BIL)
         ux = 1D0 + xi(1); lx = 1D0 - xi(1)
         uy = 1D0 + xi(2); ly = 1D0 - xi(2)
         N(1) = ux*uy/4D0
         N(2) = lx*uy/4D0
         N(3) = lx*ly/4D0
         N(4) = ux*ly/4D0

         Nxi(1,1) =  uy/4D0
         Nxi(2,1) =  ux/4D0
         Nxi(1,2) = -uy/4D0
         Nxi(2,2) =  lx/4D0
         Nxi(1,3) = -ly/4D0
         Nxi(2,3) = -lx/4D0
         Nxi(1,4) =  ly/4D0
         Nxi(2,4) = -ux/4D0

      case (eType_BIQ)
         ux = 1D0 + xi(1); mx = xi(1); lx = 1D0 - xi(1)
         uy = 1D0 + xi(2); my = xi(2); ly = 1D0 - xi(2)
         N(1) =  mx*lx*my*ly/4D0
         N(2) = -mx*ux*my*ly/4D0
         N(3) =  mx*ux*my*uy/4D0
         N(4) = -mx*lx*my*uy/4D0
         N(5) = -lx*ux*my*ly/2D0
         N(6) =  mx*ux*ly*uy/2D0
         N(7) =  lx*ux*my*uy/2D0
         N(8) = -mx*lx*ly*uy/2D0
         N(9) =  lx*ux*ly*uy

         Nxi(1,1) =  (lx - mx)*my*ly/4D0
         Nxi(2,1) =  (ly - my)*mx*lx/4D0
         Nxi(1,2) = -(ux + mx)*my*ly/4D0
         Nxi(2,2) = -(ly - my)*mx*ux/4D0
         Nxi(1,3) =  (ux + mx)*my*uy/4D0
         Nxi(2,3) =  (uy + my)*mx*ux/4D0
         Nxi(1,4) = -(lx - mx)*my*uy/4D0
         Nxi(2,4) = -(uy + my)*mx*lx/4D0
         Nxi(1,5) = -(lx - ux)*my*ly/2D0
         Nxi(2,5) = -(ly - my)*lx*ux/2D0
         Nxi(1,6) =  (ux + mx)*ly*uy/2D0
         Nxi(2,6) =  (ly - uy)*mx*ux/2D0
         Nxi(1,7) =  (lx - ux)*my*uy/2D0
         Nxi(2,7) =  (uy + my)*lx*ux/2D0
         Nxi(1,8) = -(lx - mx)*ly*uy/2D0
         Nxi(2,8) = -(ly - uy)*mx*lx/2D0
         Nxi(1,9) =  (lx - ux)*ly*uy
         Nxi(2,9) =  (ly - uy)*lx*ux

!     1D elements
      case (eType_LIN)
         N(1) = (1D0 - xi(1))/2D0
         N(2) = (1D0 + xi(1))/2D0

         Nxi(1,1) = -5D-1
         Nxi(1,2) =  5D-1
      case (eType_QUD)
         N(1) = -xi(1)*(1D0 - xi(1))/2D0
         N(2) =  xi(1)*(1D0 + xi(1))/2D0
         N(3) = (1D0 - xi(1))*(1D0 + xi(1))

         Nxi(1,1) = -5D-1 + xi(1)
         Nxi(1,2) =  5D-1 + xi(1)
         Nxi(1,3) = -2D0*xi(1)
      end select

      return
      end subroutine getShpF

!**************************************************

      subroutine calcFaceNormals(lM, lFa)
      use variables
      implicit none
      type(mshType), intent(in) :: lM
      type(faceType), intent(inout) :: lFa

      integer :: e, a, Ac, g
      real(kind=8) :: rtmp, n(nsd)
      real(kind=8), allocatable :: sV(:,:)

      allocate(sV(nsd,lM%nNo))
      sV(:,:) = 0D0
      do e=1, lFa%nEl
         do g=1, lFa%nG
            call GNNB(lFa, lM, e, g, n)
            do a=1, lFa%eNoN
               Ac = lFa%IEN(a,e) + 1
               Ac = lFa%gN(Ac)
               sV(:,Ac) = sV(:,Ac) + lFa%N(a,g)*lFa%w(g)*n(:)
            end do
         end do
      end do

      lFa%nV = 0D0
      do a=1, lFa%nNo
         Ac = lFa%gN(a)
         rtmp = SQRT(SUM(sV(:,Ac)**2))
         rtmp = MAX(rtmp, leps)
         lFa%nV(:,a) = sV(:,Ac)/rtmp

         rtmp = SQRT(SUM(lFa%nV(:,a)**2))
      end do

      return
      end subroutine calcFaceNormals

!**************************************************

      subroutine GNN(eNoN, Nxi, x, Nx, Jac, ks)
      use variables, only: nsd
      implicit none

      integer, intent(in) :: eNoN
      real(kind=8), intent(in) :: Nxi(nsd,eNoN), x(nsd,eNoN)
      real(kind=8), intent(out) :: Nx(nsd,eNoN), Jac, ks(nsd,nsd)

      integer :: a
      real(kind=8) :: xXi(nsd,nsd), xiX(nsd,nsd)

      Nx  = 0D0
      xXi = 0D0
      if (nsd .eq. 2) then
         do a=1, eNoN
            xXi(:,1) = xXi(:,1) + x(:,a)*Nxi(1,a)
            xXi(:,2) = xXi(:,2) + x(:,a)*Nxi(2,a)
         end do

         Jac = xXi(1,1)*xXi(2,2) - xXi(1,2)*xXi(2,1)

         xiX(1,1) =  xXi(2,2)/Jac
         xiX(1,2) = -xXi(1,2)/Jac
         xiX(2,1) = -xXi(2,1)/Jac
         xiX(2,2) =  xXi(1,1)/Jac

         ks(1,1) = xiX(1,1)*xiX(1,1) + xiX(2,1)*xiX(2,1)
         ks(1,2) = xiX(1,1)*xiX(1,2) + xiX(2,1)*xiX(2,2)
         ks(2,2) = xiX(1,2)*xiX(1,2) + xiX(2,2)*xiX(2,2)
         ks(2,1) = ks(1,2)

         do a=1, eNoN
            Nx(1,a) = Nx(1,a)+ Nxi(1,a)*xiX(1,1) + Nxi(2,a)*xiX(2,1)
            Nx(2,a) = Nx(2,a)+ Nxi(1,a)*xiX(1,2) + Nxi(2,a)*xiX(2,2)
         end do
      else
         do a=1, eNoN
            xXi(:,1) = xXi(:,1) + x(:,a)*Nxi(1,a)
            xXi(:,2) = xXi(:,2) + x(:,a)*Nxi(2,a)
            xXi(:,3) = xXi(:,3) + x(:,a)*Nxi(3,a)
         end do

         Jac = xXi(1,1)*xXi(2,2)*xXi(3,3) + &
               xXi(1,2)*xXi(2,3)*xXi(3,1) + &
               xXi(1,3)*xXi(2,1)*xXi(3,2) - &
               xXi(1,1)*xXi(2,3)*xXi(3,2) - &
               xXi(1,2)*xXi(2,1)*xXi(3,3) - &
               xXi(1,3)*xXi(2,2)*xXi(3,1)

         xiX(1,1) = (xXi(2,2)*xXi(3,3) - xXi(2,3)*xXi(3,2))/Jac
         xiX(1,2) = (xXi(3,2)*xXi(1,3) - xXi(3,3)*xXi(1,2))/Jac
         xiX(1,3) = (xXi(1,2)*xXi(2,3) - xXi(1,3)*xXi(2,2))/Jac
         xiX(2,1) = (xXi(2,3)*xXi(3,1) - xXi(2,1)*xXi(3,3))/Jac
         xiX(2,2) = (xXi(3,3)*xXi(1,1) - xXi(3,1)*xXi(1,3))/Jac
         xiX(2,3) = (xXi(1,3)*xXi(2,1) - xXi(1,1)*xXi(2,3))/Jac
         xiX(3,1) = (xXi(2,1)*xXi(3,2) - xXi(2,2)*xXi(3,1))/Jac
         xiX(3,2) = (xXi(3,1)*xXi(1,2) - xXi(3,2)*xXi(1,1))/Jac
         xiX(3,3) = (xXi(1,1)*xXi(2,2) - xXi(1,2)*xXi(2,1))/Jac

         ks(1,1) = xiX(1,1)*xiX(1,1)+xiX(2,1)*xiX(2,1)+xiX(3,1)*xiX(3,1)
         ks(1,2) = xiX(1,2)*xiX(1,1)+xiX(2,2)*xiX(2,1)+xiX(3,2)*xiX(3,1)
         ks(1,3) = xiX(1,3)*xiX(1,1)+xiX(2,3)*xiX(2,1)+xiX(3,3)*xiX(3,1)
         ks(2,2) = xiX(1,2)*xiX(1,2)+xiX(2,2)*xiX(2,2)+xiX(3,2)*xiX(3,2)
         ks(2,3) = xiX(1,2)*xiX(1,3)+xiX(2,2)*xiX(2,3)+xiX(3,2)*xiX(3,3)
         ks(3,3) = xiX(1,3)*xiX(1,3)+xiX(2,3)*xiX(2,3)+xiX(3,3)*xiX(3,3)
         ks(2,1) = ks(1,2)
         ks(3,1) = ks(1,3)
         ks(3,2) = ks(2,3)

         do a=1, eNoN
            Nx(1,a) = Nx(1,a) + Nxi(1,a)*xiX(1,1) + &
                                Nxi(2,a)*xiX(2,1) + &
                                Nxi(3,a)*xiX(3,1)

            Nx(2,a) = Nx(2,a) + Nxi(1,a)*xiX(1,2) + &
                                Nxi(2,a)*xiX(2,2) + &
                                Nxi(3,a)*xiX(3,2)

            Nx(3,a) = Nx(3,a) + Nxi(1,a)*xiX(1,3) + &
                                Nxi(2,a)*xiX(2,3) + &
                                Nxi(3,a)*xiX(3,3)
         end do
      end if

      return
      end subroutine GNN

!**************************************************

      subroutine GNNB(lFa, lM, e, g, n)
      use variables
      implicit none
      integer, intent(in) :: e, g
      real(kind=8), intent(out) :: n(nsd)
      type(faceType), intent(in) :: lFa
      type(mshType), intent(in) :: lM

      integer :: i, a, Ac, b, Bc, Ec, eNoN, insd
      real(kind=8) :: v(nsd), rtmp

      logical, allocatable :: setIt(:)
      integer, allocatable :: ptr(:)
      real(kind=8), allocatable :: lX(:,:), xXi(:,:)

      Ec   = lFa%gE(e)
      eNoN = lM%eNoN
      insd = nsd - 1

      allocate(lX(nsd,eNoN), ptr(eNoN), setIt(eNoN))
      setIt = .TRUE.
      do a=1, lFa%eNoN
         Ac = lFa%IEN(a,e) + 1
         Ac = lFa%gN(Ac)
         do b=1, eNoN
            if (setIt(b)) then
               Bc = lM%IEN(b,Ec) + 1
               if (Bc .EQ. Ac) exit
            end if
         end do
         ptr(a)   = b
         setIt(b) = .FALSE.
      end do
      a = lFa%eNoN
      do b=1, eNoN
         if (setIt(b)) then
            a      = a + 1
            ptr(a) = b
         end if
      end do

!     Correcting the position vector if mesh is moving
      do a=1, eNoN
         Ac = lM%IEN(a,Ec) + 1
         lX(:,a) = lM%x(:,Ac) + lM%dx(:,Ac)
      end do

!     Calculating surface deflation
      allocate(xXi(nsd,insd))
      xXi = 0D0
      do a=1, lFa%eNoN
         b = ptr(a)
         do i=1, insd
            xXi(:,i) = xXi(:,i) + lFa%Nx(i,a,g)*lX(:,b)
         end do
      end do

      if (nsd .eq. 2) then
         n(1) =  xXi(2,1)
         n(2) = -xXi(1,1)
      else if (nsd .eq. 3) then
         n(1) = xXi(2,1)*xXi(3,2) - xXi(3,1)*xXi(2,2)
         n(2) = xXi(3,1)*xXi(1,2) - xXi(1,1)*xXi(3,2)
         n(3) = xXi(1,1)*xXi(2,2) - xXi(2,1)*xXi(1,2)
      end if

!     Changing the sign if neccessary. a locates on the face and b
!     outside of the face, in the parent element
      a = ptr(1)
      b = ptr(lFa%eNoN+1)
      v = lX(:,a) - lX(:,b)
      rtmp = 0D0
      do i=1, nsd
         rtmp = rtmp + n(i)*v(i)
      end do
      if (rtmp .lt. 0D0) n = -n

      deallocate(lX, xXi)

      return
      end subroutine GNNB

!**************************************************

      subroutine calcWSSG(lM, lFa)
      use variables
      implicit none
      type(mshType), intent(in) :: lM
      type(faceType), intent(inout) :: lFa

      integer :: i, j, a, e, g, Ac, Ec, eNoN
      real(kind=8) :: w, Jac, wssg(nsd,nsd), tmp(nsd,nsd), rtmp, &
         nV(nsd), Tdn(nsd), ndTdn
      real(kind=8), allocatable :: xl(:,:), wl(:,:), N(:), Nx(:,:), &
         gnV(:,:), lnV(:,:), sA(:), sF(:)

      eNoN  = lM%eNoN
      call calcFaceNormals(lM, lFa)

      allocate(gnV(nsd,lM%nNo), lnV(nsd,eNoN))
      gnV = 0D0
      do a=1, lFa%nNo
         Ac = lFa%gN(a)
         gnV(:,Ac) = lFa%nV(:,a)
      end do

      allocate(xl(nsd,eNoN), wl(nsd,eNoN), N(eNoN), Nx(nsd,eNoN), &
         sA(lM%nNo), sF(lM%nNo))
      sA = 0D0
      sF = 0D0
      do e=1, lFa%nEl
         Ec = lFa%gE(e)

         nV  = 0D0
         lnV = 0D0
         do a=1, eNoN
            Ac = lM%IEN(a,Ec) + 1
            lnV(:,a) = gnV(:,Ac)
            nV(:)    = nV(:) + lnV(:,a)
            xl(:,a)  = lM%x(:,Ac) + lM%dx(:,Ac)
            wl(:,a)  = lM%wss(:,Ac)
         end do
         nV(:) = nV(:)/real(lFa%eNoN,kind=8)
         do a=1, eNoN
            if (SQRT(SUM(lnV(:,a)**2)) .lt. leps) lnV(:,a) = nV(:)
         end do

         do g=1, lM%nG
            if (g .eq. 1) call GNN(eNoN, lM%Nx(:,:,g), xl, Nx, Jac, tmp)
            w    = lM%w(g)*Jac
            N(:) = lM%N(:,g)

            wssg = 0D0
            nV   = 0D0
            do a=1, eNoN
               nV(:) = nV(:) + N(a)*lnV(:,a)
               do i=1, nsd
                  do j=1, nsd
                     wssg(i,j) = wssg(i,j) + Nx(i,a)*wl(j,a)
                  end do
               end do
            end do

            Tdn   = 0D0
            ndTdn = 0D0
            do i=1, nsd
               do j=1, nsd
                  Tdn(i) = Tdn(i) + wssg(i,j)*nV(j)
               end do
               ndTdn = ndTdn + Tdn(i)*nV(i)
            end do
            Tdn(:) = Tdn(:) - ndTdn*nV(:)

            rtmp = SQRT(SUM(Tdn(:)**2))
            do a=1, eNoN
               Ac = lM%IEN(a,Ec) + 1
               sA(Ac) = sA(Ac) + w*N(a)
               sF(Ac) = sF(Ac) + w*N(a)*rtmp
            end do
         end do ! g
      end do ! e

      do a=1, lFa%nNo
         Ac = lFa%gN(a)
         lFa%wssg(a) = sF(Ac)/sA(Ac)
      end do

      return
      end subroutine calcWSSG

!**************************************************

      subroutine writeVTP(lFa, fname)
      use variables
      implicit none
      type(faceType), intent(in) :: lFa
      character(len=strL), intent(in) :: fname

      type(vtkXMLType) :: vtp
      integer :: istat, vtkType

      real(kind=8), allocatable :: tmpX(:)

      select case (lFa%eNoN)
      case(2)
         vtkType = 3   ! Line !
      case(3)
         vtkType = 21   ! Quadr-Edge !
      end select

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

         call putVTK_pointData(vtp,"FS_TAWSSG", tawssg, istat)
         if (istat .lt. 0) exit

         call putVTK_pointData(vtp,"FS_minWSS", mwss(1,:), istat)
         if (istat .lt. 0) exit

         call putVTK_pointData(vtp,"FS_maxWSS", mwss(2,:), istat)
         if (istat .lt. 0) exit

         call putVTK_pointData(vtp,"GlobalNodeID", lFa%gN, istat)
         if (istat .lt. 0) exit

         allocate(tmpX(lFa%nNo))
         tmpX = (mwss(2,:)-mwss(1,:))
         call putVTK_pointData(vtp,"FS_FWSS", tmpX, istat)
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


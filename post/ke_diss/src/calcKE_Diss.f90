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
      use vtkLegacyMod
      use params
      implicit none
      
      type meshType
         integer :: nNo
         integer :: nEl
         integer :: eNoN
         integer :: nG
         integer :: eType = eType_NA
         integer :: vtkType
         integer, allocatable :: IEN(:,:)
         real(kind=8) :: vol
         real(kind=8), allocatable :: w(:)
         real(kind=8), allocatable :: xi(:,:)
         real(kind=8), allocatable :: N(:,:)
         real(kind=8), allocatable :: Nx(:,:,:)
         real(kind=8), allocatable :: x(:,:)
         real(kind=8), allocatable :: dx(:,:)
         real(kind=8), allocatable :: u(:,:)
      end type meshType
      
      integer :: nsd
      integer :: nstart, nend, nfreq
      real(kind=8) :: ake, adiss, rho, mu
      
      real(kind=8), allocatable :: ke(:), phi(:)
      
      type(meshType) :: msh
      
      contains
!-------------------------------------------------
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
!-------------------------------------------------
         subroutine destroyMesh(lM)
         implicit none
         type(meshType), intent(inout) :: lM
         
         if (allocated(lM%IEN)) deallocate(lM%IEN)
         if (allocated(lM%w))   deallocate(lM%w)
         if (allocated(lM%xi))  deallocate(lM%xi)
         if (allocated(lM%N))   deallocate(lM%N)
         if (allocated(lM%Nx))  deallocate(lM%Nx)
         if (allocated(lM%x))   deallocate(lM%x)
         if (allocated(lM%dx))  deallocate(lM%dx)
         if (allocated(lM%u))   deallocate(lM%u)
         
         return
         end subroutine destroyMesh
!-------------------------------------------------
      end module variables
      
!**************************************************

      program ke_diss_calc
      use variables
      implicit none
      integer :: i, a, cnt, ntot, fid, ntime, ike, wtf
      real(kind=8) :: rtmp, dt, time
      character(len=strL) :: fname
      character(len=strL) :: prefix, postfix
      
      interface
         subroutine readVTU(lM, fname, ipass)
         use variables
         implicit none
         type(meshType), intent(inout) :: lM
         character(len=strL), intent(in) :: fname
         integer, intent(in), optional :: ipass
         end subroutine readVTU
         
         subroutine readVTK(lM, fname, ipass)
         use variables
         implicit none
         type(meshType), intent(inout) :: lM
         character(len=strL), intent(in) :: fname
         integer, intent(in), optional :: ipass
         end subroutine readVTK
      end interface
      
      fid = 10
      i = iargc()
      if (i .ne. 1) then
         write(stdout,ftab4) "Error: at least one input argument needed"
         STOP
      end if
      call getarg(1,fname)
      
      wtf = 0
      ike = 0
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
      read(fid,'(A)') postfix
      read(fid,*)
      read(fid,*) ike, wtf
      read(fid,*)
      read(fid,*) rho, mu
      close(fid)
      
      if (nfreq .ne. 0) then
         ntot = (nend - nstart)/nfreq + 1
      else
         ntot = 1
      end if

      if (ike .gt. 0) then
         write(fname,'(A)') "ke_diss_"//trim(STR(nstart))//"_"// &
         trim(STR(nend))//".dat"
         open(101,file=trim(fname))
         write(101,'(A)') 'Variables=t, Vol, KE, DISS'
      end if
      
      write(stdout,'(A)') "############################################"
      write(stdout,ftab1)
      write(stdout,ftab1) "Computing KE and Dissipaton "// &
         "for time range: "//trim(STR(nstart))//" to "//trim(STR(nend))
      write(stdout,ftab1)
      
      do cnt=1, ntot
         ntime = nstart + (cnt-1)*nfreq
         time  = DBLE(ntime)*dt
         call getfName(fName, prefix, postfix, ntime)
         write(stdout,ftab2) "Reading file "//trim(fname)// &
            ", nTime: "//trim(STR(ntime))
         if (cnt .eq. 1) then
            select case (trim(postfix))
            case ("vtu")
               call readVTU(msh, fName, 1)
            case ("vtk")
               call readVTK(msh, fName, 1)
            end select
         else
            select case (trim(postfix))
            case ("vtu")
               call readVTU(msh, fName)
            case ("vtk")
               call readVTK(msh, fName)
            end select
         end if
         call calcKE_Diss(msh)
         write(101,'(2X,F9.4,3(2X,1pE18.6))') time, msh%vol, ake, adiss
         
         if (wtf .ne. 0) then
            write(fname,'(A)') "ke_diss_"//trim(STR(ntime))//".vtu"
            write(stdout,ftab2) "Writing to file "//trim(fname)
            call writeVTU(msh, fName)
         end if
      end do
      write(stdout,ftab1)
      write(stdout,'(A)') "############################################"
      close(101)
      
      end program ke_diss_calc
      
!**************************************************
      
      subroutine getfName(fName, prefix, postfix, ntime)
      use stdParams
      use genUtils
      implicit none
      character(len=strL), intent(in) :: prefix, postfix
      character(len=strL), intent(out) :: fName
      integer, intent(in) :: nTime
      
      if (nTime .lt. 1000) then
         write(fname,'(A,I3.3,A,A)') trim(prefix), ntime, ".", &
            trim(postfix)
      else
         write(fname,'(A)') trim(prefix)//trim(STR(ntime))//"."// &
            trim(postfix)
      end if
      
      return
      end subroutine getfName
      
!**************************************************
      
      subroutine readVTU(lM, fname, ipass)
      use variables
      implicit none
      type(meshType), intent(inout) :: lM
      character(len=strL), intent(in) :: fname
      integer, intent(in), optional :: ipass
      
      type(vtkXMLType) :: vtu
      integer :: istat, iopt
      integer :: e, nNo, eNoN, nEl
      real(kind=8), allocatable :: tmpX(:,:)
      
      iopt = 0
      if (present(ipass)) iopt = 1
      
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
         
         if (iopt .eq. 1) then
            lM%nNo = nNo
            lM%nEl = nEl
            lM%eNon = eNoN
            allocate(lM%ien(eNoN,nEl), lM%x(nsd,nNo), lM%dx(nsd,nNo), &
               lM%u(nsd,nNo))
            call selectele(lM)
         else
            if (nNo.ne.lM%nNo .or. nEl.ne.lM%nEl .or. &
               eNoN.ne.lM%eNoN) then
               call destroyMesh(lM)
               lM%nNo = nNo
               lM%nEl = nEl
               lM%eNon = eNoN
               allocate(lM%ien(eNoN,nEl), lM%x(nsd,nNo), &
                  lM%dx(nsd,nNo), lM%u(nsd,nNo))
               call selectele(lM)
            end if
         end if
         
         allocate(tmpX(maxNSD,nNo))
         call getVTK_pointCoords(vtu, tmpX, istat)
         if (istat .lt. 0) exit
         lM%x(:,:) = tmpX(1:nsd,:)
         
         call getVTK_elemIEN(vtu, lM%ien, istat)
         if (istat .lt. 0) exit
         lM%IEN(:,:) = lM%IEN(:,:) + 1
         
         call getVTK_pointData(vtu, 'FS_Velocity', tmpX, istat)
         if (istat .lt. 0) exit
         lM%u(:,:) = tmpX(1:nsd,:)
         
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
      
      subroutine readVTK(lM, fname, ipass)
      use variables
      implicit none
      type(meshType), intent(inout) :: lM
      character(len=strL), intent(in) :: fname
      integer, intent(in), optional :: ipass
      
      type(vtkUnstrucGridType) :: vtk
      integer :: istat, iopt
      integer :: e, nNo, eNoN, nEl
      integer :: i, j, k, ii
      real(kind=8), allocatable :: tmpX(:,:)
      
      iopt = 0
      if (present(ipass)) iopt = 1
      
      istat = 0
      do while (istat .eq. 0)
         
         call loadLegacyVTK(vtk, trim(fname), istat)
         if (istat .lt. 0) exit
         
         nNo  = vtk%nNo
         nEl  = vtk%nel
         eNoN = vtk%eNoN
         
         if (iopt .eq. 1) then
            lM%nNo = nNo
            lM%nEl = nEl
            lM%eNon = eNoN
            allocate(lM%ien(eNoN,nEl), lM%x(nsd,nNo), lM%dx(nsd,nNo), &
               lM%u(nsd,nNo))
            call selectele(lM)
         else
            if (nNo.ne.lM%nNo .or. nEl.ne.lM%nEl .or. &
               eNoN.ne.lM%eNoN) then
               call destroyMesh(lM)
               lM%nNo = nNo
               lM%nEl = nEl
               lM%eNon = eNoN
               allocate(lM%ien(eNoN,nEl), lM%x(nsd,nNo), &
                  lM%dx(nsd,nNo), lM%u(nsd,nNo))
               call selectele(lM)
            end if
         end if
         
         lM%x(:,:) = vtk%x(1:nsd,:)
         lM%IEN(:,:) = vtk%ien(:,:) + 1
         lM%dx = 0D0
         
         do k=1, vtk%ptData%nVar
            select case (trim(vtk%ptData%varName(k)))
            case ("FS_Velocity")
               allocate(tmpX(maxNSD,nNo))
               ii = vtk%ptData%ioff(k)
               do j=1, nNo
                  do i=1, maxNSD
                     ii = ii + 1
                     tmpX(i,j) = real(vtk%ptData%arr(ii), kind=8)
                  end do
               end do
               lM%u(:,:) = tmpX(1:nsd,:)
               deallocate(tmpX)
            case ("MS_Displacement")
               allocate(tmpX(maxNSD,nNo))
               ii = vtk%ptData%ioff(k)
               do j=1, nNo
                  do i=1, maxNSD
                     ii = ii + 1
                     tmpX(i,j) = real(vtk%ptData%arr(ii), kind=8)
                  end do
               end do
               lM%dx(:,:) = tmpX(1:nsd,:)
               deallocate(tmpX)
            end select
         end do
         
         call flushLegacyVTK(vtk)
         exit
      end do
      
      if (istat .lt. 0) then
         write(stdout,ftab4) "Error: vtk file read error"
         STOP
      end if
      
      return
      end subroutine readVTK
      
!**************************************************
      
      subroutine selectele(lM)
      use variables
      implicit none
      type(meshType), intent(inout) :: lM
      
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
      
      subroutine calcKE_Diss(lM)
      use variables
      implicit none
      type(meshType), intent(inout) :: lM
      
      integer :: a, e, Ac, i, g, eNoN, nEl, nNo
      real(kind=8) :: w, vol, Jac, tmp(nsd,nsd), rtmp, lKE, lPhi
      real(kind=8), allocatable :: xl(:,:), ul(:,:), sA(:), N(:), &
         Nx(:,:)
      
      nNo  = lM%nNo
      nEl  = lM%nEl
      eNoN = lM%eNoN
      
      if (allocated(ke)) deallocate(ke)
      if (allocated(phi)) deallocate(phi)
      allocate(ke(nNo), phi(nNo), sA(nNo))
      
      ke  = 0D0
      phi = 0D0
      do a=1, nNo
         rtmp = 0D0
         do i=1, nsd
            rtmp = rtmp + lM%u(i,a)**2
         end do
         ke(a) = 5D-1*rho*rtmp
      end do
      
      allocate(xl(nsd,eNoN), ul(nsd,eNoN), N(eNoN), Nx(nsd,eNoN))
      
      vol   = 0D0
      ake   = 0D0
      adiss = 0d0
      sA    = 0D0
      do e=1, nEl
         do a=1, eNoN
            Ac = lM%IEN(a,e)
            xl(:,a) = lM%x(:,Ac) + lM%dx(:,Ac)
            ul(:,a) = lM%u(:,Ac)
         end do
         
         do g=1, lM%nG
            if (g .eq. 1) then
               call GNN(eNoN, lM%Nx(:,:,g), xl, Nx, Jac, tmp)
            end if
            w    = lM%w(g)*Jac
            N(:) = lM%N(:,g)
            
            if (nsd .eq. 2) then
               call diss2D(eNoN, Nx, ul, rtmp)
            else
               call diss3D(eNoN, Nx, ul, rtmp)
            end if
            
            lKE  = 0D0
            lPhi = 0D0
            do a=1, eNoN
               Ac = lM%IEN(a,e)
               sA(Ac)  = sA(Ac)  + w*N(a)
               lKE  = lKE + w*N(a)*ke(Ac)
               lPhi = lPhi + w*N(a)*rtmp
               phi(Ac) = phi(Ac) + w*N(a)*rtmp
            end do
            vol   = vol + w
            ake   = ake + lKE
            adiss = adiss + lPhi
         end do
      end do
      ake    = ake / vol
      adiss  = adiss / vol
      lM%vol = vol
      
      do a=1, nNo
         phi(a) = phi(a)/sA(a)
      end do
      
      return
      end subroutine calcKE_Diss
      
!**************************************************

      subroutine diss2D(eNoN, Nx, ul, lPhi)
      use variables
      implicit none
      integer, intent(in) :: eNoN
      real(kind=8), intent(in) :: Nx(nsd,eNoN), ul(nsd,eNoN)
      real(kind=8), intent(inout) :: lPhi
      
      integer :: a
      real(kind=8) :: ux(nsd,nsd)
      
      ux = 0D0
      do a=1, eNoN
         ux(1,1) = ux(1,1) + Nx(1,a)*ul(1,a)
         ux(2,1) = ux(2,1) + Nx(2,a)*ul(1,a)
         ux(1,2) = ux(1,2) + Nx(1,a)*ul(2,a)
         ux(2,2) = ux(2,2) + Nx(2,a)*ul(2,a)
      end do
      
      lPhi = 2D0*(ux(1,1)**2 + ux(2,2)**2) + (ux(1,2) + ux(2,1))**2
      lPhi = mu * lPhi
      
      return
      end subroutine diss2D
      
!**************************************************
      
      subroutine diss3D(eNoN, Nx, ul, lPhi)
      use variables
      implicit none
      integer, intent(in) :: eNoN
      real(kind=8), intent(in) :: Nx(nsd,eNoN), ul(nsd,eNoN)
      real(kind=8), intent(out) :: lPhi
      
      integer :: a
      real(kind=8) :: ux(nsd,nsd), T1
      
      ux = 0D0
      do a=1, eNoN
         ux(1,1) = ux(1,1) + Nx(1,a)*ul(1,a)
         ux(2,1) = ux(2,1) + Nx(2,a)*ul(1,a)
         ux(3,1) = ux(3,1) + Nx(3,a)*ul(1,a)
         ux(1,2) = ux(1,2) + Nx(1,a)*ul(2,a)
         ux(2,2) = ux(2,2) + Nx(2,a)*ul(2,a)
         ux(3,2) = ux(3,2) + Nx(3,a)*ul(2,a)
         ux(1,3) = ux(1,3) + Nx(1,a)*ul(3,a)
         ux(2,3) = ux(2,3) + Nx(2,a)*ul(3,a)
         ux(3,3) = ux(3,3) + Nx(3,a)*ul(3,a)
      end do
      
      lPhi = 2D0 * ( ux(1,1)**2 + ux(2,2)**2 + ux(3,3)**2 )
      lPhi = lPhi + ( ux(1,2) + ux(2,1) )**2 &
                  + ( ux(2,3) + ux(3,2) )**2 &
                  + ( ux(3,1) + ux(1,3) )**2
      
      lPhi = mu * lPhi
      
      return
      end subroutine diss3D
      
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
      
      subroutine writeVTU(lM, fname)
      use variables
      implicit none
      type(meshType), intent(inout) :: lM
      character(len=strL), intent(in) :: fname
      
      type(vtkXMLType) :: vtu
      integer :: istat
      
      istat = 0
      do while (istat .eq. 0)
         call vtkInitWriter(vtu, fname, istat)
         if (istat .lt. 0) exit
         
         call putVTK_pointCoords(vtu, lM%x, istat)
         if (istat .lt. 0) exit
         
         lM%IEN = lM%IEN - 1
         call putVTK_elemIEN(vtu, lM%IEN, lM%vtkType, istat)
         if (istat .lt. 0) exit
         lM%IEN = lM%IEN + 1
         
         call putVTK_pointData(vtu,"FS_KE", ke, istat)
         if (istat .lt. 0) exit
         
         call putVTK_pointData(vtu,"FS_DISS", phi, istat)
         if (istat .lt. 0) exit
         
         call putVTK_pointData(vtu,"MS_Displacement", lM%dx, istat)
         if (istat .lt. 0) exit
         
         call vtkWriteToFile(vtu, istat)
         if (istat .lt. 0) exit
         
         call flushVTK(vtu)
         exit
      end do
      
      if (istat .lt. 0) then
         write(stdout,ftab4) "Error: vtk write error"
         STOP
      end if
      
      return
      end subroutine writeVTU

!**************************************************

!**************************************************

      module params

      integer, parameter :: eType_NA = 100, eType_LIN = 101, &
         eType_TRI = 102, eType_TET = 103, eType_BIL = 104, &
         eType_QUD = 105, eType_BIQ = 106, eType_BRK = 107, &
         eType_NRB = 108, eType_WDG = 109, eType_QTR = 110, &
         eType_QTE = 111

      end module params

!**************************************************

      module variables
      use genUtils
      use stdParams
      use vtkXMLMod
      use params
      implicit none

      type fsType
         integer :: eNoN
         integer :: nG
         integer :: eType = eType_NA
         real(kind=8), allocatable :: w(:)
         real(kind=8), allocatable :: xi(:,:)
         real(kind=8), allocatable :: N(:,:)
         real(kind=8), allocatable :: Nx(:,:,:)
      end type fsType

      type dataArrType
         integer :: nv
         integer :: m, n
         integer, allocatable :: s(:)
         integer, allocatable :: e(:)
         real(kind=8), allocatable :: va(:,:)
         character(len=strL), allocatable :: vname(:)
      end type dataArrType

      type meshType
         logical :: qShpF
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
      end type meshType

      integer :: nsd, nstd

      integer :: nfiles
      character(len=strL), allocatable :: vList(:)
      character(len=strL), allocatable :: fnList(:)
      character(len=strL), allocatable :: feList(:)
      real(kind=8), allocatable :: hList(:)

      type(dataArrType) :: dn, de
      type(meshType) :: mshN, mshE

      interface norm
         module procedure norms, normv
      end interface norm

      interface destroy
         module procedure destroyMesh, destroyData
      end interface destroy

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
         function norms(u) result(res)
         implicit none
         real(kind=8), intent(in) :: u(:)
         real(kind=8) res
         integer i, m

         m   = size(u)
         res = 0.0D0
         do i=1, m
            res = res + u(i)*u(i)
         end do

         return
         end function norms
!-------------------------------------------------
         function normv(u, v) result(res)
         implicit none
         real(kind=8), intent(in) :: u(:), v(:)
         real(kind=8) res
         integer i, m

         if (size(u) .ne. size(v)) then
            write(stdout,ftab4) "ERROR: incompatible vectors to "// &
     &         "compute norm"
            stop
         end if
         m   = size(u)
         res = 0.0D0
         do i=1, m
            res = res + (u(i)-v(i))*(u(i)-v(i))
         end do

         return
         end function normv
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

         return
         end subroutine destroyMesh
!-------------------------------------------------
         subroutine destroyData(d)
         implicit none
         type(dataArrType), intent(inout) :: d

         if (allocated(d%va))    deallocate(d%va)
         if (allocated(d%vname)) deallocate(d%vname)

         return
         end subroutine destroyData
!-------------------------------------------------
         function integ(lM, u, flag) result(res)
         implicit none
         type(meshType), intent(inout) :: lM
         real(kind=8), intent(in) :: u(:)
         logical, intent(in), optional :: flag
         real(kind=8) :: res

         logical :: pflag
         integer :: a, e, g, Ac, eNoN, nEl, nNo
         real(kind=8) :: Jac, tmp(nsd,nsd), sHat
         type(fsType) :: fs

         real(kind=8), allocatable :: Nx(:,:), xl(:,:), ul(:)

         if (size(u) .ne. lM%nNo) then
            write(stdout,ftab4) "ERROR: inconsistent size of input data"
            stop
         end if

         pflag = .false.
         if (present(flag)) pflag = flag

         if (.not.pflag) then
            fs%eNoN  = lM%eNoN
            fs%nG    = lM%nG
            fs%eType = lM%eType
            allocate(fs%w(fs%nG), fs%xi(nsd,fs%nG), fs%N(fs%eNoN,fs%nG),&
     &         fs%Nx(nsd,fs%eNoN,fs%nG))
            do g=1, fs%nG
               fs%w  = lM%w
               fs%xi = lM%xi
               fs%N  = lM%N
               fs%Nx = lM%Nx
            end do
         else
            call setReducedIntegFS(fs, lM%eType)
            allocate(fs%w(fs%nG), fs%xi(nsd,fs%nG), fs%N(fs%eNoN,fs%nG),&
     &         fs%Nx(nsd,fs%eNoN,fs%nG))
            call getGP(nsd, fs%eType, fs%nG, fs%w, fs%xi)

            do g=1, fs%nG
               call getShpF(nsd, fs%eType, fs%eNoN, fs%xi(:,g), &
     &            fs%N(:,g), fs%Nx(:,:,g))
            end do
         end if

         nNo  = lM%nNo
         nEl  = lM%nEl
         eNoN = fs%eNoN
         allocate(Nx(nsd,eNoN), xl(nsd,eNoN), ul(eNoN))

         res = 0D0
         do e=1, nEl
            xl = 0D0
            ul = 0D0
            do a=1, eNoN
               Ac = lM%IEN(a,e)
               xl(:,a) = lM%x(:,Ac)
               ul(a)   = u(Ac)
            end do

            do g=1, fs%nG
               call GNN(eNoN, fs%Nx(:,:,g), xl, Nx, Jac, tmp)

               sHat = 0D0
               do a=1, eNoN
                  sHat = sHat + (ul(a)*fs%N(a,g))
               end do
               res = res + (fs%w(g)*sHat*Jac)
            end do
         end do
         deallocate(Nx, xl, ul)

         return
         end function integ
!-------------------------------------------------
      end module variables

!**************************************************

      program L2Norm
      use variables
      implicit none

      logical :: pflag, qShpF
      character(len=strL) :: fname, stmp
      integer i, j, a, s, e, fid, nvar, slen
      real(kind=8) r, p, ofc

      real(kind=8), allocatable :: l2e1(:), l2e2(:,:), l2n(:,:), &
     &   sa(:)

      fid = 1889
      i = iargc()
      if (i .ne. 1) then
         write(stdout,ftab4) "Error: at least one input argument needed"
         STOP
      end if
      call getarg(1,fname)

      write(stdout,ftab1) repeat('=', 48)
      write(stdout,ftab1,advance='no') &
     &   "Number of spatial dimensions (nsd): "
      read(*,*) nsd

      if (nsd.ne.2 .and. nsd.ne.3) then
         write(stdout,ftab4) "Error: incorrect spatial dimensions"
         STOP
      end if

      if (nsd .eq. 2) then
         nstd = 3
      else
         nstd = 6
      end if

      call readInputs(fname)

      nvar = dn%nv
      allocate(l2e1(nvar), l2e2(nfiles,nvar), l2n(nfiles,nvar))
      l2e1  = 0D0
      l2e2  = 0D0
      l2n   = 0D0
      qShpF = .false.
      do i=1, nfiles
         write(stdout,ftab1) repeat('=', 48)

         write(*,ftab1) "Reading numerical solution <--- "//trim(fnList(i))
         call readVTU(mshN, dn, fnList(i))
         if (i .eq. 1) qShpF = mshN%qShpF

         write(*,ftab1) "Reading exact solution <--- "//trim(feList(i))
         call readVTU(mshE, de, feList(i))

         do j=1, nvar
            s = dn%s(j)
            e = dn%e(j)

            pflag = .false.
            if (mshN%qShpF) then
               stmp = TO_LOWER(dn%vname(j))
               slen = len(trim(stmp))
               if (stmp(1:slen) .eq. 'pressure' .or. &
     &             stmp(4:slen) .eq. 'pressure') pflag = .true.
            end if

            allocate(sa(dn%n))
            if (i .eq. 1) then
               do a=1, dn%n
                  sa(a) = norm(de%va(s:e,a))
               end do
               l2e1(j) = sqrt(integ(mshE, sa, pflag))
            end if

            do a=1, dn%n
               sa(a) = norm(dn%va(s:e,a), de%va(s:e,a))
            end do
            l2n(i,j) = sqrt(integ(mshN, sa, pflag))/l2e1(j)
            deallocate(sa)
         end do

         call destroy(mshN)
         call destroy(mshE)
      end do

      do i=1, nfiles
         if (i .eq. 1) then
            l2e2(i,:) = l2n(i,:)
         else
            r = LOG(hList(i)/hList(i-1))
            do j=1, nvar
               p = 1.0D0
               if (qShpF) then
                  stmp = TO_LOWER(dn%vname(j))
                  slen = len(trim(stmp))
                  if (stmp(1:slen) .eq. 'pressure' .or. &
     &                stmp(4:slen) .eq. 'pressure') then
                     p = 1.0D0
                  else
                     p = 2.0D0
                  end if
               end if
               l2e2(i,j) = l2e2(i-1,j)*EXP(2.0D0*p*r)
            end do
         end if
      end do

      write(stdout,ftab1) repeat('=', 48)
      do j=1, dn%nv
         ofc = 0D0
         do i=2, nfiles
            r  = LOG(hList(i)/hList(i-1))
            ofc = ofc + (LOG(l2n(i,j)/l2n(i-1,j)) / r)
         end do
         ofc = ofc / real(nfiles-1, kind=8)

         select case (dn%vname(j))
         case ("Pressure", "ST_Pressure", "NS_Pressure", "FS_Pressure",&
     &       "SS_Pressure")
            write(fname,'(A)') "l2norm_p.dat"

         case ("Displacement", "ST_Displacement", "MS_Displacement")
            write(fname,'(A)') "l2norm_u.dat"

         case ("Velocity", "ST_Velocity", "NS_Velocity", "FS_Velocity",&
     &       "SS_Velocity")
            write(fname,'(A)') "l2norm_v.dat"

         case ("Stress", "ST_Stress", "FS_Stress", "SS_Stress")
            write(fname,'(A)') "l2norm_sg.dat"

         case ("Deformation_gradient", "ST_Deformation_gradient", &
     &      "FS_Deformation_gradient")
            write(fname,'(A)') "l2norm_F.dat"

         end select

         write(*,ftab1) "Mean order of convergence for "// &
     &      TRIM(dn%vname(j))//" is: "//TRIM(STR(ofc))

         open(fid,file=trim(fname))
         write(fid,'(A)') "variables=h, rel, err2, uref"
         do i=1, nfiles
            write(fid,'(X,F10.6,4(X,1pE16.6))') hList(i), l2n(i,j), &
     &          l2e2(i,j), l2e1(j)
         end do
         close(fid)
      end do
      write(stdout,ftab1) repeat('=', 48)

      call destroy(dn)
      call destroy(de)

      deallocate(l2e1, l2e2, l2n)

      end program L2Norm

!**************************************************

      subroutine readInputs(fname)
      use variables
      implicit none
      character(len=strL), intent(in) :: fname
      integer i, fid, ntok, slen, dof, istat
      character(len=strL) :: rLine, stmp
      character(len=strL), dimension(maxToks) :: tokL

      fid = 1265
      istat = 0
      open(fid,file=trim(fname))
      OUTER_LOOP: do
         call findKwrd(fid, "numVar", istat)
         if (istat .ne. 0) exit
         read(fid,*) dn%nv

         call findKwrd(fid, "numFiles", istat)
         if (istat .ne. 0) exit
         read(fid,*) nfiles

         allocate(dn%vname(dn%nv), dn%s(dn%nv), dn%e(dn%nv))
         allocate(fnList(nfiles), feList(nfiles), hList(nfiles))

         call findKwrd(fid, "variableNameList", istat)
         if (istat .ne. 0) exit
         dn%s(1) = 1
         dn%m = 0
         do i=1, dn%nv
            read(fid,'(A)') rLine
            call parseString(rLine, tokL, ntok)
            if (ntok .eq. 1) then
               write(stdout,ftab4) "ERROR: variable dof not provided"
               istat = -1
               exit OUTER_LOOP
            end if
            dn%vname(i) = tokL(1)
            stmp = tokL(2)
            slen = len(trim(stmp))
            read(stmp(1:slen),*) dof
            if (i .gt. 1) dn%s(i) = dn%e(i-1) + 1
            dn%e(i) = dn%s(i) + dof - 1
            dn%m = dn%m + dof
         end do

         call findKwrd(fid, "exactSolutionFileList", istat)
         if (istat .ne. 0) exit
         do i=1, nfiles
            read(fid,'(A)') feList(i)
         end do

         call findKwrd(fid, "numericalSolutionFileList", istat)
         if (istat .ne. 0) exit
         do i=1, nfiles
            read(fid,'(A)') fnList(i)
         end do

         call findKwrd(fid, "gridSpacingList", istat)
         if (istat .ne. 0) exit
         do i=1, nfiles
            read(fid,*) hList(i)
         end do

         exit OUTER_LOOP
      end do OUTER_LOOP
      close(fid)

      if (istat .ne. 0) then
         write(stdout,ftab4) "ERROR: reading inputs"
         stop
      end if

      de%nv = dn%nv
      de%m  = dn%m
      allocate(de%vname(de%nv), de%s(dn%nv), de%e(dn%nv))
      de%vname(:) = dn%vname(:)
      de%s(:) = dn%s(:)
      de%e(:) = dn%e(:)

      return
      contains
         !==========================================
         subroutine findKwrd(fileId, sKwrd, istat)
         implicit none
         integer, intent(in) :: fileId
         integer, intent(inout) :: istat
         character(len=*), intent(in) :: sKwrd

         integer :: kwrdL, slen
         character(len=strL) :: sLine

         istat = 0
         kwrdL = len(trim(sKwrd))
         do
            read(fileId,'(A)',end=001) sLine
            if (sLine(1:1) .eq. '#') then
               slen  = len(trim(sLine))
               sLine = sLine(2:slen)
               sLine = adjustl(sLine)
               if (sLine(1:kwrdL) .eq. trim(sKwrd)) return
            end if
         end do

 001     write(stdout,ftab4) "ERROR: EOF reached while finding "// &
         "keyword <"//trim(sKwrd)//">"
         istat = -1
         return

         end subroutine findKwrd
         !==========================================
      end subroutine readInputs

!**************************************************

      subroutine readVTU(lM, d, fname)
      use variables
      implicit none
      type(meshType), intent(inout) :: lM
      type(dataArrType), intent(inout) :: d
      character(len=strL), intent(in) :: fname

      type(vtkXMLType) :: vtu
      integer :: istat, i, s, e, l
      character(len=strL) :: sname
      real(kind=8), allocatable :: tmpX1(:), tmpX2(:,:)

      istat = 0
      do while (istat .eq. 0)

         call loadVTK(vtu, trim(fname), istat)
         if (istat .lt. 0) exit

         call getVTK_numPoints(vtu, lM%nNo, istat)
         if (istat .lt. 0) exit

         call getVTK_numElems(vtu, lM%nEl, istat)
         if (istat .lt. 0) exit

         call getVTK_nodesPerElem(vtu, lM%eNoN, istat)
         if (istat .lt. 0) exit

         allocate(lM%ien(lM%eNoN,lM%nEl), lM%x(nsd,lM%nNo))
         call selectele(lM)

         allocate(tmpX2(maxNSD,lM%nNo))
         call getVTK_pointCoords(vtu, tmpX2, istat)
         if (istat .lt. 0) exit
         lM%x(:,:) = tmpX2(1:nsd,:)
         deallocate(tmpX2)

         call getVTK_elemIEN(vtu, lM%ien, istat)
         if (istat .lt. 0) exit
         lM%IEN(:,:) = lM%IEN(:,:) + 1

         d%n = lM%nNo
         if (allocated(d%va)) deallocate(d%va)
         allocate(d%va(d%m, d%n))
         d%va = 0.0D0
         do i=1, d%nv
            s = d%s(i)
            e = d%e(i)
            l = e - s + 1
            sname = d%vname(i)

            select case(trim(sname))
            case ("Pressure", "ST_Pressure", "NS_Pressure", &
     &         "FS_Pressure", "SS_Pressure")
               if (l .ne. 1) then
                  write(stdout,ftab4) "ERROR: inconsistent dof for "// &
     &               "var: "//trim(sname)
                  istat = -1
                  exit
               end if
               allocate(tmpX1(lM%nNo))
               call getVTK_pointData(vtu, sname, tmpX1, istat)
               if (istat .lt. 0) exit

            case ("Stress", "ST_Stress", "FS_Stress", "SS_Stress")
               if (l .ne. nstd) then
                  write(stdout,ftab4) "ERROR: inconsistent dof for "// &
     &               "var: "//trim(sname)
                  istat = -1
                  exit
               end if
               allocate(tmpX2(nstd,lM%nNo))
               call getVTK_pointData(vtu, sname, tmpX2, istat)
               if (istat .lt. 0) exit

            case ("Deformation_gradient", "ST_Deformation_gradient", &
     &         "FS_Deformation_gradient")
               if (l .ne. nsd*nsd) then
                  write(stdout,ftab4) "ERROR: inconsistent dof for "// &
     &               "var: "//trim(sname)
                  istat = -1
                  exit
               end if
               allocate(tmpX2(nsd*nsd,lM%nNo))
               call getVTK_pointData(vtu, sname, tmpX2, istat)
               if (istat .lt. 0) exit

            case default
               if (l .ne. nsd) then
                  write(stdout,ftab4) "ERROR: inconsistent dof for "// &
     &               "var: "//trim(sname)
                  istat = -1
                  exit
               end if
               allocate(tmpX2(maxNSD,lM%nNo))
               call getVTK_pointData(vtu, sname, tmpX2, istat)
               if (istat .lt. 0) exit

            end select

            if (l .eq. 1) then
               d%va(s,:) = tmpX1(:)
               deallocate(tmpX1)
            else
               d%va(s:e,:) = tmpX2(1:l,:)
               deallocate(tmpX2)
            end if
         end do
         if (istat .lt. 0) exit

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
            lM%qShpF   = .false.
         case (6)
            lM%eType   = eType_WDG
            lM%nG      = 6
            lM%vtkType = 13
            lM%qShpF   = .false.
         case (4)
            lM%eType   = eType_TET
            lM%nG      = 4
            lM%vtkType = 10
            lM%qShpF   = .false.
         case (10)
            lM%eType   = eType_QTE
            lM%nG      = 15
            lM%vtkType = 24
            lM%qShpF   = .true.
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
            lM%qShpF   = .false.
         case (4)
            lM%eType   = eType_BIL
            lM%nG      = 4
            lM%vtkType = 9
            lM%qShpF   = .false.
         case (6)
            lM%eType   = eType_QTR
            lM%nG      = 7
            lM%vtkType = 22
            lM%qShpF   = .true.
         case (9)
            lM%eType   = eType_BIQ
            lM%nG      = 9
            lM%vtkType = 28
            lM%qShpF   = .true.
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

      subroutine setReducedIntegFS(fs, eType)
      use variables
      implicit none
      type(fsType), intent(inout) :: fs
      integer, intent(in) :: eType

      select case (eType)
      case (eType_QTE)
         fs%eType   = eType_TET
         fs%nG      = 4
         fs%eNoN    = 4
      case (eType_QTR)
         fs%eType   = eType_TRI
         fs%nG      = 3
         fs%eNoN    = 3
      case (eType_BIQ)
         fs%eType   = eType_BIL
         fs%nG      = 4
         fs%eNoN    = 4
      case default
         write(stdout,ftab4) "ERROR: unknown higher order element type"
         stop
      end select

      return
      end subroutine setReducedIntegFS

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
      case (eType_QTE)
         w(1)     = 0.030283678097089D0
         w(2:5)   = 0.006026785714286D0
         w(6:9)   = 0.011645249086029D0
         w(10:15) = 0.010949141561386D0

         s = 0.25D0
         xi(1,1) = s; xi(2,1) = s; xi(3,1) = s

         s = 0.333333333333333D0
         t = 0.0D0
         xi(1,2) = t; xi(2,2) = s; xi(3,2) = s
         xi(1,3) = s; xi(2,3) = t; xi(3,3) = s
         xi(1,4) = s; xi(2,4) = s; xi(3,4) = t
         xi(1,5) = s; xi(2,5) = s; xi(3,5) = s

         s = 0.090909090909091D0
         t = 0.727272727272727D0
         xi(1,6) = t; xi(2,6) = s; xi(3,6) = s
         xi(1,7) = s; xi(2,7) = t; xi(3,7) = s
         xi(1,8) = s; xi(2,8) = s; xi(3,8) = t
         xi(1,9) = s; xi(2,9) = s; xi(3,9) = s

         s = 0.066550153573664D0
         t = 0.433449846426336D0
         xi(1,10) = s; xi(2,10) = s; xi(3,10) = t
         xi(1,11) = s; xi(2,11) = t; xi(3,11) = s
         xi(1,12) = s; xi(2,12) = t; xi(3,12) = t
         xi(1,13) = t; xi(2,13) = t; xi(3,13) = s
         xi(1,14) = t; xi(2,14) = s; xi(3,14) = t
         xi(1,15) = t; xi(2,15) = s; xi(3,15) = s

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
      CASE(eType_QTR)
         w(1)   = 0.225000000000000D0 * 5D-1
         w(2:4) = 0.125939180544827D0 * 5D-1
         w(5:7) = 0.132394152788506D0 * 5D-1

         s = 0.333333333333333D0
         xi(1,1) = s; xi(2,1) = s

         s = 0.797426985353087D0
         t = 0.101286507323456D0
         xi(1,2) = s; xi(2,2) = t
         xi(1,3) = t; xi(2,3) = s
         xi(1,4) = t; xi(2,4) = t

         s = 0.059715871789770D0
         t = 0.470142064105115D0
         xi(1,5) = s; xi(2,5) = t
         xi(1,6) = t; xi(2,6) = s
         xi(1,7) = t; xi(2,7) = t

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

      case (eType_QTE)
         s     = 1.0D0 - xi(1) - xi(2) - xi(3)
         N(1)  = xi(1)*(2.0D0*xi(1) - 1.0D0)
         N(2)  = xi(2)*(2.0D0*xi(2) - 1.0D0)
         N(3)  = xi(3)*(2.0D0*xi(3) - 1.0D0)
         N(4)  = s    *(2.0D0*s     - 1.0D0)
         N(5)  = 4.0D0*xi(1)*xi(2)
         N(6)  = 4.0D0*xi(2)*xi(3)
         N(7)  = 4.0D0*xi(1)*xi(3)
         N(8)  = 4.0D0*xi(1)*s
         N(9)  = 4.0D0*xi(2)*s
         N(10) = 4.0D0*xi(3)*s

         Nxi(1,1)  =  4.0D0*xi(1) - 1.0D0
         Nxi(2,1)  =  0.0D0
         Nxi(3,1)  =  0.0D0
         Nxi(1,2)  =  0.0D0
         Nxi(2,2)  =  4.0D0*xi(2) - 1.0D0
         Nxi(3,2)  =  0.0D0
         Nxi(1,3)  =  0.0D0
         Nxi(2,3)  =  0.0D0
         Nxi(3,3)  =  4.0D0*xi(3) - 1.0D0
         Nxi(1,4)  =  1.0D0 - 4.0D0*s
         Nxi(2,4)  =  1.0D0 - 4.0D0*s
         Nxi(3,4)  =  1.0D0 - 4.0D0*s
         Nxi(1,5)  =  4.0D0*xi(2)
         Nxi(2,5)  =  4.0D0*xi(1)
         Nxi(3,5)  =  0.0D0
         Nxi(1,6)  =  0.0D0
         Nxi(2,6)  =  4.0D0*xi(3)
         Nxi(3,6)  =  4.0D0*xi(2)
         Nxi(1,7)  =  4.0D0*xi(3)
         Nxi(2,7)  =  0.0D0
         Nxi(3,7)  =  4.0D0*xi(1)
         Nxi(1,8)  =  4.0D0*( s - xi(1) )
         Nxi(2,8)  = -4.0D0*xi(1)
         Nxi(3,8)  = -4.0D0*xi(1)
         Nxi(1,9)  = -4.0D0*xi(2)
         Nxi(2,9)  =  4.0D0*( s - xi(2) )
         Nxi(3,9)  = -4.0D0*xi(2)
         Nxi(1,10) = -4.0D0*xi(3)
         Nxi(2,10) = -4.0D0*xi(3)
         Nxi(3,10) =  4.0D0*( s - xi(3) )

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

      case (eType_QTR)
         s    = 1.0D0 - xi(1) - xi(2)
         N(1) = xi(1)*( 2.0D0*xi(1) - 1.0D0 )
         N(2) = xi(2)*( 2.0D0*xi(2) - 1.0D0 )
         N(3) = s    *( 2.0D0*s     - 1.0D0 )
         N(4) = 4.0D0*xi(1)*xi(2)
         N(5) = 4.0D0*xi(2)*s
         N(6) = 4.0D0*xi(1)*s

         Nxi(1,1) =  4.0D0*xi(1) - 1.0D0
         Nxi(2,1) =  0.0D0
         Nxi(1,2) =  0.0D0
         Nxi(2,2) =  4.0D0*xi(2) - 1.0D0
         Nxi(1,3) =  1.0D0 - 4.0D0*s
         Nxi(2,3) =  1.0D0 - 4.0D0*s
         Nxi(1,4) =  4.0D0*xi(2)
         Nxi(2,4) =  4.0D0*xi(1)
         Nxi(1,5) = -4.0D0*xi(2)
         Nxi(2,5) =  4.0D0*( s - xi(2) )
         Nxi(1,6) =  4.0D0*( s - xi(1) )
         Nxi(2,6) = -4.0D0*xi(1)

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

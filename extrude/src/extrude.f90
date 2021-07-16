!**************************************************

      module var
      use stdParams
      use genUtils
      use vtkXMLMod
      
      integer, parameter :: nsd=3
      integer :: numFaces, numEdges, numMesh
      integer :: nStart, nEnd, nFreq, nTime
      double precision :: dt, time
      character(1) :: vtkFrmt
      logical :: isBinary
      
      type :: faceType
         character(len=stdl) :: naem
         integer :: fid, nNo, nEl, eNoN
         integer, dimension(:,:), allocatable :: ien
         real(kind=8), dimension(:,:), allocatable :: x
         integer, dimension(:), allocatable :: gN, gE
      end type faceType

      type extrType
         integer :: dirID(nsd), nlev, flatEdg, fwrite
         double precision :: stepSize, dir(nsd), refPt(nsd)
         type(faceType) :: eface
      end type extrtype

      type edgeType
         integer :: eid, nNo, nEl, eNoN
         integer :: fid(2), iExtr
         integer, allocatable, dimension(:) :: gN
         integer, allocatable, dimension(:,:) :: ien
         type(extrType) :: extr
      end type edgeType
      
      type meshType
         integer :: nNo, nEl, eNoN
         integer, dimension(:,:), allocatable :: ien
         real(kind=8), dimension(:,:), allocatable :: x
         real(kind=8), dimension(:,:), allocatable :: norm
         type(edgeType), dimension(:), allocatable :: edges
         type(faceType), dimension(:), allocatable :: faces
         real(kind=8) :: vol, sarea
         integer, dimension(:), allocatable :: iFixed
      end type meshType
      
      type smoothType
         integer :: iter
         integer :: icorr
         double precision :: alpha
         double precision :: tol
      end type smoothType
      
      type(smoothType) :: smoother
      type(meshType), allocatable, dimension(:) :: meshList
      type(edgeType), allocatable, dimension(:) :: edgeList
      type(faceType), allocatable, dimension(:) :: faceList
      
      double precision, dimension(:), allocatable :: volDes
      
      end module var

!**************************************************

      program extrudeWall
      use var
      implicit none
      integer :: iM, iFa, jFa, iEdg
      integer :: i, Ac, istat, fid
      double precision, allocatable :: dVdt(:)
      character(len=stdL) :: fName, strng
      
      interface
         subroutine smoothMsh(lM, vol0)
         use var
         implicit none
         type(meshType), intent(inout) :: lM
         double precision, intent(in), optional :: vol0
         end subroutine smoothMsh
      end interface
            
      numFaces    = 0
      numEdges    = 0
      call readInputs
      
      if (vtkFrmt .eq. 'B') then
         isBinary = .true.
      else
         isBinary = .false.
      end if
      
      do iFa=1, numFaces
         call loadVTP(faceList(iFa))
      end do
      
      write(stdout,'(A)') "======================================"
      write(stdout,'(A)')
      write(stdout,ftab1) "Template faces loaded.."
      write(stdout,'(A)')
      do iFa=1, numFaces
         write(stdout,ftab2) "Face: "//trim(faceList(iFa)%naem)
         write(stdout,ftab3) "nPoints =  "//trim(STR(faceList(iFa)%nNo))
         write(stdout,ftab3) "nElems =  "//trim(STR(faceList(iFa)%nEl))
      end do
      
      do iEdg=1, numEdges
         do i=1, numFaces
            if (edgeList(iEdg)%fid(1) .eq. faceList(i)%fid) iFa = i
            if (edgeList(iEdg)%fid(2) .eq. faceList(i)%fid) jFa = i
         end do
         call calcBoundEdge(faceList(iFa),faceList(jFa),edgeList(iEdg))
      end do
      
      write(stdout,'(A)')
      do iM=1, numMesh
      
         write(stdout,'(A)') "======================================"
         write(stdout,'(A)')
         nTime = nStart + (iM-1)*nFreq
         write(stdout,ftab1) "Processing mesh at nTime="// &
            trim(STR(nTime))
         allocate(meshList(iM)%faces(numFaces))
         write(stdout,ftab2) "Initializing faces.."
         do iFa=1, numFaces
            call initMeshFaces(meshList(iM)%faces(iFa), faceList(iFa))

            write(fName,'(A,I2.2,A,I2.2,A)') "01-morphed-mshsrf_",iFa, &
               "/", nTime, "_registered.vtk"
            strng = meshList(iM)%faces(iFa)%naem
            i = len(trim(strng))
            write(stdout,ftab3) "Reading mesh face <"//strng(1:i-4)// &
               "> from file <"//trim(fName)//">"
            call readFaceData(meshList(iM)%faces(iFa), fName)
         end do
         
         allocate(meshList(iM)%edges(numEdges))
         do iEdg=1, numEdges
            write(stdout,ftab2) "Initializing edge "// &
               trim(STR(edgeList(iEdg)%eid))
            call initMeshEdges(meshList(iM)%edges(iEdg), edgeList(iEdg))
            
            if (meshList(iM)%edges(iEdg)%iExtr .gt. 0) then
               do i=1, numFaces
                  if (meshList(iM)%faces(i)%fid .eq. &
                      meshList(iM)%edges(iEdg)%fid(1)) iFa = i
                  if (meshList(iM)%faces(i)%fid .eq. &
                      meshList(iM)%edges(iEdg)%fid(2)) jFa = i
               end do

               write(stdout,ftab2) "Extruding edge "// &
                  trim(STR(meshList(iM)%edges(iEdg)%eid))
               call extrude( meshList(iM)%faces(iFa),  &
                             meshList(iM)%faces(jFa),  &
                             meshList(iM)%edges(iEdg) )
            end if
         end do
         
         write(stdout,ftab2) &
            'Combining faces to create a mesh surface..'
         call createMsh(meshList(iM))
         
         if ((smoother%iter.gt.0) .or. (smoother%icorr.eq.2)) then
            if (smoother%icorr .eq. 2) then
               write(stdout,ftab2) &
                  "Correcting volume to the desired value.."
               call smoothMsh(meshList(iM), volDes(iM))
            else
               write(stdout,ftab2) "Smoothing.."
               call smoothMsh(meshList(iM))
            end if
         end if
         
         write(stdout,ftab2) "Writing mesh and face data to file.."
         write(fName,'("02-marker-data/marker_",I2.2,".vtk")') nTime
         call writeVTK(meshList(iM)%nNo, meshList(iM)%eNoN, &
            meshList(iM)%nEl, meshList(iM)%x, meshList(iM)%ien, fName)
         do iFa=1, numFaces
            i = len(trim(meshList(iM)%faces(iFa)%naem))
            strng = meshList(iM)%faces(iFa)%naem(1:i-4)
            write(fName,'(A,I2.2,A,I2.2,A)') "01-morphed-mshsrf_", iFa, &
               "/"//trim(strng)//"_", nTime, ".vtk"
            call writeVTK(meshList(iM)%faces(iFa)%nNo, &
               meshList(iM)%faces(iFa)%eNoN, meshList(iM)%faces(iFa)%nEl, &
               meshList(iM)%faces(iFa)%x, meshList(iM)%faces(iFa)%ien, fName)
         end do
         
      end do
      
      allocate(dVdt(numMesh)); dVdt = 0D0
      fid = 100
      write(fName, '("volFlowRaw.dat")')
      open(fid, file=trim(fName))
      write(fid, '(A)') "Variables=t, V, Q"
      do iM=1, numMesh
         nTime = nStart + (iM-1)*nFreq
         time = dble(nTime-nStart) * dt
         if (iM .gt. 1) then
            dVdt(iM) = (meshList(iM)%vol - meshList(iM-1)%vol) / &
               (dble(nFreq)*dt)
         end if
         write(fid,'(A)') STR(time)//"   "// &
            STR(meshList(iM)%vol)//"   "//STR(dVdt(iM))
      end do
      close(fid)
      
      contains
      
         !==========================================
         
         subroutine loadVTU(lM)
         implicit none
         integer :: i,j
         type(meshType), intent(inout) :: lM
         type(vtkXMLType) :: vtu
         
         call loadVTK(vtu,"data/mesh-complete.mesh.vtu",istat)
            if ( istat.lt.0 ) then
               write(stdout,ftab4) "ERROR: VTK file read error.."
               STOP
            end if
         
         call getVTK_numPoints(vtu,lM%nNo,istat)
            if ( istat.lt.0 ) then
               write(stdout,ftab4) "ERROR: VTK file read error.."
               STOP
            end if
         
         call getVTK_numElems(vtu,lM%nEl,istat)
            if ( istat.lt.0 ) then
               write(stdout,ftab4) "ERROR: VTK file read error.."
               STOP
            end if
         
         call getVTK_nodesPerElem(vtu,lM%eNoN,istat)
            if ( istat.lt.0 ) then
               write(stdout,ftab4) "ERROR: VTK file read error.."
               STOP
            end if
         
         allocate(lM%ien(lM%eNoN,lM%nEl))
         allocate(lM%x(nsd,lM%nNo))
         
         call getVTK_pointCoords(vtu,lM%x,istat)
            if ( istat.lt.0 ) then
               write(stdout,ftab4) "ERROR: VTK file read error.."
               STOP
            end if
         
         call getVTK_elemIEN(vtu,lM%ien,istat)
            if ( istat.lt.0 ) then
               write(stdout,ftab4) "ERROR: VTK file read error.."
               STOP
            end if
         
         call flushVTK(vtu)
         
         end subroutine loadVTU
         
         !==========================================
            
         subroutine loadVTP(fa)
         implicit none
         integer :: i,j
         type(faceType), intent(inout) :: fa
         type(vtkXMLType) :: vtp
         
         call loadVTK(vtp,fa%naem,istat)
            if ( istat.lt.0 ) then
               write(stdout,ftab4) "ERROR: VTK file read error.."
               STOP
            end if
         
         call getVTK_numPoints(vtp,fa%nNo,istat)
            if ( istat.lt.0 ) then
               write(stdout,ftab4) "ERROR: VTK file read error.."
               STOP
            end if
         
         call getVTK_numElems(vtp,fa%nEl,istat)
            if ( istat.lt.0 ) then
               write(stdout,ftab4) "ERROR: VTK file read error.."
               STOP
            end if
         
         call getVTK_nodesPerElem(vtp,fa%eNoN,istat)
            if ( istat.lt.0 ) then
               write(stdout,ftab4) "ERROR: VTK file read error.."
               STOP
            end if
         
         allocate(fa%ien(fa%eNoN,fa%nEl))
         allocate(fa%x(nsd,fa%nNo))
         allocate(fa%gN(fa%nNo))
         allocate(fa%gE(fa%nEl))
         
         call getVTK_pointCoords(vtp,fa%x,istat)
            if ( istat.lt.0 ) then
               write(stdout,ftab4) "ERROR: VTK file read error.."
               STOP
            end if
         
         call getVTK_elemIEN(vtp,fa%ien,istat)
            if ( istat.lt.0 ) then
               write(stdout,ftab4) "ERROR: VTK file read error.."
               STOP
            end if
         
         call getVTK_pointData(vtp,'GlobalNodeID',fa%gN,istat)
            if ( istat.lt.0 ) then
               write(stdout,ftab4) "WARNING: GlobalNodeID not found"
               deallocate(fa%gN)
            end if
         
         call getVTK_elemData(vtp,'GlobalElementID',fa%gE,istat)
            if ( istat.lt.0 ) then
               write(stdout,ftab4) "WARNING: GlobalElementID not found"
               deallocate(fa%gE)
            end if
         
         call flushVTK(vtp)
         
         end subroutine loadVTP
         
         !==========================================
            
      end program extrudeWall
      
!**************************************************

      subroutine readInputs
      use var
      implicit none
      integer :: iM, iFa, iEdg, itmp
      double precision :: rtmp
      
      open(10,file='in_params.dat')
      read(10,*)
      
      ! stitch parameters !
      read(10,*)
      read(10,*)
      read(10,*) nStart, nEnd, nFreq, dt
      read(10,*)
      
      ! face parameters !
      read(10,*)
      read(10,*)
      read(10,*) numFaces
      if (numFaces .gt. 0) then
         allocate(faceList(numFaces))
      end if
      read(10,*)
      do iFa=1, numFaces
         read(10,*) faceList(iFa)%fid, faceList(iFa)%naem
      end do
      read(10,*)

      ! edge parameters !
      read(10,*)
      read(10,*) 
      read(10,*) numEdges
      if (numEdges .gt. 0) then
         allocate(edgeList(numEdges))
         do iEdg=1, numEdges
            edgeList(iEdg)%eid    = 0
            edgeList(iEdg)%eNoN   = 2
            edgeList(iEdg)%iExtr  = 0
            edgeList(iEdg)%extr%nlev     = 0
            edgeList(iEdg)%extr%flatEdg  = 0
            edgeList(iEdg)%extr%fwrite   = 0
            edgeList(iEdg)%extr%stepSize = 0D0
            edgeList(iEdg)%extr%dirID(:) = 0
            edgeList(iEdg)%extr%dir(:)   = 0D0
            edgeList(iEdg)%extr%refPt(:) = 0D0
         end do
      end if
      
      read(10,*)
      do iEdg=1, numEdges
         read(10,*) edgeList(iEdg)%eid, edgeList(iEdg)%fid(:), edgeList(iEdg)%iExtr
      end do
      read(10,*)
      
      ! extrude options !
      read(10,*)
      read(10,*)
      do iEdg=1, numEdges
         if (edgeList(iEdg)%iExtr .gt. 0) then
            read(10,*) itmp, edgeList(iEdg)%extr%dirID(:), &
               edgeList(iEdg)%extr%stepSize, edgeList(iEdg)%extr%nlev, &
               edgeList(iEdg)%extr%flatEdg, edgeList(iEdg)%extr%fwrite
         else
            read(10,*)
         end if
      end do
      read(10,*)
      
      ! smoothing parameters !
      read(10,*)
      read(10,*)
      read(10,*) smoother%iter, smoother%alpha, smoother%icorr, &
         smoother%tol
      read(10,*)
      
      ! output options !
      read(10,*)
      read(10,*) vtkFrmt
      
      !read(10,*)
      !read(10,*) useDiffWall
      !read(10,*)
      !read(10,*) ntExtr
      close(10)

      numMesh = ((nEnd - nStart) / nFreq) + 1
      allocate(meshList(numMesh))
      allocate(volDes(numMesh)); volDes = 0D0
      
      if (smoother%icorr .eq. 2) then
         open(10,file='volFlowDesired.dat')
         read(10,*)
         do iM=1, numMesh
            read(10,*) rtmp, volDes(iM), rtmp
         end do
         close(10)
      end if
      end subroutine readInputs

!**************************************************

      subroutine calcBoundEdge(pFa, qFa, lEdg)
      use var
      implicit none
      type(faceType), intent(in) :: pFa, qFa
      type(edgeType), intent(inout) :: lEdg
      integer :: nNo, nEl
      integer :: i, j, k, e, a, Ac, Ac1, Ac2
      integer, allocatable :: tmpI(:)
      double precision :: u(nsd), v(nsd), w(nsd), sgn
      
      write(stdout,ftab1) "Determining edge between template face <"// &
         trim(pFa%naem)// "> and template face <"//trim(qFa%naem)//">"
      
      call calcDirExtr(qFa, lEdg)
      
      call calcRefPt(qFa, lEdg)
      
      nNo = min(pFa%nNo, size(pFa%gN))
      allocate(tmpI(min(nNo, qFa%nNo)))
      tmpI = 0
      lEdg%nNo = 0
      do i=1, nNo
         Ac1 = pFa%gN(i)
         do j=1, qFa%nNo
            Ac2 = qFa%gN(j)
            if (Ac1 .eq. Ac2) then
               lEdg%nNo = lEdg%nNo + 1
               tmpI(lEdg%nNo) = i
            end if
         end do
      end do
      
      write(stdout,ftab2) "Number of edge nodes = "//trim(STR(lEdg%nNo))
      write(stdout,ftab2) "Determining edge connectivity.."
      lEdg%nEl = lEdg%nNo
      allocate(lEdg%ien(2,lEdg%nEl))
      lEdg%ien = 0
      k = 0
      do e=1, pFa%nEl
         j = 0
         do a=1, pFa%eNoN
            Ac = pFa%ien(a,e) + 1
            do i=1, lEdg%nNo
               if (Ac .eq. tmpI(i)) then
                  j = j+1
                  if (j .eq. 1) Ac1 = Ac
                  if (j .eq. 2) Ac2 = Ac
                  exit
               end if
            end do ! i
         end do ! a
         if (j .eq. 2) then
            u(:) = pFa%x(:,Ac1) - lEdg%extr%refPt(:)
            v(:) = pFa%x(:,Ac2) - lEdg%extr%refPt(:)
            w(:) = CROSS(u,v)
            sgn = sum( lEdg%extr%dir(:) * w(:) )
            k = k+1
            if ( sgn .ge. 0.0d0) then
               lEdg%ien(1,k) = Ac1
               lEdg%ien(2,k) = Ac2
            else
               lEdg%ien(1,k) = Ac2
               lEdg%ien(2,k) = Ac1
            end if
         end if
      end do ! e
      deallocate(tmpI)
      
      allocate(lEdg%gN(lEdg%nNo))
      lEdg%gN = 0
      lEdg%gN(1) = lEdg%ien(1,1)
      i=1; j=1
      do
         Ac = lEdg%ien(2,j)
         do j=1, lEdg%nNo
            if (Ac .eq. lEdg%ien(1,j)) then
               i=i+1
               lEdg%gN(i) = Ac
               exit
            end if
         end do
         if (i .eq. lEdg%nNo) exit
      end do
      
      end subroutine calcBoundEdge
      
!**************************************************

      subroutine calcDirExtr(lFa, lEdg)
      use var
      implicit none
      type(faceType), intent(in) :: lFa
      type(edgeType), intent(inout) :: lEdg
      
      integer :: i, e
      double precision :: mag, a(nsd), b(nsd)
      
      if (lEdg%iExtr .gt. 0) then
         a(:) = lFa%x(:,lEdg%extr%dirID(2)) - &
            lFa%x(:,lEdg%extr%dirID(1))
         b(:) = lFa%x(:,lEdg%extr%dirID(3)) - &
            lFa%x(:,lEdg%extr%dirID(1))
      else
         e = lFa%nEl/2
         a(:) = lFa%x(:,lFa%ien(2,e)) - &
            lFa%x(:,lFa%ien(1,e))
         b(:) = lFa%x(:,lFa%ien(3,e)) - &
            lFa%x(:,lFa%ien(1,e))
      end if
      
      lEdg%extr%dir(:) = CROSS(a,b)
      mag = dsqrt( sum(lEdg%extr%dir(:)**2) )
      lEdg%extr%dir(:) = lEdg%extr%dir(:)/mag
      
      end subroutine calcDirExtr
      
!**************************************************

      subroutine calcRefPt(lFa, lEdg)
      use var
      implicit none
      type(faceType), intent(in) :: lFa
      type(edgeType), intent(inout) :: lEdg
      integer :: i
      
      lEdg%extr%refPt(:) = 0D0
      do i=1, lFa%nNo
         lEdg%extr%refPt(:) = lEdg%extr%refPt(:) + lFa%x(:,i)
      end do
      lEdg%extr%refPt(:) = lEdg%extr%refPt(:)/dble(lFa%nNo)
      
      end subroutine calcRefPt
      
!**************************************************

      subroutine initMeshFaces(lFa, pFa)
      use var
      implicit none
      type(faceType), intent(inout) :: lFa, pFa
      
      lFa%fid   =  pFa%fid
      lFa%naem  =  pFa%naem
      lFa%nNo   =  pFa%nNo
      lFa%nEl   =  pFa%nEl
      lFa%eNoN  =  pFa%eNoN
      allocate(lFa%x(nsd,lFa%nNo))
      allocate(lFa%gN(lFa%nNo))
      allocate(lFa%ien(lFa%eNoN,lFa%nEl))
      lFa%gN(:) = pFa%gN(:)
      lFa%ien(:,:) = pFa%ien(:,:)
      
      end subroutine initMeshFaces

!**************************************************

      subroutine readFaceData(lFa, fName)
      use var
      implicit none
      type(faceType), intent(inout) :: lFa
      character(len=*), intent(in)  :: fName
      integer :: fid, Ac
      
      fid = 10
      open(fid, file=trim(fName))
      read(fid,*)
      read(fid,*)
      read(fid,*)
      read(fid,*)
      read(fid,*)
      do Ac=1, lFa%nNo
         read(fid,*) lFa%x(:,Ac)
      end do
      close(fid)
      
      end subroutine readFaceData
      
!**************************************************

      subroutine initMeshEdges(lEdg, pEdg)
      use var
      implicit none
      type(edgeType), intent(inout) :: lEdg, pEdg
      
      lEdg%eid    =  pEdg%eid
      lEdg%nNo    =  pEdg%nNo
      lEdg%nEl    =  pEdg%nEl
      lEdg%eNoN   =  pEdg%eNoN
      lEdg%fid(:) =  pEdg%fid(:)
      
      allocate(lEdg%gN(lEdg%nNo))
      allocate(lEdg%ien(lEdg%eNoN, lEdg%nEl))
      lEdg%gN(:)    =  pEdg%gN(:)
      lEdg%ien(:,:) = pEdg%ien(:,:)
      
      lEdg%iExtr         =   pEdg%iExtr
      lEdg%extr%nlev     =   pEdg%extr%nlev
      lEdg%extr%flatEdg  =   pEdg%extr%flatEdg
      lEdg%extr%fwrite   =   pEdg%extr%fwrite
      lEdg%extr%stepSize =   pEdg%extr%stepSize
      lEdg%extr%dirID(:) =   pEdg%extr%dirID(:)
      lEdg%extr%dir(:)   =   pEdg%extr%dir(:)
      lEdg%extr%refPt(:) =   pEdg%extr%refPt(:)
      
      end subroutine initMeshEdges
      
!**************************************************
      
      subroutine extrude(pFa, qFa, lEdg)
      use var
      implicit none
      type(faceType), intent(inout) :: pFa, qFa
      type(edgeType), intent(inout) :: lEdg
      
      integer :: i, j, k, Ac, nNo, nEl
      double precision :: ds
      integer, allocatable :: tmpI(:,:)
      double precision, allocatable :: disp(:,:), tmpX(:,:)
      
      lEdg%extr%eface%nNo = lEdg%extr%nlev * lEdg%nNo
      allocate(lEdg%extr%eface%x(nsd,lEdg%extr%eface%nNo))
      allocate(lEdg%extr%eface%gN(lEdg%extr%eface%nNo))
      allocate(disp(nsd,lEdg%nNo))
      
      call calcRefPt(qFa, lEdg)
      
!     transform the connected face (jFa) first
      do i=1, qFa%nNo
         if (lEdg%extr%flatEdg .gt. 0) then
            ds = sum( (qFa%x(:,i) - lEdg%extr%refPt(:)) &
               * lEdg%extr%dir(:) )
         else
            ds = 0.0d0
         end if
         ds = dble(lEdg%extr%nlev) * lEdg%extr%stepSize - ds
         qFa%x(:,i) = qFa%x(:,i) + ds*lEdg%extr%dir(:)
      end do
      
!     compute edge displacements
      disp = 0.0d0
      do i=1, lEdg%nNo
         Ac = lEdg%gN(i)
         do j=1, qFa%nNo
            if (pFa%gN(Ac) .eq. qFa%gN(j)) then
               disp(:,i) = qFa%x(:,j) - pFa%x(:,Ac)
               exit
            end if
         end do
      end do
      
!     extrude
      j = 0
      do k=1, lEdg%extr%nlev
         do i=1, lEdg%nNo
            j = j+1
            Ac = lEdg%gN(i)
            ds = dble(k) / ( dble(lEdg%extr%nLev) + epsilon(ds) )
            lEdg%extr%eface%x(:,j) = pFa%x(:,Ac) + ds*disp(:,i)
            lEdg%extr%eface%gN(j)  = pFa%nNo + j
         end do
      end do
      deallocate(disp)
      
      lEdg%extr%eface%nEl = (lEdg%extr%nlev-1) * lEdg%nNo * 2
      lEdg%extr%eface%eNoN = pFa%eNoN
      allocate(lEdg%extr%eface%ien(lEdg%extr%eFace%eNoN, lEdg%extr%eface%nEl))
      j = 0
      do k=1, lEdg%extr%nlev-1
         do i=1, lEdg%nNo
            j = j+1
            lEdg%extr%eface%ien(1,j) = (k-1)*lEdg%nNo + i
            lEdg%extr%eface%ien(2,j) = (k-1)*lEdg%nNo + i + 1
            lEdg%extr%eface%ien(3,j) = (k  )*lEdg%nNo + i
            if (i .eq. lEdg%nNo) lEdg%extr%eface%ien(2,j) = &
               (k-1)*lEdg%nNo + 1
            
            j = j+1
            lEdg%extr%eface%ien(1,j) = (k-1)*lEdg%nNo + i + 1
            lEdg%extr%eface%ien(2,j) = (k  )*lEdg%nNo + i + 1
            lEdg%extr%eface%ien(3,j) = (k  )*lEdg%nNo + i
            if (i .eq. lEdg%nNo) then
               lEdg%extr%eface%ien(1,j) = &
                  (k-1)*lEdg%nNo + 1
               lEdg%extr%eface%ien(2,j) = &
                  (k  )*lEdg%nNo + 1
            end if
         end do
      end do
      lEdg%extr%eface%ien(:,:) = lEdg%extr%eface%ien(:,:) - 1
      
      ! stitch extruded face with parent face !
      nNo = pFa%nNo
      nEl = pFa%nEl
      allocate(tmpX(nsd,nNo))
      allocate(tmpI(pFa%eNoN,nEl))
      tmpX = pFa%x
      tmpI = pFa%ien
      
      deallocate(pFa%x, pFa%ien)
      pFa%nNo = nNo + lEdg%extr%eface%nNo
      pFa%nEl = nEl + lEdg%extr%eface%nEl + 2*lEdg%nNo
      allocate(pFa%x(nsd,pFa%nNo))
      allocate(pFa%ien(pFa%eNoN,pFa%nEl))
      
      pFa%x(:,1:nNo) = tmpX(:,1:nNo)
      pFa%ien(:,1:nEl) = tmpI(:,1:nEl)
      deallocate(tmpX, tmpI)
      
      do i=1, lEdg%extr%eface%nNo
         j = nNo + i
         pFa%x(:,j) = lEdg%extr%eface%x(:,i)
      end do
      
      j = nEl
      do i=1, lEdg%nNo
         j = j + 1
         pFa%ien(1,j) = lEdg%gN(i  ) - 1
         if (i .eq. lEdg%nNo) then
            pFa%ien(2,j) = lEdg%gN(1) - 1
         else
            pFa%ien(2,j) = lEdg%gN(i+1) - 1
         end if
         pFa%ien(3,j) = nNo + i - 1
         
         j = j + 1
         if (i .eq. lEdg%nNo) then
            pFa%ien(1,j) = lEdg%gN(1) - 1
            pFa%ien(2,j) = nNo
         else
            pFa%ien(1,j) = lEdg%gN(i+1) - 1
            pFa%ien(2,j) = nNo + i
         end if
         pFa%ien(3,j) = nNo + i - 1
      end do
      
      j = nEl + 2*lEdg%nNo
      do i=1, lEdg%extr%eface%nEl
         pFa%ien(:,j+i) = nNo + lEdg%extr%eface%ien(:,i)
      end do
      
      if (lEdg%extr%fwrite .gt. 0) &
         call writeExtrOutput(pFa, qFa, lEdg)
      
      end subroutine extrude
      
!**************************************************

      subroutine writeExtrOutput(pFa, qFa, lEdg)
      use var
      implicit none
      type(faceType), intent(in) :: pFa, qFa
      type(edgeType), intent(in) :: lEdg
      character(len=stdl) :: fname, s1, s2
      integer :: i, j, Ac, fid
      
      s1 = ""
      s2 = ""
      fid = 10
      s1 = pFa%naem(1:len(trim(pFa%naem))-4)
      s2 = qFa%naem(1:len(trim(qFa%naem))-4)

!      write(stdout,ftab2) "Writing boundary data to file.."
!      write(fname,'(A)') &
!         "bedge_"//trim(s1)//"_"//trim(s2)//".dat"
!      open(fid,file=trim(fname))
!      write(fid,'(A)') 'variables=node,x,y,z'
!      do i=1, lEdg%nNo
!         Ac = lEdg%gN(i)
!         write(fid,'(A)',advance='no') trim(STR(Ac))
!         do j=1, nsd
!            write(fid,'(A)',advance='no') " "//STR(pFa%x(j,Ac))
!         end do
!         write(fid,'(A)')
!      end do
!      close(fid)
      
      write(fname,'(A)') &
         "extr_"//trim(s1)//"_"//trim(s2)//".vtk"
      
      if ( isBinary ) then
         write(stdout,ftab3) &
            "Writing extruded face data to file (binary).."
         open(fid,file=trim(fName),status='replace',access='stream',&
         form='unformatted',convert='big_endian')
         write(fid) '# vtk DataFile Version 3.0'//newl
         write(fid) 'Simulation Results'//newl
         write(fid) 'BINARY'//newl
         write(fid) 'DATASET UNSTRUCTURED_GRID'//newl
         write(fid) 'POINTS '//trim(STR(pFa%nNo))//' double'//newl
         do i=1, pFa%nNo
         write(fid) pFa%x(:,i)
       end do
       write(fid) 'CELLS '//trim(STR(pFa%nEl))//' '// &
          trim(STR(4*pFa%nEl))//newl
       do i=1, pFa%nEl
          write(fid) pFa%eNoN, pFa%ien(:,i)
       end do
       write(fid) 'CELL_TYPES '//trim(STR(pFa%nEl))//newl
       do i=1, pFa%nEl
          write(fid) 5
       end do
       close(fid)
      else
         write(stdout,ftab3) &
            "Writing extruded face data to file (ascii).."
         open(fid,file=trim(fName),status='replace')
         write(fid,'(A)') '# vtk DataFile Version 3.0'
         write(fid,'(A)') 'Simulation Results'
         write(fid,'(A)') 'ASCII'
         write(fid,'(A)') 'DATASET UNSTRUCTURED_GRID'
         write(fid,'(A)') 'POINTS '//trim(STR(pFa%nNo))//' double'
         do i=1, pFa%nNo
            do j=1, nsd
               write(fid,'(A)',advance='no') STR(pFa%x(j,i))//' '
            end do
            write(fid,'(A)')
         end do
         write(fid,'(A)') &
            'CELLS '//trim(STR(pFa%nEl))//' '//trim(STR(4*pFa%nEl))
         do i=1, pFa%nEl
            write(fid,'(A)',advance='no') '3'
            do j=1, pFa%eNoN
               write(fid,'(A)',advance='no') ' '// &
                  trim(STR(pFa%ien(j,i)))
            end do
            write(fid,'(A)')
         end do
         write(fid,'(A)') 'CELL_TYPES '//trim(STR(pFa%nEl))
         do i=1, pFa%nEl
            write(fid,'(A)') '5'
         end do
         close(fid)
      end if
      
      end subroutine writeExtrOutput
      
!**************************************************

      subroutine createMsh(lM)
      use var
      implicit none
      type(meshType), intent(inout) :: lM
      integer :: iFa, jFa, iEdg, nNo, nEl
      integer :: i, j, k, Ac, e, a
      integer, dimension(:), allocatable :: inclFace, mapIEN
      logical :: bflag
      
      allocate(inclFace(numFaces))
      
      inclFace = 0
      do iEdg=1, numEdges
         do iFa=1, numFaces
            if ( lM%faces(iFa)%fid .eq. lM%edges(iEdg)%fid(1) ) exit
         end do
         do jFa=1, numFaces
            if ( lM%faces(jFa)%fid .eq. lM%edges(iEdg)%fid(2) ) exit
         end do
         inclFace(iFa) = 1
         inclFace(jfa) = 1
      end do
      
      lM%nNo = 0
      lM%nEl = 0
      do iFa=1, numFaces
         if (inclFace(iFa) .eq. 1) then
            lM%nNo = lM%nNo + lM%faces(iFa)%nNo
            lM%nEl = lM%nEl + lM%faces(iFa)%nEl
         end if
      end do
      lM%nNo = lM%nNo - sum(lM%edges(:)%nNo)
      lM%eNoN = lM%faces(1)%eNoN
      
      allocate(lM%x(nsd,lM%nNo))
      allocate(lM%ien(lM%eNoN,lM%nEl))
      allocate(lM%norm(nsd, lM%nEl))
      allocate(lM%iFixed(lM%nNo))
      lM%iFixed(:) = 0
      
      nNo = 0; nEl = 0
      do iEdg=1, numedges
         do iFa=1, numFaces
            if ( lM%faces(iFa)%fid .eq. lM%edges(iEdg)%fid(1) ) exit
         end do
         
         do jFa=1, numFaces
            if ( lM%faces(jFa)%fid .eq. lM%edges(iEdg)%fid(2) ) exit
         end do
         
         if (inclFace(iFa) .eq. 1) then
            lM%x(:,nNo+1:nNo+lM%faces(iFa)%nNo) = &
               lM%faces(iFa)%x(:,:)
            lM%ien(:,nEl+1:nEl+lM%faces(iFa)%nEl) = &
               lM%faces(iFa)%ien(:,:)
         end if
         
         nNo = nNo + inclFace(iFa)*lM%faces(iFa)%nNo
         nEl = nEl + inclFace(iFa)*lM%faces(iFa)%nEl
         if (inclFace(jFa) .eq. 1) then
            allocate(mapIEN(lM%faces(jFa)%nNo))
            mapIEN = 0
            k = (lM%edges(iEdg)%extr%nlev-1) * lM%edges(iEdg)%nNo
            do j=1, lM%faces(jFa)%nNo
               bflag = .false.
               do i=1, lM%edges(iEdg)%nNo
                  Ac = lM%edges(iEdg)%gN(i)
                  if (lM%faces(iFa)%gN(Ac) .eq. lM%faces(jFa)%gN(j)) then
                     bflag = .true.
                     exit
                  end if
               end do
               if (.not.bflag) then
                  nNo = nNo + 1
                  lM%x(:,nNo) = lM%faces(jFa)%x(:,j)
                  mapIEN(j) = nNo
               else
                  if (lM%edges(iEdg)%iExtr .eq. 1) then
                     mapIEN(j) = lM%edges(iEdg)%extr%eface%gN(k+i)
                  else
                     mapIEN(j) = lM%edges(iEdg)%gN(i)
                  end if
                  if (mapIEN(j) .gt. lM%nNo) then
                     write(stdout,ftab4) &
                        "Error: out of bounds while mapping edges.."
                     STOP
                  end if
                  lM%iFixed(mapIEN(j)) = 1
               end if
            end do
            
            do e=1, lM%faces(jFa)%nEl
               nEl = nEl + 1
               do a=1, lM%faces(jFa)%eNoN
                  Ac = lM%faces(jFa)%IEN(a,e) + 1
                  lM%ien(a,nEl) = mapIEN(Ac) - 1
               end do
            end do
            deallocate(mapIEN)
         end if
         inclFace(iFa) = 0
         inclFace(jFa) = 0
      end do
      
      deallocate(inclFace)

      end subroutine createMsh
      
!**************************************************

      subroutine smoothMsh(lM, volD)
      use var
      implicit none
      type(meshType), intent(inout) :: lM
      double precision, intent(in), optional :: volD
      
      logical :: skip
      integer :: iter, i, a, e, Ac, iel
      double precision :: vol, vol0, ds, ds0, dir(nsd)
      
      integer, dimension(:), allocatable :: tmpI
      integer, dimension(:,:), allocatable :: eList
      double precision, dimension(:,:), allocatable :: tmpX
      logical, dimension(:), allocatable :: eFixed
      
      allocate(tmpX(nsd,lM%nNo))
      allocate(tmpI(lM%nNo))
      allocate(eFixed(lM%nEl))
      eFixed = .false.
      
      lM%ien(:,:) = lM%ien(:,:) + 1
      
      call calcVolume(lM)
      write(stdout,ftab3) "Initial volume: "//STR(lM%vol)
      if (present(volD)) then
         vol0 = volD
         write(stdout, ftab3) "Desired volume: "//STR(vol0)
      else
         vol0 = lM%vol
      end if
      
      do e=1, lM%nEl
         do a=1, lM%eNoN
            Ac = lM%ien(a,e)
            if (lM%iFixed(Ac) .eq. 1) eFixed(e) = .true.
         end do
      end do
      
      tmpI = 0
      do e=1, lM%nEl
         do a=1, lM%eNoN
            Ac = lM%ien(a,e)
            tmpI(Ac) = tmpI(Ac) + 1
         end do
      end do
      
      allocate(eList(maxval(tmpI),lM%nNo))
      tmpI = 0
      do e=1, lM%nEl
         do a=1, lM%eNoN
            Ac = lM%ien(a,e)
            tmpI(Ac) = tmpI(Ac) + 1
            eList(tmpI(Ac), Ac) = e
         end do
      end do
      
      if (smoother%icorr .eq. 2) then
         write(stdout, ftab3) "Deforming surface to the desired volume"
         ds0 = dsqrt(lM%sarea/dble(lM%nEl))
         tmpX = lM%x
         do
            call calcVolume(lM)
            vol = lM%vol
            if (dabs(vol0-vol) .lt. smoother%tol) exit
            ds = ds0*(vol0-vol) / (3D0*vol0)
            do i=1, lM%nNo
               skip = .false.
               do e=1, tmpI(i)
                  iel = eList(e,i)
                  if (eFixed(iel) .or. i.gt.lM%faces(1)%nNo) &
                     skip = .true.
               end do
               if (skip) then
                  tmpX(:,i) = lM%x(:,i)
                  cycle
               end if
               dir(:) = 0D0
               do e=1, tmpI(i)
                  iel = eList(e,i)
                  if (eFixed(iel)) cycle
                  dir(:) = dir(:) + lM%norm(:,iel)
               end do
               dir(:) = dir(:) / dble(tmpI(i))
               tmpX(:,i) = lM%x(:,i) + ds*dir(:)
            end do
            lM%x = tmpX
         end do
         call calcVolume(lM)
         vol0 = lM%vol
      end if
      
      do iter=1, smoother%iter
         write(stdout,ftab3) 'Smoothing iteration '//trim(STR(iter))
         tmpX = 0.0d0
         do i=1, lM%nNo
            skip = .false.
            do e=1, tmpI(i)
               iel = eList(e,i)
               if (eFixed(iel)) skip = .true.
            end do
            if (skip) then
               tmpX(:,i) = lM%x(:,i)
               cycle
            end if
            do e=1, tmpI(i)
               iel = eList(e,i)
               do a=1, lM%eNoN
                  Ac = lM%ien(a,iel)
                  if (Ac .ne. i) then
                     tmpX(:,i) = tmpX(:,i) + lM%x(:,Ac)
                  end if
               end do
            end do
            tmpX(:,i) = tmpX(:,i) / dble(2*tmpI(i))
         end do

         do Ac=1, lM%nNo
            lM%x(:,Ac) = smoother%alpha*tmpX(:,Ac) + &
               (1.0d0 - smoother%alpha)*lM%x(:,Ac)
         end do
      end do
      
      if (smoother%icorr .gt. 0) then
         ds0 = dsqrt(lM%sarea/dble(lM%nEl))
         write(stdout,ftab3) "Correction for reduction in volume "// &
            "during smoothing.."
         tmpX = lM%x
         do
            call calcVolume(lM)
            vol = lM%vol
            if (dabs(vol0-vol) .lt. smoother%tol) exit
            ds = ds0*(vol0-vol) / (3D0*vol0)
            do i=1, lM%nNo
               skip = .false.
               do e=1, tmpI(i)
                  iel = eList(e,i)
                  if (eFixed(iel) .or. i.gt.lM%faces(1)%nNo) &
                     skip = .true.
               end do
               if (skip) then
                  tmpX(:,i) = lM%x(:,i)
                  cycle
               end if
               dir(:) = 0D0
               do e=1, tmpI(i)
                  iel = eList(e,i)
                  if (eFixed(iel)) cycle
                  dir(:) = dir(:) + lM%norm(:,iel)
               end do
               dir(:) = dir(:) / dble(tmpI(i))
               tmpX(:,i) = lM%x(:,i) + ds*dir(:)
            end do
            lM%x = tmpX
         end do
      end if
      
      call calcVolume(lM)
      write(stdout,ftab3) "Final volume: "//STR(lM%vol)
      lM%ien(:,:) = lM%ien(:,:) - 1
      
      deallocate(tmpX,tmpI,eList,eFixed)
      
      end subroutine smoothMsh
      
!**************************************************

      subroutine calcVolume(lM)
      use var
      implicit none
      type(meshType), intent(inout) :: lM
      
      integer :: Ac, a, e
      double precision :: elemArea, magN
      double precision :: u(nsd), v(nsd), elemCent(nsd)
      
      call calcNormalsDir(lM)
      
      lM%sarea = 0D0
      lM%vol   = 0D0
      do e=1, lM%nEl
         u(:) = lM%x(:,lM%ien(2,e)) - lM%x(:,lM%ien(1,e))
         v(:) = lM%x(:,lM%ien(3,e)) - lM%x(:,lM%ien(1,e))
         lM%norm(:,e) = CROSS(u, v) 
         magN = dsqrt( sum( (lM%norm(:,e))**2 ) )
         lM%norm(:,e) = lM%norm(:,e) / magN
         
         elemCent(:) = 0D0
         do a=1, lM%eNoN
            Ac = lM%ien(a,e)
            elemCent(:) = elemCent(:) + lM%x(:,Ac)
         end do
         elemCent(:) = elemCent(:) / dble(lM%eNoN)
         
         lM%sarea = lM%sarea + magN/2D0
         lM%vol   = lM%vol + (sum(elemCent(:)*lM%norm(:,e))*magN/2D0)
      end do
      lM%vol = lM%vol/3D0
      
      end subroutine calcVolume

!**************************************************
   
      subroutine calcNormalsDir(lM)
      use var
      implicit none
      type(meshType), intent(inout) :: lM
      
      integer :: e, a, Ac, iel, itmp
      double precision :: s, sMin
      double precision, dimension(nsd) :: u, v, w, refV
      
      refV(:) = -1D1
      sMin = huge(sMin)
      do a=1, lM%nNo
         s = dsqrt( sum( (lM%x(:,a)-refV(:))**2 ) )
         if (s .lt. sMin) then
            sMin = s
            Ac = a
         end if
      end do
      
      do e=1, lM%nEl
          if ( (lM%ien(1,e) .eq. Ac) .or. &
               (lM%ien(2,e) .eq. Ac) .or. &
               (lM%ien(3,e) .eq. Ac) ) then
             iel = e
             exit
          end if
      end do
      
      u(:) = lM%x(:,lM%ien(2,iel)) - lM%x(:,lM%ien(1,iel))
      v(:) = lM%x(:,lM%ien(3,iel)) - lM%x(:,lM%ien(1,iel))
      w = CROSS(u,v)
      
      u(:) = 0D0
      do a=1, lM%eNoN
         u(:) = u(:) + lM%x(:,lM%ien(a,e))
      end do
      u(:) = u(:) / dble(lM%eNoN)
      s = sum( (u(:) - refV(:))*w(:) )

      if(s .ge. 0D0) then
         do e=1, lM%nEl
            itmp = lM%ien(2,e)
            lM%ien(2,e) = lM%ien(1,e)
            lM%ien(1,e) = itmp
         end do
      end if
      
      end subroutine calcNormalsDir

!**************************************************

      subroutine writeVTK(nNo, eNoN, nEl, x, ien, fName)
      use var
      implicit none
      integer, intent(in) :: nNo, eNoN, nEl, ien(eNoN, nEl)
      double precision, intent(in) :: x(nsd, nNo)
      character(len=*), intent(in) :: fName
      integer :: i, Ac, a, e, fid

      fid = 100
      if ( isBinary ) then
         open(fid,file=trim(fName),status='replace',access='stream',&
         form='unformatted',convert='big_endian')
         write(fid) '# vtk DataFile Version 3.0'//newl
         write(fid) 'Simulation Results'//newl
         write(fid) 'BINARY'//newl
         write(fid) 'DATASET UNSTRUCTURED_GRID'//newl
         write(fid) 'POINTS '//trim(STR(nNo))//' double'//newl
         do Ac=1, nNo
            write(fid) x(:,Ac)
         end do
         write(fid) 'CELLS '//trim(STR(nEl))//' '// &
            trim(STR(4*nEl))//newl
         do e=1, nEl
            write(fid) eNoN, ien(:,e)
         end do
         write(fid) 'CELL_TYPES '//trim(STR(nEl))//newl
         do e=1, nEl
            write(fid) 5
         end do
         close(fid)
      else
         open(fid,file=trim(fName),status='replace')
         write(fid,'(A)') '# vtk DataFile Version 3.0'
         write(fid,'(A)') 'Simulation Results'
         write(fid,'(A)') 'ASCII'
         write(fid,'(A)') 'DATASET UNSTRUCTURED_GRID'
         write(fid,'(A)') 'POINTS '//trim(STR(nNo))//' double'
         do Ac=1, nNo
            do i=1, nsd
               write(fid,'(A)',advance='no') STR(x(i,Ac))//' '
            end do
            write(fid,'(A)')
         end do
         write(fid,'(A)') &
            'CELLS '//trim(STR(nEl))//' '//trim(STR(4*nEl))
         do e=1, nEl
            write(fid,'(A)',advance='no') '3'
            do a=1, eNoN
               write(fid,'(A)',advance='no') ' '//trim(STR(ien(a,e)))
            end do
            write(fid,'(A)')
         end do
         write(fid,'(A)') 'CELL_TYPES '//trim(STR(nEl))
         do e=1, nEl
            write(fid,'(A)') '5'
         end do
         close(fid)
      end if

      end subroutine writeVTK
      
!**************************************************

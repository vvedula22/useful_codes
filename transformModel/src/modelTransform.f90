!**************************************************

      module variables
      use vtkXMLMod
      use vtkLegacyMod
      implicit none

      integer, parameter :: nsd=3

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
      end type faceType

      type imageType
         integer :: iLims(2,3)
         integer :: pLims(2,3)
         real(kind=8) :: origin(3)
         real(kind=8) :: spacng(3)
         integer(IK), dimension(:), allocatable :: ival
      end type imageType

      interface destroy
         module procedure destroyMesh, destroyFace, destroyImage
      end interface destroy

      real(kind=8) :: sX2R(4,4), tX2R(4,4), tR2X(4,4)

      contains
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

      program modelTrans
      use MATFUN
      use variables
      implicit none

      type(faceType)  :: fa
      integer :: fid, i, j
      character(len=strl) :: fname, stmp
      real(kind=8) :: T(4,4)

      fid = 101
      i = iargc()
      if (i .ne. 1) then
         write(stdout,ftab4) "Error: at least one input argument needed"
         STOP
      end if
      call getarg(1,fname)

      open(fid,file=trim(fname))
      read(fid,*)
      read(fid,*) stmp
      read(fid,*)
      read(fid,*) sX2R
      read(fid,*)
      read(fid,*) tX2R
      close(fid)

      call readVTP(fa, stmp)

      sX2R = TRANSPOSE(sX2R)
      write(*,ftab1) "Transform: XYZ (source) --> RAS"
      write(*,ftab2) "xyz2ras (source) matrix: "
      do i=1, 4
         write(*,ftab3,advance='no')
         do j=1, 4
            write(*,'(A)',advance='no') " "//STR(sX2R(i,j))
         end do
         write(*,'(A)')
      end do

      call transform(fa, sX2R)

      tX2R = TRANSPOSE(tX2R)
      tR2X = MAT_INV(tX2R, 4)
      write(*,ftab2) "ras2xyz matrix: "
      do i=1, 4
         write(*,ftab3,advance='no')
         do j=1, 4
            write(*,'(A)',advance='no') " "//STR(tR2X(i,j))
         end do
         write(*,'(A)')
      end do

      write(*,ftab1) "Transform: RAS --> XYZ (target)"
      call transform(fa, tR2X)

      i = len(trim(stmp))
      fname = stmp(1:i-4)
      write(fname,'(A)') trim(fname)//"_transformed.vtp"
      call writeVTP(fa, fname)

      call destroy(fa)

      return
      end program modelTrans

!**************************************************

      subroutine readVTP(lFa, fname)
      use variables
      implicit none
      type(faceType), intent(inout) :: lFa
      character(len=strL), intent(in) :: fname

      type(vtkXMLType) :: vtp
      integer :: istat
      integer, allocatable :: tmpI(:)

      write(*,ftab1) '-------------------------------------------------'
      write(*,ftab1) "Reading model from <"//trim(fname)//">"
      istat = 0
      do while (istat .eq. 0)
         call loadVTK(vtp, trim(fname), istat)
         if (istat .lt. 0) exit

         call getVTK_numPoints(vtp, lFa%nNo, istat)
         if (istat .lt. 0) exit

         call getVTK_numElems(vtp, lFa%nEl, istat)
         if (istat .lt. 0) exit

         call getVTK_nodesPerElem(vtp, lFa%eNoN, istat)
         if (istat .lt. 0) exit

         allocate(lFa%x(nsd,lFa%nNo),  lFa%IEN(lFa%eNoN,lFa%nEl))

         call getVTK_pointCoords(vtp, lFa%x, istat)
         if (istat .lt. 0) exit

         call getVTK_elemIEN(vtp, lFa%ien, istat)
         if (istat .lt. 0) exit

         allocate(tmpI(lFa%nNo))
         call getVTK_pointData(vtp, 'GlobalNodeID', tmpI, istat)
         if (istat .eq. 0) then
            allocate(lFa%gN(lFa%nNo))
            lFa%gN = tmpI
         else
            istat = 0
         end if
         deallocate(tmpI)

         allocate(tmpI(lFa%nEl))
         call getVTK_elemData(vtp, 'GlobalElementID', tmpI, istat)
         if (istat .eq. 0) then
            allocate(lFa%gE(lFa%nEl))
            lFa%gE = tmpI
         else
            istat = 0
         end if
         deallocate(tmpI)

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

      subroutine writeVTP(lFa, fname)
      use variables
      implicit none
      type(faceType), intent(inout) :: lFa
      character(len=strL), intent(in) :: fname

      type(vtkXMLType) :: vtp
      integer :: istat

      write(*,ftab1) "Writing transformed model to <"//trim(fname)//">"
      write(*,ftab1) '-------------------------------------------------'
      if (nsd .EQ. 3) then
         select case (lFa%eNoN)
         case(3)
            lFa%vtkType = 5   ! Tri !
         case(4)
            lFa%vtkType = 9   ! Quad !
         end select
      else
         select case (lFa%eNoN)
         case(2)
            lFa%vtkType = 3   ! Line !
         case(3)
            lFa%vtkType = 21   ! Quadr-Edge !
         end select
      end if

      istat = 0
      do while (istat .eq. 0)
         call vtkInitWriter(vtp, fname, istat)
         if (istat .lt. 0) exit

         call putVTK_pointCoords(vtp, lFa%x, istat)
         if (istat .lt. 0) exit

         call putVTK_elemIEN(vtp, lFa%IEN, lFa%vtkType, istat)
         if (istat .lt. 0) exit

         if (allocated(lFa%gN)) then
            call putVTK_pointData(vtp,"GlobalNodeID",lFa%gN,istat)
            if (istat .lt. 0) exit
         end if

         if (allocated(lFa%gE)) then
            call putVTK_elemData(vtp,"GlobalElementID",lFa%gE,istat)
            if (istat .lt. 0) exit
         end if

         call vtkWriteToFile(vtp, istat)

         call flushVTK(vtp)
         exit
      end do

      if (istat .lt. 0) then
         write(*,ftab4) "ERROR: VTP file write error"
         STOP
      end if

      return
      end subroutine writeVTP

!**************************************************

      subroutine transform(lFa, T)
      use variables
      implicit none
      type(faceType), intent(inout) :: lFa
      real(kind=8), intent(in) :: T(4,4)

      integer :: a, nNo
      real(kind=8), allocatable :: x(:,:)

      nNo = lFa%nNo
      allocate(x(nsd,nNo))
      x = lFa%x
      do a=1, lFa%nNo
         lFa%x(1,a) = T(1,1)*x(1,a) + T(1,2)*x(2,a) + T(1,3)*x(3,a) + T(1,4)
         lFa%x(2,a) = T(2,1)*x(1,a) + T(2,2)*x(2,a) + T(2,3)*x(3,a) + T(2,4)
         lFa%x(3,a) = T(3,1)*x(1,a) + T(3,2)*x(2,a) + T(3,3)*x(3,a) + T(3,4)
      end do
      deallocate(x)

      return
      end subroutine transform

!**************************************************


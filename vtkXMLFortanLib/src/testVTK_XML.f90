!**************************************************

      module variables
      use vtkXMLMod
      use vtkLegacyMod
      implicit none

!      integer(IK), parameter :: nsd=3
      integeR(IK) :: nsd

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

      integer(IK) :: istat

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

      program testVTK_XML
      use variables
      implicit none

      type(mshType)   :: msh
      type(faceType)  :: fa
      type(imageType) :: im

      !call readVTK(msh)

      !call writeVTK(msh) ! not yet implemented !

      nsd = 3
      call readVTU(msh)

      call writeVTU(msh)

      call readVTP(fa)

      call writeVTP(fa)

      call readVTI(im)

      call writeVTI(im)

      call destroy(fa)

      nsd = 2
      call readVTP2D(fa)
      call writeVTP2D(fa)

      call destroy(msh)
      call destroy(fa)
      call destroy(im)

      return
      end program testVTK_XML

!**************************************************

      subroutine readVTK
      use variables
      implicit none
      type(vtkUnstrucGridType) :: vtk

      call loadLegacyVTK(vtk,"in.vtk",istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file read error.."
            STOP
         end if

      return
      end subroutine readVTK

!**************************************************

      subroutine readVTU(lM)
      use variables
      implicit none
      type(mshType), intent(inout) :: lM
      type(vtkXMLType) :: vtu

      write(stdout,'(A)') '-----------------------'
      write(stdout,ftab1) 'Reading file <in.vtu>'

      call loadVTK(vtu,"in.vtu",istat)
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

      allocate(lM%IEN(lM%eNoN,lM%nEl))
      allocate(lM%x(nsd,lM%nNo))

      call getVTK_pointCoords(vtu,lM%x,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file read error.."
            STOP
         end if

      call getVTK_elemIEN(vtu,lM%IEN,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file read error.."
            STOP
         end if

      call flushVTK(vtu)

      return
      end subroutine readVTU

!**************************************************

      subroutine writeVTU(lM)
      use variables
      implicit none
      type(mshType), intent(inout) :: lM
      type(vtkXMLType) :: vtu

      write(stdout,ftab1) 'Writing file <out.vtu>'
      if (nsd .EQ. 3) then
         select case (lM%eNoN)
         case(4)
            lM%vtkType = 10  ! Tet !
         case(6)
            lM%vtkType = 13  ! Wedge !
         case(8)
            lM%vtkType = 12  ! Brick/Hex !
         end select
      else
         select case (lM%eNoN)
         case(3)
            lM%vtkType = 5   ! Tri !
         case(4)
            lM%vtkType = 9   ! Bilinear/Quad !
         case(9)
            lM%vtkType = 28  ! Biquad/Quadr-Quad !
         end select
      end if

      call vtkInitWriter(vtu,"out.vtu",istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file write error.."
            STOP
         end if

      call putVTK_pointCoords(vtu,lM%x,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file write error.."
            STOP
         end if

      call putVTK_elemIEN(vtu,lM%IEN,lM%vtkType,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file write error.."
            STOP
         end if

      call vtkWriteToFile(vtu,istat)

      call flushVTK(vtu)
      write(stdout,'(A)') '-----------------------'

      return
      end subroutine writeVTU

!**************************************************

      subroutine readVTP(lFa)
      use variables
      implicit none
      type(faceType), intent(inout) :: lFa
      type(vtkXMLType) :: vtp

      write(stdout,'(A)') '-----------------------'
      write(stdout,ftab1) 'Reading file <in.vtp>'

      call loadVTK(vtp,"in.vtp",istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file read error.."
            STOP
         end if

      call getVTK_numPoints(vtp,lFa%nNo,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file read error.."
            STOP
         end if

      call getVTK_numElems(vtp,lFa%nEl,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file read error.."
            STOP
         end if

      call getVTK_nodesPerElem(vtp,lFa%eNoN,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file read error.."
            STOP
         end if

      allocate(lFa%IEN(lFa%eNoN,lFa%nEl))
      allocate(lFa%x(nsd,lFa%nNo))
      allocate(lFa%gN(lFa%nNo), lFa%gE(lFa%nEl))

      call getVTK_pointCoords(vtp,lFa%x,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file read error.."
            STOP
         end if

      call getVTK_elemIEN(vtp,lFa%ien,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file read error.."
            STOP
         end if

      call getVTK_pointData(vtp,'GlobalNodeID',lFa%gN,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file read error.."
            STOP
         end if

      call getVTK_elemData(vtp,'GlobalElementID',lFa%gE,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file read error.."
            STOP
         end if

      call flushVTK(vtp)

      return
      end subroutine readVTP

!**************************************************

      subroutine writeVTP(lFa)
      use variables
      implicit none
      type(faceType), intent(inout) :: lFa
      type(vtkXMLType) :: vtp

      write(stdout,ftab1) 'Writing file <out.vtp>'

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

      call vtkInitWriter(vtp,"out.vtp",istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file write error.."
            STOP
         end if

      call putVTK_pointCoords(vtp,lFa%x,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file write error.."
            STOP
         end if

      call putVTK_elemIEN(vtp,lFa%IEN,lFa%vtkType,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file write error.."
            STOP
         end if

      call putVTK_pointData(vtp,"GlobalNodeID",lFa%gN,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file write error.."
            STOP
         end if

      call putVTK_elemData(vtp,"GlobalElementID",lFa%gE,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file write error.."
            STOP
         end if

      call vtkWriteToFile(vtp,istat)

      call flushVTK(vtp)
      write(stdout,'(A)') '-----------------------'

      return
      end subroutine writeVTP

!**************************************************

      subroutine readVTP2D(lFa)
      use variables
      implicit none
      type(faceType), intent(inout) :: lFa
      type(vtkXMLType) :: vtp

      real(kind=8), allocatable :: tmpX(:,:)

      write(stdout,'(A)') '-----------------------'
      write(stdout,ftab1) 'Reading file <in2D.vtp>'

      call loadVTK(vtp,"in2D.vtp",istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file read error.."
            STOP
         end if

      call getVTK_numPoints(vtp,lFa%nNo,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file read error.."
            STOP
         end if

      call getVTK_numElems(vtp,lFa%nEl,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file read error.."
            STOP
         end if

      call getVTK_nodesPerElem(vtp,lFa%eNoN,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file read error.."
            STOP
         end if

      allocate(lFa%IEN(lFa%eNoN,lFa%nEl))
      allocate(lFa%x(nsd,lFa%nNo), tmpX(maxNSD,lFa%nNo))
      allocate(lFa%gN(lFa%nNo), lFa%gE(lFa%nEl))

      call getVTK_pointCoords(vtp,tmpX,istat)
      lFa%x(:,:) = tmpX(1:nsd,:)
      deallocate(tmpX)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file read error.. (pcoord)"
            STOP
         end if

      call getVTK_elemIEN(vtp,lFa%ien,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file read error.. (IEN)"
            STOP
         end if

      call getVTK_pointData(vtp,'GlobalNodeID',lFa%gN,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file read error.. (NodeID)"
            STOP
         end if

      call getVTK_elemData(vtp,'GlobalElementID',lFa%gE,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file read error.. (CellID)"
            STOP
         end if

      call flushVTK(vtp)

      return
      end subroutine readVTP2D

!**************************************************

      subroutine writeVTP2D(lFa)
      use variables
      implicit none
      type(faceType), intent(inout) :: lFa
      type(vtkXMLType) :: vtp

      write(stdout,ftab1) 'Writing file <out2D.vtp>'

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

      call vtkInitWriter(vtp,"out2D.vtp",istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file write error.."
            STOP
         end if

      call putVTK_pointCoords(vtp,lFa%x,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file write error.."
            STOP
         end if

      call putVTK_elemIEN(vtp,lFa%IEN,lFa%vtkType,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file write error.."
            STOP
         end if

      call putVTK_pointData(vtp,"GlobalNodeID",lFa%gN,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file write error.."
            STOP
         end if

      call putVTK_elemData(vtp,"GlobalElementID",lFa%gE,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file write error.."
            STOP
         end if

      call vtkWriteToFile(vtp,istat)

      call flushVTK(vtp)
      write(stdout,'(A)') '-----------------------'

      return
      end subroutine writeVTP2D

!**************************************************

      subroutine readVTI(lIm)
      use variables
      implicit none
      type(imageType), intent(inout) :: lIm
      type(vtkXMLType) :: vti
      integer :: n

      write(stdout,'(A)') '-----------------------'
      write(stdout,ftab1) 'Reading file <in.vti>'

      call loadVTK(vti,"in.vti",istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file read error.."
            STOP
         end if

      call getVTK_imageExtent(vti,lIm%iLims,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file read error.."
            STOP
         end if

      call getVTK_imageOrigin(vti,lIm%origin,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file read error.."
            STOP
         end if

      call getVTK_imageSpacing(vti,lIm%spacng,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file read error.."
            STOP
         end if

      call getVTK_pieceExtent(vti,lIm%pLims,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file read error.."
            STOP
         end if

      n = (lIm%pLims(2,1)-lIm%pLims(1,1)+1) * &
          (lIm%pLims(2,2)-lIm%pLims(1,2)+1) * &
          (lIm%pLims(2,3)-lIm%pLims(1,3)+1)
      allocate(lIm%ival(n))
      call getVTK_pointData(vti,'ImageScalars',lIm%ival,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file read error.."
            STOP
         end if

      call flushVTK(vti)

      return
      end subroutine readVTI

!**************************************************

      subroutine writeVTI(lIm)
      use variables
      implicit none
      type(imageType), intent(inout) :: lIm
      type(vtkXMLType) :: vti

      write(stdout,ftab1) 'Writing file <out.vti>'

      call vtkInitWriter(vti,"out.vti",istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file write error.."
            STOP
         end if

      call putVTK_imageExtent(vti,lIm%iLims,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file write error.."
            STOP
         end if

      call putVTK_imageOrigin(vti,lIm%origin,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file write error.."
            STOP
         end if

      call putVTK_imageSpacing(vti,lIm%spacng,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file write error.."
            STOP
         end if

      call putVTK_pieceExtent(vti,lIm%pLims,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file write error.."
            STOP
         end if

      call putVTK_pointData(vti,'ImageScalars',lIm%ival,istat)
         if ( istat.lt.0 ) then
            write(stdout,ftab4) "ERROR: VTK file write error.."
            STOP
         end if

      call vtkWriteToFile(vti,istat)

      call flushVTK(vti)
      write(stdout,'(A)') '-----------------------'

      return
      end subroutine writeVTI

!**************************************************


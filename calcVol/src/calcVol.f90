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
      
      type meshType
         integer :: nNo, nEl, eNoN
         integer, dimension(:,:), allocatable :: ien
         real(kind=8), dimension(:,:), allocatable :: x
         real(kind=8), dimension(:,:), allocatable :: norm
         real(kind=8) :: vol, sarea
      end type meshType
      
      type(meshType) :: msh
      
      contains

        function CROSS(u,v) result(w)
		implicit none
		double precision, intent(in) :: u(:),v(:)
		double precision :: w(size(u,1))
		integer :: m
        
		m = size(u,1)
		if (m .ne. size(v,1)) then
			write(stdout,'(14X,A)') &
		   "Error: inconsistent vector dimensions for cross product"
		end if
		 
		w(1) = u(2)*v(3) - u(3)*v(2)
		w(2) = u(3)*v(1) - u(1)*v(3)
		w(3) = u(1)*v(2) - u(2)*v(1)
		 
		end function CROSS
        
      end module var

!**************************************************

      program vol
      use var
      implicit none
      integer :: i, Ac, istat, fid
      character(len=strL) :: fName
      
      i = IARGC()
      if (i .eq. 0) then
         write(stdout,ftab4) "ERROR: Input file name not specified"
         STOP
      else if (i .gt. 1) then
         write(stdout,ftab4) "ERROR: Too many arguments"
         STOP
      end if
      
      call getarg(1,fName)
      
      write(stdout,ftab1) "Loading mesh... <"//TRIM(fName)//">"
      call loadVTP(msh)
      
      write(stdout,ftab1) "Computing mesh volume..."
      call calcVolume(msh)
      write(stdout,ftab2) "Volume: "//TRIM(STR(msh%vol))
      
      contains
      
         !==========================================
         
         subroutine loadVTP(lM)
         implicit none
         integer :: i,j
         type(meshType), intent(inout) :: lM
         type(vtkXMLType) :: vtp
         
         call loadVTK(vtp,fName,istat)
            if ( istat.lt.0 ) then
               write(stdout,ftab4) "ERROR: VTK file read error.."
               STOP
            end if
         
         call getVTK_numPoints(vtp,lM%nNo,istat)
            if ( istat.lt.0 ) then
               write(stdout,ftab4) "ERROR: VTK file read error.."
               STOP
            end if
         
         call getVTK_numElems(vtp,lM%nEl,istat)
            if ( istat.lt.0 ) then
               write(stdout,ftab4) "ERROR: VTK file read error.."
               STOP
            end if
         
         call getVTK_nodesPerElem(vtp,lM%eNoN,istat)
            if ( istat.lt.0 ) then
               write(stdout,ftab4) "ERROR: VTK file read error.."
               STOP
            end if
         
         allocate(lM%ien(lM%eNoN,lM%nEl))
         allocate(lM%x(nsd,lM%nNo))
         
         call getVTK_pointCoords(vtp,lM%x,istat)
            if ( istat.lt.0 ) then
               write(stdout,ftab4) "ERROR: VTK file read error.."
               STOP
            end if
         
         call getVTK_elemIEN(vtp,lM%ien,istat)
            if ( istat.lt.0 ) then
               write(stdout,ftab4) "ERROR: VTK file read error.."
               STOP
            end if
         lM%ien = lM%ien + 1
         
         call flushVTK(vtp)
         
         end subroutine loadVTP
         
      end program vol
      
!**************************************************

      subroutine calcVolume(lM)
      use var
      implicit none
      type(meshType), intent(inout) :: lM
      
      integer :: Ac, a, e
      double precision :: elemArea, magN
      double precision :: u(nsd), v(nsd), elemCent(nsd)
      
      if (.not.allocated(lM%norm)) allocate(lM%norm(nsd,lM%nEl))
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

!**************************************************

      module variables
      integer, parameter :: nDim=3
      integer, parameter :: vtkCellType=10  ! 10-Tet; 5-tri !
      character(len=4), parameter :: outType='VTKB'

      integer, parameter :: stdl=250
      integer, parameter :: maxDataType = 10
      character, parameter :: eol=char(10)

      integer, parameter :: maxPtDOF=24
      integer, parameter :: maxElDOF=6
      integer, parameter :: maxNSD=3
      
      type outPtData
         integer :: nVar,tDOF
         integer, dimension(maxPtDOF) :: ndof,ioff,ikind
         character(len=stdl), dimension(maxPtDOF) :: varName
         character(len=stdl), dimension(maxPtDOF) :: varType
         real*4, dimension(:), allocatable :: arr
      end type outPtData
      
      type outElData
         integer :: nVar,tDOF
         integer, dimension(maxElDOF) :: ndof,ioff,ikind
         character(len=stdl), dimension(maxElDOF) :: varName
         character(len=stdl), dimension(maxElDOF) :: varType
         real*4, dimension(:), allocatable :: arr
      end type outElData

      type mesh
         integer :: nodes,elems,nodesPerElem,nGauss
         integer, dimension(:,:), allocatable :: connec
         double precision, dimension(:,:), allocatable :: coords
         
         double precision :: tMshInt
         double precision, dimension(:), allocatable :: wG
         double precision, dimension(:), allocatable :: J,skewness,AR
         double precision, dimension(:,:), allocatable :: xi,N
         double precision, dimension(:,:,:), allocatable :: Nx
         
         logical :: isZPresent
      end type mesh
      
      type(mesh) :: iMesh
      type(outPtData) :: ptData
      type(outElData) :: elData
      
      integer :: scalarDOF,vectorDOF

      end module variables
      
!**************************************************

      module allFun
      implicit none
      
      interface STR
         module procedure :: ITSTR, RTSTR, DTSTR, NDTSTR
      end interface STR
      
      contains

         double precision function Jacobian(nsd,eNoN,x,Nxi)
         implicit none
         integer, intent(in) :: nsd,eNoN
         double precision, intent(in) :: x(nsd,eNoN),Nxi(nsd,eNoN)
         
         integer :: a
         double precision :: xXi(nsd,nsd)

         xXi = 0.0d0
         if ( nsd.eq.2 ) then
            do a=1, eNoN
               xXi(:,1) = xXi(:,1) + x(:,a)*Nxi(1,a)
               xXi(:,2) = xXi(:,2) + x(:,a)*Nxi(2,a)
            end do
         else if ( nsd.eq.3) then
            do a=1, eNoN
               xXi(:,1) = xXi(:,1) + x(:,a)*Nxi(1,a)
               xXi(:,2) = xXi(:,2) + x(:,a)*Nxi(2,a)
               xXi(:,3) = xXi(:,3) + x(:,a)*Nxi(3,a)
            end do
         end if
         Jacobian = det(xXi,nsd)
         
         return
         end function Jacobian
         
!**************************************************
         
         recursive double precision function det(A,nd) result(detr)
         implicit none
         integer, intent(in) :: nd
         double precision, intent(in) :: A(nd,nd)
         integer :: i,j,m,n
         double precision :: Ared(nd-1,nd-1)
         
         detr = 0.0d0
         if ( nd.eq.2 ) then
            detr = A(1,1)*A(2,2) - A(1,2)*A(2,1)
            return
         else
            do i=1, nd
               n=0
               do j=1, nd
                  if ( j.eq.i ) then
                     cycle
                  else
                     n = n+1
                     Ared(:,n) = A(2:nd,j)
                  end if
               end do
               detr = detr + (-1.0d0)**dble(1+i) * &
                     A(1,i) * det(Ared,nd-1)
            end do
         end if
         
         return
         end function det

!**************************************************

         pure function ITSTR(iVal) result(str)
         implicit none
         integer, intent(in) :: iVal
         integer :: n,ist,j,k,slen
         character(len=80) str
         
         str = ''
         if ( iVal.lt.0 ) then
            str(1:1) = '-'
            ist = 2
         else
            ist = 1
         end if
         
         slen = 2 - max(0,sign(1,iVal)) + &
            maxval( (/ (min(abs(iVal)/10**n,1)*n,n=1, 9) /) )
         
         n = iVal
         do j=slen, ist, -1
            k = modulo(abs(n),10) + 1
            str(j:j) = '0123456789'(k:k)
            n = n/10
         end do
         
         return
         end function ITSTR
         
!**************************************************
         
         pure function RTSTR(rVal) result(str)
         implicit none
         integer, parameter :: l=8
         real*4, intent(in) :: rVal
         character(len=l) :: str
         
         str = NDTSTR(dble(rVal),l)
         
         return
         end function RTSTR
         
!**************************************************
         
         pure function DTSTR(dVal) result(str)
         implicit none
         integer, parameter :: l=8
         double precision, intent(in) :: dVal
         character(len=l) :: str
         
         str = NDTSTR(dVal,l)
         
         return
         end function DTSTR
         
!**************************************************
         
         pure function NDTSTR(dVal,l) result(str)
         implicit none
         integer, intent(in) :: l
         double precision, intent(in) :: dVal
         character(len=l) :: str
         
         integer :: i,j,k,ipos,cnt,ex,abex,nex
         double precision :: absd
         
         ! check NaN !
         if ( dVal.ne.dVal ) then
            if ( l.ge.3) then
               str = ''
               str(l-2:l) = 'NaN'
            else
               str = 'NaN'(1:l)
            end if
            return
         end if
         
         ! check infinity !
         absd = abs(dVal)
         if ( absd.gt.huge(absd) ) then
            if ( l.ge.8 ) then
               str = ''
               str(l-7:l) = 'Infinity'
            else
               str = 'Infinity'(1:l)
            end if
            return
         end if
         
         ! check zero !
         if ( absd.lt.tiny(absd) ) then
            ipos = 1
            if ( l.ge.3 ) then
               str(1:3) = '0.0'
               ipos = 4
            end if
            do i=ipos, l
               str(i:i) = '0'
            end do
            return
         end if
         
         ! count exp and number !
         ex = 0; nex = 0
         if ( absd.ne.0d0 ) ex = floor(log10(absd)) ! exponent !
         abex = abs(ex)
         
         ! no of digits in exponent !
         if ( ex.ne.0 ) nex = floor(log10(real(abex,8))) + 1
         
         cnt = nex+1  ! atleast 1 number before exponent !
         if ( dVal.lt.0d0 ) cnt = cnt+1 ! for - sign of number !
         if ( ex.lt.0 ) cnt = cnt+1 ! for - sign of exponent !
         if ( ex.ne.0 ) cnt = cnt+1 ! for letter E !
         
         if ( l.lt.cnt ) then   ! check total sum !
            do i=1, l
               str(i:i) = '*'
            end do
            return
         end if
         
         if ( ex.ne.0 ) then
            do ipos=l, l-nex+1, -1
               k = modulo(abex,10) + 1
               str(ipos:ipos) = '0123456789'(k:k)
               abex = abex/10
            end do
            if ( ex.lt.0 ) then
               str(ipos:ipos) = '-'
               ipos = ipos-1
            end if
            str(ipos:ipos) = 'E'
            ipos = ipos-1
         else
            ipos = l
         end if
         
         if ( l.gt.cnt ) then
            absd = absd*(1D1**(l-cnt-1-ex))
            do i=ipos, ipos-(l-cnt)+2, -1
               k = modulo(absd,1D1) + 1
               str(i:i) = '0123456789'(k:k)
               absd = absd/1D1
            end do
            ipos = i
            str(ipos:ipos) = '.'
            ipos = ipos-1
            k = modulo(absd,1D1) + 1
            str(ipos:ipos) = '0123456789'(k:k)
         else ! l.eq.cnt
            absd = absd*(1D1**(l-cnt-ex))
            k = floor(modulo(absd,1D1)) + 1
            str(ipos:ipos) = '0123456789'(k:k)
         end if
         if ( dVal.lt.0.0d0 ) str(1:1) = '-'

         return
         end function NDTSTR

!**************************************************

         pure function vInteg(nsd,eNoN,mat) result(vint)
         implicit none
         integer, intent(in) :: nsd,eNoN
         double precision, intent(in) :: mat(nsd,eNoN-1)
         double precision :: vint
         
         if ( nsd.eq.2 ) then
            vint =  mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1)
            vint =  vint/2.0d0
         else
            vint =  mat(1,1)*( mat(2,2)*mat(3,3) - mat(3,2)*mat(2,3) ) + &
                  mat(2,1)*( mat(3,2)*mat(1,3) - mat(1,2)*mat(3,3) ) + &
                  mat(3,1)*( mat(1,2)*mat(2,3) - mat(2,2)*mat(1,3) )
            vint =  vint/6.0d0
         end if
         
         end function vInteg
         
      end module allFun

!**************************************************

      program vtkConv
      use variables
      use allFun
      implicit none
      character(len=stdl) :: stmp,fName
      integer :: ist,iend,ifreq
      integer :: i,ntime,nframe
      
      open(10,file="inputs.dat")
      read(10,*)
      read(10,*) ist, iend, ifreq
      read(10,*)
      read(10,*) stmp
      close(10)
      
      nframe = (iend-ist)/ifreq + 1
      
      if ( vtkCellType.eq.5 ) then ! tri !
         iMesh%nodesPerElem = 3
         iMesh%nGauss = 3
      else if (vtkCellType.eq.10 ) then ! tet !
         iMesh%nodesPerElem = 4
         iMesh%nGauss = 4
      else
         write(*,'(2X,A)') 'Unknown VTK cell type..'
         STOP
      end if
      
      iMesh%isZpresent = .true.
      scalarDOF = 1
      vectorDOF = maxNSD

      do i=1, nframe
         ntime = ist + (i-1)*ifreq
         
         write(fname,'(A)') trim(stmp)//trim(STR(ntime))//".vtk"
         write(*,*)
         write(*,*) '*********************************************'
         write(*,*)
         write(*,'(2X,A)') "Reading file <--- "//trim(fName)
         write(*,*)
         call readVTK(fname)
         
         allocate(iMesh%wG(iMesh%nGauss))
         allocate(iMesh%xi(nDim,iMesh%nodesPerElem))
         allocate(iMesh%N (iMesh%nodesPerElem,iMesh%nGauss))
         allocate(iMesh%Nx(nDim,iMesh%nodesPerElem,iMesh%nGauss))
         allocate(iMesh%J(iMesh%elems))
         allocate(iMesh%skewness(iMesh%elems))
         allocate(iMesh%AR(iMesh%elems))
         
         call xigauss(nDim,iMesh%nodesPerElem,iMesh%nGauss, &
            iMesh%wG,iMesh%xi,iMesh%N,iMesh%Nx)
         
         call calcMeshProps
         
         write(fname,'(A)') trim(stmp)//"_b_"//trim(STR(ntime))//".vtk"
         write(*,*)
         write(*,'(2X,A)') "Writing file ---> "//trim(fName)
         write(*,*)
         call writeVTK(fname)
         
         call freeMem
      end do
      
      end program vtkConv

!**************************************************

      subroutine readVTK(fname)
      use variables
      use allFun
      implicit none
      character(len=*), intent(in) :: fname
      integer :: i,j,k,fid,iPos
      logical :: isBinary,flag
      character(len=stdl) :: rLine
      character :: c
      
      integer :: nn,ne
      integer, allocatable, dimension(:,:) :: connec
      real*4,  allocatable, dimension(:,:) :: coords
      
      integer :: ivar,ist,iend,vecl,dtype
      integer :: itmp
      real*4  :: rtmp
      integer, dimension(:), allocatable :: tmpI
      
      write(*,'(2X,A)') 'Reading vtk file..'
      isBinary = .false.
      fid = 10
      
      inquire(file=trim(fName), exist=flag)
      if ( .not.flag ) then
         write(*,'(2X,3A)') 'File ',trim(fName),' does not exist'
         STOP
      end if
      
      open(fid,file=trim(fName),status='old')
      read(fid,*) ! # vtk DataFile Version 3.0
      read(fid,*) ! Simulation Results < dummy hdr >
      read(fid,'(A)') rLine ! ASCII or BINARY
      
      rLine = adjustl(rLine)
      if ( rLine(1:6) .eq. 'BINARY' ) then
         isBinary = .true.
         close(fid)
         open(fid,file=trim(fName),status='unknown',access='stream',&
         form='unformatted',convert='big_endian')
         iPos = 1
      end if
      rewind(fid)

      write(*,'(4X,A)') 'VTK properties:'
      if ( isBinary ) then
         write(*,'(8X,3A)') &
         'VTK file ',trim(fName),' is in BINARY format'
      else
         write(*,'(8X,3A)') &
         'VTK file ',trim(fName),' is in ASCII Format..'
      end if
      
      ! points !
      call findKwrd(fid,'POINTS',iPos,rLine)
      rLine = adjustl(rLine(7:))
      i = len(trim(rLine))
      rLine = rLine(1:i-6)  ! skipping 'float'
      read(rLine,*) nn
      
      write(*,'(8X,2A)') 'No of nodes: ',trim(STR(nn))
      allocate(coords(nDim,nn)); coords = 0.0
      
      if ( isBinary ) then
         do i=1, nn
            read(fid,pos=iPos) coords(:,i)
            if ( iMesh%isZPresent .and. nDim.eq.2 ) then
               iPos = iPos + maxNSD*kind(rtmp)
            else
               iPos = iPos + nDim*kind(rtmp)
            end if
         end do
      else
         do i=1, nn
            read(fid,*) coords(:,i)
         end do
      end if
      
      ! cell connectivity !
      call findKwrd(fid,'CELLS',iPos,rLine)
      rLine = adjustl(rLine(6:))
      read(rLine,*) ne, i
      if ( i.ne. ne*(iMesh%nodesPerElem+1) ) then
         write(*,'(8X,A)') &
         'Error: inconsistent element type and connectivity..'
         write(*,'(12X,2A)') 'No of elements: ',trim(STR(ne))
         write(*,'(12X,2A)') 'Given nodes per element: ', &
         trim(STR(iMesh%nodesPerElem))
         write(*,'(12X,2A)') 'No of values to be read: ', &
         trim(STR(ne*(iMesh%nodesPerElem+1)))
         write(*,'(12X,2A)') 'No of values (VTK): ', &
         trim(STR(i))
         STOP
      end if
      
      write(*,'(8X,2A)') 'No of elems: ',trim(STR(ne))
      allocate(connec(iMesh%nodesPerElem,ne)); connec = 0
      
      if ( isBinary ) then
         do i=1, ne
            read(fid,pos=iPos) itmp, connec(:,i)
            iPos = iPos + (iMesh%nodesPerElem+1)*kind(itmp)
         end do
      else
         do i=1, ne
            read(fid,*) itmp, connec(:,i)
         end do
      end if
      
      ! reset connectivity order from 1 !
      do i=1, ne
         connec(:,i) = connec(:,i) + 1
      end do
      
      call findKwrd(fid,'CELL_TYPES',iPos,rLine)
      if ( isBinary ) then
         read(fid,pos=iPos) itmp
      else
         read(fid,*) itmp
      end if
      if ( itmp.ne.vtkCellType ) then
         write(*,'(8X,A)') 'Error: inconsistent cell type'
         write(*,'(8X,2A)') 'Given type: ',trim(STR(vtkCellType))
         write(*,'(8X,2A)') 'VTK cell type in file: ', &
         trim(STR(itmp))
         STOP
      end if
      
      ! data section !
      ptData%nvar = 0
      ptData%tDOF = 0
      ptData%varName(:) = ''
      ptData%varType(:) = ''
      elData%nvar = 0
      elData%tDOF = 0
      elData%varName(:) = ''
      elData%varType(:) = ''
      dType = 0
      do
         rLine = ''
         if ( isBinary ) then
            do i=1, stdl
               read(fid,pos=iPos,end=001) c
               iPos = iPos + 1
               if ( c.eq.eol ) exit
               rLine(i:i) = c
            end do
         else
            read(fid,'(A)',end=001) rLine
         end if
         
         rLine = adjustl(rLine)
         if ( rLine(1:10).eq.'POINT_DATA' ) then
            allocate(ptData%arr(nn*maxPtDOF))
            ptData%arr = 0.0
            dType = 1
            cycle
         else if ( rLine(1:9).eq.'CELL_DATA' ) then
            allocate(elData%arr(ne*maxElDOF))
            elData%arr = 0.0
            dType = 2
            cycle
         end if
         
         if ( dType.eq.1 ) then
            if ( rLine(1:7).eq.'SCALARS' ) then
               ptData%nvar = ptData%nvar + 1
               ptData%ndof(ptData%nvar) = scalarDOF
               ptData%tdof = ptData%tdof + ptData%ndof(ptData%nvar)
               rLine = rLine(9:)
               i = len(trim(rLine))
               do j=1, i
                  if ( rLine(j:j).eq.' ' ) exit
               end do
               ptData%varName(ptData%nvar) = rLine(1:j-1)
               ptData%varType(ptData%nvar) = rLine(j+1:i)
               select case ( trim(ptData%varType(ptData%nvar)) )
               case ('int')
                  ptData%ikind(ptData%nvar) = kind(itmp)
               case ('float')
                  ptData%ikind(ptData%nvar) = kind(rtmp)
               end select
               write(*,'(8X,5A)') 'Reading ',&
               trim(ptData%varName(ptData%nvar)), ' scalars (', &
               trim(ptData%varType(ptData%nvar)), ')'
               
               ist  = (ptData%tdof-ptData%ndof(ptData%nvar))*nn + 1
               vecl =  ptData%ndof(ptData%nvar)*nn
               iend =  ist + vecl - 1
               ptData%ioff(ptData%nvar) = ist-1
               if ( isBinary ) then
                  do i=1, stdl
                     read(fid,pos=iPos) c
                     iPos = iPos+1
                     if(c.eq.eol) exit
                  end do
                  if ( trim(ptData%varType(ptData%nvar)).eq.'int' ) then
                     if ( allocated(tmpI) ) deallocate(tmpI)
                     allocate(tmpI(vecl))
                     read(fid,pos=iPos) tmpI(1:vecl)
                     ptData%arr(ist:iend) = real(tmpI(1:vecl))
                  else if ( trim(ptData%varType(ptData%nvar)).eq.'float' ) then
                     read(fid,pos=iPos) ptData%arr(ist:iend)
                  else
                     write(*,'(8X,A)') 'Unknown data type..'
                     STOP
                  end if
                  iPos = iPos + vecl*ptData%ikind(ptData%nvar)
               else
                  read(fid,'(A)') rLine
                  do i=1, nn
                     read(fid,*) &
                     (ptData%arr(ist-1+(scalarDOF*(i-1)+j)),j=1,scalarDOF)
                  end do
               end if
            else if ( rLine(1:7).eq.'VECTORS' ) then
               ptData%nvar = ptData%nvar + 1
               ptData%ndof(ptData%nvar) = vectorDOF
               ptData%tdof = ptData%tdof + ptData%ndof(ptData%nvar)
               rLine = rLine(9:)
               i = len(trim(rLine))
               do j=1, i
                  if ( rLine(j:j).eq.' ' ) exit
               end do
               ptData%varName(ptData%nvar) = rLine(1:j-1)
               ptData%varType(ptData%nvar) = rLine(j+1:i)
               select case ( trim(ptData%varType(ptData%nvar)) )
               case ('int')
                  ptData%ikind(ptData%nvar) = kind(itmp)
               case ('float')
                  ptData%ikind(ptData%nvar) = kind(rtmp)
               end select
               write(*,'(8X,5A)') 'Reading ',&
               trim(ptData%varName(ptData%nvar)), ' vectors (', &
               trim(ptData%varType(ptData%nvar)), ')'

               ist  = (ptData%tdof-ptData%ndof(ptData%nvar))*nn + 1
               vecl =  ptData%ndof(ptData%nvar)*nn
               iend =  ist + vecl - 1
               ptData%ioff(ptData%nvar) = ist-1
               if ( isBinary ) then
                  if ( trim(ptData%varType(ptData%nvar)).eq.'int' ) then
                     if ( allocated(tmpI) ) deallocate(tmpI)
                     allocate(tmpI(vecl))
                     read(fid,pos=iPos) tmpI(1:vecl)
                     ptData%arr(ist:iend) = real(tmpI(1:vecl))
                  else if ( trim(ptData%varType(ptData%nvar)).eq.'float' ) then
                     read(fid,pos=iPos) ptData%arr(ist:iend)
                  else
                     write(*,'(8X,A)') 'Unknown data type..'
                     STOP
                  end if
                  iPos = iPos + vecl*ptData%ikind(ptData%nvar)
               else
                  do i=1, nn
                     read(fid,*) &
                     (ptData%arr(ist-1+(vectorDOF*(i-1)+j)),j=1,vectorDOF)
                  end do ! i
               end if ! isBinary

            end if ! rLine: scalars/vectors
            
         else if ( dType.eq.2 ) then
            if ( rLine(1:7).eq.'SCALARS' ) then
               elData%nvar = elData%nvar + 1
               elData%ndof(elData%nvar) = scalarDOF
               elData%tdof = elData%tdof + elData%ndof(elData%nvar)
               rLine = rLine(9:)
               i = len(trim(rLine))
               do j=1, i
                  if ( rLine(j:j).eq.' ' ) exit
               end do
               elData%varName(elData%nvar) = rLine(1:j-1)
               elData%varType(elData%nvar) = rLine(j+1:i)
               select case ( trim(elData%varType(elData%nvar)) )
               case ('int')
                  elData%ikind(elData%nvar) = kind(itmp)
               case ('float')
                  elData%ikind(elData%nvar) = kind(rtmp)
               end select
               write(*,'(8X,5A)') 'Reading ',&
               trim(elData%varName(elData%nvar)), ' scalars (', &
               trim(elData%varType(elData%nvar)), ')'
               
               ist  = (elData%tdof-elData%ndof(elData%nvar))*ne + 1
               vecl =  elData%ndof(elData%nvar)*ne
               iend =  ist + vecl - 1
               elData%ioff(elData%nvar) = ist-1
               if ( isBinary ) then
                  do i=1, stdl
                     read(fid,pos=iPos) c
                     iPos = iPos+1
                     if(c.eq.eol) exit
                  end do
                  if ( trim(elData%varType(elData%nvar)).eq.'int' ) then
                     if ( allocated(tmpI) ) deallocate(tmpI)
                     allocate(tmpI(vecl))
                     read(fid,pos=iPos) tmpI(1:vecl)
                     elData%arr(ist:iend) = real(tmpI(1:vecl))
                  else if ( trim(elData%varType(elData%nvar)).eq.'float' ) then
                     read(fid,pos=iPos) elData%arr(ist:iend)
                  else
                     write(*,'(8X,A)') 'Unknown data type..'
                     STOP
                  end if
                  iPos = iPos + vecl*elData%ikind(elData%nvar)
               else
                  read(fid,'(A)') rLine
                  do i=1, ne
                     read(fid,*) &
                     (elData%arr(ist-1+(scalarDOF*(i-1)+j)),j=1,scalarDOF)
                  end do
               end if
            else if ( rLine(1:7).eq.'VECTORS' ) then
               elData%nvar = elData%nvar + 1
               elData%ndof(elData%nvar) = vectorDOF
               elData%tdof = elData%tdof + elData%ndof(elData%nvar)
               rLine = rLine(9:)
               i = len(trim(rLine))
               do j=1, i
                  if ( rLine(j:j).eq.' ' ) exit
               end do
               elData%varName(elData%nvar) = rLine(1:j-1)
               elData%varType(elData%nvar) = rLine(j+1:i)
               select case ( trim(elData%varType(elData%nvar)) )
               case ('int')
                  elData%ikind(elData%nvar) = kind(itmp)
               case ('float')
                  elData%ikind(elData%nvar) = kind(rtmp)
               end select
               write(*,'(8X,5A)') 'Reading ',&
               trim(elData%varName(elData%nvar)), ' vectors (', &
               trim(elData%varType(elData%nvar)), ')'

               ist  = (elData%tdof-elData%ndof(elData%nvar))*ne + 1
               vecl =  elData%ndof(elData%nvar)*ne
               iend =  ist + vecl - 1
               elData%ioff(elData%nvar) = ist-1
               if ( isBinary ) then
                  if ( trim(elData%varType(elData%nvar)).eq.'int' ) then
                     if ( allocated(tmpI) ) deallocate(tmpI)
                     allocate(tmpI(vecl))
                     read(fid,pos=iPos) tmpI(1:vecl)
                     elData%arr(ist:iend) = real(tmpI(1:vecl))
                  else if ( trim(elData%varType(elData%nvar)).eq.'float' ) then
                     read(fid,pos=iPos) elData%arr(ist:iend)
                  else
                     write(*,'(8X,A)') 'Unknown data type..'
                     STOP
                  end if
                  iPos = iPos + vecl*elData%ikind(elData%nvar)
               else
                  do i=1, ne
                     read(fid,*) &
                     (elData%arr(ist-1+(vectorDOF*(i-1)+j)),j=1,vectorDOF)
                  end do ! i
               end if ! isBinary
            end if ! rLine: scalars/vectors
         end if ! dType
      end do
         
 001   close(fid)
      
      do ivar=1, ptData%nvar
         ! calculate actual node position (MS disp) !
         if ( trim(ptData%varName(ivar)).eq.'MS_Displacement' ) then
            write(*,'(12X,2A)') 'Adjusting mesh coords using ',&
            trim(ptData%varName(ivar))
            ist = ptData%ioff(ivar) + 1 !(tdof-ptData%ndof(ivar))*nn+1
            do i=1, nn
               do j=1, nDim
                  k = maxNSD*(i-1)+j
                  coords(j,i) = coords(j,i) + &
                  ptData%arr(ist-1+k)
               end do ! j
            end do ! i
         end if
      end do ! ivar

      iMesh%nodes = nn
      iMesh%elems = ne
      allocate(iMesh%connec(nDim,iMesh%elems))
      allocate(iMesh%coords(nDim,iMesh%nodes))
      
      iMesh%connec = connec
      iMesh%coords = dble(coords)
      
      deallocate(connec)
      deallocate(coords)
      if (allocated(tmpI)) deallocate(tmpI)
      
      contains
      
         subroutine findKwrd(fileId,sKwrd,iPos,str)
         implicit none
         integer, intent(in) :: fileId
         integer, intent(inout) :: iPos
         character(len=*), intent(in) :: sKwrd
         character(len=stdl), intent(out) :: str
         
         integer :: i,kwrdL
         character :: c
         
         kwrdL = len(trim(sKwrd))
         do
            if ( isBinary ) then
               str = ''
               do i=1, stdl
                  read(fileId,pos=iPos) c
                  iPos = iPos + 1
                  if(c.eq.eol) exit
                  str(i:i) = c
               end do
            else
               read(fileId,'(A)') str
            end if
            str = adjustl(str)
            if ( str(1:kwrdL).eq.trim(sKwrd) ) return
         end do

         end subroutine findKwrd

      end subroutine readVTK
      
!**************************************************

      subroutine xigauss(nsd,eNoN,nG,wG,xi,N,Nx)
      implicit none
      integer, intent(in) :: nsd,eNoN,nG
      double precision, intent(out) :: wG(nG)
      double precision, intent(out) :: xi(nsd,eNoN)
      double precision, intent(out) :: N(eNoN,nG),Nx(nsd,eNoN,nG)
      
      integer :: g
      double precision :: w,r,s
      
      write(*,'(2X,A)') 'Setting shape function parameters..'
      if ( eNon.eq.3 ) then ! tri elem !
         w = 1.0d0/6.0d0
         r = 2.0d0/3.0d0
         s = 1.0d0/6.0d0
         xi(1,1) = r; xi(2,1) = s
         xi(1,2) = s; xi(2,2) = r
         xi(1,3) = s; xi(2,3) = s
         do g=1, nG
            wG(g)  = w
            N(1,g) = xi(1,g)
            N(2,g) = xi(2,g)
            N(3,g) = 1D0 - xi(1,g) - xi(2,g)
            
            Nx(1,1,g) =  1D0
            Nx(2,1,g) =  0D0
            Nx(1,2,g) =  0D0
            Nx(2,2,g) =  1D0
            Nx(1,3,g) = -1D0
            Nx(2,3,g) = -1D0
         end do

      else if ( eNoN.eq.4 ) then ! tet elem !
         w = 1.0d0/24.0d0
         r = (5.0d0 + 3.0d0*DSQRT(5.0d0))/20.0d0
         s = (5.0d0 -       DSQRT(5.0d0))/20.0d0
         xi(1,1) = r; xi(2,1) = s; xi(3,1) = s
         xi(1,2) = s; xi(2,2) = r; xi(3,2) = s
         xi(1,3) = s; xi(2,3) = s; xi(3,3) = r
         xi(1,4) = s; xi(2,4) = s; xi(3,4) = s
         do g=1, nG
            wG(g)  = w
            N(1,g) = xi(1,g)
            N(2,g) = xi(2,g)
            N(3,g) = xi(3,g)
            N(4,g) = 1D0 - xi(1,g) - xi(2,g) - xi(3,g)
            
            Nx(1,1,g) =  1D0
            Nx(2,1,g) =  0D0
            Nx(3,1,g) =  0D0
            Nx(1,2,g) =  0D0
            Nx(2,2,g) =  1D0
            Nx(3,2,g) =  0D0
            Nx(1,3,g) =  0D0
            Nx(2,3,g) =  0D0
            Nx(3,3,g) =  1D0
            Nx(1,4,g) = -1D0
            Nx(2,4,g) = -1D0
            Nx(3,4,g) = -1D0
         end do
      else
         write(*,'(4X,A)') 'Error: unknown element type..'
         STOP
      end if
      
      end subroutine xigauss
      
!**************************************************

      subroutine calcMeshProps
      use variables
      implicit none
      integer :: ist,iend,vecl
      
      write(*,'(2X,A)') 'Computing mesh properties..'

      call calcJacobian
      
      call calcSkewness
      
      call calcAR
      
      if ( .not.allocated(elData%arr) ) &
         allocate(elData%arr(3*iMesh%elems))
      elData%nvar = elData%nvar + 1
      elData%ndof(elData%nvar) = 1
      elData%tdof = elData%tdof + elData%ndof(elData%nvar)
      elData%varName(elData%nvar) = 'EL_Jacobian'
      elData%varType(elData%nvar) = 'float'
      ist  = (elData%tdof-elData%ndof(elData%nvar))*iMesh%elems + 1
      vecl =  elData%ndof(elData%nvar)*iMesh%elems
      iend =  ist + vecl - 1
      elData%ioff(elData%nvar) = ist-1
      elData%arr(ist:iend) = real(iMesh%J(1:vecl))

      elData%nvar = elData%nvar + 1
      elData%ndof(elData%nvar) = 1
      elData%tdof = elData%tdof + elData%ndof(elData%nvar)
      elData%varName(elData%nvar) = 'EL_Skewness'
      elData%varType(elData%nvar) = 'float'
      ist  = (elData%tdof-elData%ndof(elData%nvar))*iMesh%elems + 1
      vecl =  elData%ndof(elData%nvar)*iMesh%elems
      iend =  ist + vecl - 1
      elData%ioff(elData%nvar) = ist-1
      elData%arr(ist:iend) = real(iMesh%skewness(1:vecl))

      elData%nvar = elData%nvar + 1
      elData%ndof(elData%nvar) = 1
      elData%tdof = elData%tdof + elData%ndof(elData%nvar)
      elData%varName(elData%nvar) = 'EL_AspectRatio'
      elData%varType(elData%nvar) = 'float'
      ist  = (elData%tdof-elData%ndof(elData%nvar))*iMesh%elems + 1
      vecl =  elData%ndof(elData%nvar)*iMesh%elems
      iend =  ist + vecl - 1
      elData%ioff(elData%nvar) = ist-1
      elData%arr(ist:iend) = real(iMesh%AR(1:vecl))

      end subroutine calcMeshProps

!**************************************************

      subroutine calcJacobian
      use variables
      use allFun
      implicit none
      integer :: i,e,a,g,Ac,cnt
      integer :: ist,iend,vecl
      double precision :: sHat,maxJ
      double precision, allocatable, dimension(:,:) :: xl
      logical :: sgn
      
      write(*,'(8X,A)') 'Computing normalized element Jacobian..'
      allocate(xl(nDim,iMesh%nodesPerElem))
      cnt = 0
      sgn = .false.
      do e=1, iMesh%elems
         do a=1, iMesh%nodesPerElem
            Ac = iMesh%connec(a,e)
            xl(:,a) = iMesh%coords(:,Ac)
         end do
         iMesh%J(e) = Jacobian(nDim,iMesh%nodesPerElem,&
         xl,iMesh%Nx(:,:,1))
         if ( iMesh%J(1).lt.0.0d0 ) sgn = .true.
         if ( sgn ) iMesh%J(e) = -iMesh%J(e)
         if(iMesh%J(e).lt.0.0d0) then
            cnt = cnt+1
            !write(*,'(8X,4A)') &
            !'WARNING: Element Jacobian < 0; iel= ',trim(STR(e)),&
            !' J=',STR(iMesh%J(e))
         end if
      end do ! iel
      
      maxJ = dabs(maxval(iMesh%J))
      iMesh%J = iMesh%J / maxJ
      write(*,'(10X,4A)') 'Min. normalized Jacobian: ', &
      STR(minval(iMesh%J)), ' at el: ', &
      trim(STR(minloc(iMesh%J,DIM=1)))
      
      if ( minval(iMesh%J).lt.0.0d0 ) then
         write(*,'(12X,2A)') 'Warning: No of elems with J < 0: ', &
         trim(STR(cnt))
      end if

      iMesh%tMshInt = 0.0d0
      do e=1, iMesh%elems
         do g=1, iMesh%nGauss
            sHat = 0.0d0
            do a=1, iMesh%nodesPerElem
               sHat = sHat + iMesh%N(a,g)
            end do
            iMesh%tMshInt = iMesh%tMshInt + &
               iMesh%wG(g)*iMesh%J(e)*maxJ*sHat
         end do
      end do
      
      write(*,'(10X,2A)') 'Mesh Integral (Jacobian): ', &
      STR(iMesh%tMshInt)
      
      deallocate(xl)
      
      end subroutine calcJacobian
      
!**************************************************

      subroutine calcSkewness
      use variables
      use allFun
      implicit none
      integer :: i,j,e,a,Ac,cnt
      double precision :: s,Rc,integ_eq,integ_el
      double precision, allocatable, dimension(:) :: detD,mult
      double precision, allocatable, dimension(:,:) :: Dmat,Dsub

      write(*,'(A)')
      write(*,'(8X,A)') 'Computing element Skewness..'
      allocate(mult(iMesh%nodesPerElem+1))
      if ( vtkCellType.eq.5 ) then
         mult = (/1.0d0, -1.0d0, 1.0d0, -1.0d0/)
      else if ( vtkCellType.eq.10 ) then
         mult = (/1.0d0, 1.0d0, -1.0d0, 1.0d0, 1.0d0/)
      else
         write(*,'(10X,A)') 'Unknown element type..'
         return
      end if
      
      allocate(Dmat(iMesh%nodesPerElem,iMesh%nodesPerElem+1))
      allocate(Dsub(iMesh%nodesPerElem,iMesh%nodesPerElem))
      allocate(detD(iMesh%nodesPerElem+1))
      
      Dmat = 1.0d0
      Dsub = 0.0d0
      s = 0.0d0
      iMesh%tMshInt = 0.0d0
      do e=1, iMesh%elems
         do a=1, iMesh%nodesPerElem
            Ac = iMesh%connec(a,e)
            Dmat(a,2:nDim+1) = iMesh%coords(:,Ac)
            s = sum(iMesh%coords(:,Ac)**2)
            Dmat(a,1) = s
         end do
         
         do j=1, iMesh%nodesPerElem+1
            cnt = 0
            do i=1, iMesh%nodesPerElem+1
               if ( i.eq.j ) cycle
               cnt = cnt+1
               Dsub(:,cnt) = Dmat(:,i)
            end do
            detD(j) = mult(j)*det(Dsub,iMesh%nodesPerElem)
         end do
         Rc = dsqrt( sum(detD(2:iMesh%nodesPerElem)**2) - &
                  4.0d0*detD(1)*detD(iMesh%nodesPerElem+1) ) / &
               (2.0d0*dabs(detD(1)))
         
         if ( vtkCellType.eq.5 ) then
            integ_eq = (dsqrt(27.0d0)*Rc**2)/4.0d0
            integ_el = dabs(detD(1))/2.0d0
         else if ( vtkCellType.eq.10 ) then
            integ_eq = (8.0d0*Rc**3)/dsqrt(243.0d0)
            integ_el = dabs(detD(1))/6.0d0
         end if

         iMesh%skewness(e) = dabs(integ_eq - integ_el)/integ_eq
         
         iMesh%tMshInt = iMesh%tMshInt + integ_el
      end do

      write(*,'(10X,10A)') 'Skewness (min, max): ', &
      '(',STR(minval(iMesh%skewness)),', ', &
      STR(maxval(iMesh%skewness)),') at el: (', &
      trim(STR(minloc(iMesh%skewness, DIM=1))),', ', &
      trim(STR(maxloc(iMesh%skewness, DIM=1))),  ')'

      write(*,'(10X,2A)') 'Mesh Integral (Skewness): ', &
      STR(iMesh%tMshInt)

      deallocate(mult,Dmat,Dsub,detD)

      end subroutine calcSkewness
      
!**************************************************

      subroutine calcAR
      use variables
      use allFun
      implicit none
      integer :: i,j,e,a,b,Ac,Ac1,Ac2,ap,cnt
      integer :: irow,icol
      double precision :: minS,maxS
      double precision, dimension(:), allocatable :: s,dt
      double precision, dimension(:,:), allocatable :: xl,Dsub
      double precision, dimension(:,:), allocatable :: rowMat,colMat
      
      allocate(s(iMesh%nodesPerElem))
      
      if ( vtkCellType.eq.10 ) then
         allocate(xl(iMesh%nodesPerElem,nDim))
         allocate(rowMat(iMesh%nodesPerElem-1,iMesh%nodesPerElem))
         allocate(colMat(nDim-1,nDim))
         allocate(Dsub(nDim,nDim),dt(nDim))
         rowMat = reshape( (/ 1, 2, 3, 1, 2, 4, 1, 3, 4, 2, 3, 4 /),&
               shape(rowMat) )
         colMat = reshape( (/ 1, 2, 2, 3, 1, 3 /), shape(colMat))
      end if

      write(*,'(A)')
      write(*,'(8X,A)') 'Computing element Aspect Ratio (AR)..'
      s = 0.0d0
      do e=1, iMesh%elems
         if ( vtkCellType.eq.5 ) then
            do a=1, iMesh%nodesPerElem
               Ac1 = iMesh%connec(a,e)
               ap = a+1
               if ( a.eq.iMesh%nodesPerElem ) ap = 1
               Ac2 = iMesh%connec(ap,e)
               s(a) = dsqrt( sum( (iMesh%coords(:,Ac1) - &
                              iMesh%coords(:,Ac2))**2 ) )
            end do
         else if ( vtkCellType.eq.10 ) then
            do a=1, iMesh%nodesPerElem
               Ac = iMesh%connec(a,e)
               xl(a,:) = iMesh%coords(:,Ac)
            end do
            
            do a=1, iMesh%nodesPerElem
               do b=1, nDim
                  Dsub = 1.0d0
                  do i=1, iMesh%nodesPerElem-1
                     irow = rowMat(i,a)
                     do j=1, nDim-1
                        icol = colMat(j,b)
                        Dsub(i,j) = xl(irow,icol)
                     end do ! j
                  end do ! i
                  dt(b) = det(Dsub,nDim)
               end do ! b
               s(a) = 0.5d0*dsqrt( sum( dt(:)**2 ) )
            end do ! a
         else
            write(*,'(10X,A)') 'Unknown element type..'
            return
         end if ! vtkCellType
         
         minS = minval(s)
         maxS = maxval(s)
         iMesh%AR(e) = maxS/minS
      end do ! e

      write(*,'(10X,4A)') 'Max. Aspect Ratio: ', &
      STR(maxval(iMesh%AR)), ' at el: ',    &
      trim(STR(maxloc(iMesh%AR, DIM=1)))
      
      deallocate(s,dt,xl,Dsub,colMat,rowMat)
      
      end subroutine calcAR
      
!**************************************************

      subroutine writeVTK(fName)
      use variables
      use allFun
      implicit none
      character(len=*), intent(in) :: fName
      integer :: i,j,k,fid
      logical :: isBinary
      
      integer :: ivar,ist,iend,vecl

      write(*,'(2X,A)') 'Writing vtk file..'
      isBinary = .false.
      if ( outType .eq. 'VTKB' ) isBinary = .true.
      fid = 10

      do ivar=1, ptData%nvar
         if ( trim(ptData%varName(ivar)).eq.'MS_Displacement' ) then
            ist = ptData%ioff(ivar) + 1
            do i=1, iMesh%nodes
               do j=1, nDim
                  k = maxNSD*(i-1)+j
                  iMesh%coords(j,i) = iMesh%coords(j,i) - &
                  ptData%arr(ist-1+k)
               end do ! j
            end do ! i
         end if
      end do ! ivar
      
      if ( isBinary ) then
         open(fid,file=trim(fName),status='replace',access='stream',&
         form='unformatted',convert='big_endian')
         write(fid) '# vtk DataFile Version 3.0'//eol
         write(fid) 'Simulation Results'//eol
         write(fid) 'BINARY'//eol
         write(fid) 'DATASET UNSTRUCTURED_GRID'//eol
         write(fid) 'POINTS '//trim(STR(iMesh%nodes))//' float'//eol
      else
         open(fid,file=trim(fName),status='replace')
         write(fid,'(A)') '# vtk DataFile Version 3.0'
         write(fid,'(A)') 'Simulation Results'
         write(fid,'(A)') 'ASCII'
         write(fid,'(A)') 'DATASET UNSTRUCTURED_GRID'
         write(fid,'(A)') 'POINTS '//trim(STR(iMesh%nodes))//' float'
      end if
      
      ! points !
      write(*,'(8X,A)') 'Writing field POINTS to file..'
      call writeData(fid,isBinary,'','float',nDim,nDim*iMesh%nodes,&
      real(iMesh%coords(:,:)))
      
      ! cell connectivity and cell_types !
      if( isBinary ) then
         write(*,'(8X,A)') 'Writing field CELLS to file..'
         write(fid) 'CELLS '//trim(STR(iMesh%elems))//' '// &
         trim(STR( iMesh%elems*(iMesh%nodesPerElem+1) ))//eol
         do i=1, iMesh%elems
            write(fid) iMesh%nodesPerElem,iMesh%connec(:,i)-1
         end do
         
         write(*,'(8X,A)') 'Writing field CELL_TYPES to file..'
         write(fid) 'CELL_TYPES '//trim(STR(iMesh%elems))//eol
         do i=1, iMesh%elems
            write(fid) vtkCellType
         end do
      else
         write(*,'(8X,A)') 'Writing field CELLS to file..'
         write(fid,'(A)') 'CELLS '//trim(STR(iMesh%elems))//' '// &
         trim(STR( iMesh%elems*(iMesh%nodesPerElem+1) ))
         do i=1, iMesh%elems
            write(fid,'(A)',advance='no') &
            trim(STR(iMesh%nodesPerElem))
            do j=1, iMesh%nodesPerElem
               write(fid,'(A)',advance='no') &
               ' '//trim(STR(iMesh%connec(j,i)-1))
            end do
            write(fid,'(A)')
         end do
         
         write(*,'(8X,A)') 'Writing field CELL_TYPES to file..'
         write(fid,'(A)') 'CELL_TYPES '//trim(STR(iMesh%elems))
         do i=1, iMesh%elems
            write(fid,'(A)') trim(STR(vtkCellType))
         end do
      end if
      
      ! point_data !
      if ( allocated(ptData%arr) ) then
         if ( isBinary ) then
            write(fid) 'POINT_DATA '//trim(STR(iMesh%nodes))//eol
         else
            write(fid,'(A)') 'POINT_DATA '//trim(STR(iMesh%nodes))
         end if
         do ivar=1, ptData%nvar
            ist  = ptData%ioff(ivar) + 1
            vecl =  ptData%ndof(ivar)*iMesh%nodes
            iend =  ist + vecl - 1
            call writeData(fid,isBinary,trim(ptData%varName(ivar)),&
            trim(ptData%varType(ivar)),ptData%ndof(ivar),vecl, &
            ptData%arr(ist:iend))
         end do
      else
         write(*,'(8X,A)') 'Warning: no point data found to write..'
      end if

      ! cell_data !
      if ( allocated(elData%arr) ) then
         if ( isBinary ) then
            write(fid) 'CELL_DATA '//trim(STR(iMesh%elems))//eol
         else
            write(fid,'(A)') 'CELL_DATA '//trim(STR(iMesh%elems))
         end if
         do ivar=1, elData%nvar
            ist  = elData%ioff(ivar) + 1
            vecl =  elData%ndof(ivar)*iMesh%elems
            iend =  ist + vecl - 1
            call writeData(fid,isBinary,trim(elData%varName(ivar)),&
            trim(elData%varType(ivar)),elData%ndof(ivar),vecl, &
            elData%arr(ist:iend))
         end do
      else
         write(*,'(8X,A)') 'Warning: no cell data found to write..'
      end if

      close(fid)
      
      end subroutine writeVTK
      
!**************************************************

      subroutine writeData(fileId,isBin,sName,sType,ndof,n,arr)
      use variables
      use allFun
      implicit none
      logical, intent(in) :: isBin
      integer, intent(in) :: fileId,ndof,n
      character(len=*), intent(in) :: sName,sType
      real*4, intent(in) :: arr(n)
      
      integer :: i
      
      if ( trim(sName).ne.'' .and. trim(sType).ne.'' ) then
         write(*,'(8X,A)') &
         'Writing field '//trim(sName)//' to file..'
         if ( isBin ) then
            if ( ndof.eq.1 ) then
               write(fileId) &
               'SCALARS '//trim(sName)//' '//trim(sType)//eol
               write(fileId) 'LOOKUP_TABLE default'//eol
            else
               write(fileId) &
               'VECTORS '//trim(sName)//' '//trim(sType)//eol
            end if
         else
            if ( ndof.eq.1 ) then
               write(fileId,'(A)') &
               'SCALARS '//trim(sName)//' '//trim(sType)
               write(fileId,'(A)') &
               'LOOKUP_TABLE default'
            else
               write(fileId,'(A)') &
               'VECTORS '//trim(sName)//' '//trim(sType)
            end if
         end if
      end if
      
      if ( isBin ) then
         if ( ndof.eq.2 ) then
            select case ( trim(sType) )
            case ('int')
               do i=1, n/2
                  write(fileId) &
                  int(arr(2*i-1)),  int(arr(2*i)),  0
               end do
            case ('float')
               do i=1, n/2
                  write(fileId) &
                  arr(2*i-1), arr(2*i), 0.0
               end do
            end select
         else
            select case ( trim(sType) )
            case ('int')
               write(fileId) int(arr(:))
            case ('float')
               write(fileId) arr(1:n)
            end select
         end if
      else
         do i=1, n/ndof
            select case (ndof)
            case (1)
               if ( trim(sType).eq.'int' ) then
                  write(fileId,'(A)') trim(STR(int(arr(i))))
               else if ( trim(sType).eq.'float' ) then
                  write(fileId,'(A)') STR(arr(i))
               end if
            case (2)
               if ( trim(sType).eq.'int' ) then
                  write(fileId,'(A)') trim(STR(int(arr(2*i-1))))// &
                  ' '//trim(STR(int(arr(2*i))))//' 0'
               else if ( trim(sType).eq.'float' ) then
                  write(fileId,'(A)') STR(arr(2*i-1))//' '// &
                  STR(arr(2*i))//' 0.0'
               end if
            case (3)
               if ( trim(sType).eq.'int' ) then
                  write(fileId,'(A)') trim(STR(int(arr(3*i-2))))// &
                  ' '//trim(STR(int(arr(3*i-1))))//' '// &
                  trim(STR(int(arr(3*i))))
               else if ( trim(sType).eq.'float' ) then
                  write(fileId,'(A)') STR(arr(3*i-2))//' '// &
                  STR(arr(3*i-1))//' '//STR(arr(3*i))
               end if
            end select
         end do
      end if

      end subroutine writeData
      
!**************************************************

      subroutine freeMem
      use variables
      use allFun
      implicit none
      
      deallocate(iMesh%connec)
      deallocate(iMesh%coords)
      deallocate(iMesh%wG)
      deallocate(iMesh%J)
      deallocate(iMesh%skewness)
      deallocate(iMesh%AR)
      deallocate(iMesh%xi)
      deallocate(iMesh%N)
      deallocate(iMesh%Nx)
      deallocate(ptData%arr)
      deallocate(elData%arr)
      
      end subroutine freeMem
      
 !**************************************************



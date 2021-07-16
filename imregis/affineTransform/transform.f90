    module variables
    type vectors 
      integer :: nodes,elems
      integer, allocatable, dimension(:,:) :: conn
      double precision,allocatable,dimension(:,:) :: coord     
    end type
    type(vectors) :: so,tar
    integer :: nDim
    integer :: iSmooth,smoothIter
    double precision :: scaling
    double precision, allocatable, dimension(:) :: trans
    double precision, allocatable, dimension(:,:) :: rotMat
    character *80 :: fName
    end module variables
    
!**************************************************

    program affine
    use variables
    implicit none
    integer :: iFile,i,j

    call readInputs

    allocate(so%coord(3,so%nodes),so%conn(3,so%elems))

    open(20,file=trim(fName))
    read(20,*)
    read(20,*)
    read(20,*)
    do i=1, so%nodes
      read(20,*) (so%coord(j,i),j=1,3)
    end do
    do i=1, so%elems
      read(20,*) (so%conn(j,i), j=1,3)
    end do
    close(20)

    call transform

    if(iSmooth.eq.1) call smoother

    call writeOutput
      
    deallocate(so%coord,so%conn,tar%coord,tar%conn)

    end program affine

!**************************************************

    subroutine readInputs
    use variables
    implicit none
    integer :: i,j

    open(10,file='in-transform.dat')
    read(10,*)
    read(10,*) nDim
    read(10,*)
    read(10,*) fName
    read(10,*)
    read(10,*) so%nodes, so%elems

    allocate(rotMat(nDim,nDim),trans(nDim))
    read(10,*)
    read(10,*)
    do j=1, nDim
      read(10,*) (rotMat(j,i),i=1, nDim)
    end do
    read(10,*)
    read(10,*) scaling
    read(10,*)
    read(10,*) (trans(i),i=1, nDim)
    read(10,*)
    read(10,*) iSmooth, smoothIter
    close(10)

    if(nDim.ne.2 .and. nDim.ne.3) then
      write(*,*) 'ERROR: Wrong choice of nDim. Can be only 2 or 3'
      STOP
    end if
    
    end subroutine readInputs

!**************************************************

    subroutine transform
    use variables
    implicit none
    integer :: i

    tar%nodes = so%nodes
    tar%elems = so%elems
    allocate(tar%coord(3,tar%nodes),tar%conn(3,tar%elems))
    tar%conn = so%conn
    
    do i=1, tar%nodes
      tar%coord(:,i) = scaling*matmul(rotMat,so%coord(:,i)) + trans
    end do
    
    end subroutine transform

!**************************************************
    
    subroutine smoother
    use variables
    implicit none
    integer :: i,j,nn,iter
    integer :: nAdj(tar%nodes),nList(tar%nodes)
    integer,dimension(20,tar%nodes) :: adjNodes
    double precision :: sum
    double precision, dimension(3,tar%nodes) :: pts
    logical, dimension(tar%nodes) :: iFixed

    iFixed = .false.

    do i=1, tar%nodes
      if(.not.iFixed(i)) then
        call getAdjNodes(tar%conn,tar%elems,i,nList,nn)
        nAdj(i) = nn
        do j=1, nAdj(i)
          adjNodes(j,i) = nList(j)
        end do
      end if
    end do ! i
   
    do iter=1, smoothIter
      write(*,*) 'Smoothing iteration..', iter

      do i=1, tar%nodes
        pts(:,i) = tar%coord(:,i)
      end do

      do i=1, tar%nodes
        if(.not.iFixed(i)) then         
          sum = 0.0d0
          do j=1, nAdj(i)
            sum = sum + tar%coord(1,adjNodes(j,i))
          end do 
          pts(1,i) = sum/nAdj(i)

          sum = 0.0d0
          do j=1, nAdj(i)
            sum = sum + tar%coord(2,adjNodes(j,i))
          end do
          pts(2,i) = sum/nAdj(i)

          sum = 0.0d0
          do j=1, nAdj(i)
            sum = sum + tar%coord(3,adjNodes(j,i))
          end do 
          pts(3,i) = sum/nAdj(i)
        end if
      end do ! i
  
      do i=1, tar%nodes
        tar%coord(:,i) = pts(:,i)
      end do
    end do
    
    end subroutine smoother

!**************************************************

    subroutine getAdjNodes(connec,nel,n,nList,cnt)
    use variables
    implicit none
    integer, intent(in) :: nel,n,connec(3,nel)
    integer, intent(inout) :: cnt,nList(cnt)
    integer :: i
    
    cnt=0    
    do i=1, nel
      if( (connec(1,i).eq.n) .or. &
          (connec(2,i).eq.n) .or. &
          (connec(3,i).eq.n) ) then
          nList(cnt+1) = connec(1,i)
          nList(cnt+2) = connec(2,i)
          nList(cnt+3) = connec(3,i)
          cnt = cnt+3
      end if
    end do  ! i

    call editList(nList,cnt)

    end subroutine getAdjNodes

!**************************************************

    subroutine editList(data, n)
    implicit none
    integer, intent(inout) :: n, data(n)
    integer :: i, j, cnt, newdata(n)
       
    do i=1, n
     do j=i, n
      if(data(i).ge.data(j)) then
       cnt=data(i)
       data(i) = data(j)
       data(j) = cnt
      end if 
     end do 
    end do
       
    cnt=1
    newdata(cnt) = data(cnt)
    do i=2, n
     if(data(i-1).ne.data(i)) then
      cnt=cnt+1
      newdata(cnt) = data(i)
     end if
    end do

    data(1:cnt) = newdata(1:cnt)
    n=cnt

    end subroutine editlist

!**************************************************

    subroutine writeOutput
    use variables
    implicit none
    integer :: i,j,dum

    dum = len(trim(fName))
    dum = dum-4
    write(fName,'(A,"-transformed.dat")') fName(1:dum)
    open(50,file=trim(fName))
    write(50,*) 'VARIABLES = "x", "y", "z"'
    write(50,*) 'ZONE N=', tar%nodes, ', E=', tar%elems, ', ZONETYPE=FETRIANGLE'
    write(50,*) 'DATAPACKING=POINT'
    do i=1, tar%nodes
      write(50,'(3(2X,1p,E15.6))') (tar%coord(j,i), j=1,3)
    end do
    do i=1, tar%elems
      write(50,'(3(2X,I6))') (tar%conn(j,i), j=1,3)
    end do
    close(50)

    end subroutine writeOutput

!**************************************************
    

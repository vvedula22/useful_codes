    module variables
    type vectors 
      integer :: nodes,elems
      integer, allocatable, dimension(:) :: pt
      integer, allocatable, dimension(:,:) :: conn
      double precision,allocatable,dimension(:,:) :: coord     
      character *80 :: fName     
    end type
    type(vectors) :: so,tar,newTar
    end module variables
    
!**************************************************

    module affineVar
    integer :: nDim,nPts
    integer :: LDA,LDU,LDVT,LWORK,INFO
    double precision :: xVar,yVar,det,scaling
    double precision, dimension(:), allocatable :: xMean,yMean,WORK
    double precision, dimension(:), allocatable :: dMat,trans
    double precision, dimension(:,:), allocatable :: xCoord,yCoord,coVar,tmp
    double precision, dimension(:,:), allocatable :: sMat,uMat,vtMat,RMat       
    end module affineVar
    
!**************************************************

    program affine
    use variables
    use affineVar
    implicit none
    integer :: i,j

    call readFiles

    do j=1, nDim
      do i=1, nPts
        xCoord(i,j) = so%coord(j,so%pt(i))
      end do
    end do
    
    do j=1, nDim
      do i=1, nPts
        yCoord(i,j) = tar%coord(j,tar%pt(i))
      end do
    end do

    open(1000,file='affine-out.log')
    write(1000,*) 'source:', trim(so%fName)
    write(1000,*) 'target:', trim(tar%fName)
    call calAffine

    call transform
    close(1000)

    end program affine

!**************************************************

    subroutine readFiles
    use variables
    use affineVar
    implicit none
    integer :: i,j

    open(10,file='affine-in.dat')
    read(10,*)
    read(10,*) nDim, nPts
    allocate(so%pt(nPts),tar%pt(nPts))
    read(10,*)
    read(10,*) so%fName
    read(10,*)
    read(10,*) so%nodes, so%elems
    read(10,*)
    read(10,*) (so%pt(i), i=1,nPts)
    read(10,*)
    read(10,*) tar%fName
    read(10,*) 
    read(10,*) tar%nodes, tar%elems
    read(10,*)
    read(10,*) (tar%pt(i), i=1, nPts)
    close(10)

    if(nDim.ne.2 .and. nDim.ne.3) then
      write(*,*) 'ERROR: Wrong choice of nDim. Can be only 2 or 3'
      STOP
    end if

    allocate(so%coord(3,so%nodes),so%conn(3,so%elems))
    allocate(tar%coord(3,tar%nodes),tar%conn(3,tar%elems))
    call allocateVar
    
    open(20,file=trim(so%fName))
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

    open(30,file=trim(tar%fName))
    read(30,*)
    read(30,*)
    read(30,*) 
    do i=1, tar%nodes
      read(30,*) (tar%coord(j,i), j=1,3)
    end do
    do i=1, tar%elems
      read(30,*) (tar%conn(j,i), j=1,3)
    end do      
    close(30)

    end subroutine readFiles

!**************************************************

    subroutine allocateVar
    use affineVar
    implicit none
 
    allocate(xMean(nDim),yMean(nDim))
    allocate(xCoord(nPts,nDim),yCoord(nPts,nDim))
    allocate(coVar(nDim,nDim),tmp(nDim,nDim))
    allocate(uMat(nDim,nDim),vtMat(nDim,nDim))
    allocate(dMat(nDim),sMat(nDim,nDim))
    allocate(RMat(nDim,nDim),trans(nDim))
    
    end subroutine allocateVar

!**************************************************    

    subroutine calAffine
    use affineVar
    implicit none
    integer :: i,j,k

    write(1000,*) 'Mean:'
    do i=1, nDim
      xMean(i) = sum(xCoord(:,i))
      yMean(i) = sum(yCoord(:,i))
    end do
    xMean = xMean/dble(nPts)
    yMean = yMean/dble(nPts)
    write(1000,*) 'X'
    write(1000,*) (xMean(i), i=1,nDim)
    write(1000,*) 'Y'
    write(1000,*) (yMean(i), i=1,nDim)
    write(1000,*)
    
    write(1000,*) 'Variance'
    xVar = 0.0d0
    yVar = 0.0d0
    do j=1, nPts
      do i=1, nDim
        xVar = xVar + (xCoord(j,i)-xMean(i))**2
        yVar = yVar + (yCoord(j,i)-yMean(i))**2
      end do
    end do
    xVar = xVar/dble(nPts)
    yVar = yVar/dble(nPts)
    write(1000,*) 'X'
    write(1000,*) xVar
    write(1000,*) 'Y'
    write(1000,*) yVar
    write(1000,*)

    write(1000,*) 'Covariance'
    coVar = 0.0d0
    do i=1, nDim
      do j=1, nDim
        do k=1, nPts
          coVar(i,j) = coVar(i,j) + (yCoord(k,i)-yMean(i)) * (xCoord(k,j)-xMean(j))
        end do
        coVar(i,j) = coVar(i,j)/dble(nPts)
      end do
    end do
    do j=1, nDim
      write(1000,*) (coVar(j,i),i=1,3)
    end do
    write(1000,*)

    LDA=nDim; LDU=nDim; LDVT=nDim
    LWORK = 100*nDim
    allocate(WORK(LWORK))
    call DGESVD('A','A',nDim,nDim,coVar,LDA,dMat,uMat,LDU,vtMat,LDVT,WORK,LWORK,INFO)
    if(INFO.ne.0) then
      print*, 'ERROR: Singular Value Decomposition failed'
      STOP
    end if
    
    select case(nDim)
    case(2)
     det = coVar(1,1)*coVar(2,2) - coVar(1,2)*coVar(2,1)
    case(3)
     det = coVar(1,1)*coVar(2,2)*coVar(3,3) + coVar(1,2)*coVar(2,3)*coVar(3,1) + coVar(1,3)*coVar(2,1)*coVar(3,2) &
         - coVar(1,3)*coVar(2,2)*coVar(3,1) - coVar(1,2)*coVar(2,1)*coVar(3,3) - coVar(1,1)*coVar(2,3)*coVar(3,2)
    end select
    write(1000,*) 'Det:', det    
    if(abs(det).lt.1E-12) &
      print*, 'Warning: determinant near zero'

    sMat = 0.0d0
    do i=1, nDim
      sMat(i,i) = 1.0d0
    end do
    if(det.lt.0) then
     sMat(nDim,nDim) = -1.0d0
     write(1000,*) 'Sign of the last diagonal element changed..'
    endif

    write(1000,*) 'Affine Transformation:'
    write(1000,*) 'Rotation:'
    tmp = matmul(uMat,sMat)
    RMat = matmul(tmp,vtMat)
    do j=1, nDim
      write(1000,*) (RMat(j,i), i=1,nDim)
    end do
    write(1000,*)

    write(1000,*) 'Scaling:'
    scaling = 0.0d0
    do i=1, nDim
      scaling = scaling + dMat(i)*sMat(i,i)
    end do
    scaling = scaling/xVar
    write(1000,*) scaling
    write(1000,*)

    write(1000,*) 'Translation:'
    trans = yMean - scaling*matmul(RMat,xMean)
    write(1000,*) trans
    write(1000,*)
    close(1000)
    
    end subroutine calAffine

!**************************************************

    subroutine transform
    use variables
    use affineVar
    implicit none
    integer :: i,j,dum

    newTar%nodes = so%nodes
    newTar%elems = so%elems
    allocate(newTar%coord(3,newTar%nodes),newTar%conn(3,newTar%elems))
    newTar%conn = so%conn
    
    do i=1, newTar%nodes
      newTar%coord(:,i) = scaling*matmul(RMat,so%coord(:,i)) + trans
    end do

    dum = len(trim(so%fName))
    dum = dum-4
    write(newTar%fName,'(A,"-new.dat")') so%fName(1:dum)

    open(50,file=trim(newTar%fName))
    write(50,*) 'VARIABLES = "x", "y", "z"'
    write(50,*) 'ZONE N=', newTar%nodes, ', E=', newTar%elems, ', ZONETYPE=FETRIANGLE'
    write(50,*) 'DATAPACKING=POINT'
    do i=1, newTar%nodes
      write(50,'(3(2X,1p,E15.6))') (newTar%coord(j,i), j=1,3)
    end do
    do i=1, newTar%elems
      write(50,'(3(2X,I6))') (newTar%conn(j,i), j=1,3)
    end do
    close(50)

    end subroutine transform

!**************************************************
    

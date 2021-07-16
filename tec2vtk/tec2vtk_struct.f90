

    program dat2VTK
    implicit none
    integer, parameter :: nvar=6,iout=2
    integer, parameter :: nx=256,ny=256,nz=148
    integer :: i,j,k
    double precision :: x(nx),y(ny),z(nz)
    double precision, dimension(nx,ny,nz) :: u,v,w
    character, parameter :: eol=char(10)
    
    ! read tecplot dat file !
    write(*,*) 'Reading tecplot dat file..'
    open(10,file='in.dat')
    read(10,*)
    do i=1, nvar
        read(10,*)
    end do
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    do k=1, nz
        do j=1, ny
            do i=1, nx
                read(10,*) x(i),y(j),z(k), &
                u(i,j,k),v(i,j,k),w(i,j,k)
            end do
        end do
    end do
    close(10)
    
    select case(iout)
    case(1)
        ! write ascii vtk file !
        write(*,*) 'Writing ascii vtk file..'
        open(10,file='out_ascii.vtk')
        write(10,'(A)') '# vtk DataFile Version 3.0'
        write(10,'(A)') 'Simulation Results'
        write(10,'(A)') 'ASCII'
        write(10,'(A)') 'DATASET RECTILINEAR_GRID'
        write(10,'(A,3(X,I4))') 'DIMENSIONS', nx,ny,nz
        write(10,'(A,X,I4,X,A)') 'X_COORDINATES',nx,'double'
        do i=1, nx
            write(10,'(ES14.6E3)',advance='no') x(i)
        end do
        write(10,*)
        write(10,'(A,X,I4,X,A)') 'Y_COORDINATES',ny,'double'
        do j=1, ny
            write(10,'(ES14.6E3)',advance='no') y(j)
        end do
        write(10,*)
        write(10,'(A,X,I4,X,A)') 'Z_COORDINATES',nz,'double'
        do k=1, nz
            write(10,'(ES14.6E3)',advance='no') z(k)
        end do
        write(10,*)
        write(10,'(A,X,I10)') 'POINT_DATA', nx*ny*nz
        write(10,'(A)') 'VECTORS Velocity double'
        do k=1, nz
            do j=1, ny
                do i=1, nx
                    write(10,'(3(ES14.6E3,X))') &
                    u(i,j,k),v(i,j,k),w(i,j,k)
                end do
            end do
        end do
        close(10)
    case(2)
        ! write binary vtk file !
        write(*,*) 'Writing binary vtk file..'
        open(10,file='out_binary.vtk',status='unknown',access='stream',&
        form='unformatted',convert='big_endian')
        write(10) "# vtk DataFile Version 3.0"//eol
        write(10) "Simulation Results"//eol
        write(10) "BINARY"//eol
        write(10) "DATASET RECTILINEAR_GRID"//eol
        write(10) "DIMENSIONS "//trim(int2str(nx))//" "//&
        trim(int2str(ny))//" "//trim(int2str(nz))//eol
        write(10) "X_COORDINATES "//trim(int2str(nx))//" double"//eol
        write(10) x
        write(10) "Y_COORDINATES "//trim(int2str(ny))//" double"//eol
        write(10) y
        write(10) "Z_COORDINATES "//trim(int2str(nz))//" double"//eol
        write(10) z
        write(10) "POINT_DATA "//trim(int2str(nx*ny*nz))//eol
        write(10) "VECTORS Velocity double"//eol
        do k=1, nz
            do j=1, ny
                do i=1, nx
                    write(10) u(i,j,k),v(i,j,k),w(i,j,k)
                end do
            end do
        end do
        close(10)
    end select
    
    contains
        character(len=20) function int2str(i)
        integer, intent(in) :: i
        write(int2str,*) i
        int2str = adjustl(int2str)
        end function int2str
    end program

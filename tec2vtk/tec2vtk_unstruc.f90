

    program dat2VTK
    implicit none
    integer, parameter :: nvar=6,iout=2
    integer, parameter :: nodes=10228, elems=20421
    integer :: i,j
    double precision, dimension(nodes) :: x,y,z
    double precision, dimension(nodes) :: u,v,w
    integer, dimension(3,elems) :: conn
    character, parameter :: eol=char(10)
    
    ! read tecplot dat file !
    write(*,*) 'Reading tecplot dat file..'
    open(10,file='in.dat')
    read(10,*)
    read(10,*)
    read(10,*)
    do i=1, nodes
        read(10,*) x(i),y(i),z(i)!,u(i),v(i),w(i)
    end do
    do i=1, elems
        read(10,*) (conn(j,i),j=1, 3)
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
        write(10,'(A)') 'DATASET UNSTRUCTURED_GRID'
        write(10,'(A,X,I10,X,A)') 'POINTS', nodes, 'double'
        do i=1, nodes
            write(10,'(3(ES14.6E3,X))') x(i),y(i),z(i)
        end do
        write(10,'(A,2(X,I10))') 'CELLS', elems, elems*4
        do i=1, elems
            write(10,'(I1,X,3(I10,X))') 3,(conn(j,i)-1,j=1, 3)
        end do
        write(10,'(A,X,I10)') 'CELL_TYPES', elems
        do i=1, elems
            write(10,'(A)') '5'
        end do
        !write(10,'(A,X,I10)') 'POINT_DATA', nodes
        !write(10,'(A)') 'VECTORS Velocity double'
        !do i=1, nodes
        !   write(10,'(3(ES14.6E3,X))') u(i),v(i),w(i)
        !end do
        close(10)
    case(2)
        ! write binary vtk file !
        write(*,*) 'Writing binary vtk file..'
        open(10,file='out_binary.vtk',status='unknown',access='stream',&
        form='unformatted',convert='big_endian')
        write(10) "# vtk DataFile Version 3.0"//eol
        write(10) "Simulation Results"//eol
        write(10) "BINARY"//eol
        write(10) "DATASET UNSTRUCTURED_GRID"//eol
        write(10) "POINTS "//trim(int2str(nodes))//" double"//eol
        do i=1, nodes
            write(10) x(i),y(i),z(i)
        end do
        write(10) "CELLS "//trim(int2str(elems))//" "//&
        trim(int2str(elems*4))//eol
        do i=1, elems
            write(10) 3,(conn(j,i)-1,j=1, 3)
        end do
        write(10) "CELL_TYPES "//trim(int2str(elems))//eol
        do i=1, elems
            write(10) 5
        end do
        !write(10) "POINT_DATA "//trim(int2str(nodes))//eol
        !write(10) "VECTORS Velocity double"//eol
        !do i=1, nodes
        !   write(10) u(i),v(i),w(i)
        !end do
        close(10)
    end select
    
    contains
        character(len=20) function int2str(i)
        integer, intent(in) :: i
        write(int2str,*) i
        int2str = adjustl(int2str)
        end function int2str
    end program

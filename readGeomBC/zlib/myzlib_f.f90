!**************************************************

	program zlibTemp
	implicit none
	integer, parameter :: stdl=256
	integer, parameter :: np=10,nsd=3
	character, parameter :: eol=char(0)
	character, parameter :: newl = char(10)
	integer :: fid,i,j,nsize
	character(len=256) :: fName,mode,str
	double precision, dimension(3,np) :: coords, coords1
	
	fName = ''
	mode = ''
	
	open(10,file='ascii.dat')
	do i=1, np
		read(10,*) (coords(j,i),j=1, nsd)
	end do
	close(10)
	
	write(fName,'(A)') "bindata.dat"//eol
	write(mode,'(A)') "write"//eol
	call openFile(trim(fName),trim(mode),fid)
	if ( fid.eq.0 ) then
		write(*,'(2A)') 'Error opening file ', trim(fName)
		STOP
	end if
	
	str = "Hello World"//newl//eol
	call writeString(fid,str)
	
	nsize = np*nsd
	call writeData(fid,coords,nsize)
	call closeFile(fid)


	write(fName,'(A)') "bindata.dat"//eol
	write(mode,'(A)') "read"//eol
	call openFile(trim(fName),trim(mode),fid)
	if ( fid.eq.0 ) then
		write(*,'(2A)') 'Error opening file ', trim(fName)
		STOP
	end if
	
	call readString(fid,str)
	
	write(*,'(A)') trim(str)
	
	nsize = np*nsd
	call readData(fid,coords1,nsize)
	call closeFile(fid)
	
	do i=1, np
		write(*,*) (coords(j,i),j=1, nsd)
	end do

	write(*,*)
	do i=1, np
		write(*,*) (coords1(j,i),j=1, nsd)
	end do

	end program zlibTemp
	
!**************************************************
	

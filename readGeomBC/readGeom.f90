!**************************************************

		module params
		integer, parameter :: stdl=256
		integer, parameter :: svMagicNo=362436
		logical, parameter :: isBinary=.true.
		character, parameter :: eol=achar(0)
		character, parameter :: newl=achar(10)
		end module params
		
!**************************************************

		module var
		use params
		use, intrinsic :: ISO_C_BINDING
		integer :: nArr
		integer(kind=c_int), allocatable, dimension(:), target :: intArr
		real(kind=c_double), allocatable, dimension(:), target :: dbleArr
		end module var
		
!**************************************************

		module geomData
		use params
		integer, parameter :: nsd=3
		integer, parameter :: nelType=8
		integer :: nprocs,nvar,np,nedges,nfaces
		integer :: nmodes,nel,nelb,nen,nenb
		integer :: nelblk,nelblkb
		integer :: ndirbc,nshf
		integer :: bnen,bpoly,bnsh,bnshlb,bnenbl,blcsyst
		integer :: numNBC,numEBC
		integer :: nxadj,nadjncy,nvwgt
		integer, dimension(:), allocatable :: belIDs,bcMap,bcCodes,iPer
		integer, dimension(:), allocatable :: xadj,adjncy,vwgt
		integer, dimension(:,:), allocatable :: connec,connecb
		double precision, dimension(:,:), allocatable :: coords,bcInp
		
		! data blocks !
		integer :: elType(nelType,nsd)
		data elType	/	2,	2,	2,	2,	3,	3,	3,	3,	&
						4,	3,	3,	4,	9,	6,	6,	9,	&
						8,	4,	6,	6,	27,	10,	18,	18	/
		
		end module geomData
		
!**************************************************

		module allfun
		use params
		implicit none
		
		interface STR
			module procedure :: ITSTR, RTSTR, DTSTR, NDTSTR
		end interface STR

		contains
			character(len=stdl) function strtok(str, dlms)
			implicit none
			character(len=*), intent(in) :: str
			character(len=*), intent(in) :: dlms
			
			integer :: ist, iend

			integer, save :: ist0,slen
			character*255, save :: str0
			
			if ( str(1:1).ne.eol ) then
				ist0 = 1
				str0 = str
				slen = len(trim(str0))
			end if
			
			ist = ist0	
			do
				if ( (ist.le.slen) .and. (index(dlms,str0(ist:ist)).ne.0) ) then
					ist = ist+1
				else
					exit
				end if
			end do
			
			if ( ist.gt.slen ) then
				strtok = eol
				return
			end if
			
			iend = ist
			do
				if ( (iend.le.slen) .and. (index(dlms,str0(iend:iend)).eq.0) ) then
					iend = iend+1
				else
					exit
				end if
			end do
			
			strtok = str0(ist:iend-1)
			ist0 = iend+1
			return
			
			end function strtok

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
		
		end module allfun
		
!**************************************************

		program readGeomBC
		use modGzip
		use var
		use geomData
		use allfun
		implicit none
		integer(c_int) :: fid,rt
		character(len=stdl) :: fName,mode,hdr
		
		integer :: i,j,cnt,ijunk,nLines,nBytesRead
		
		fName = ''
		mode  = ''
		hdr   = ''
		nBytesRead = 0
		nLines = 0
		
		write(fName,'(A)') "geombc.dat.1"//eol
		write(mode, '(A)') "rb"//eol
		fid = gzOpen(fName,mode)

		write(hdr,'(A)') "number of processors"
		nArr = 1
		call resetArr
		call readHeader(fid,hdr,"integer")
		nprocs = intArr(1)
		nLines = nLines+1

		write(hdr,'(A)') "number of variables"
		call readHeader(fid,hdr,"integer")
		nvar = intArr(1)
		nLines = nLines+1

		write(hdr,'(A)') "number of nodes"
		call readHeader(fid,hdr,"integer")
		np = intArr(1)
		nLines = nLines+1

		write(hdr,'(A)') "number of nodes in the mesh"
		call readHeader(fid,hdr,"integer")
		nLines = nLines+1

		write(hdr,'(A)') "number of edges in the mesh"
		call readHeader(fid,hdr,"integer")
		nedges = intArr(1)
		nLines = nLines+1

		write(hdr,'(A)') "number of faces in the mesh"
		call readHeader(fid,hdr,"integer")
		nfaces = intArr(1)
		nLines = nLines+1
		
		write(hdr,'(A)') "number of modes"
		call readHeader(fid,hdr,"integer")
		nmodes = intArr(1)
		nLines = nLines+1
		
		write(hdr,'(A)') "number of shapefunctions solved on processor"
		call readHeader(fid,hdr,"integer")
		nshf = intArr(1)
		nLines = nLines+1
		
		write(hdr,'(A)') "number of interior elements"
		call readHeader(fid,hdr,"integer")
		nel = intArr(1)
		nLines = nLines+1
		
		write(hdr,'(A)') "number of boundary elements"
		call readHeader(fid,hdr,"integer")
		nelb = intArr(1)
		nLines = nLines+1
		
		write(hdr,'(A)') "maximum number of element nodes"
		call readHeader(fid,hdr,"integer")
		nen = intArr(1)
		nLines = nLines+1
		
		write(hdr,'(A)') "number of interior tpblocks"
		call readHeader(fid,hdr,"integer")
		nelblk = intArr(1)
		nLines = nLines+1
		
		write(hdr,'(A)') "number of boundary tpblocks"
		call readHeader(fid,hdr,"integer")
		nelblkb = intArr(1)
		nLines = nLines+1
		
		write(hdr,'(A)') "number of nodes with Dirichlet BCs"
		call readHeader(fid,hdr,"integer")
		ndirbc = intArr(1)
		nLines = nLines+1
		
		! find boundary nen based on nen and elemType catalog !
		nenb = 0
		do i=1, nelType
			if ( nen.eq.elType(i,nsd) ) then
				nenb = max( elType(i,nsd-1), nenb )
			end if
		end do
		
		if ( nenb.eq.0 ) then
			write(*,'(A)') "Error: wrong elem type (nen: ",nen,")"
			STOP
		end if
		
		! data blocks !
		write(*,'(A)')
		write(*,'(A)') "***** Data Blocks *****"
		! coordinates data !
		write(hdr,'(A)') "co-ordinates"
		nArr = 2
		call resetArr
		call readHeader(fid,hdr,"integer")
		nLines = nLines+1
		if ( intArr(1).ne.np .or. intArr(2).ne.nsd ) then
			write(*,'(A)') "Error: inconsistent no of nodes and nsd.."
			write(*,'(2A)') "Command: ",trim(hdr)
			write(*,'(8A)') "nodes nsd: ", &
			trim(STR(np))," ",trim(STR(nsd))," ", &
			trim(STR(intArr(1)))," ",trim(STR(intArr(2)))
			STOP
		end if
		
		nArr = np*nsd
		call resetArr
		call readData(fid,"double")
		allocate(coords(nsd,np)); coords=0.0d0
		cnt = 0
		do j=1, nsd
			do i=1, np
				cnt = cnt+1
				coords(j,i) = dbleArr(cnt)
			end do
		end do
		nBytesRead = nBytesRead + nArr*kind(dbleArr)
		
		! connectivity interior data !
		write(hdr,'(A)') "connectivity interior linear tetrahedron"
		nArr = 7
		call resetArr
		call readHeader(fid,hdr,"integer")
		nLines = nLines+1
		nelblk	= intArr(1)
		bnen 	= intArr(2)
		bpoly	= intArr(3)
		bnsh	= intArr(4)
		bnshlb	= intArr(5)
		bnenbl	= intArr(6)
		blcsyst	= intArr(7)
		
		if ( nel.ne.intArr(1) .or. nen.ne.intArr(4) ) then
			write(*,'(A)') "Error: inconsistent no of elems and nen.."
			write(*,'(2A)') "Command: ",trim(hdr)
			write(*,'(8A)') "nel nen: ", &
			trim(STR(nel))," ",trim(STR(nen))," ", &
			trim(STR(intArr(1)))," ",trim(STR(intArr(4)))
			STOP
		end if
		
		nArr = nel*nen
		call resetArr
		call readData(fid,"integer")
		allocate(connec(nen,nel)); connec=0
		cnt = 0
		do j=1, nen
			do i=1, nel
				cnt = cnt+1
				connec(j,i) = intArr(cnt)
			end do
		end do
		nBytesRead = nBytesRead + nArr*kind(intArr)
		
		! connectivity boundary data !
		write(hdr,'(A)') "connectivity boundary linear tetrahedron"
		nArr = 8
		call resetArr
		call readHeader(fid,hdr,"integer")
		nLines = nLines+1
		nelblkb	= intArr(1)
		numNBC	= intArr(8)

		if ( nelb.ne.intArr(1) .or. nenb.ne.intArr(5) ) then
			write(*,'(A)') &
			"Error: inconsistent no of boundary elems and nenb.."
			write(*,'(2A)') "Command: ",trim(hdr)
			write(*,'(8A)') "nelb nenb: ", &
			trim(STR(nelb))," ",trim(STR(nenb))," ", &
			trim(STR(intArr(1)))," ",trim(STR(intArr(5)))
			STOP
		end if
		
		nArr = nelb*nenb
		call resetArr
		call readData(fid,"integer")
		allocate(connecb(nenb,nelb)); connecb=0
		cnt = 0
		do j=1, nenb
			do i=1, nelb
				cnt = cnt+1
				connecb(j,i) = intArr(cnt)
			end do
		end do
		nBytesRead = nBytesRead + nArr*kind(intArr)
		
		! boundary element ids !
		write(hdr,'(A)') "ienb to sms linear tetrahedron"
		nArr = 1
		call resetArr
		call readHeader(fid,hdr,"integer")
		nLines = nLines+1

		nArr = nelb
		call resetArr
		call readData(fid,"integer")
		allocate(belIDs(nelb)); belIDs=0
		belIDs(:) = intArr(:)
		nBytesRead = nBytesRead + nArr*kind(intArr)
		
		! bc mapping array !
		write(hdr,'(A)') "bc mapping array"
		nArr = 1
		call resetArr
		call readHeader(fid,hdr,"integer")
		nLines = nLines+1
		
		nArr = np
		call resetArr
		call readData(fid,"integer")
		allocate(bcMap(np)); bcMap=0
		bcMap(:) = intArr(:)
		nBytesRead = nBytesRead + nArr*kind(intArr)
		
		! bc codes array !
		write(hdr,'(A)') "bc codes array"
		nArr = 1
		call resetArr
		call readHeader(fid,hdr,"integer")
		numEBC = intArr(1)
		nLines = nLines+1

		nArr = numEBC
		call resetArr
		call readData(fid,"integer")
		allocate(bcCodes(numEBC)); bcCodes=0
		bcCodes(:) = intArr(:)
		nBytesRead = nBytesRead + nArr*kind(intArr)
		
		! boundary condition array !
		write(hdr,'(A)') "boundary condition array"
		nArr = 1
		call resetArr
		call readHeader(fid,hdr,"integer")
		nLines = nLines+1
		
		nArr = numEBC*12
		call resetArr
		call readData(fid,"double")
		allocate(bcInp(12,numEBC)); bcInp=0.0d0
		do j=1, 12
			do i=1, numEBC
				cnt = cnt+1
				bcInp(j,i) = dbleArr(cnt)
			end do
		end do
		nBytesRead = nBytesRead + nArr*kind(dbleArr)
		
		! periodic masters array !
		write(hdr,'(A)') "periodic masters array"
		nArr = 1
		call resetArr
		call readHeader(fid,hdr,"integer")
		nLines = nLines+1
		
		nArr = np
		call resetArr
		call readData(fid,"integer")
		allocate(iPer(np)); iPer=0
		iPer(:) = intArr(:)
		nBytesRead = nBytesRead + nArr*kind(intArr)
		
		! xadj !
		write(hdr,'(A)') "keyword xadj"
		nArr = 2
		call resetArr
		call readHeader(fid,hdr,"integer")
		nLines = nLines+1
		nxadj = intArr(1) + 1
		
		nArr = nxadj
		call resetArr
		call readData(fid,"integer")
		allocate(xadj(nxadj)); xadj=0
		xadj(:) = intArr(:)
		nBytesRead = nBytesRead + nArr*kind(intArr)
		
		! adjncy !
		write(hdr,'(A)') "keyword adjncy"
		nArr = 1
		call resetArr
		call readHeader(fid,hdr,"integer")
		nLines = nLines+1
		nadjncy = (nfaces-nenb)*2
		
		nArr = nadjncy
		call resetArr
		call readData(fid,"integer")
		allocate(adjncy(nadjncy)); adjncy=0
		adjncy(:) = intArr(:)
		nBytesRead = nBytesRead + nArr*kind(intArr)
		
		! vwgt !
		write(hdr,'(A)') "keyword vwgt"
		nArr = 1
		call resetArr
		call readHeader(fid,hdr,"integer")
		nLines = nLines+1
		nvwgt = nen
		
		nArr = nvwgt
		call resetArr
		call readData(fid,"integer")
		allocate(vwgt(nvwgt)); vwgt=0
		vwgt(:) = intArr(:)
		nBytesRead = nBytesRead + nArr*kind(intArr)
		
		! close file !
		rt = gzClose(fid)

		write(*,'(A)')
		write(*,'(A)') "Summary:"
		write(*,'(2A)') "No of lines read: ", trim(STR(nLines))
		write(*,'(3A)') "No of bytes read: ", trim(STR(nBytesRead))
		write(*,'(A)')

		end program readGeomBC
		
!**************************************************

		subroutine resetArr
		use var
		implicit none
		
		if ( allocated(intArr) ) deallocate(intArr)
		if ( allocated(dbleArr) ) deallocate(dbleArr)
		
		allocate(intArr(nArr))
		allocate(dbleArr(nArr))
		
		intArr = 0
		dbleArr = 0.0d0
		
		end subroutine resetArr

!**************************************************

		subroutine readHeader(fid,phrase,dataType)
		use modGzip
		use var
		use allfun
		implicit none
		integer, intent(in) :: fid
		character(len=*), intent(in) :: phrase,dataType
		
		character(len=stdl) :: token,dlms
		character(len=stdl), target :: rLine
		integer :: cnt, skip, slen, fsweep
		integer(kind=c_int) :: rt
		
		integer, target :: magicNo
		character, target :: c
		type(c_ptr) :: c_ptr1
		
		dlms = ''
		token = ''
		fsweep = 1
		do
			do
				c_ptr1 = gzGets(fid,c_loc(rLine(1:1)),stdl)
				if ( .not.c_associated(c_ptr1) ) exit
				
				if( rLine(1:1).eq."#" .or. len(trim(rLine)).le.1 ) cycle
				
				dlms = ':'
				token = strtok(trim(rLine), trim(dlms))
				token = adjustl(token)
				if ( trim(phrase).eq.trim(token) ) then
					write(*,'(2A)',advance='no') trim(token)," : "
					dlms = " ;<>"
					token = strtok(eol,dlms)
					slen = len(trim(adjustl(token)))
					read(token(1:slen),'(I10)') skip
					write(*,'(3A)',advance='no') "< ", trim(STR(skip))," > "
					cnt = 0
					do while( token.ne.eol .and. cnt.lt.nArr )
						cnt = cnt+1
						token = strtok(eol,dlms)
						slen = len(trim(adjustl(token)))
						if ( trim(dataType).eq."integer" ) then
							read(token(1:slen),*) intArr(cnt)
							write(*,'(2A)',advance='no') &
							trim(STR(intArr(cnt)))," "
						else if ( trim(dataType).eq."double" ) then
							read(token(1:slen),*) dbleArr(cnt)
							write(*,'(2A)',advance='no') &
							STR(dbleArr(cnt))
						else
							write(*,'(A)')
							write(*,'(4X,A)') "Unknown data type.."
							STOP
						end if
					end do
					write(*,'(A)')
					return
				end if
				
				if ( trim(token).eq."byteorder magic number" .and. &
					fsweep.eq.1 ) then
					write(*,'(3A)',advance='no') trim(token)," : "
					rt = gzRead(fid,c_loc(magicNo),kind(magicNo))
					write(*,'(A)') trim(STR(magicNo))
					if ( magicNo.ne.svMagicNo ) then
						write(*,'(A)') "Error: file magic num does not &
						match with solver magic num.."
						STOP
					end if
					!write(*,'(3X,A)') "Magic num matches with file =>  &
					!correct endian.."
					rt = gzRead(fid,c_loc(c),kind(c)) ! new line char !
				end if
			end do
			
			if ( fsweep.eq.1 ) then
				write(*,'(3A)') "Warning: could not find key phrase &
				<", trim(phrase),"> in first sweep"
				write(*,'(A)') "Rewinding file again.."
				rt = gzRewind(fid)
			end if
			fsweep = fsweep + 1
			if ( fsweep.gt.2) exit
		end do
		
		write(*,'(3A)') "Error: key phrase <", trim(phrase), &
		"> does not exist"
		STOP
		
		end subroutine readHeader
		
!**************************************************

		subroutine readData(fid,dataType)
		use modGzip
		use var
		use allfun
		implicit none
		integer, intent(in) :: fid
		character(len=*) :: dataType
		character, target :: c
		integer(kind=c_int) :: rt
		
		if ( trim(dataType).eq."integer" ) then
			rt = gzRead(fid,c_loc(intArr(1)),kind(intArr)*nArr)
		else if ( trim(dataType).eq."double" ) then
			rt = gzRead(fid,c_loc(dbleArr(1)),kind(dbleArr)*nArr)
		else
			write(*,'(A)') "Error: unknown type"
			STOP
		end if
		
		rt = gzRead(fid,c_loc(c),kind(c)) ! new line char !

		end subroutine readData
		
!**************************************************


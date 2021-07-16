!**************************************************
        
        module typeParams
        integer, parameter :: IK1 = selected_int_kind(2)                  ! integer 8 bits
        integer, parameter :: IK2 = selected_int_kind(4)                  ! integer 16 bits
        integer, parameter :: IK4 = selected_int_kind(9)                  ! integer 32 bits
        integer, parameter :: IK8 = selected_int_kind(18)                 ! integer 64 bits
        integer, parameter :: IK = IK4                                  ! default integer type

        integer, parameter :: RK4 = selected_real_kind(6,37)            ! real 32 bits (single)
        integer, parameter :: RK8 = selected_real_kind(15,307)          ! real 64 bits (double)
        integer, parameter :: RK16 = selected_real_kind(33,4931)         ! real 128 bits (long double)
        integer, parameter :: RK = RK8                                  ! default real type (double)
        
        ! interface procedure to transfer bits (type cast) 
        ! from any type to default type
        interface transferBits
            module procedure ::    trBitsIK1, trBitsIK1A, &
                                trbitsIK2, trBitsIK2A, &
                                trBitsIK4, trBitsIK4A, &
                                trBitsIK8, trBitsIK8A, &
                                trBitsRK4, trBitsRK4A, &
                                trBitsRK8, trBitsRK8A
        end interface transferBits
        
        contains
        
            !==========================================
            
            subroutine trBitsIK1(p1,np1,ikind,res)
            implicit none
            integer(IK4), intent(in) :: np1
            integer(IK1), intent(in) :: p1(np1)
            integer(IK1), intent(in) :: ikind
            integer(IK4), intent(out) :: res
            integer(IK1) :: intK1
            
            intK1 = transfer(p1,intK1)
            res = int(intK1, kind=IK4)
            end subroutine trBitsIK1
            
            !==========================================
            
            subroutine trBitsIK2(p1,np1,ikind,res)
            implicit none
            integer(IK4), intent(in) :: np1
            integer(IK1), intent(in) :: p1(np1)
            integer(IK2), intent(in) :: ikind
            integer(IK4), intent(out) :: res
            integer(IK2) :: intK2
            
            intK2 = transfer(p1,intK2)
            res = int(intK2, kind=IK4)
            end subroutine trBitsIK2
            
            !==========================================
            
            subroutine trBitsIK4(p1,np1,ikind,res)
            implicit none
            integer(IK4), intent(in) :: np1
            integer(IK1), intent(in) :: p1(np1)
            integer(IK4), intent(in) :: ikind
            integer(IK4), intent(out) :: res
            
            res = transfer(p1,res)
            end subroutine trBitsIK4
            
            !==========================================
            
            subroutine trBitsIK8(p1,np1,ikind,res)
            implicit none
            integer(IK4), intent(in) :: np1
            integer(IK1), intent(in) :: p1(np1)
            integer(IK8), intent(in) :: ikind
            integer(IK4), intent(out) :: res
            integer(IK8) :: intK8
            
            intK8 = transfer(p1,intK8)
            res = int(intK8,kind=IK8)
            end subroutine trBitsIK8
            
            !==========================================
            
            subroutine trBitsRK4(p1,np1,ikind,res)
            implicit none
            integer(IK4), intent(in) :: np1
            integer(IK1), intent(in) :: p1(np1)
            integer(IK4), intent(in) :: ikind
            real(RK8), intent(out) :: res
            real(RK4) :: realK4
            
            realK4 = transfer(p1,realK4)
            res = real(realK4,kind=RK8)
            end subroutine trBitsRK4
            
            !==========================================
            
            subroutine trBitsRK8(p1,np1,ikind,res)
            implicit none
            integer(IK4), intent(in) :: np1
            integer(IK1), intent(in) :: p1(np1)
            integer(IK8), intent(in) :: ikind
            real(RK8), intent(out) :: res
            
            res = transfer(p1,res)
            end subroutine trBitsRK8
            
            !==========================================
            
            subroutine trBitsIK1A(p1,np1,ikind,res,nout)
            implicit none
            integer(IK4), intent(in) :: np1
            integer(IK1), intent(in) :: p1(np1)
            integer(IK1), intent(in) :: ikind
            integer(IK4), intent(in) :: nout
            integer(IK4), intent(out) :: res(nout)
            integer(IK1), dimension(:), allocatable :: intK1
            
            allocate(intK1(nout)); intK1=0_IK1
            intK1 = transfer(p1,intK1)
            res = int(intK1,kind=IK1)
            deallocate(intK1)
            end subroutine trBitsIK1A
            
            !==========================================
            
            subroutine trBitsIK2A(p1,np1,ikind,res,nout)
            implicit none
            integer(IK4), intent(in) :: np1
            integer(IK1), intent(in) :: p1(np1)
            integer(IK2), intent(in) :: ikind
            integer(IK4), intent(in) :: nout
            integer(IK4), intent(out) :: res(nout)
            integer(IK2), dimension(:), allocatable :: intK2
            
            allocate(intK2(nout)); intK2=0_IK2
            intK2 = transfer(p1,intK2)
            res = int(intK2,kind=IK2)
            deallocate(intK2)
            end subroutine trBitsIK2A
            
            !==========================================
            
            subroutine trBitsIK4A(p1,np1,ikind,res,nout)
            implicit none
            integer(IK4), intent(in) :: np1
            integer(IK1), intent(in) :: p1(np1)
            integer(IK4), intent(in) :: ikind
            integer(IK4), intent(in) :: nout
            integer(IK4), intent(out) :: res(nout)
            
            res=0_IK4
            res = transfer(p1,res)
            end subroutine trBitsIK4A
            
            !==========================================
            
            subroutine trBitsIK8A(p1,np1,ikind,res,nout)
            implicit none
            integer(IK4), intent(in) :: np1
            integer(IK1), intent(in) :: p1(np1)
            integer(IK8), intent(in) :: ikind
            integer(IK4), intent(in) :: nout
            integer(IK4), intent(out) :: res(nout)
            integer(IK8), dimension(:), allocatable :: intK8
            
            allocate(intK8(nout)); intK8=0_IK8
            intK8 = transfer(p1,intK8)
            res = int(intK8,kind=IK8)
            deallocate(intK8)
            end subroutine trBitsIK8A
            
            !==========================================
            
            subroutine trBitsRK4A(p1,np1,ikind,res,nout)
            implicit none
            integer(IK4), intent(in) :: np1
            integer(IK1), intent(in) :: p1(np1)
            integer(IK4), intent(in) :: ikind
            integer(IK4), intent(in) :: nout
            real(RK8), intent(out) :: res(nout)
            real(RK4), dimension(:), allocatable :: realK4
            
            allocate(realK4(nout)); realK4=0._RK4
            realK4 = transfer(p1,realK4)
            res = real(realK4,kind=RK8)
            deallocate(realK4)
            end subroutine trBitsRK4A
            
            !==========================================
            
            subroutine trBitsRK8A(p1,np1,ikind,res,nout)
            implicit none
            integer(IK4), intent(in) :: np1
            integer(IK1), intent(in) :: p1(np1)
            integer(IK8), intent(in) :: ikind
            integer(IK4), intent(in) :: nout
            real(RK8), intent(out) :: res(nout)
            
            res=0._RK8
            res = transfer(p1,res)
            end subroutine trBitsRK8A
            
            !==========================================
        
        end module typeParams
        
!**************************************************
        
        module stdParams
        use typeParams
        character, parameter :: eol=achar(0)                            ! end ofline char
        character, parameter :: newl=achar(10)                            ! new line char
        character(len=64), parameter :: b64List = &                        ! base 64 char list
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"
        
        character(len=8), parameter :: ftab1="(4X,A)"
        character(len=8), parameter :: ftab2="(8X,A)"
        character(len=9), parameter :: ftab3="(12X,A)"
        character(len=9), parameter :: ftab4="(14X,A)"

        integer(IK), parameter :: stdl=400                                ! string length
        integer(IK), parameter :: stdout=6                                ! output to screen
        integer(IK), parameter :: maxToks=30                            ! parameter to set max. string tokens
        logical, parameter :: debug=.false.                                ! run in debug mode if .true.
        end module stdParams
        
!**************************************************

        module genUtils
        use typeParams
        use stdParams, only : stdl,eol,stdout,maxToks
        implicit none
        
        interface STR
            module procedure :: ITSTR, RTSTR, DTSTR, NDTSTR
        end interface STR

        contains

            !==========================================

            subroutine parse(str,toks,ntoks)
            implicit none
            character(len=*), intent(in) :: str
            character(len=*), dimension(maxToks), intent(out) :: toks
            integer(IK), intent(out) :: ntoks
            
            character(len=stdl) :: dlms,token
            
            dlms = ''
            token = ''
            
            dlms = '< =">'
            ntoks = 1
            toks(1) = strtok(trim(str), trim(dlms))
            toks(1) = adjustl(toks(1))
            do
                token = strtok(eol,dlms)
                if ( token.ne.eol ) then
                    ntoks = ntoks+1
                    toks(ntoks) = adjustl(token)
                else
                    exit
                end if
            end do
            
            end subroutine parse
        
            !==========================================

            character(len=stdl) function strtok(str, dlms)
            implicit none
            character(len=*), intent(in) :: str
            character(len=*), intent(in) :: dlms
            
            integer(IK) :: ist, iend

            integer(IK), save :: ist0,slen
            character(len=stdl), save :: str0
            
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

            !==========================================

            function getToken(toks,ntoks,kwrd) result (res)
            implicit none
            integer, intent(in) :: ntoks
            character(len=*), dimension(ntoks), intent(in) :: toks
            character(len=*), intent(in) :: kwrd
            character(len=stdl) :: res
            
            integer :: i
            
            res=''
            do i=1, ntoks
                if (trim(toks(i)).eq.trim(kwrd) ) then
                    res = toks(i)
                    return
                end if
            end do
            
            return
            
            end function getToken
            
            !==========================================
            
            function getTokenValue(toks,ntoks,kwrd) result (res)
            implicit none
            integer, intent(in) :: ntoks
            character(len=*), dimension(ntoks), intent(in) :: toks
            character(len=*), intent(in) :: kwrd
            character(len=stdl) :: res
            
            integer :: i
            
            res=''
            do i=1, ntoks
                if (trim(toks(i)).eq.trim(kwrd) ) then
                    res = toks(i+1)
                    return
                end if
            end do
            
            return
            
            end function getTokenValue

            !==========================================

            pure function ITSTR(iVal) result(str)
            implicit none
            integer(IK), intent(in) :: iVal
            integer(IK) :: n,ist,j,k,slen
            character(len=stdl) str
            
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
            
            !==========================================
            
            pure function RTSTR(rVal) result(str)
            implicit none
            integer(IK), parameter :: l=8
            real(RK4), intent(in) :: rVal
            character(len=l) :: str
            
            str = NDTSTR(dble(rVal),l)
            
            return
            end function RTSTR
            
            !==========================================
            
            pure function DTSTR(dVal) result(str)
            implicit none
            integer, parameter :: l=8
            real(RK), intent(in) :: dVal
            character(len=l) :: str
            
            str = NDTSTR(dVal,l)
            
            return
            end function DTSTR
            
            !==========================================
            
            pure function NDTSTR(dVal,l) result(str)
            implicit none
            integer(IK), intent(in) :: l
            real(RK), intent(in) :: dVal
            character(len=l) :: str
            
            integer(IK) :: i,j,k,ipos,cnt,ex,abex,nex
            real(RK) :: absd
            
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

			!==========================================

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

            !==========================================
        
        end module genUtils
    
!**************************************************

        module vtkXMLlib
        use typeParams
        use stdParams, only : stdl
        
        integer(IK), parameter :: nVTKElms=5
        integer(IK), parameter :: nPieceTyps=5
        integer(IK), parameter :: nPieceElms=6
        integer(IK), parameter :: nPieceData=5
        integer(IK), parameter :: nPieceAtts=9
        integer(IK), parameter :: nDataElms=5
        integer(IK), parameter :: nDataTyps=10
        integer(IK), parameter :: nDataFrmt=3
        integer(IK), parameter :: nDataEnc=2
        
        character(len=stdl), dimension(nVTKElms)   :: libVTKElms
        character(len=stdl), dimension(nPieceTyps) :: libVTKPcTyps
        character(len=stdl), dimension(nPieceElms) :: libPieceElms
        character(len=stdl), dimension(nPieceData) :: libPcPtClData
        character(len=stdl), dimension(nPieceAtts) :: libPieceAtts
        character(len=stdl), dimension(nDataElms)  :: libDataElms
        character(len=stdl), dimension(nDataTyps)  :: libDataTyps
        character(len=stdl), dimension(nDataFrmt)  :: libDataFrmt
        character(len=stdl), dimension(nDataFrmt)  :: libDataEnc

        contains
        
            subroutine initVTKXMLlib
            implicit none
            
            libVTKElms=""
            libVTKPcTyps=""
            libPieceElms=""
            libPcPtClData=""
            libPieceAtts=""
            libDataElms=""
            libDataTyps=""
            libDataFrmt=""
            
            ! libVTKElms(5) !
            libVTKElms(1) = "type"
            libVTKElms(2) = "version"
            libVTKElms(3) = "byte_order"
            libVTKElms(4) = "header_type"
            libVTKElms(5) = "compressor"
            
            ! libVTKPcTyps(5) !
            libVTKPcTyps(1) = "ImageData"
            libVTKPcTyps(2) = "RectilinearGrid"
            libVTKPcTyps(3) = "StructuredGrid"
            libVTKPcTyps(4) = "PolyData"
            libVTKPcTyps(5) = "UnstructuredGrid"
            
            ! libPieceElms(6) !
            libPieceElms(1) = "NumberOfPoints"
            libPieceElms(2) = "NumberOfCells"
            libPieceElms(3) = "NumberOfVerts"
            libPieceElms(4) = "NumberOfLines"
            libPieceElms(5) = "NumberOfStrips"
            libPieceElms(6) = "NumberOfPolys"

            ! libPieceData(5) !
            libPcPtClData(1)  = "Scalars"
            libPcPtClData(2)  = "Vectors"
            libPcPtClData(3)  = "Normals"
            libPcPtClData(4)  = "Tensors"
            libPcPtClData(5)  = "Tcoords"
            
            ! libPieceAtts(9) !
            libPieceAtts(1) = "PointData"
            libPieceAtts(2) = "CellData"
            libPieceAtts(3) = "Points"
            libPieceAtts(4) = "Coords"
            libPieceAtts(5) = "Verts"
            libPieceAtts(6) = "Lines"
            libPieceAtts(7) = "Strips"
            libPieceAtts(8) = "Polys"
            libPieceAtts(9) = "Cells"

            ! libDataElms(5) !
            libDataElms(1) = "type"
            libDataElms(2) = "Name"
            libDataElms(3) = "NumberOfComponents"
            libDataElms(4) = "format"
            libDataElms(5) = "offset"
            !libDataElms(6) = "RangeMin"
            !libDataElms(7) = "RangeMax"

            ! libDataTyps(10) !
            libDataTyps(1)  = "Int8"
            libDataTyps(2)  = "UInt8"
            libDataTyps(3)  = "Int16"
            libDataTyps(4)  = "UInt16"
            libDataTyps(5)  = "Int32"
            libDataTyps(6)  = "UInt32"
            libDataTyps(7)  = "Int64"
            libDataTyps(8)  = "UInt64"
            libDataTyps(9)  = "Float32"
            libDataTyps(10) = "Float64"

            ! libDataFrmt(3) !
            libDataFrmt(1) = "ascii"
            libDataFrmt(2) = "binary"
            libDataFrmt(3) = "appended"
            
            ! libDataEnc(2) !
            libDataEnc(1) = "raw"
            libDataEnc(2) = "base64"

            end subroutine initVTKXMLlib
        
        end module vtkXMLlib

!**************************************************

        module vtkXMLMod
        use stdParams
        use typeParams
        use genUtils
        use vtkXMLlib
        implicit none
        
        private
        public :: vtkXMLType
        public :: loadVTK,flushVTK
        public :: getVTK_nodesPerElem
        public :: getVTK_numElems
        public :: getVTK_numPoints
        public :: getVTK_elemIEN
        public :: getVTK_pointCoords
        public :: getVTK_pointData
        public :: getVTK_elemData
        
        interface getVTK_pointData
           module procedure getVTK_pointDataIntS, getVTK_pointDataRealS, &
                            getVTK_pointDataIntV, getVTK_pointDataRealV
        end interface

        interface getVTK_elemData
           module procedure getVTK_elemDataIntS, getVTK_elemDataRealS, &
                            getVTK_elemDataIntV, getVTK_elemDataRealV
        end interface
        
        logical :: flag
        character(len=stdl) :: rLine,stmp
        character(len=stdl) :: startKwrd,stopKwrd
        character(len=stdl), dimension(maxToks) :: tokenList
        character :: c
        integer(IK) :: itok,ntoks,slen,iPos
        integer(IK) :: rank,maxRank
        integer(IK) :: pcStPos,pcEndPos
        
        type dataArrType
            private
            character(len=stdl), dimension(nDataElms)  :: dElms
            character(len=stdl) :: dType,dName,dFrmt,hdrType
            logical :: isInt
            integer(IK) :: hdrKind
            integer(IK) :: iOffst,appRank,ikind,rank
            integer(IK) :: stPos,endPos,nBytes
            integer(IK) :: nElms,nComps,nVals
            integer(IK4), dimension(:), allocatable :: iarr
            real(RK8), dimension(:), allocatable :: darr
        end type dataArrType
        
        type pieceAttType
            private
            integer(IK) :: n,stPos,endPos
            character(len=stdl) :: pName
            character(len=stdl) :: ptClField,ptClFieldName
            type(dataArrType), dimension(:), allocatable :: dataArr
        end type pieceAttType
        
        type vtkXMLType
            private
            logical :: isBinApp
            character(len=stdl) :: fileName
            character(len=stdl) :: dataFormat
            character(len=stdl) :: dataEncdng
            character(len=stdl) :: vtkPcType
            integer(IK) :: fid
            integer(IK) :: stAppendPos
            integer(IK) :: endAppendPos
            integer(IK) :: offsets(100)

            character(len=stdl), dimension(nVTKElms)   :: vtkElms
            integer(IK), dimension(nPieceElms) :: pieceElms
            type(pieceAttType), dimension(nPieceAtts) :: pcAtt
        end type vtkXMLType
        
        contains

            !==========================================
            
            subroutine loadVTK(vtk,fName,istat)
            implicit none
            type(vtkXMLType), intent(inout) :: vtk
            character(len=*), intent(in) :: fName
            integer(IK) :: istat
            
            istat = 0
            inquire(file=trim(fName), exist=flag)
            if ( .not.flag ) then
                write(stdout,ftab4) &
                    "ERROR: File "//trim(fName)//" does not exist"
                istat=-1; return
            end if
            
            slen = len(trim(fName))
            select case (fName(slen-2:slen))
            case ("vtu","vtp")
            case default
                write(stdout,ftab4) &
                    "ERROR: unknown file extension &
                    (can only be vtu or vtp)"
                istat=-1; return
            end select
            
            call initVTKXMLlib
            
            call initVTKXMLstruct(vtk,fName)
            
            call readHeader(vtk,istat)
            if ( istat.lt.0 ) return
            
            call parseVTKKernel(vtk,istat)
            if ( istat.lt.0 ) return
            
            call vtkDataLoader(vtk,istat)
            if ( istat.lt.0 ) return
            
            return
            
            end subroutine loadVTK

            !==========================================
            
            subroutine flushVTK(vtk)
            implicit none
            type(vtkXMLType), intent(inout) :: vtk
            integer(IK) :: iatt,i
            
            do iatt=1, nPieceAtts
                if ( vtk%pcAtt(iatt)%n.lt.1 .or. &
                     .not.allocated(vtk%pcAtt(iatt)%dataArr) ) cycle
                do i=1, vtk%pcAtt(iatt)%n
                    if ( allocated(vtk%pcAtt(iatt)%dataArr(i)%iarr) ) &
                        deallocate(vtk%pcAtt(iatt)%dataArr(i)%iarr)
                    if ( allocated(vtk%pcAtt(iatt)%dataArr(i)%darr) ) &
                        deallocate(vtk%pcAtt(iatt)%dataArr(i)%darr)
                    vtk%pcAtt(iatt)%dataArr(i)%dElms(:) = ""
                    vtk%pcAtt(iatt)%dataArr(i)%dType = ""
                    vtk%pcAtt(iatt)%dataArr(i)%dName = ""
                    vtk%pcAtt(iatt)%dataArr(i)%dFrmt = ""
                    vtk%pcAtt(iatt)%dataArr(i)%hdrType = ""
                end do
            end do
            close(vtk%fid)
            
            write(stdout,ftab1) "Flushed VTK object.."
            
            end subroutine flushVTK
            
            !==========================================

            subroutine initVTKXMLstruct(vtk,fName)
            implicit none
            type(vtkXMLType), intent(inout) :: vtk
            character(len=*), intent(in) :: fName
            integer :: fid
            
            vtk%fileName         = trim(fName)
            vtk%isBinApp         = .false.
            vtk%dataFormat         = ""
            vtk%dataEncdng         = ""
            vtk%vtkPcType         = ""
            vtk%stAppendPos     = 0
            vtk%endAppendPos    = 0
            vtk%offsets(:)      = 0
            vtk%vtkElms(:)         = ""
            vtk%pieceElms(:)     = 0
            vtk%pcAtt(:)%pName     = ""
            
            do fid=11, 1024
                inquire(unit=fid, opened=flag)
                if ( .not.flag ) exit
            end do
            
            vtk%fid = fid
            
            end subroutine initVTKXMLstruct

            !==========================================

            subroutine readHeader(vtk,istat)
            implicit none
            type(vtkXMLType), intent(inout) :: vtk
            integer(IK), intent(inout) :: istat
            
            rLine = ""
            stmp = ""
            tokenList(:) = ""
            iPos = 0
            
            write(stdout,ftab1) &
                "Reading file: <"//trim(vtk%fileName)//">"

            open(vtk%fid,file=trim(vtk%fileName),status='unknown')
            read(vtk%fid,'(A)') rLine
            rLine = adjustl(rLine)
            if ( rLine(2:5).eq."?xml" ) then
                read(vtk%fid,'(A)') rLine
            end if
            
            ! get data format and encoding !
            call findKwrdXML(vtk,"<DataArray",iPos,rLine,istat)
            if ( istat.lt.0 ) return
            iPos = 0
            call parse(rLine,tokenList,ntoks)
            stmp = getTokenValue(tokenList,ntoks,"format")
            if ( trim(stmp).eq."binary" .or. &
                 trim(stmp).eq."appended" ) then
                close(vtk%fid)
                open(vtk%fid,file=trim(vtk%fileName),form="unformatted", &
                access="stream",convert="big_endian")
                vtk%isBinApp = .true.
                iPos = 1
            end if
            vtk%dataFormat = trim(stmp)
            write(stdout,ftab2) &
                "Data format: <"//trim(vtk%dataFormat)//">"
            
            vtk%dataEncdng = "base64"
            if ( vtk%dataFormat.eq."appended" ) then
                call findKwrdXML(vtk,"<AppendedData",iPos,rLine,istat)
                if ( istat.lt.0 ) return
                do
                    read(vtk%fid,pos=iPos,end=001) c
                    iPos = iPos + 1
                    if ( c.eq."_" ) exit
                end do
                vtk%stAppendPos = iPos
                call parse(rLine,tokenList,ntoks)
                stmp = getTokenValue(tokenList,ntoks,"encoding")
                vtk%dataEncdng = trim(stmp)
                call findKwrdXML(vtk,"</AppendedData",iPos,rLine,istat)
                if ( istat.lt.0 ) return
                vtk%endAppendPos = iPos - len(trim(rLine))
            end if
            write(stdout,ftab2) &
                "Data encoding: <"//trim(vtk%dataEncdng)//">"
                
            if ( debug ) then
                if ( vtk%dataFormat.eq."appended" ) then
                    write(stdout,ftab2) &
                        "Data begins at pos: "//trim(STR(vtk%stAppendPos))
                    write(stdout,ftab2) &
                        "Data ends at pos: "// trim(STR(vtk%endAppendPos))
                end if
            end if
            
            iPos = 0
            if ( vtk%isBinApp ) iPos = 1
            rewind(vtk%fid)
            return
            
 001        write(stdout,ftab4) "ERROR: end of file reached.."
            istat = -1; return
            
            end subroutine readHeader
            
            !==========================================

            subroutine parseVTKKernel(vtk,istat)
            implicit none
            type(vtkXMLType), intent(inout) :: vtk
            integer(IK), intent(inout) :: istat
                        
            rLine = ""
            stmp = ""
            tokenList(:) = ""
            
            call findKwrdXML(vtk,"<VTKFile",iPos,rLine,istat)
            if ( istat.lt.0 ) return
            if ( debug ) write(stdout,ftab2) trim(rLine)
            call parse(rLine,tokenList,ntoks)
            
            do itok=1, nVTKElms
                stmp = getTokenValue(tokenList,ntoks,libVTKElms(itok))
                vtk%vtkElms(itok) = trim(stmp)
            end do
            
            do itok=1, nPieceTyps
                if ( trim(vtk%vtkElms(1)).eq.trim(libVTKPcTyps(itok)) ) exit
            end do
            if ( itok.le.3 ) then
                write(stdout,ftab4) &
                    "ERROR: unknown piece type"
                write(stdout,ftab4) &
                    "Piece <"//trim(tokenList(1))//'>'
                write(stdout,ftab4) &
                    "Piece can be <UnstructuredGrid> or <PolyData> only"
                istat=-1; return
            end if
            vtk%vtkPcType = libVTKPcTyps(itok)
            if ( len(trim(vtk%vtkElms(5))).gt.0 ) &
                write(stdout,ftab2) &
                    "Data compression: "//trim(vtk%vtkElms(5))
                    
            ! Piece elements parser !
            call findKwrdXML(vtk,"<Piece",iPos,rLine,istat)
            if ( istat.lt.0 ) return
            if ( debug ) write(stdout,ftab2) trim(rLine)
            call parse(rLine,tokenList,ntoks)
            
            do itok=1, nPieceElms
                stmp = getTokenValue(tokenList,ntoks,libPieceElms(itok))
                slen = len(trim(stmp))
                if ( slen.gt.0 ) then
                    read(stmp(1:slen),*) vtk%pieceElms(itok)
                end if
            end do
            
            pcStPos = iPos! - len(trim(rLine))
            call findKwrdXML(vtk,"</Piece",iPos,rLine,istat)
            if ( istat.lt.0 ) return
            if ( vtk%isBinApp ) then
                pcEndPos = iPos - len(trim(rLine))
            else
                pcEndPos = iPos - 1
            end if
            
            call readPieceAttributes(vtk,istat)
            
            return
            end subroutine parseVTKKernel
            
            !==========================================

            subroutine readPieceAttributes(vtk,istat)
            implicit none
            type(vtkXMLType), intent(inout) :: vtk
            integer(IK), intent(inout) :: istat
            integer :: cntr,iatt,i
            
            rLine = ""
            stmp = ""
            tokenList(:) = ""
            cntr = 0
            
            do iatt=1, nPieceAtts
                call resetFilePos(vtk%fid,vtk%fileName,pcStPos,vtk%isBinApp)
                iPos = pcStPos
                call findKwrdXML(vtk,"<"//libPieceAtts(iatt),iPos,rLine,istat,pcEndPos)
                if ( istat.lt.0 ) return
                if ( len(trim(rLine)).lt.1 ) cycle
                if ( debug ) write(stdout,ftab2) trim(rLine)

                call parse(rLine,tokenList,ntoks)
                stmp = getToken(tokenList,ntoks,libPieceAtts(iatt))
                vtk%pcAtt(iatt)%pName = trim(stmp)
                
                if ( iatt.le.2 ) then
                    do itok=1, nPieceData
                        stmp = getToken(tokenList,ntoks,libPcPtClData(itok))
                        vtk%pcAtt(iatt)%ptClField = trim(stmp)
                        stmp = getTokenValue(tokenList,ntoks,libPcPtClData(itok))
                        vtk%pcAtt(iatt)%ptClFieldName = trim(stmp)
                    end do
                else
                    vtk%pcAtt(iatt)%ptClField = ""
                    vtk%pcAtt(iatt)%ptClFieldName = ""
                end if

                vtk%pcAtt(iatt)%stPos = iPos! - len(trim(rLine))
                call findKwrdXML(vtk,"</"//libPieceAtts(iatt),iPos,rLine,istat,pcEndPos)
                if ( istat.lt.0 ) return
                if ( vtk%isBinApp ) then
                    vtk%pcAtt(iatt)%endPos= iPos - len(trim(rLine))
                else
                    vtk%pcAtt(iatt)%endPos= iPos - 1
                end if
                
                vtk%pcAtt(iatt)%n = 0
                iPos = vtk%pcAtt(iatt)%stPos
                if ( .not.vtk%isBinApp ) &
                    call resetFilePos(vtk%fid,vtk%fileName,iPos,vtk%isBinApp)
                
                ! Count DataArray elements within a Piece attribute !
                if ( debug ) then
                    write(stdout,ftab3) &
                        "Counting <DataArray> elements in <"//&
                        trim(vtk%pcAtt(iatt)%pName)//">"
                end if
                do
                    call findKwrdXML(vtk,"<DataArray",iPos,rLine,istat,vtk%pcAtt(iatt)%endPos)
                    if ( istat.lt.0 ) return
                    if ( iPos.ge.vtk%pcAtt(iatt)%endPos ) exit
                    if ( debug ) write(stdout,ftab3) trim(rLine)
                    call parse(rLine,tokenList,ntoks)
                    if ( trim(tokenList(1)).eq."DataArray" ) &
                        vtk%pcAtt(iatt)%n = vtk%pcAtt(iatt)%n + 1
                end do ! inner loop over DataArray elms !
                if ( vtk%pcAtt(iatt)%n.eq.0 ) then
                    if ( debug ) write(stdout,ftab3) "None found.."
                    cycle
                end if

                if ( debug ) then
                    write(stdout,ftab3) "No. of <DataArray> elems: "//&
                        trim(STR(vtk%pcAtt(iatt)%n))
                end if
                
                allocate(vtk%pcAtt(iatt)%dataArr(vtk%pcAtt(iatt)%n))
                iPos = vtk%pcAtt(iatt)%stPos
                call resetFilePos(vtk%fid,vtk%fileName,iPos,vtk%isBinApp)
                if ( debug ) then
                    write(stdout,ftab2) &
                        "Reading <DataArray> elements in <"//&
                        trim(vtk%pcAtt(iatt)%pName)//">"
                end if
                do i=1, vtk%pcAtt(iatt)%n
                    call findKwrdXML(vtk,"<DataArray",iPos,rLine,istat)
                    if ( istat.lt.0 ) return
                    vtk%pcAtt(iatt)%dataArr(i)%stPos = iPos+1
                    
                    if ( debug ) write(stdout,ftab3) trim(rLine)
                    
                    call parse(rLine,tokenList,ntoks)
                    do itok=1, nDataElms
                        stmp = getTokenValue(tokenList,ntoks,libDataElms(itok))
                        vtk%pcAtt(iatt)%dataArr(i)%dElms(itok) = trim(stmp)
                    end do
                    
                    if ( vtk%dataFormat.ne."appended" ) then
                        call findKwrdXML(vtk,"</DataArray>",iPos,rLine,istat)
                        if ( istat.lt.0 ) return
                        if ( vtk%isBinApp ) then
                            vtk%pcAtt(iatt)%dataArr(i)%endPos = &
                                iPos - len(trim(rLine))
                        else
                            vtk%pcAtt(iatt)%dataArr(i)%endPos = &
                                iPos
                        end if
                    end if
                    
                    cntr = cntr+1
                    vtk%pcAtt(iatt)%dataArr(i)%rank = cntr
                    
                    vtk%pcAtt(iatt)%dataArr(i)%dType = &
                        trim(vtk%pcAtt(iatt)%dataArr(i)%dElms(1))
                    vtk%pcAtt(iatt)%dataArr(i)%dName = &
                        trim(vtk%pcAtt(iatt)%dataArr(i)%dElms(2))
                    
                    slen = len(trim(vtk%pcAtt(iatt)%dataArr(i)%dElms(3)))
                    if ( slen.gt.0 ) then
                        read(vtk%pcAtt(iatt)%dataArr(i)%dElms(3)(1:slen),*) &
                            vtk%pcAtt(iatt)%dataArr(i)%nComps
                    else
                        vtk%pcAtt(iatt)%dataArr(i)%nComps = 0 !1
                    end if
                    
                    vtk%pcAtt(iatt)%dataArr(i)%dFrmt = &
                        trim(vtk%pcAtt(iatt)%dataArr(i)%dElms(4))
                    
                    slen = len(trim(vtk%pcAtt(iatt)%dataArr(i)%dElms(5)))
                    if ( slen.gt.0 ) then
                        read(vtk%pcAtt(iatt)%dataArr(i)%dElms(5)(1:slen),*) &
                            vtk%pcAtt(iatt)%dataArr(i)%iOffst
                    else
                        vtk%pcAtt(iatt)%dataArr(i)%iOffst = 0
                    end if
                    
                    if ( vtk%isBinApp ) then
                        if ( vtk%dataFormat.eq."binary" ) then
                            call adjustdataArray(vtk%fid, &
                            vtk%pcAtt(iatt)%dataArr(i)%stPos, &
                            vtk%pcAtt(iatt)%dataArr(i)%endPos )

                            iPos = vtk%pcAtt(iatt)%dataArr(i)%stPos
                            vtk%pcAtt(iatt)%dataArr(i)%nbytes = &
                                vtk%pcAtt(iatt)%dataArr(i)%endPos - &
                                vtk%pcAtt(iatt)%dataArr(i)%stPos
                        else if ( vtk%dataFormat.eq."appended" ) then
                            vtk%pcAtt(iatt)%dataArr(i)%stPos = &
                                vtk%pcAtt(iatt)%dataArr(i)%iOffst + &
                                vtk%stAppendPos
                            vtk%offsets(vtk%pcAtt(iatt)%dataArr(i)%rank) = &
                                vtk%pcAtt(iatt)%dataArr(i)%stPos
                        end if
                    else
                        vtk%pcAtt(iatt)%dataArr(i)%nbytes = &
                            vtk%pcAtt(iatt)%dataArr(i)%endPos - &
                            vtk%pcAtt(iatt)%dataArr(i)%stPos
                    end if
                    
                end do ! idata
            end do ! outer loop over Piece atts !
            maxRank = cntr
            
            return
            end subroutine readPieceAttributes
            
            !==========================================

            subroutine vtkDataLoader(vtk,istat)
            implicit none
            type(vtkXMLType), intent(inout) :: vtk
            integer(IK), intent(inout) :: istat
            integer(IK) :: iatt,i,rank

            if ( .not.debug ) write(stdout,ftab2) &
            "Piece Type: "//trim(vtk%vtkPcType)
            do iatt=1, nPieceAtts
                if ( vtk%pcAtt(iatt)%n.lt.1 .or. &
					.not.allocated(vtk%pcAtt(iatt)%dataArr) ) cycle
                if ( .not.debug ) write(stdout,ftab3) &
                    "Piece Attribute: "//trim(vtk%pcAtt(iatt)%pName)
                do i=1, vtk%pcAtt(iatt)%n
                    if ( vtk%dataFormat.eq."appended" ) then
                        rank = vtk%pcAtt(iatt)%dataArr(i)%rank
                        do
                            rank = rank+1
                            if ( rank.gt.maxRank ) then
                                ipos = vtk%endAppendPos-1
                                do
                                    read(vtk%fid,pos=iPos,end=001) c
                                    if ( c.eq.' ' .or. c.eq.eol .or. &
                                        c.eq.'    ') then
                                        iPos = iPos-1 
                                        cycle
                                    else
                                        exit
                                    end if
                                end do
                                vtk%pcAtt(iatt)%dataArr(i)%endPos = iPos
                                exit
                            end if
                                
                            if ( vtk%offsets(rank).gt. &
                                 vtk%pcAtt(iatt)%dataArr(i)%stPos ) then
                                 vtk%pcAtt(iatt)%dataArr(i)%endPos = &
                                    vtk%offsets(rank)
                                exit
                            end if
                        end do
                        vtk%pcAtt(iatt)%dataArr(i)%nbytes = &
                            vtk%pcAtt(iatt)%dataArr(i)%endPos - &
                                vtk%pcAtt(iatt)%dataArr(i)%stPos
                    end if ! appended
                    
                    if ( vtk%dataFormat.eq."ascii" ) then
                        call resetFilePos(vtk%fid,vtk%fileName,&
                            vtk%pcAtt(iatt)%dataArr(i)%stPos-1,vtk%isBinApp)
                    end if
                    
                    if ( .not.debug ) write(stdout,ftab4) &
                        "DataArray name: "// &
                        trim(vtk%pcAtt(iatt)%dataArr(i)%dName)

                    call vtkXMLDataParser(vtk,vtk%pcAtt(iatt)%dataArr(i),iatt,i,istat)
                    if ( istat.lt.0 ) return
                end do
            end do
            
            return
            
 001        write(stdout,ftab4) "ERROR: end of file reached.."
            istat = -1; return
            
            end subroutine vtkDataLoader
            
            !==========================================
            
            subroutine vtkXMLDataParser(vtk,dataArr,iatt,idata,istat)
            implicit none
            type(vtkXMLType), intent(inout) :: vtk
            type(dataArrType), intent(inout), target :: dataArr
            integer(IK), intent(in) :: iatt,idata
            integer(IK), intent(inout) :: istat
            type(dataArrType), pointer :: dPtr
            
            dPtr => dataArr
            
            if ( len(trim(vtk%vtkElms(4))).gt.0 ) then
                dPtr%hdrType = trim(vtk%vtkElms(4))
            else
                dPtr%hdrType = "UInt32"
            end if
            call selectDataType(dPtr%hdrType,dPtr%hdrKind,dPtr%isInt,istat)
            if ( istat.lt.0 ) return
            
            if ( debug ) then
                write(stdout,ftab1) repeat("*",76)
                write(stdout,ftab3) &
                    "Piece Type: "//trim(vtk%vtkPcType)
                write(stdout,ftab3) &
                    "Piece Attribute: "//trim(vtk%pcAtt(iatt)%pName)
                write(stdout,ftab3) &
                    "Header type: "//trim(dPtr%hdrType)
                write(stdout,ftab3) &
                    "DataArray type: "//trim(dPtr%dType)
                write(stdout,ftab3) &
                    "DataArray name: "//trim(dPtr%dName)
                write(stdout,ftab3) &
                    "DataArray frmt: "//trim(dPtr%dFrmt)
                write(stdout,ftab3) &
                    "DataArray nComps: "//trim(STR(dPtr%nComps))
                write(stdout,ftab3) &
                    "DataArray offset: "//trim(STR(dPtr%ioffst))
                write(stdout,ftab3) &
                    "DataArray start pos: "//trim(STR(dPtr%stPos))
                write(stdout,ftab3) &
                    "DataArray end pos:   "//trim(STR(dPtr%endPos))
                write(stdout,ftab3) &
                    "DataArray nbytes:    "//trim(STR(dPtr%nbytes))
            end if
            
            select case (trim(vtk%vtkPcType))
            case ("UnstructuredGrid")
                call vtkParseUnstrucGrid(vtk,dataArr,iatt,idata,istat)
                if ( istat.lt.0 ) return
            
            case ("PolyData")
                call vtkParsePolyData(vtk,dataArr,iatt,idata,istat)
                if ( istat.lt.0 ) return
            
            end select
            
            if ( vtk%isBinApp ) then
                if ( vtk%vtkElms(5).eq."vtkZLibDataCompressor" ) then
                    call readZlibBinaryData(vtk,dPtr,istat)
                else
                    if ( dPtr%nElms.eq.0 ) &
                        return
                    call readBinaryData(vtk,dPtr,istat)
                end if
            else
                if ( dPtr%nElms.eq.0 ) &
                    return
                call readAsciiData(vtk,dPtr,istat)
            end if

            if ( istat.lt.0 ) return
            
            end subroutine vtkXMLDataParser
                
            !==========================================
                
            subroutine vtkParseUnstrucGrid(vtk,dataArr,iatt,idata,istat)
            implicit none
            type(vtkXMLType), intent(inout) :: vtk
            type(dataArrType), intent(inout), target :: dataArr
            integer(IK), intent(in) :: iatt,idata
            integer(IK), intent(inout) :: istat
            type(dataArrType), pointer :: dPtr
            integer(IK) :: nPoints,nCells
            
            dPtr => dataArr

            nPoints = vtk%pieceElms(1)
            nCells  = vtk%pieceElms(2)
            if ( nPoints.lt.1 ) then
                write(stdout,ftab4) &
                "ERROR: VTK Piece element NumberOfPoints not defined.."
                istat=-1; return
            end if
            
            if ( nCells.lt.1 ) then
                write(stdout,ftab4) &
                "ERROR: VTK Piece element NumberOfCells not defined.."
                istat=-1; return
            end if

            dPtr%isInt = .false.
            call selectDataType(dPtr%dType,dPtr%ikind,dPtr%isInt,istat)
            if ( istat.lt.0 ) return
            
            select case (trim(vtk%pcAtt(iatt)%pName))
            case ("PointData")
                if ( nPoints.gt.0 .and. dPtr%nComps.eq.0 ) then
                    dPtr%nComps = getNumComps(vtk,dataArr,nPoints)
                    if ( .not.vtk%isBinApp ) dPtr%nComps = 1
                end if
                dPtr%nVals = nPoints
                dPtr%nElms = nPoints * dPtr%nComps
                if ( debug ) then
                    write(stdout,ftab4) &
                        "nPoints "// trim(STR(dPtr%nVals)) //&
                        "; nComps "// trim(STR(dPtr%nComps)) //&
                        "; nElems "// trim(STR(dPtr%nElms))
                end if
                
            case ("CellData")
                if ( nCells.gt.0 .and. dPtr%nComps.eq.0 ) then
                    dPtr%nComps = getNumComps(vtk,dataArr,nCells)
                    if ( .not.vtk%isBinApp ) dPtr%nComps = 1
                end if
                dPtr%nVals = nCells
                dPtr%nElms = nCells * dPtr%nComps
                if ( debug ) then
                    write(stdout,ftab4) &
                        "nCells "// trim(STR(dPtr%nVals)) //&
                        "; nComps "// trim(STR(dPtr%nComps)) //&
                        "; nElems "// trim(STR(dPtr%nElms))
                end if
            
            case ("Points")
                if ( dPtr%nComps.ne.3 ) then
                    write(stdout,ftab4) &
                    "WARNING: Element <NumberOfComponents> in <Points> &
                     attribute < 3"
                end if
                if ( nPoints.gt.0 .and. dPtr%nComps.eq.0 ) &
                    dPtr%nComps = getNumComps(vtk,dataArr,nPoints)
                dPtr%nVals = nPoints
                dPtr%nElms = nPoints * dPtr%nComps
                if ( debug ) then
                    write(stdout,ftab4) &
                        "nPoints "// trim(STR(dPtr%nVals)) //&
                        "; nComps "// trim(STR(dPtr%nComps)) //&
                        "; nElems "// trim(STR(dPtr%nElms))
                end if
                
            case ("Cells")
                if ( nCells.gt.0 .and. dPtr%nComps.eq.0 ) then
                    dPtr%nComps = getNumComps(vtk,dataArr,nCells)
                    if ( .not.vtk%isBinApp ) dPtr%nComps = 1
                end if
                dPtr%nVals = nCells
                dPtr%nElms = nCells * dPtr%nComps
                if ( debug ) then
                    write(stdout,ftab4) &
                        "nCells "// trim(STR(dPtr%nVals)) //&
                        "; nComps "// trim(STR(dPtr%nComps)) //&
                        "; nElems "// trim(STR(dPtr%nElms))
                end if
            
            end select
            
            end subroutine vtkParseUnstrucGrid
            
            !==========================================
            
            subroutine vtkParsePolyData(vtk,dataArr,iatt,idata,istat)
            implicit none
            type(vtkXMLType), intent(inout) :: vtk
            type(dataArrType), intent(inout), target :: dataArr
            integer(IK), intent(in) :: iatt,idata
            integer(IK), intent(inout) :: istat
            type(dataArrType), pointer :: dPtr
            integer(IK) :: nPoints,nVerts,nLines,nStrips,nPolys
            
            dPtr => dataArr
            
            nPoints = vtk%pieceElms(1)
            nVerts  = vtk%pieceElms(3)
            nLines  = vtk%pieceElms(4)
            nStrips = vtk%pieceElms(5)
            nPolys  = vtk%pieceElms(6)
            
            if ( nPoints.lt.1 ) then
                write(stdout,ftab4) &
                "ERROR: VTK Piece element NumberOfPoints not defined.."
                istat=-1; return
            end if
            
            if ( nPolys.lt.1 ) then
                write(stdout,ftab4) &
                "ERROR: VTK Piece element NumberOfPolys not defined.."
                istat=-1; return
            end if

            dPtr%isInt = .false.
            call selectDataType(dPtr%dType,dPtr%ikind,dPtr%isInt,istat)
            
            select case (trim(vtk%pcAtt(iatt)%pName))
            case ("PointData")
                if ( nPoints.gt.0 .and. dPtr%nComps.eq.0 ) then
                    dPtr%nComps = getNumComps(vtk,dataArr,nPoints)
                    if ( .not.vtk%isBinApp ) dPtr%nComps = 1
                end if
                dPtr%nVals = nPoints
                dPtr%nElms = nPoints * dPtr%nComps
                if ( debug ) then
                    write(stdout,ftab4) &
                        "nPoints "// trim(STR(dPtr%nVals)) //&
                        "; nComps "// trim(STR(dPtr%nComps)) //&
                        "; nElems "// trim(STR(dPtr%nElms))
                end if
                
            case ("CellData")
                if ( nPolys.gt.0 .and. dPtr%nComps.eq.0 ) then
                    dPtr%nComps = getNumComps(vtk,dataArr,nPolys)
                    if ( .not.vtk%isBinApp ) dPtr%nComps = 1
                end if
                dPtr%nVals = nPolys
                dPtr%nElms = nPolys * dPtr%nComps
                if ( debug ) then
                    write(stdout,ftab4) &
                        "nPolys "// trim(STR(dPtr%nVals)) //&
                        "; nComps "// trim(STR(dPtr%nComps)) //&
                        "; nElems "// trim(STR(dPtr%nElms))
                end if
            
            case ("Points")
                if ( nPoints.gt.0 .and. dPtr%nComps.eq.0 ) then
                    dPtr%nComps = getNumComps(vtk,dataArr,nPoints)
                    if ( .not.vtk%isBinApp ) dPtr%nComps = 1
                end if
                if ( dPtr%nComps.ne.3 ) then
                    write(stdout,ftab4) &
                    "WARNING: Element <NumberOfComponents> in <Points> &
                     attribute < 3"
                end if
                dPtr%nVals = nPoints
                dPtr%nElms = nPoints * dPtr%nComps
                if ( debug ) then
                    write(stdout,ftab4) &
                        "nPoints "// trim(STR(dPtr%nVals)) //&
                        "; nComps "// trim(STR(dPtr%nComps)) //&
                        "; nElems "// trim(STR(dPtr%nElms))
                end if
            
            case ("Verts")
                if ( nVerts.gt.0 .and. dPtr%nComps.eq.0 ) then
                    dPtr%nComps = getNumComps(vtk,dataArr,nVerts)
                    if ( .not.vtk%isBinApp ) dPtr%nComps = 1
                end if
                dPtr%nVals = nVerts
                dPtr%nElms = nVerts * dPtr%nComps
                if ( debug ) then
                    write(stdout,ftab4) &
                        "nVerts "// trim(STR(dPtr%nVals)) //&
                        "; nComps "// trim(STR(dPtr%nComps)) //&
                        "; nElems "// trim(STR(dPtr%nElms))
                end if
            
            case ("Lines")
                if ( nLines.gt.0 .and. dPtr%nComps.eq.0 ) then
                    dPtr%nComps = getNumComps(vtk,dataArr,nLines)
                    if ( .not.vtk%isBinApp ) dPtr%nComps = 1
                end if
                dPtr%nVals = nLines
                dPtr%nElms = nLines * dPtr%nComps
                if ( debug ) then
                    write(stdout,ftab4) &
                        "nLines "// trim(STR(dPtr%nVals)) //&
                        "; nComps "// trim(STR(dPtr%nComps)) //&
                        "; nElems "// trim(STR(dPtr%nElms))
                end if
            
            case ("Strips")
                if ( nStrips.gt.0 .and. dPtr%nComps.eq.0 ) then
                    dPtr%nComps = getNumComps(vtk,dataArr,nStrips)
                    if ( .not.vtk%isBinApp ) dPtr%nComps = 1
                end if
                dPtr%nVals = nStrips
                dPtr%nElms = nStrips * dPtr%nComps
                if ( debug ) then
                    write(stdout,ftab4) &
                        "nStrips "// trim(STR(dPtr%nVals)) //&
                        "; nComps "// trim(STR(dPtr%nComps)) //&
                        "; nElems "// trim(STR(dPtr%nElms))
                end if
            
            case ("Polys")
                if ( nPolys.gt.0 .and. dPtr%nComps.eq.0 ) then
                    dPtr%nComps = getNumComps(vtk,dataArr,nPolys)
                    if ( .not.vtk%isBinApp ) dPtr%nComps = 1
                end if
                dPtr%nVals = nPolys
                dPtr%nElms = nPolys * dPtr%nComps
                if ( debug ) then
                    write(stdout,ftab4) &
                        "nPolys "// trim(STR(dPtr%nVals)) //&
                        "; nComps "// trim(STR(dPtr%nComps)) //&
                        "; nElems "// trim(STR(dPtr%nElms))
                end if
            end select
            
            end subroutine vtkParsePolyData
            
            !==========================================
            
            function getNumComps(vtk,dArr,m) result(n)
            implicit none
            type(vtkXMLType), intent(in) :: vtk
            type(dataArrType), intent(inout) :: dArr
            integer(IK), intent(in) :: m
            integer(IK) :: n,npadd
            
            tokenList(:) = ""
            rLine = ""
            if ( vtk%isBinApp ) then
                n = dArr%nbytes
                npadd = 0
                if ( vtk%dataEncdng.eq."base64" ) then
                    n = n*3_IK/4_IK
                    if ( mod(n,3_IK).gt.0 ) &
                        npadd = 3_IK-mod(n,3_IK)
                end if
                n = n - npadd - dArr%hdrKind
                n = n / dArr%ikind / m
            else
                read(vtk%fid,'(A)') rLine
                call parse(trim(rLine),tokenList,ntoks)
                if ( ntoks.ne.0 ) n = ntoks
                call resetFilePos(vtk%fid,vtk%fileName,dArr%stPos-1,vtk%isBinApp)
            end if
            
            end function getNumComps

            !==========================================
            
            subroutine readBinaryData(vtk,dPtr,istat)
            implicit none
            type(vtkXMLType), intent(inout) :: vtk
            integer(IK), intent(inout) :: istat
            type(dataArrType), pointer :: dPtr
            
            character(len=:), allocatable :: code
            integer(IK) :: i,j,hdr,ikind
            integer(IK) :: np,np1,np2,npadd
            integer(IK1), allocatable, dimension(:) :: pIK1,p1,p2
            
            ikind = dPtr%ikind
            
            np  = dPtr%nBytes
            np1 = 1_IK*dPtr%hdrKind
            np2 = dPtr%nelms * ikind
            npadd = 0_IK
            
            allocate(pIK1(np),p1(np1),p2(np2))
            pIK1 = 0_IK1; p1 = 0_IK1; p2 = 0_IK1
            
            if ( vtk%dataEncdng.eq."base64" ) then
                code = repeat(" ",np)
                j = 0_IK
                do i=dPtr%stPos, dPtr%endPos-1
                    j = j+1_IK
                    read(vtk%fid,pos=i,end=001) code(j:j)
                end do

                np = np*3_IK/4_IK
                if ( mod(np1+np2,3_IK).ne.0 ) &
                    npadd = 3_IK - mod(np1+np2,3_IK)
            end if
            
            if ( np.ne.(np1+np2+npadd) ) then
                write(stdout,ftab4) &
                    "ERROR: inconsistent data array dimension.."
                write(stdout,ftab4) &
                    "Expected size: "//trim(STR(np1+np2+npadd))
                write(stdout,ftab4) &
                    "Actual size in file: "//trim(STR(np))
                istat=-1; return
            end if
            
            if ( vtk%dataEncdng.eq."base64" ) then
                call decode_bits(code,pIK1)
            else
                read(vtk%fid,pos=dPtr%stPos,end=001) pIK1(1:np)
            end if
            
            p1(1:np1) = pIK1(1:np1)
            p2(1:np2) = pIK1(np1+1:np1+np2)
            
            select case(dPtr%hdrkind)
            case(IK1)
                call transferBits(p1,np1,int(ikind,kind=IK1),hdr)
            case(IK2)
                call transferBits(p1,np1,int(ikind,kind=IK2),hdr)
            case(IK4)
                call transferBits(p1,np1,int(ikind,kind=IK4),hdr)
            case(IK8)
                call transferBits(p1,np1,int(ikind,kind=IK8),hdr)
            end select
            
            if ( dPtr%isInt ) then
                allocate(dPtr%iarr(dPtr%nElms))
                select case(ikind)
                case(IK1)
                    call transferBits(p2,np2,int(ikind,kind=IK1),dPtr%iarr,dPtr%nelms)
                case(IK2)
                    call transferBits(p2,np2,int(ikind,kind=IK2),dPtr%iarr,dPtr%nelms)
                case(IK4)
                    call transferBits(p2,np2,int(ikind,kind=IK4),dPtr%iarr,dPtr%nelms)
                case(IK8)
                    call transferBits(p2,np2,int(ikind,kind=IK8),dPtr%iarr,dPtr%nelms)
                    write(stdout,ftab4) &
                        "WARNING: Typecasting from INT64 to INT32 &
                        could lead to errors.."
                end select
            else
                allocate(dPtr%darr(dPtr%nelms))
                select case(ikind)
                case(RK4)
                    call transferBits(p2,np2,int(ikind,kind=IK4),dPtr%darr,dPtr%nelms)
                case(RK8)
                    call transferBits(p2,np2,int(ikind,kind=IK8),dPtr%darr,dPtr%nelms)
                end select
            end if
            
            if ( debug ) then
                if ( vtk%dataEncdng.eq."base64" ) &
                    write(stdout,ftab4) "Encoded data: "//trim(code)
                write(stdout,ftab4) "Decoded data: "
                write(stdout,'(13X,A)',advance="no")
                if ( dPtr%isInt ) then
                    do i=1, dPtr%nElms
                        write(stdout,'(A)',advance="no") " "//&
                            trim(STR(dPtr%iarr(i)))
                    end do
                else
                    do i=1, dPtr%nElms
                        write(stdout,'(A)',advance="no") " "//&
                            STR(dPtr%darr(i))
                    end do
                end if
                write(stdout,'(A)')
            end if
            
            return
            
 001        write(stdout,ftab4) "ERROR: end of file reached.."
            istat=-1; return
            
            end subroutine readBinaryData
            
            !==========================================
            
            subroutine readZlibBinaryData(vtk,dPtr,istat)
            implicit none
            type(vtkXMLType), intent(inout) :: vtk
            integer(IK), intent(inout) :: istat
            type(dataArrType), pointer :: dPtr
            
            character(len=:), allocatable :: code
            integer(IK) :: hdr(3),nlen,ikind
            integer(IK) :: i,j,ist,iend,iBlk
            integer(IK) :: np,np1,np2,npadd
            integer(IK), dimension(:), allocatable :: szBlk
            integer(IK1), allocatable, dimension(:) :: pIK1,p1,p2
            
            ikind = dPtr%ikind
            
            iPos = dPtr%stPos
            np1 = 3_IK * dPtr%hdrKind
            allocate(p1(np1)); p1=0_IK1
            npadd = 0_IK
            
            if ( vtk%dataEncdng.eq."base64" ) then
                nlen = 4_IK * dPTr%hdrKind
                code = repeat(" ",nlen)
                j = 0_IK
                do i=iPos, iPos+nlen-1_IK
                    j = j+1_IK
                    read(vtk%fid,pos=i,end=001) code(j:j)
                end do
                call decode_bits(code,p1)
                iPos = iPos + nlen
            else        
                read(vtk%fid,pos=iPos,end=001) p1(1:np1)
                iPos = iPos + np1
            end if

            select case ( dPtr%hdrKind )
            case (IK1)
                call transferBits(p1,np1,int(dPtr%hdrKind,kind=IK1),hdr,3_IK)
            case (IK2)
                call transferBits(p1,np1,int(dPtr%hdrKind,kind=IK2),hdr,3_IK)
            case (IK4)
                call transferBits(p1,np1,int(dPtr%hdrKind,kind=IK4),hdr,3_IK)
            case (IK8)
                call transferBits(p1,np1,int(dPtr%hdrKind,kind=IK8),hdr,3_IK)
            end select
            
            if ( debug ) then
                write(stdout,ftab4) &
                    "num blocks "//trim(STR(hdr(1)))
                write(stdout,ftab4) &
                    "size of block (before compr.) "// &
                    trim(STR(hdr(2)))
                write(stdout,ftab4) &
                    "size of last block (before compr.) "// &
                    trim(STR(hdr(3)))
            end if
            
            if ( hdr(1).lt.1_IK ) return

            np1 = (hdr(2)*(hdr(1)-1) + hdr(3))/ikind
            if ( np1.ne.dPtr%nElms ) then
                dPtr%nElms = np1
                dPtr%nComps = dPtr%nElms / dPtr%nVals
                if ( debug ) then
                    write(stdout,ftab4) "Reset : nComps "// &
                        trim(STR(dPtr%nComps))//" "//"nElms "//&
                        trim(STR(dPtr%nElms))
                end if
            end if
            deallocate(p1)
            
            np1 = hdr(1) * dPtr%hdrKind
            allocate(p1(np1)); p1(:) = 0_IK1
            allocate(szBlk(hdr(1))); szBlk(:) = 0_IK
            
            npadd = 0_IK
            if ( vtk%dataEncdng.eq."base64" ) then
                if ( mod(np1,3_IK).gt.0_IK ) &
                    npadd = 3_IK-mod(np1,3_IK)
                nlen = (np1+npadd)*4_IK/3_IK
                j = 0_IK
                code = repeat(" ",nlen)
                do i=iPos, iPos+nlen-1_IK
                    j = j+1_IK
                    read(vtk%fid,pos=i,end=001) code(j:j)
                end do
                call decode_bits(code,p1)
                iPos = iPos + nlen
            else                
                read(vtk%fid,pos=iPos,end=001) p1(1:np1)
                iPos = iPos + np1
            end if
            
            select case ( dPtr%hdrKind )
            case(IK1)
                call transferBits(p1,np1,int(dPtr%hdrKind,kind=IK1),szBlk,hdr(1))
            case(IK2)
                call transferBits(p1,np1,int(dPtr%hdrKind,kind=IK2),szBlk,hdr(1))
            case(IK4)
                call transferBits(p1,np1,int(dPtr%hdrKind,kind=IK4),szBlk,hdr(1))
            case(IK8)
                call transferBits(p1,np1,int(dPtr%hdrKind,kind=IK8),szBlk,hdr(1))
            end select
            
            if ( debug ) then
                do iBlk=1_IK, hdr(1)
                    write(stdout,ftab4) &
                        "size of block#"//trim(STR(iBlk))// &
                        " (after compr.) "//trim(STR(szBlk(iBlk)))
                end do
            end if
            
            if ( dPtr%isInt ) then
                allocate(dPtr%iarr(dPtr%nElms))
            else
                allocate(dPtr%darr(dPtr%nElms))
            end if
            
            if ( vtk%dataEncdng.eq."base64" ) then
                npadd = 0_IK
                np = sum(szBlk)
                if ( allocated(pIK1) ) deallocate(pIK1)
                allocate(pIK1(np)); pIK1 = 0_IK1
                if ( mod(np,3_IK).gt.0_IK ) &
                    npadd = 3_IK-mod(np,3_IK)
                nlen = (np + npadd)*4_IK/3_IK
                code = repeat(" ",nlen)
                j = 0_IK
                do i=iPos, iPos+nlen-1_IK
                    j = j+1_IK
                    read(vtk%fid,pos=i,end=001) code(j:j)
                end do
                call decode_bits(code,pIK1)
                iPos = iPos + nlen
            end if
            
            ist = 1_IK; np = 0_IK
            do iBlk=1_IK, hdr(1)
                np1 = szBlk(iBlk)
                np2 = hdr(2)
                if ( iBlk.eq.hdr(1) ) np2 = hdr(3)

                if ( allocated(p1) ) deallocate(p1)
                if ( allocated(p2) ) deallocate(p2)
                allocate(p1(np1)); p1(:) = 0_IK1
                allocate(p2(np2)); p2(:) = 0_IK1
                
                if ( vtk%dataEncdng.eq."base64" ) then
                    p1(1:np1) = pIK1(np+1:np+np1)
                    np = np+np1
                else
                    read(vtk%fid,pos=ipos) p1(1:np1)
                    iPos = iPos + np1
                end if
                
                call infZlibData(p1,np1,p2,np2,istat)
                if ( istat.lt.0 ) return
                
                nlen = np2/ikind
                iend = ist-1 + nlen

                if ( dPtr%isInt ) then
                    select case ( ikind )
                    case (IK1)
                        call transferBits(p2,np2,int(ikind,kind=IK1),dPtr%iarr(ist:iend),nlen)
                    case (IK2)
                        call transferBits(p2,np2,int(ikind,kind=IK2),dPtr%iarr(ist:iend),nlen)
                    case (IK4)
                        call transferBits(p2,np2,int(ikind,kind=IK4),dPtr%iarr(ist:iend),nlen)
                    case (IK8)
                        call transferBits(p2,np2,int(ikind,kind=IK8),dPtr%iarr(ist:iend),nlen)
                    end select
                    
                    if ( debug ) then
                        write(stdout,'(13X,A)',advance="no")
                        do i=ist,iend
                            write(stdout,'(A)',advance="no") &
                                trim(STR(dPtr%iarr(i)))//" "
                        end do
                        write(stdout,'(A)')
                    end if

                else
                    select case ( ikind )
                    case (RK4)
                        call transferBits(p2,np2,int(ikind,kind=IK4),dPtr%darr(ist:iend),nlen)
                    case (RK8)
                        call transferBits(p2,np2,int(ikind,kind=IK8),dPtr%darr(ist:iend),nlen)
                    end select
                    
                    if ( debug ) then
                        write(stdout,'(13X,A)',advance="no")
                        do i=ist,iend
                            write(stdout,'(A)',advance="no") &
                                STR(dPtr%darr(i))//" "
                        end do
                        write(stdout,'(A)')
                    end if

                    ist = iend+1
                end if
                ist = iend+1
            end do
            
            return
            
 001        write(stdout,ftab4) "ERROR: end of file reached.."
            istat=-1; return
            
            end subroutine readZlibBinaryData
            
            !==========================================
            
            subroutine decode_bits(code,bits)
            implicit none
            character(len=*), intent(in) :: code
            integer(IK1), intent(out) :: bits(:)

            integer(IK1) :: sixb(1:4)
            integer(IK) :: c,e,Nb,i
            
            Nb = size(bits,dim=1,kind=IK)
            e = 1_IK
            do c=1_IK, len(code), 4_IK
                sixb(:) = 0_IK1
                sixb(1) = index(b64List,code(c  :c  ),kind=IK1) - 1_IK1
                sixb(2) = index(b64List,code(c+1:c+1),kind=IK1) - 1_IK1
                sixb(3) = index(b64List,code(c+2:c+2),kind=IK1) - 1_IK1
                sixb(4) = index(b64List,code(c+3:c+3),kind=IK1) - 1_IK1
                
                call mvbits(sixb(1),0,6,bits(e),2)
                call mvbits(sixb(2),4,2,bits(e),0)
                
                if ( e+1.le.Nb ) then
                    call mvbits(sixb(2),0,4,bits(e+1),4)
                    call mvbits(sixb(3),2,4,bits(e+1),0)
                end if
                
                if ( e+2.le.Nb ) then
                    call mvbits(sixb(3),0,2,bits(e+2),6)
                    call mvbits(sixb(4),0,6,bits(e+2),0)
                end if
                e = e+3_IK
            end do
            
            end subroutine decode_bits

            !==========================================
            
            subroutine readAsciiData(vtk,dPtr,istat)
            implicit none
            type(vtkXMLType), intent(inout) :: vtk
            integer(IK), intent(inout) :: istat
            type(dataArrType), pointer :: dPtr
            integer(IK) :: i
            
            if ( dPtr%isInt ) then
                allocate(dPtr%iarr(dPtr%nElms))
                read(vtk%fid,*,end=001) (dPtr%iarr(i),i=1, dPtr%nElms)
                if ( debug ) then
                    write(stdout,'(13X,A)',advance="no")
                    do i=1, dPtr%nElms
                        write(stdout,'(A)',advance="no") " "//&
                            trim(STR(dPtr%iarr(i)))
                    end do
                    write(stdout,'(A)')
                end if
            else
                allocate(dPtr%darr(dPtr%nElms))
                read(vtk%fid,*,end=001) (dPtr%darr(i),i=1, dPtr%nElms)
                if ( debug ) then
                    write(stdout,'(13X,A)',advance="no")
                    do i=1, dPtr%nElms
                        write(stdout,'(A)',advance="no") " "//&
                            STR(dPtr%darr(i))
                    end do
                    write(stdout,'(A)')
                end if
            end if
            
            return
            
 001        write(stdout,ftab4) "ERROR: end of file reached.."
            istat=-1; return
            
            end subroutine readAsciiData

            !==========================================

            subroutine findKwrdXML(vtk,sKwrd,iPos,strng,istat,ePos)
            implicit none
            type(vtkXMLType), intent(in) :: vtk
            integer(IK), intent(inout) :: iPos
            character(len=*), intent(in) :: sKwrd
            character(len=stdl), intent(out) :: strng
            integer(IK), intent(in), optional :: ePos
            integer(IK), intent(inout) :: istat
            
            integer(IK) :: i,kwrdL,cnt
            character :: c
            
            kwrdL = len(trim(sKwrd))
            do
                strng = ''
                if ( vtk%isBinApp ) then
                    do i=1, stdl
                        read(vtk%fid,pos=iPos,end=001) c
                        iPos = iPos + 1
                        if ( c.eq. '<' ) then
                            cnt = 1
                            strng(cnt:cnt) = c
                            do
                                read(vtk%fid,pos=iPos,end=001) c
                                iPos = iPos + 1
                                if ( present(ePos) .and. iPos.ge.ePos) &
                                then
                                    strng = ''
                                    return
                                end if
                                cnt = cnt+1
                                strng(cnt:cnt) = c
                                if ( cnt.eq.kwrdL ) then
                                    if (strng(1:cnt).eq.trim(skwrd)) then
                                        cycle
                                    else
                                        strng = ''
                                        exit
                                    end if
                                end if
                                if ( c.eq.'>' ) exit
                            end do
                            strng = adjustl(strng)
                            if ( strng(1:kwrdL).eq.trim(sKwrd) ) return
                        end if
                        if(c.eq.eol) exit
                    end do
                else
                    strng = ''
                    read(vtk%fid,'(A)',end=001) strng
                    iPos = iPos+1
                    if ( present(ePos) .and. iPos.ge.ePos ) then
                        strng = ''
                        return
                    end if
                end if
                strng = adjustl(strng)
                if ( strng(1:kwrdL).eq.trim(sKwrd) ) return
            end do
            
            return
            
 001        write(stdout,ftab4) "ERROR: end of file reached.."
            istat = -1; return

            end subroutine findKwrdXML
        
            !==========================================

            subroutine adjustDataArray(fid,spos,epos)
            implicit none
            integer(IK), intent(in) :: fid
            integer(IK), intent(inout) :: spos,epos
            integer(IK) :: npos
            character :: c
            
            write(stdout,'(A)') trim(STR(spos))//" "//trim(STR(epos))
            npos = spos
            do
                read(fid,pos=npos) c
                if ( c.eq.' ' .or. c.eq.eol .or. c.eq.'    ') then
                    npos = npos+1
                    cycle
                else
                    exit
                end if
            end do
            spos = npos
            
            do
                read(fid,pos=npos) c
                write(stdout,'(A)') trim(STR(npos))//" '"// &
                c//"' "//trim(STR(ichar(c)))
                if ( c.eq.' ' .or. c.eq.eol .or. c.eq.'    ' ) exit
                npos = npos+1
            end do
            epos = npos-1
            
            end subroutine adjustDataArray
                
            !==========================================
            
            subroutine resetFilePos(fileId,fName,fPos,isBin)
            implicit none
            integer, intent(in) :: fileId,fPos
            logical, intent(in) :: isBin
            character(len=*), intent(in) :: fName
            integer(IK), parameter :: SEEK_SET=0
            integer(IK) :: i,ierr
            
            if ( isBin ) then
                call fseek(fileId,fPos,SEEK_SET,ierr)
            else
                close(fileId)
                open(fileId,file=trim(fName),status='old')
                do i=1, fPos
                    read(fileId,*)
                end do
            end if
            
            end subroutine resetFilePos
            
            !==========================================
            
            subroutine selectDataType(varType,ikind,isint,istat)
            implicit none
            character(len=*), intent(in) :: varType
            integer(IK), intent(inout) :: istat
            integer(IK), intent(out) :: ikind
            logical, intent(out) :: isint
            
            isint = .true.
            ikind = 0
            select case (trim(varType))
            case ("Int8", "UInt8")
                ikind = IK1
            case ("Int16", "UInt16")
                ikind = IK2
            case ("Int32", "UInt32")
                ikind = IK4
            case ("Int64", "UInt64")
                ikind = IK8
            case ("Float32")
                ikind = RK4
                isint = .false.
            case ("Float64")
                ikind = RK8
                isint = .false.
            case default
                write(stdout,ftab4) "ERROR: unknown data type.."
                istat=-1; return
            end select
            
            end subroutine selectDataType
            
            !==========================================

            subroutine getVTK_numPoints(vtk,nn,istat)
            implicit none
            type(vtkXMLType), intent(in) :: vtk
            integer(IK), intent(inout) :: istat
            integer(IK), intent(out) :: nn
            integer(IK) :: i
            
            nn = vtk%pieceElms(1)
            if ( nn.eq.0 ) then
                write(stdout,ftab4) &
                    "ERROR: could not find POINTS attribute"
                istat=-1; return
            end if
            
            end subroutine getVTK_numPoints

            !==========================================

            subroutine getVTK_numElems(vtk,ne,istat)
            implicit none
            type(vtkXMLType), intent(in) :: vtk
            integer(IK), intent(inout) :: istat
            integer(IK), intent(out) :: ne
            integer(IK) :: i
            
            if ( vtk%pieceElms(2).gt.0 ) then
                ne = vtk%pieceElms(2)
            else if ( vtk%pieceElms(6).gt.0 ) then
                ne = vtk%pieceElms(6)
            end if
            
            if ( ne.eq.0 ) then
                write(stdout,ftab4) &
                    "ERROR: could not find CELLS or POLYS attributes"
                istat=-1; return
            end if
            
            end subroutine getVTK_numElems

            !==========================================

            subroutine getVTK_nodesPerElem(vtk,eNoN,istat)
            implicit none
            type(vtkXMLType), intent(in) :: vtk
            integer(IK), intent(inout) :: istat
            integer(IK), intent(out) :: eNoN
            integer(IK) :: iatt,i,itmp
            
            if ( vtk%pieceElms(2).gt.0 ) then    ! nCells !
                iatt = 9
            else if (vtk%pieceElms(6).gt.0 ) then ! nPolys !
                iatt = 8
            else
                write(stdout,ftab4) &
                    "ERROR: could not find CELLS or POLYS attributes"
                istat=-1; return
            end if
            
            eNoN = 0; itmp = 0
            do i=1, vtk%pcAtt(iatt)%n
                if ( trim(vtk%pcAtt(iatt)%dataArr(i)%dName) .eq. &
                    "connectivity" ) then
                    eNoN = vtk%pcAtt(iatt)%dataArr(i)%nComps
                end if
                if ( trim(vtk%pcAtt(iatt)%dataArr(i)%dName) .eq. &
                    "offsets" ) then
                    itmp = vtk%pcAtt(iatt)%dataArr(i)%iarr(2) - &
                           vtk%pcAtt(iatt)%dataArr(i)%iarr(1)
                end if
            end do
            
            eNoN = max(eNoN, itmp)
            if ( eNoN.eq.0 ) then
                write(stdout,ftab4) &
                    "ERROR: unexpected VTK behavior (eNoN)"
                istat=-1; return
            end if
            
            end subroutine getVTK_nodesPerElem

            !==========================================

            subroutine getVTK_elemIEN(vtk,ien,istat)
            implicit none
            type(vtkXMLType), intent(inout) :: vtk
            integer(IK), intent(out) :: ien(:,:)
            integer(IK), intent(inout) :: istat
            integer(IK) :: iatt,ne,eNoN,i,j,k,l
            
            eNoN = size(ien,1)
            ne   = size(ien,2)
            if ( vtk%pieceElms(2).gt.0 ) then    ! nCells !
                iatt = 9
            else if (vtk%pieceElms(6).gt.0 ) then ! nPolys !
                iatt = 8
            else
                write(stdout,ftab4) &
                    "ERROR: could not find CELLS or POLYS attributes"
                istat=-1; return
            end if
            
            do i=1, vtk%pcAtt(iatt)%n
                if ( trim(vtk%pcAtt(iatt)%dataArr(i)%dName) .eq. &
                    "connectivity" ) then
                    if ( vtk%pcAtt(iatt)%dataArr(i)%nElms .ne. &
                         (eNoN*ne) ) then
                        write(stdout,ftab4) &
                            "ERROR: mismatch in IEN params.."
                        write(stdout,ftab4) &
                            trim(STR(vtk%pcAtt(iatt)%dataArr(i)%nElms))
                        write(stdout,ftab4) &
                            trim(STR(eNoN))//" "//trim(STR(ne))
                        istat=-1; return
                    end if
                    if ( .not.allocated(vtk%pcAtt(iatt)%dataArr(i)%iarr) ) then
                        write(stdout,ftab4) &
                            "ERROR: connectivity array unallocated.."
                        istat=-1; return
                    end if
                    
                    l = 0
                    do j=1, ne
                        do k=1, eNoN
                            l = l+1
                            ien(k,j) = &
                                vtk%pcAtt(iatt)%dataArr(i)%iarr(l)
                        end do
                    end do
                    exit
                end if
            end do

            if ( i.gt.vtk%pcAtt(iatt)%n ) then
                write(stdout,ftab4) &
                    "ERROR: could not find connectivity in Cells or &
                    Polys attributes"
                istat=-1; return
            end if

            end subroutine getVTK_elemIEN
            
            !==========================================

            subroutine getVTK_pointCoords(vtk,x,istat)
            implicit none
            type(vtkXMLType), intent(inout) :: vtk
            integer(IK), intent(inout) :: istat
            real(RK), intent(out) :: x(:,:)
            integer(IK) :: iatt,i,j,k,l,nd,nn
            
            nd = size(x,1)
            nn = size(x,2)
            if ( vtk%pieceElms(1).eq.nn ) then    ! nPoints !
                iatt = 3
            else
                write(stdout,ftab4) &
                    "ERROR: could not find POINTS attribute"
                istat=-1; return
            end if
            
            i=1
            if ( trim(vtk%pcAtt(iatt)%dataArr(i)%dName) .ne. &
                "Points" ) then
                write(stdout,ftab4) &
                    "WARNING: <DataArray> name does not match Points.."
            end if
            
            if ( vtk%pcAtt(iatt)%dataArr(i)%nElms .ne. (nd*nn) ) then
                write(stdout,ftab4) &
                    "ERROR: mismatch in POINTS params.."
                write(stdout,ftab4) &
                    trim(STR(vtk%pcAtt(iatt)%dataArr(i)%nElms))
                write(stdout,ftab4) &
                    trim(STR(nd))//" "//trim(STR(nn))
                istat=-1; return
            end if
            
            if ( .not.allocated(vtk%pcAtt(iatt)%dataArr(i)%darr) ) then
                write(stdout,ftab4) &
                    "ERROR: Points array unallocated.."
                istat=-1; return
            end if
            
            l = 0
            do j=1, nn
                do k=1, nd
                    l = l+1
                    x(k,j) = vtk%pcAtt(iatt)%dataArr(i)%darr(l)
                end do
            end do
            
            end subroutine getVTK_pointCoords
            
            !==========================================

            subroutine getVTK_pointDataIntS(vtk,kwrd,u,istat)
            implicit none
            type(vtkXMLType), intent(inout) :: vtk
            character(len=*), intent(in) :: kwrd
            integer(IK), intent(inout) :: u(:)
            integer(IK), intent(inout) :: istat
            integer(IK) :: iatt,i,n
            logical :: flag
            
            n = size(u)
            if ( vtk%pieceElms(1).gt.0 ) then    ! nPoints !
                iatt = 1
            else
                write(stdout,ftab4) &
                    "ERROR: could not find POINTS attribute to read &
                    PointData"
                istat=-1; return
            end if
            
            flag = .false.
            do i=1, vtk%pcAtt(iatt)%n
                if ( trim(vtk%pcAtt(iatt)%dataArr(i)%dName) .eq. &
                    trim(kwrd) ) then
                    flag = .true.
                    if ( vtk%pcAtt(iatt)%dataArr(i)%nElms .ne. n ) then
                        write(stdout,ftab4) &
                            "ERROR: mismatch in array size and numNodes"
                        write(stdout,ftab4) "Actual size: "//&
                            trim(STR(vtk%pcAtt(iatt)%dataArr(i)%nElms))
                        write(stdout,ftab4) "Input size: "// &
                            trim(STR(n))
                        istat=-1; return
                    end if
                    if ( .not.allocated(vtk%pcAtt(iatt)%dataArr(i)%iarr) ) then
                        write(stdout,ftab4) &
                            "ERROR: data not found.."
                        istat=-1; return
                    end if
                    
                    u(:) = vtk%pcAtt(iatt)%dataArr(i)%iarr(:)
                    exit
                end if
            end do
            
            if ( .not.flag ) then
                write(stdout,ftab4) &
                    "ERROR: could not find <"//trim(kwrd)//"> &
                    in PointData attribute"
                istat=-1; return
            end if
            
            end subroutine getVTK_pointDataIntS
            
            !==========================================

            subroutine getVTK_pointDataIntV(vtk,kwrd,u,istat)
            implicit none
            type(vtkXMLType), intent(inout) :: vtk
            character(len=*), intent(in) :: kwrd
            integer(IK), intent(inout) :: u(:,:)
            integer(IK), intent(inout) :: istat
            integer(IK) :: iatt,i,j,k,m,n
            logical :: flag
            
            m = size(u,1)
            n = size(u,2)
            if ( vtk%pieceElms(1).gt.0 ) then    ! nPoints !
                iatt = 1
            else
                write(stdout,ftab4) &
                    "ERROR: could not find POINTS attribute to read &
                    PointData"
                istat=-1; return
            end if
            
            flag = .false.
            do i=1, vtk%pcAtt(iatt)%n
                if ( trim(vtk%pcAtt(iatt)%dataArr(i)%dName) .eq. &
                    trim(kwrd) ) then
                    flag = .true.
                    if ( (vtk%pcAtt(iatt)%dataArr(i)%nVals  .ne. n) .or. &
                         (vtk%pcAtt(iatt)%dataArr(i)%nComps .ne. m) ) then
                        write(stdout,ftab4) &
                            "ERROR: mismatch in array sizes"
                        write(stdout,ftab4) "Actual size: "//&
                            trim(STR(vtk%pcAtt(iatt)%dataArr(i)%nComps))//", "// &
                            trim(STR(vtk%pcAtt(iatt)%dataArr(i)%nVals))
                        write(stdout,ftab4) "Input size: "// &
                            trim(STR(m))//", "//trim(STR(n))
                        istat=-1; return
                    end if
                    if ( .not.allocated(vtk%pcAtt(iatt)%dataArr(i)%iarr) ) then
                        write(stdout,ftab4) &
                            "ERROR: data not found.."
                        istat=-1; return
                    end if
                    
                    do j=1, n
                        do k=1, m
                            u(k,j) = &
                              vtk%pcAtt(iatt)%dataArr(i)%iarr((j-1)*m+k)
                        end do
                    end do
                    exit
                end if
            end do
            
            if ( .not.flag ) then
                write(stdout,ftab4) &
                    "ERROR: could not find <"//trim(kwrd)//"> &
                    in PointData attribute"
                istat=-1; return
            end if
            
            end subroutine getVTK_pointDataIntV
            
            !==========================================

            subroutine getVTK_pointDataRealS(vtk,kwrd,u,istat)
            implicit none
            type(vtkXMLType), intent(inout) :: vtk
            character(len=*), intent(in) :: kwrd
            REAL(RK), intent(inout) :: u(:)
            integer(IK), intent(inout) :: istat
            integer(IK) :: iatt,i,n
            logical :: flag
            
            n = size(u)
            if ( vtk%pieceElms(1).gt.0 ) then    ! nPoints !
                iatt = 1
            else
                write(stdout,ftab4) &
                    "ERROR: could not find POINTS attribute to read &
                    PointData"
                istat=-1; return
            end if
            
            flag = .false.
            do i=1, vtk%pcAtt(iatt)%n
                if ( trim(vtk%pcAtt(iatt)%dataArr(i)%dName) .eq. &
                    trim(kwrd) ) then
                    flag = .true.
                    if ( vtk%pcAtt(iatt)%dataArr(i)%nElms .ne. n ) then
                        write(stdout,ftab4) &
                            "ERROR: mismatch in array size and numNodes"
                        write(stdout,ftab4) &
                            trim(STR(vtk%pcAtt(iatt)%dataArr(i)%nElms))
                        write(stdout,ftab4) trim(STR(n))
                        istat=-1; return
                    end if
                    if ( .not.allocated(vtk%pcAtt(iatt)%dataArr(i)%darr) ) then
                        write(stdout,ftab4) &
                            "ERROR: data not found.."
                        istat=-1; return
                    end if
                    
                    u(:) = vtk%pcAtt(iatt)%dataArr(i)%darr(:)
                    exit
                end if
            end do
            
            if ( .not.flag ) then
                write(stdout,ftab4) &
                    "ERROR: could not find <"//trim(kwrd)//"> &
                    in PointData attribute"
                istat=-1; return
            end if
            
            end subroutine getVTK_pointDataRealS
            
            !==========================================

            subroutine getVTK_pointDataRealV(vtk,kwrd,u,istat)
            implicit none
            type(vtkXMLType), intent(inout) :: vtk
            character(len=*), intent(in) :: kwrd
            real(RK), intent(inout) :: u(:,:)
            integer(IK), intent(inout) :: istat
            integer(IK) :: iatt,i,j,k,m,n
            logical :: flag
            
            m = size(u,1)
            n = size(u,2)
            if ( vtk%pieceElms(1).gt.0 ) then    ! nPoints !
                iatt = 1
            else
                write(stdout,ftab4) &
                    "ERROR: could not find POINTS attribute to read &
                    PointData"
                istat=-1; return
            end if
            
            flag = .false.
            do i=1, vtk%pcAtt(iatt)%n
                if ( trim(vtk%pcAtt(iatt)%dataArr(i)%dName) .eq. &
                    trim(kwrd) ) then
                    flag = .true.
                    if ( (vtk%pcAtt(iatt)%dataArr(i)%nVals  .ne. n) .or. &
                         (vtk%pcAtt(iatt)%dataArr(i)%nComps .ne. m) ) then
                        write(stdout,ftab4) &
                            "ERROR: mismatch in array sizes"
                        write(stdout,ftab4) "Actual size: "//&
                            trim(STR(vtk%pcAtt(iatt)%dataArr(i)%nComps))//", "// &
                            trim(STR(vtk%pcAtt(iatt)%dataArr(i)%nVals))
                        write(stdout,ftab4) "Input size: "// &
                            trim(STR(m))//", "//trim(STR(n))
                        istat=-1; return
                    end if
                    if ( .not.allocated(vtk%pcAtt(iatt)%dataArr(i)%darr) ) then
                        write(stdout,ftab4) &
                            "ERROR: data not found.."
                        istat=-1; return
                    end if
                    
                    do j=1, n
                        do k=1, m
                            u(k,j) = &
                              vtk%pcAtt(iatt)%dataArr(i)%darr((j-1)*m+k)
                        end do
                    end do
                    exit
                end if
            end do
            
            if ( .not.flag ) then
                write(stdout,ftab4) &
                    "ERROR: could not find <"//trim(kwrd)//"> &
                    in PointData attribute"
                istat=-1; return
            end if
            
            end subroutine getVTK_pointDataRealV
            
            !==========================================

            subroutine getVTK_elemDataIntS(vtk,kwrd,u,istat)
            implicit none
            type(vtkXMLType), intent(inout) :: vtk
            character(len=*), intent(in) :: kwrd
            integer(IK), intent(inout) :: u(:)
            integer(IK), intent(inout) :: istat
            integer(IK) :: iatt,i,n
            logical :: flag
            
            n = size(u)
            if ( vtk%pieceElms(2).gt.0 .or. vtk%pieceElms(6).gt.0 ) then
                iatt = 2
            else
                write(stdout,ftab4) &
                    "ERROR: could not find either CELLS or POLYS &
                    attributes to read CellData"
                istat=-1; return
            end if
            
            do i=1, vtk%pcAtt(iatt)%n
                if ( trim(vtk%pcAtt(iatt)%dataArr(i)%dName) .eq. &
                     trim(kwrd) ) then
                    flag = .true.
                    if ( vtk%pcAtt(iatt)%dataArr(i)%nElms .ne. n ) then
                        write(stdout,ftab4) &
                            "ERROR: mismatch in array size and numElems"
                        write(stdout,ftab4) "Actual size: "//&
                            trim(STR(vtk%pcAtt(iatt)%dataArr(i)%nElms))
                        write(stdout,ftab4) "Input size: "//&
                            trim(STR(n))
                        istat=-1; return
                    end if
                    if ( .not.allocated(vtk%pcAtt(iatt)%dataArr(i)%iarr) ) then
                        write(stdout,ftab4) &
                            "ERROR: data not found.."
                        istat=-1; return
                    end if
                    
                    u(:) = vtk%pcAtt(iatt)%dataArr(i)%iarr(:)
                    exit
                end if
            end do

            if ( .not.flag ) then
                write(stdout,ftab4) &
                    "ERROR: could not find <"//trim(kwrd)//"> &
                    in CellData attribute"
                istat=-1; return
            end if
            
            end subroutine getVTK_elemDataIntS
            
             !==========================================

            subroutine getVTK_elemDataIntV(vtk,kwrd,u,istat)
            implicit none
            type(vtkXMLType), intent(inout) :: vtk
            character(len=*), intent(in) :: kwrd
            integer(IK), intent(inout) :: u(:,:)
            integer(IK), intent(inout) :: istat
            integer(IK) :: iatt,i,j,k,m,n
            logical :: flag
            
            m = size(u,1)
            n = size(u,2)
            if ( vtk%pieceElms(2).gt.0 .or. vtk%pieceElms(6).gt.0 ) then
                iatt = 2
            else
                write(stdout,ftab4) &
                    "ERROR: could not find either CELLS or POLYS &
                    attributes to read CellData"
                istat=-1; return
            end if
            
            do i=1, vtk%pcAtt(iatt)%n
                if ( trim(vtk%pcAtt(iatt)%dataArr(i)%dName) .eq. &
                     trim(kwrd) ) then
                    flag = .true.
                    if ( (vtk%pcAtt(iatt)%dataArr(i)%nVals  .ne. n) .or. &
                         (vtk%pcAtt(iatt)%dataArr(i)%nComps .ne. m) ) then
                        write(stdout,ftab4) &
                            "ERROR: mismatch in array sizes"
                        write(stdout,ftab4) "Actual size: "//&
                            trim(STR(vtk%pcAtt(iatt)%dataArr(i)%nComps))//", "// &
                            trim(STR(vtk%pcAtt(iatt)%dataArr(i)%nVals))
                        write(stdout,ftab4) "Input size: "// &
                            trim(STR(m))//", "//trim(STR(n))
                        istat=-1; return
                    end if
                    if ( .not.allocated(vtk%pcAtt(iatt)%dataArr(i)%iarr) ) then
                        write(stdout,ftab4) &
                            "ERROR: data not found.."
                        istat=-1; return
                    end if
                    
                    do j=1, n
                        do k=1, m
                            u(k,j) = &
                              vtk%pcAtt(iatt)%dataArr(i)%iarr((j-1)*m+k)
                        end do
                    end do
                    exit
                end if
            end do

            if ( .not.flag ) then
                write(stdout,ftab4) &
                    "ERROR: could not find <"//trim(kwrd)//"> &
                    in CellData attribute"
                istat=-1; return
            end if
            
            end subroutine getVTK_elemDataIntV
            
            !==========================================

            subroutine getVTK_elemDataRealS(vtk,kwrd,u,istat)
            implicit none
            type(vtkXMLType), intent(inout) :: vtk
            character(len=*), intent(in) :: kwrd
            real(RK), intent(inout) :: u(:)
            integer(IK), intent(inout) :: istat
            integer(IK) :: iatt,i,n
            logical :: flag
            
            n = size(u)
            if ( vtk%pieceElms(2).gt.0 .or. vtk%pieceElms(6).gt.0 ) then
                iatt = 2
            else
                write(stdout,ftab4) &
                    "ERROR: could not find either CELLS or POLYS &
                    attributes to read CellData"
                istat=-1; return
            end if
            
            do i=1, vtk%pcAtt(iatt)%n
                if ( trim(vtk%pcAtt(iatt)%dataArr(i)%dName) .eq. &
                     trim(kwrd) ) then
                    flag = .true.
                    if ( vtk%pcAtt(iatt)%dataArr(i)%nElms .ne. n ) then
                        write(stdout,ftab4) &
                            "ERROR: mismatch in array size and numElems"
                        write(stdout,ftab4) "Actual size: "//&
                            trim(STR(vtk%pcAtt(iatt)%dataArr(i)%nElms))
                        write(stdout,ftab4) "Input size: "//&
                            trim(STR(n))
                        istat=-1; return
                    end if
                    if ( .not.allocated(vtk%pcAtt(iatt)%dataArr(i)%darr) ) then
                        write(stdout,ftab4) &
                            "ERROR: data not found.."
                        istat=-1; return
                    end if
                    
                    u(:) = vtk%pcAtt(iatt)%dataArr(i)%darr(:)
                    exit
                end if
            end do

            if ( .not.flag ) then
                write(stdout,ftab4) &
                    "ERROR: could not find <"//trim(kwrd)//"> &
                    in CellData attribute"
                istat=-1; return
            end if
            
            end subroutine getVTK_elemDataRealS
            
             !==========================================

            subroutine getVTK_elemDataRealV(vtk,kwrd,u,istat)
            implicit none
            type(vtkXMLType), intent(inout) :: vtk
            character(len=*), intent(in) :: kwrd
            real(RK), intent(inout) :: u(:,:)
            integer(IK), intent(inout) :: istat
            integer(IK) :: iatt,i,j,k,m,n
            logical :: flag
            
            m = size(u,1)
            n = size(u,2)
            if ( vtk%pieceElms(2).gt.0 .or. vtk%pieceElms(6).gt.0 ) then
                iatt = 2
            else
                write(stdout,ftab4) &
                    "ERROR: could not find either CELLS or POLYS &
                    attributes to read CellData"
                istat=-1; return
            end if
            
            do i=1, vtk%pcAtt(iatt)%n
                if ( trim(vtk%pcAtt(iatt)%dataArr(i)%dName) .eq. &
                     trim(kwrd) ) then
                    flag = .true.
                    if ( (vtk%pcAtt(iatt)%dataArr(i)%nVals  .ne. n) .or. &
                         (vtk%pcAtt(iatt)%dataArr(i)%nComps .ne. m) ) then
                        write(stdout,ftab4) &
                            "ERROR: mismatch in array sizes"
                        write(stdout,ftab4) "Actual size: "//&
                            trim(STR(vtk%pcAtt(iatt)%dataArr(i)%nComps))//", "// &
                            trim(STR(vtk%pcAtt(iatt)%dataArr(i)%nVals))
                        write(stdout,ftab4) "Input size: "// &
                            trim(STR(m))//", "//trim(STR(n))
                        istat=-1; return
                    end if
                    if ( .not.allocated(vtk%pcAtt(iatt)%dataArr(i)%darr) ) then
                        write(stdout,ftab4) &
                            "ERROR: data not found.."
                        istat=-1; return
                    end if
                    
                    do j=1, n
                        do k=1, m
                            u(k,j) = &
                              vtk%pcAtt(iatt)%dataArr(i)%darr((j-1)*m+k)
                        end do
                    end do
                    exit
                end if
            end do

            if ( .not.flag ) then
                write(stdout,ftab4) &
                    "ERROR: could not find <"//trim(kwrd)//"> &
                    in CellData attribute"
                istat=-1; return
            end if
            
            end subroutine getVTK_elemDataRealV
            
            !==========================================

        end module vtkXMLMod
        
!**************************************************







!**************************************************

      module variables
      use commod
      implicit none

      character(len=strL) :: fout
      type(faceType) :: sFa, pFa, tFa

      end module variables

!**************************************************

      program main_partFaces
      use variables
      implicit none

      integer i, e, a, Ac, Ec, fid
      character(len=strL) :: fname

      integer, allocatable :: incE(:), incN(:), ptr(:)

      fid = 1889
      i = iargc()
      if (i .ne. 1) then
         write(stdout,ftab4) "Error: at least one input argument needed"
         STOP
      end if
      call getarg(1,fname)

      write(stdout,ftab1) repeat('=', 48)
      write(stdout,ftab1,advance='no') &
     &   "Number of spatial dimensions (nsd): "
      read(*,*) nsd

      if (nsd.ne.2 .and. nsd.ne.3) then
         write(stdout,ftab4) "Error: incorrect spatial dimensions"
         STOP
      end if

      if (nsd .eq. 2) then
         nstd = 3
      else
         nstd = 6
      end if

      call readInputs(fname)

      write(*,ftab1) "Extracting faces..."
      i = MAXVAL(sFa%gE(:))
      allocate(incE(i))
      incE = 0

      do e=1, sFa%nEl
         Ec = sFa%gE(e)
         incE(Ec) = 1
      end do

      do e=1, pFa%nEl
         Ec = pFa%gE(e)
         incE(Ec) = 0
      end do

      tFa%vtkType = sFa%vtkType
      tFa%eNoN    = sFa%eNoN
      tFa%nEl     = sum(incE(:))
      allocate(tFa%gE(tFa%nEl), tFa%IEN(tFa%eNoN,tFa%nEl), &
     &   incN(sFa%nNo))

      i    = 0
      incN = 0
      do e=1, sFa%nEl
         Ec = sFa%gE(e)
         if (incE(Ec) .eq. 0) cycle
         i = i + 1
         tFa%gE(i) = Ec
         tFa%IEN(:,i) = sFa%IEN(:,e)
         do a=1, sFa%eNoN
            Ac = sFa%IEN(a,e)
            incN(Ac) = 1
         end do
      end do

      tFa%nNo = sum(incN)
      allocate(tFa%x(nsd,tFa%nNo), tFa%gN(tFa%nNo), ptr(sFa%nNo))
      i   = 0
      ptr = 0
      do a=1, sFa%nNo
         if (incN(a) .eq. 0) cycle
         i = i + 1
         ptr(a) = i
         tFa%gN(i)  = sFa%gN(a)
         tFa%x(:,i) = sFa%x(:,a)
      end do

      do e=1, tFa%nEl
         do a=1, tFa%eNoN
            Ac = tFa%IEN(a,e)
            Ac = ptr(Ac)
            if (Ac .eq. 0) then
               write(*,ftab4) "ERROR: unexpected behavior"
               STOP
            end if
            tFa%IEN(a,e) = Ac - 1
         end do
      end do

      write(*,ftab2) "Num nodes in new face: "//trim(STR(tFa%nNo))
      write(*,ftab2) "Num polys in new face: "//trim(STR(tFa%nEl))

      write(*,ftab1) "Writing separated face ---> "//trim(fout)
      call writeVTP(tFa, fout)
      write(stdout,ftab1) repeat('=', 48)

      deallocate(incE, incN)
      call destroy(sFa)
      call destroy(pFa)
      call destroy(tFa)

      return
      end program main_partFaces

!**************************************************

      subroutine readInputs(fin)
      use variables
      implicit none
      character(len=strL), intent(in) :: fin

      integer fid, istat
      character(len=strL) :: fname

      fid = 1265
      istat = 0
      open(fid,file=trim(fin))
      OUTER_LOOP: do
         call findKwrd(fid, "sourceFace", istat)
         if (istat .ne. 0) exit
         read(fid,*) fname
         write(*,ftab1) "Reading source face <--- "//trim(fname)
         call readVTP(sFa, fname)

         call findKwrd(fid, "separatedFace", istat)
         if (istat .ne. 0) exit
         read(fid,*) fname
         write(*,ftab1) "Reading separated face <--- "//trim(fname)
         call readVTP(pFa, fname)

         if (sFa%vtkType .ne. pFa%vtkType) then
            write(stdout,ftab4) "ERROR: inconsistent face types"
            istat = -1
            exit
         end if

         call findKwrd(fid, "outputFace", istat)
         if (istat .ne. 0) exit
         read(fid,*) fout

         exit OUTER_LOOP
      end do OUTER_LOOP
      close(fid)

      if (istat .ne. 0) then
         write(stdout,ftab4) "ERROR: reading inputs"
         stop
      end if

      return
      contains
         !==========================================
         subroutine findKwrd(fileId, sKwrd, istat)
         implicit none
         integer, intent(in) :: fileId
         integer, intent(inout) :: istat
         character(len=*), intent(in) :: sKwrd

         integer :: kwrdL, slen
         character(len=strL) :: sLine

         istat = 0
         kwrdL = len(trim(sKwrd))
         do
            read(fileId,'(A)',end=001) sLine
            if (sLine(1:1) .eq. '#') then
               slen  = len(trim(sLine))
               sLine = sLine(2:slen)
               sLine = adjustl(sLine)
               if (sLine(1:kwrdL) .eq. trim(sKwrd)) return
            end if
         end do

 001     write(stdout,ftab4) "ERROR: EOF reached while finding "// &
         "keyword <"//trim(sKwrd)//">"
         istat = -1
         return

         end subroutine findKwrd
         !==========================================
      end subroutine readInputs

!**************************************************

      subroutine readVTP(lFa, fname)
      use commod
      implicit none
      character(len=strL), intent(in) :: fname
      type(faceType), intent(inout) :: lFa

      type(vtkXMLType) :: vtp
      integer :: istat
      real(kind=8), allocatable :: tmpX1(:), tmpX2(:,:)

      istat = 0
      do while (istat .eq. 0)

         call loadVTK(vtp, trim(fname), istat)
         if (istat .lt. 0) exit

         call getVTK_numPoints(vtp, lFa%nNo, istat)
         if (istat .lt. 0) exit

         call getVTK_numElems(vtp, lFa%nEl, istat)
         if (istat .lt. 0) exit

         call getVTK_nodesPerElem(vtp, lFa%eNoN, istat)
         if (istat .lt. 0) exit

         allocate(lFa%x(nsd,lFa%nNo), lFa%ien(lFa%eNoN,lFa%nEl), &
     &      lFa%gN(lFa%nNo), lFa%gE(lFa%nEl))
         call selecteleb(lFa)

         allocate(tmpX2(maxNSD,lFa%nNo))
         call getVTK_pointCoords(vtp, tmpX2, istat)
         if (istat .lt. 0) exit
         lFa%x(:,:) = tmpX2(1:nsd,:)
         deallocate(tmpX2)

         call getVTK_elemIEN(vtp, lFa%ien, istat)
         if (istat .lt. 0) exit
         lFa%IEN(:,:) = lFa%IEN(:,:) + 1

         call getVTK_pointData(vtp, 'GlobalNodeID', lFa%gN, istat)
         if (istat .lt. 0) exit

         call getVTK_elemData(vtp, 'GlobalElementID', lFa%gE, istat)
         if (istat .lt. 0) exit

         call flushVTK(vtp)
         exit
      end do

      if (istat .lt. 0) then
         write(stdout,ftab4) "Error: vtk file read error"
         STOP
      end if

      return
      end subroutine readVTP

!**********************************************************************

      subroutine writeVTP(lFa, fName)
      use commod
      implicit none

      type(faceType), intent(in) :: lFa
      character(len=strL), intent(in) :: fName

      type(vtkXMLtype) :: vtp
      integer :: iStat

      call vtkInitWriter(vtp, trim(fName), iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTP file write error (init)"
         stop
      end if

      call putVTK_pointCoords(vtp, lFa%x, iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTP file write error (coords)"
         stop
      end if

      call putVTK_elemIEN(vtp, lFa%IEN, lFa%vtktype, iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTP file write error (ien)"
         stop
      end if

      call putVTK_pointData(vtp, "GlobalNodeID", lFa%gN, iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTP file write error (point data)"
         stop
      end if

      call putVTK_elemData(vtp, "GlobalElementID", lFa%gE, iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTP file write error (cell data)"
         stop
      end if

      call vtkWriteToFile(vtp, iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTP file write error"
         stop
      end if

      call flushVTK(vtp)

      return
      end subroutine writeVTP

!**************************************************


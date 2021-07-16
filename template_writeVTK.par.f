c            fid = 10000+cm%tF()
c            WRITE(fName,'(A)') "res_"//STR(cm%tF())//"_"//STR(cTS)//
c     2         "_"//STR(cntr)//".vtk"
c            OPEN(fid,FILE=TRIM(fName),STATUS='UNKNOWN',ACCESS='STREAM',
c     2         FORM='UNFORMATTED',CONVERT='BIG_ENDIAN')
c            WRITE(fid) "# vtk DataFile Version 3.0"//eol
c            WRITE(fid) "Residue"//eol
c            WRITE(fid) "BINARY"//eol
c            WRITE(fid) "DATASET UNSTRUCTURED_GRID"//eol
c            WRITE(fid) "POINTS "//STR(SUM(msh(:)%nNo))//" double"//eol
c            itmp = 0
c            DO iM=1, nMsh
c               DO a=1, msh(iM)%nNo
c                  Ac = msh(iM)%gN(a)
c                  IF (nsd .EQ. 2) THEN
c                     WRITE(fid) x(1:nsd,Ac), 0D0
c                  ELSE
c                     WRITE(fid) x(1:nsd,Ac)
c                  END IF
c               END DO
c               itmp = itmp + (msh(iM)%eNoN+1)*msh(iM)%nEl
c            END DO
c            WRITE(fid) eol//"CELLS "//STR(SUM(msh(:)%nEl))//" "//
c     2         STR(itmp)//eol
c            DO iM=1, nMsh
c               itmp = 0
c               DO e=1, msh(iM)%nEl
c                  WRITE(fid) msh(iM)%eNoN, msh(iM)%IEN(:,e)+iTmp-1
c               END DO
c               itmp = itmp + msh(iM)%nNo
c            END DO
c            WRITE(fid) eol//"CELL_TYPES "//STR(SUM(msh(:)%nEl))//eol
c            DO iM=1, nMsh
c               DO e=1, msh(iM)%nEl
c                  WRITE(fid) msh(iM)%vtkType
c               END DO
c            END DO
c            WRITE(fid) eol//"POINT_DATA "//STR(SUM(msh(:)%nNo))//eol
c            WRITE(fid) "VECTORS Res_M double"//eol
c            DO iM=1, nMsh
c               DO a=1, msh(iM)%nNo
c                  Ac = msh(iM)%gN(a)
c                  IF (nsd .EQ. 2) THEN
c                     WRITE(fid) R(1:nsd,Ac), 0D0
c                  ELSE
c                     WRITE(fid) R(1:nsd,Ac)
c                  END IF
c               END DO
c            END DO
c            WRITE(fid) eol//"SCALARS Res_C double"//eol
c            WRITE(fid) "LOOKUP_TABLE default"//eol
c            DO iM=1, nMsh
c               DO a=1, msh(iM)%nNo
c                  Ac = msh(iM)%gN(a)
c                  WRITE(fid) R(nsd+1,Ac)
c               END DO
c            END DO
c            WRITE(fid) eol//"CELL_DATA "//STR(SUM(msh(:)%nEl))//eol
c            WRITE(fid) "SCALARS Proc_ID int"//eol
c            WRITE(fid) "LOOKUP_TABLE default"//eol
c            DO iM=1, nMsh
c               DO e=1, msh(iM)%nEl
c                  WRITE(fid) cm%tF()
c               END DO
c            END DO
c            CALL FLUSH(fid)
c            CLOSE(fid)
            

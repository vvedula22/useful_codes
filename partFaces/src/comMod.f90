!**************************************************

      module params

      integer, parameter :: eType_NA = 100, eType_LIN = 101, &
         eType_TRI = 102, eType_TET = 103, eType_BIL = 104, &
         eType_QUD = 105, eType_BIQ = 106, eType_BRK = 107, &
         eType_NRB = 108, eType_WDG = 109, eType_QTR = 110, &
         eType_QTE = 111

      end module params

!**************************************************

      module commod
      use genUtils
      use stdParams
      use vtkXMLMod
      use params
      implicit none

      type fsType
         integer :: eNoN
         integer :: nG
         integer :: eType = eType_NA
         real(kind=8), allocatable :: w(:)
         real(kind=8), allocatable :: xi(:,:)
         real(kind=8), allocatable :: N(:,:)
         real(kind=8), allocatable :: Nx(:,:,:)
      end type fsType

      type dataArrType
         integer :: nv
         integer :: m, n
         integer, allocatable :: s(:)
         integer, allocatable :: e(:)
         real(kind=8), allocatable :: va(:,:)
         character(len=strL), allocatable :: vname(:)
      end type dataArrType

      type meshType
         logical :: qShpF
         integer :: nNo
         integer :: nEl
         integer :: eNoN
         integer :: nG
         integer :: eType = eType_NA
         integer :: vtkType
         integer, allocatable :: IEN(:,:)
         real(kind=8), allocatable :: w(:)
         real(kind=8), allocatable :: xi(:,:)
         real(kind=8), allocatable :: N(:,:)
         real(kind=8), allocatable :: Nx(:,:,:)
         real(kind=8), allocatable :: x(:,:)
      end type meshType

      type, extends(meshType) :: faceType
         integer, allocatable :: gE(:)
         integer, allocatable :: gN(:)
      end type faceType

      integer :: nsd, nstd

      interface norm
         module procedure norms, normv
      end interface norm

      interface destroy
         module procedure destroyMesh, destroyFace, destroyData
      end interface destroy

      contains
!-------------------------------------------------
         pure function cross(U) result(V)
         implicit none
         real(kind=8), intent(in) :: U(:,:)
         real(kind=8) V(size(U,1))

         if (size(U,1) .eq. 2) then
            V(1) =  U(2,1)
            V(2) = -U(1,1)
         else
            V(1) = U(2,1)*U(3,2) - U(3,1)*U(2,2)
            V(2) = U(3,1)*U(1,2) - U(1,1)*U(3,2)
            V(3) = U(1,1)*U(2,2) - U(2,1)*U(1,2)
         end if

         return
         end function cross
!-------------------------------------------------
         function norms(u) result(res)
         implicit none
         real(kind=8), intent(in) :: u(:)
         real(kind=8) res
         integer i, m

         m   = size(u)
         res = 0.0D0
         do i=1, m
            res = res + u(i)*u(i)
         end do

         return
         end function norms
!-------------------------------------------------
         function normv(u, v) result(res)
         implicit none
         real(kind=8), intent(in) :: u(:), v(:)
         real(kind=8) res
         integer i, m

         if (size(u) .ne. size(v)) then
            write(stdout,ftab4) "ERROR: incompatible vectors to "// &
     &         "compute norm"
            stop
         end if
         m   = size(u)
         res = 0.0D0
         do i=1, m
            res = res + (u(i)-v(i))*(u(i)-v(i))
         end do

         return
         end function normv
!-------------------------------------------------
         subroutine destroyMesh(lM)
         implicit none
         type(meshType), intent(inout) :: lM

         if (allocated(lM%IEN)) deallocate(lM%IEN)
         if (allocated(lM%w))   deallocate(lM%w)
         if (allocated(lM%xi))  deallocate(lM%xi)
         if (allocated(lM%N))   deallocate(lM%N)
         if (allocated(lM%Nx))  deallocate(lM%Nx)
         if (allocated(lM%x))   deallocate(lM%x)

         return
         end subroutine destroyMesh
!-------------------------------------------------
         subroutine destroyFace(lFa)
         implicit none
         type(faceType), intent(inout) :: lFa

         if (allocated(lFa%IEN))  deallocate(lFa%IEN)
         if (allocated(lFa%w))    deallocate(lFa%w)
         if (allocated(lFa%xi))   deallocate(lFa%xi)
         if (allocated(lFa%N))    deallocate(lFa%N)
         if (allocated(lFa%Nx))   deallocate(lFa%Nx)
         if (allocated(lFa%x))    deallocate(lFa%x)
         if (allocated(lFa%gN))   deallocate(lFa%gN)
         if (allocated(lFa%gE))   deallocate(lFa%gE)

         return
         end subroutine destroyFace
!-------------------------------------------------
         subroutine destroyData(d)
         implicit none
         type(dataArrType), intent(inout) :: d

         if (allocated(d%va))    deallocate(d%va)
         if (allocated(d%vname)) deallocate(d%vname)

         return
         end subroutine destroyData
!-------------------------------------------------
         function integ(lM, u, flag) result(res)
         implicit none
         type(meshType), intent(inout) :: lM
         real(kind=8), intent(in) :: u(:)
         logical, intent(in), optional :: flag
         real(kind=8) :: res

         logical :: pflag
         integer :: a, e, g, Ac, eNoN, nEl, nNo
         real(kind=8) :: Jac, tmp(nsd,nsd), sHat
         type(fsType) :: fs

         real(kind=8), allocatable :: Nx(:,:), xl(:,:), ul(:)

         if (size(u) .ne. lM%nNo) then
            write(stdout,ftab4) "ERROR: inconsistent size of input data"
            stop
         end if

         pflag = .false.
         if (present(flag)) pflag = flag

         if (.not.pflag) then
            fs%eNoN  = lM%eNoN
            fs%nG    = lM%nG
            fs%eType = lM%eType
            allocate(fs%w(fs%nG), fs%xi(nsd,fs%nG), fs%N(fs%eNoN,fs%nG),&
     &         fs%Nx(nsd,fs%eNoN,fs%nG))
            do g=1, fs%nG
               fs%w  = lM%w
               fs%xi = lM%xi
               fs%N  = lM%N
               fs%Nx = lM%Nx
            end do
         else
            call setReducedIntegFS(fs, lM%eType)
            allocate(fs%w(fs%nG), fs%xi(nsd,fs%nG), fs%N(fs%eNoN,fs%nG),&
     &         fs%Nx(nsd,fs%eNoN,fs%nG))
            call getGP(nsd, fs%eType, fs%nG, fs%w, fs%xi)

            do g=1, fs%nG
               call getShpF(nsd, fs%eType, fs%eNoN, fs%xi(:,g), &
     &            fs%N(:,g), fs%Nx(:,:,g))
            end do
         end if

         nNo  = lM%nNo
         nEl  = lM%nEl
         eNoN = fs%eNoN
         allocate(Nx(nsd,eNoN), xl(nsd,eNoN), ul(eNoN))

         res = 0D0
         do e=1, nEl
            xl = 0D0
            ul = 0D0
            do a=1, eNoN
               Ac = lM%IEN(a,e)
               xl(:,a) = lM%x(:,Ac)
               ul(a)   = u(Ac)
            end do

            do g=1, fs%nG
               call GNN(eNoN, fs%Nx(:,:,g), xl, Nx, Jac, tmp)

               sHat = 0D0
               do a=1, eNoN
                  sHat = sHat + (ul(a)*fs%N(a,g))
               end do
               res = res + (fs%w(g)*sHat*Jac)
            end do
         end do
         deallocate(Nx, xl, ul)

         return
         end function integ
!-------------------------------------------------
      end module commod

!**************************************************

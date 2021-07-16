!**************************************************

      subroutine selectele(lM)
      use commod
      implicit none
      type(meshType), intent(inout) :: lM

      integer :: g

      if (nsd .eq. 3) then
         select case (lM%eNoN)
         case (8)
            lM%eType   = eType_BRK
            lM%nG      = 8
            lM%vtkType = 12
            lM%qShpF   = .false.
         case (6)
            lM%eType   = eType_WDG
            lM%nG      = 6
            lM%vtkType = 13
            lM%qShpF   = .false.
         case (4)
            lM%eType   = eType_TET
            lM%nG      = 4
            lM%vtkType = 10
            lM%qShpF   = .false.
         case (10)
            lM%eType   = eType_QTE
            lM%nG      = 15
            lM%vtkType = 24
            lM%qShpF   = .true.
         case default
            write(stdout,ftab4) &
               "Error: unknown combination of nsd & eNoN"
            STOP
         end select
      else
         select case (lM%eNoN)
         case (3)
            lM%eType   = eType_TRI
            lM%nG      = 3
            lM%vtkType = 5
            lM%qShpF   = .false.
         case (4)
            lM%eType   = eType_BIL
            lM%nG      = 4
            lM%vtkType = 9
            lM%qShpF   = .false.
         case (6)
            lM%eType   = eType_QTR
            lM%nG      = 7
            lM%vtkType = 22
            lM%qShpF   = .true.
         case (9)
            lM%eType   = eType_BIQ
            lM%nG      = 9
            lM%vtkType = 28
            lM%qShpF   = .true.
         case default
            write(stdout,ftab4) &
               "Error: unknown combination of nsd & eNoN"
            STOP
         end select
      end if

      allocate(lM%w(lM%nG), lM%xi(nsd,lM%nG), lM%N(lM%eNoN,lM%nG), &
         lM%Nx(nsd,lM%eNoN,lM%nG))

      call getGP(nsd, lM%eType, lM%nG, lM%w, lM%xi)

      do g=1, lM%nG
         call getShpF(nsd, lM%eType, lM%eNoN, lM%xi(:,g), lM%N(:,g), &
            lM%Nx(:,:,g))
      end do

      return
      end subroutine selectele

!**************************************************

      subroutine selecteleb(lFa)
      use commod
      implicit none
      type(faceType), intent(inout) :: lFa

      integer :: g, insd

      insd = nsd-1

      if (insd .eq. 2) then
         select case (lFa%eNoN)
         case (3)
            lFa%eType   = eType_TRI
            lFa%nG      = 3
            lFa%vtkType = 5
            lFa%qShpF   = .false.
         case (4)
            lFa%eType   = eType_BIL
            lFa%nG      = 4
            lFa%vtkType = 9
            lFa%qShpF   = .false.
         case (6)
            lFa%eType   = eType_QTR
            lFa%nG      = 7
            lFa%vtkType = 22
            lFa%qShpF   = .true.
         case (9)
            lFa%eType   = eType_BIQ
            lFa%nG      = 9
            lFa%vtkType = 28
            lFa%qShpF   = .true.
         case default
            write(stdout,ftab4) &
               "Error: unknown combination of nsd & eNoN"
            STOP
         end select
      else if (insd .eq. 1) then
         select case (lFa%eNoN)
         case (2)
            lFa%eType   = eType_LIN
            lFa%nG      = 2
            lFa%vtkType = 3
            lFa%qShpF   = .false.
         case (3)
            lFa%eType   = eType_QUD
            lFa%nG      = 3
            lFa%vtkType = 21
            lFa%qShpF   = .false.
         case default
            write(stdout,ftab4) &
               "Error: unknown combination of nsd & eNoN"
            STOP
         end select
      end if

      allocate(lFa%w(lFa%nG), lFa%xi(nsd,lFa%nG), &
         lFa%N(lFa%eNoN,lFa%nG), lFa%Nx(nsd,lFa%eNoN,lFa%nG))

      call getGP(insd, lFa%eType, lFa%nG, lFa%w, lFa%xi)

      do g=1, lFa%nG
         call getShpF(nsd, lFa%eType, lFa%eNoN, lFa%xi(:,g), lFa%N(:,g),&
     &      lFa%Nx(:,:,g))
      end do

      return
      end subroutine selecteleb

!**************************************************

      subroutine setReducedIntegFS(fs, eType)
      use variables
      implicit none
      type(fsType), intent(inout) :: fs
      integer, intent(in) :: eType

      select case (eType)
      case (eType_QTE)
         fs%eType   = eType_TET
         fs%nG      = 4
         fs%eNoN    = 4
      case (eType_QTR)
         fs%eType   = eType_TRI
         fs%nG      = 3
         fs%eNoN    = 3
      case (eType_BIQ)
         fs%eType   = eType_BIL
         fs%nG      = 4
         fs%eNoN    = 4
      case default
         write(stdout,ftab4) "ERROR: unknown higher order element type"
         stop
      end select

      return
      end subroutine setReducedIntegFS

!**************************************************

      subroutine getGP(insd, eType, nG, w, xi)
      use params
      implicit none
      integer, intent(in) :: insd, eType, nG
      real(kind=8), intent(out) :: w(nG), xi(insd,nG)

      real(kind=8) s, t, lz, uz

!     3D elements
      select case (eType)
      case (eType_BRK)
         w = 1D0
         s =  1D0/SQRT(3D0)
         t = -1D0/SQRT(3D0)
         xi(1,1) = s; xi(2,1) = s; xi(3,1) = s
         xi(1,2) = t; xi(2,2) = s; xi(3,2) = s
         xi(1,3) = t; xi(2,3) = s; xi(3,3) = t
         xi(1,4) = s; xi(2,4) = s; xi(3,4) = t
         xi(1,5) = s; xi(2,5) = t; xi(3,5) = s
         xi(1,6) = t; xi(2,6) = t; xi(3,6) = s
         xi(1,7) = t; xi(2,7) = t; xi(3,7) = t
         xi(1,8) = s; xi(2,8) = t; xi(3,8) = t
      case (eType_TET)
         w = 1D0/24D0
         s = (5D0 + 3D0*SQRT(5D0))/2D1
         t = (5D0 -     SQRT(5D0))/2D1
         xi(1,1) = s; xi(2,1) = t; xi(3,1) = t
         xi(1,2) = t; xi(2,2) = s; xi(3,2) = t
         xi(1,3) = t; xi(2,3) = t; xi(3,3) = s
         xi(1,4) = t; xi(2,4) = t; xi(3,4) = t
      case (eType_WDG)
         w  =  1D0/6D0
         s  =  2D0/3D0
         t  =  1D0/6D0
         uz =  1D0/SQRT(3D0)
         lz = -1D0/SQRT(3D0)
         xi(1,1) = s; xi(2,1) = t; xi(3,1) = lz
         xi(1,2) = t; xi(2,2) = s; xi(3,2) = lz
         xi(1,3) = t; xi(2,3) = t; xi(3,3) = lz
         xi(1,4) = s; xi(2,4) = t; xi(3,4) = uz
         xi(1,5) = t; xi(2,5) = s; xi(3,5) = uz
         xi(1,6) = t; xi(2,6) = t; xi(3,6) = uz
      case (eType_QTE)
         w(1)     = 0.030283678097089D0
         w(2:5)   = 0.006026785714286D0
         w(6:9)   = 0.011645249086029D0
         w(10:15) = 0.010949141561386D0

         s = 0.25D0
         xi(1,1) = s; xi(2,1) = s; xi(3,1) = s

         s = 0.333333333333333D0
         t = 0.0D0
         xi(1,2) = t; xi(2,2) = s; xi(3,2) = s
         xi(1,3) = s; xi(2,3) = t; xi(3,3) = s
         xi(1,4) = s; xi(2,4) = s; xi(3,4) = t
         xi(1,5) = s; xi(2,5) = s; xi(3,5) = s

         s = 0.090909090909091D0
         t = 0.727272727272727D0
         xi(1,6) = t; xi(2,6) = s; xi(3,6) = s
         xi(1,7) = s; xi(2,7) = t; xi(3,7) = s
         xi(1,8) = s; xi(2,8) = s; xi(3,8) = t
         xi(1,9) = s; xi(2,9) = s; xi(3,9) = s

         s = 0.066550153573664D0
         t = 0.433449846426336D0
         xi(1,10) = s; xi(2,10) = s; xi(3,10) = t
         xi(1,11) = s; xi(2,11) = t; xi(3,11) = s
         xi(1,12) = s; xi(2,12) = t; xi(3,12) = t
         xi(1,13) = t; xi(2,13) = t; xi(3,13) = s
         xi(1,14) = t; xi(2,14) = s; xi(3,14) = t
         xi(1,15) = t; xi(2,15) = s; xi(3,15) = s

!     2D elements
      case (eType_TRI)
         w = 1D0/6D0
         s = 2D0/3D0
         t = 1D0/6D0
         xi(1,1) = s; xi(2,1) = t
         xi(1,2) = t; xi(2,2) = s
         xi(1,3) = t; xi(2,3) = t
      case (eType_BIL)
         w = 1D0
         s =  1D0/SQRT(3D0)
         t = -1D0/SQRT(3D0)
         xi(1,1) = s; xi(2,1) = s
         xi(1,2) = t; xi(2,2) = s
         xi(1,3) = t; xi(2,3) = t
         xi(1,4) = s; xi(2,4) = t
      case (eType_BIQ)
         w(1) = 25D0/81D0; w(2) = 25D0/81D0; w(3) = 25D0/81D0
         w(4) = 25D0/81D0; w(5) = 40D0/81D0; w(6) = 40D0/81D0
         w(7) = 40D0/81D0; w(8) = 40D0/81D0; w(9) = 64D0/81D0
         s    = SQRT(6D-1)
         xi(1,1) =  -s; xi(2,1) =  -s
         xi(1,2) =   s; xi(2,2) =  -s
         xi(1,3) =   s; xi(2,3) =   s
         xi(1,4) =  -s; xi(2,4) =   s
         xi(1,5) = 0D0; xi(2,5) =  -s
         xi(1,6) =   s; xi(2,6) = 0D0
         xi(1,7) = 0D0; xi(2,7) =   s
         xi(1,8) =  -s; xi(2,8) = 0D0
         xi(1,9) = 0D0; xi(2,9) = 0D0
      CASE(eType_QTR)
         w(1)   = 0.225000000000000D0 * 5D-1
         w(2:4) = 0.125939180544827D0 * 5D-1
         w(5:7) = 0.132394152788506D0 * 5D-1

         s = 0.333333333333333D0
         xi(1,1) = s; xi(2,1) = s

         s = 0.797426985353087D0
         t = 0.101286507323456D0
         xi(1,2) = s; xi(2,2) = t
         xi(1,3) = t; xi(2,3) = s
         xi(1,4) = t; xi(2,4) = t

         s = 0.059715871789770D0
         t = 0.470142064105115D0
         xi(1,5) = s; xi(2,5) = t
         xi(1,6) = t; xi(2,6) = s
         xi(1,7) = t; xi(2,7) = t

!     1D elements
      case (eType_LIN)
         w = 1D0
         s = 1D0/SQRT(3D0)
         xi(1,1) = -s
         xi(1,2) =  s
      case (eType_QUD)
         w(1) = 5D0/9D0; w(2) = 5D0/9D0; w(3) = 8D0/9D0
         s = SQRT(6D-1)
         xi(1,1) = -s
         xi(1,2) =  s
         xi(1,3) = 0D0
      end select

      return
      end subroutine getGP

!**************************************************

      subroutine getShpF(insd, eType, eNoN, xi, N, Nxi)
      use params
      implicit none
      integer, intent(in) :: insd, eType, eNoN
      real(kind=8), intent(out) :: xi(insd), N(eNoN), Nxi(insd,eNoN)

      real(kind=8) :: s, t, mx, my, ux, uy, uz, lx, ly, lz

!     3D elements
      select case (eType)
      case (eType_BRK)
         ux = 1D0 + xi(1); lx = 1D0 - xi(1)
         uy = 1D0 + xi(2); ly = 1D0 - xi(2)
         uz = 1D0 + xi(3); lz = 1D0 - xi(3)
         N(1) = ux*uy*uz/8D0
         N(2) = lx*uy*uz/8D0
         N(3) = lx*uy*lz/8D0
         N(4) = ux*uy*lz/8D0
         N(5) = ux*ly*uz/8D0
         N(6) = lx*ly*uz/8D0
         N(7) = lx*ly*lz/8D0
         N(8) = ux*ly*lz/8D0

         Nxi(1,1) =  uy*uz/8D0
         Nxi(2,1) =  ux*uz/8D0
         Nxi(3,1) =  ux*uy/8D0
         Nxi(1,2) = -uy*uz/8D0
         Nxi(2,2) =  lx*uz/8D0
         Nxi(3,2) =  lx*uy/8D0
         Nxi(1,3) = -uy*lz/8D0
         Nxi(2,3) =  lx*lz/8D0
         Nxi(3,3) = -lx*uy/8D0
         Nxi(1,4) =  uy*lz/8D0
         Nxi(2,4) =  ux*lz/8D0
         Nxi(3,4) = -ux*uy/8D0
         Nxi(1,5) =  ly*uz/8D0
         Nxi(2,5) = -ux*uz/8D0
         Nxi(3,5) =  ux*ly/8D0
         Nxi(1,6) = -ly*uz/8D0
         Nxi(2,6) = -lx*uz/8D0
         Nxi(3,6) =  lx*ly/8D0
         Nxi(1,7) = -ly*lz/8D0
         Nxi(2,7) = -lx*lz/8D0
         Nxi(3,7) = -lx*ly/8D0
         Nxi(1,8) =  ly*lz/8D0
         Nxi(2,8) = -ux*lz/8D0
         Nxi(3,8) = -ux*ly/8D0

      case (eType_TET)
         N(1) = xi(1)
         N(2) = xi(2)
         N(3) = xi(3)
         N(4) = 1D0 - xi(1) - xi(2) - xi(3)

         Nxi(1,1) =  1D0
         Nxi(2,1) =  0D0
         Nxi(3,1) =  0D0
         Nxi(1,2) =  0D0
         Nxi(2,2) =  1D0
         Nxi(3,2) =  0D0
         Nxi(1,3) =  0D0
         Nxi(2,3) =  0D0
         Nxi(3,3) =  1D0
         Nxi(1,4) = -1D0
         Nxi(2,4) = -1D0
         Nxi(3,4) = -1D0

      case (eType_WDG)
         ux = xi(1) ; uy = xi(2) ; uz = 1D0 - ux - uy
         s = (1D0 + xi(3))/2D0; t = (1D0 - xi(3))/2D0
         N(1) = ux*t
         N(2) = uy*t
         N(3) = uz*t
         N(4) = ux*s
         N(5) = uy*s
         N(6) = uz*s

         Nxi(1,1) =  t
         Nxi(2,1) =  0D0
         Nxi(3,1) = -ux/2D0
         Nxi(1,2) =  0D0
         Nxi(2,2) =  t
         Nxi(3,2) = -uy/2D0
         Nxi(1,3) = -t
         Nxi(2,3) = -t
         Nxi(3,3) = -uz/2D0
         Nxi(1,4) =  s
         Nxi(2,4) =  0D0
         Nxi(3,4) =  ux/2D0
         Nxi(1,5) =  0D0
         Nxi(2,5) =  s
         Nxi(3,5) =  uy/2D0
         Nxi(1,6) = -s
         Nxi(2,6) = -s
         Nxi(3,6) =  uz/2D0

      case (eType_QTE)
         s     = 1.0D0 - xi(1) - xi(2) - xi(3)
         N(1)  = xi(1)*(2.0D0*xi(1) - 1.0D0)
         N(2)  = xi(2)*(2.0D0*xi(2) - 1.0D0)
         N(3)  = xi(3)*(2.0D0*xi(3) - 1.0D0)
         N(4)  = s    *(2.0D0*s     - 1.0D0)
         N(5)  = 4.0D0*xi(1)*xi(2)
         N(6)  = 4.0D0*xi(2)*xi(3)
         N(7)  = 4.0D0*xi(1)*xi(3)
         N(8)  = 4.0D0*xi(1)*s
         N(9)  = 4.0D0*xi(2)*s
         N(10) = 4.0D0*xi(3)*s

         Nxi(1,1)  =  4.0D0*xi(1) - 1.0D0
         Nxi(2,1)  =  0.0D0
         Nxi(3,1)  =  0.0D0
         Nxi(1,2)  =  0.0D0
         Nxi(2,2)  =  4.0D0*xi(2) - 1.0D0
         Nxi(3,2)  =  0.0D0
         Nxi(1,3)  =  0.0D0
         Nxi(2,3)  =  0.0D0
         Nxi(3,3)  =  4.0D0*xi(3) - 1.0D0
         Nxi(1,4)  =  1.0D0 - 4.0D0*s
         Nxi(2,4)  =  1.0D0 - 4.0D0*s
         Nxi(3,4)  =  1.0D0 - 4.0D0*s
         Nxi(1,5)  =  4.0D0*xi(2)
         Nxi(2,5)  =  4.0D0*xi(1)
         Nxi(3,5)  =  0.0D0
         Nxi(1,6)  =  0.0D0
         Nxi(2,6)  =  4.0D0*xi(3)
         Nxi(3,6)  =  4.0D0*xi(2)
         Nxi(1,7)  =  4.0D0*xi(3)
         Nxi(2,7)  =  0.0D0
         Nxi(3,7)  =  4.0D0*xi(1)
         Nxi(1,8)  =  4.0D0*( s - xi(1) )
         Nxi(2,8)  = -4.0D0*xi(1)
         Nxi(3,8)  = -4.0D0*xi(1)
         Nxi(1,9)  = -4.0D0*xi(2)
         Nxi(2,9)  =  4.0D0*( s - xi(2) )
         Nxi(3,9)  = -4.0D0*xi(2)
         Nxi(1,10) = -4.0D0*xi(3)
         Nxi(2,10) = -4.0D0*xi(3)
         Nxi(3,10) =  4.0D0*( s - xi(3) )

!     2D elements
      case (eType_TRI)
         N(1) = xi(1)
         N(2) = xi(2)
         N(3) = 1D0 - xi(1) - xi(2)

         Nxi(1,1) =  1D0
         Nxi(2,1) =  0D0
         Nxi(1,2) =  0D0
         Nxi(2,2) =  1D0
         Nxi(1,3) = -1D0
         Nxi(2,3) = -1D0

      case (eType_BIL)
         ux = 1D0 + xi(1); lx = 1D0 - xi(1)
         uy = 1D0 + xi(2); ly = 1D0 - xi(2)
         N(1) = ux*uy/4D0
         N(2) = lx*uy/4D0
         N(3) = lx*ly/4D0
         N(4) = ux*ly/4D0

         Nxi(1,1) =  uy/4D0
         Nxi(2,1) =  ux/4D0
         Nxi(1,2) = -uy/4D0
         Nxi(2,2) =  lx/4D0
         Nxi(1,3) = -ly/4D0
         Nxi(2,3) = -lx/4D0
         Nxi(1,4) =  ly/4D0
         Nxi(2,4) = -ux/4D0

      case (eType_BIQ)
         ux = 1D0 + xi(1); mx = xi(1); lx = 1D0 - xi(1)
         uy = 1D0 + xi(2); my = xi(2); ly = 1D0 - xi(2)
         N(1) =  mx*lx*my*ly/4D0
         N(2) = -mx*ux*my*ly/4D0
         N(3) =  mx*ux*my*uy/4D0
         N(4) = -mx*lx*my*uy/4D0
         N(5) = -lx*ux*my*ly/2D0
         N(6) =  mx*ux*ly*uy/2D0
         N(7) =  lx*ux*my*uy/2D0
         N(8) = -mx*lx*ly*uy/2D0
         N(9) =  lx*ux*ly*uy

         Nxi(1,1) =  (lx - mx)*my*ly/4D0
         Nxi(2,1) =  (ly - my)*mx*lx/4D0
         Nxi(1,2) = -(ux + mx)*my*ly/4D0
         Nxi(2,2) = -(ly - my)*mx*ux/4D0
         Nxi(1,3) =  (ux + mx)*my*uy/4D0
         Nxi(2,3) =  (uy + my)*mx*ux/4D0
         Nxi(1,4) = -(lx - mx)*my*uy/4D0
         Nxi(2,4) = -(uy + my)*mx*lx/4D0
         Nxi(1,5) = -(lx - ux)*my*ly/2D0
         Nxi(2,5) = -(ly - my)*lx*ux/2D0
         Nxi(1,6) =  (ux + mx)*ly*uy/2D0
         Nxi(2,6) =  (ly - uy)*mx*ux/2D0
         Nxi(1,7) =  (lx - ux)*my*uy/2D0
         Nxi(2,7) =  (uy + my)*lx*ux/2D0
         Nxi(1,8) = -(lx - mx)*ly*uy/2D0
         Nxi(2,8) = -(ly - uy)*mx*lx/2D0
         Nxi(1,9) =  (lx - ux)*ly*uy
         Nxi(2,9) =  (ly - uy)*lx*ux

      case (eType_QTR)
         s    = 1.0D0 - xi(1) - xi(2)
         N(1) = xi(1)*( 2.0D0*xi(1) - 1.0D0 )
         N(2) = xi(2)*( 2.0D0*xi(2) - 1.0D0 )
         N(3) = s    *( 2.0D0*s     - 1.0D0 )
         N(4) = 4.0D0*xi(1)*xi(2)
         N(5) = 4.0D0*xi(2)*s
         N(6) = 4.0D0*xi(1)*s

         Nxi(1,1) =  4.0D0*xi(1) - 1.0D0
         Nxi(2,1) =  0.0D0
         Nxi(1,2) =  0.0D0
         Nxi(2,2) =  4.0D0*xi(2) - 1.0D0
         Nxi(1,3) =  1.0D0 - 4.0D0*s
         Nxi(2,3) =  1.0D0 - 4.0D0*s
         Nxi(1,4) =  4.0D0*xi(2)
         Nxi(2,4) =  4.0D0*xi(1)
         Nxi(1,5) = -4.0D0*xi(2)
         Nxi(2,5) =  4.0D0*( s - xi(2) )
         Nxi(1,6) =  4.0D0*( s - xi(1) )
         Nxi(2,6) = -4.0D0*xi(1)

!     1D elements
      case (eType_LIN)
         N(1) = (1D0 - xi(1))/2D0
         N(2) = (1D0 + xi(1))/2D0

         Nxi(1,1) = -5D-1
         Nxi(1,2) =  5D-1
      case (eType_QUD)
         N(1) = -xi(1)*(1D0 - xi(1))/2D0
         N(2) =  xi(1)*(1D0 + xi(1))/2D0
         N(3) = (1D0 - xi(1))*(1D0 + xi(1))

         Nxi(1,1) = -5D-1 + xi(1)
         Nxi(1,2) =  5D-1 + xi(1)
         Nxi(1,3) = -2D0*xi(1)
      end select

      return
      end subroutine getShpF

!**************************************************

      subroutine GNN(eNoN, Nxi, x, Nx, Jac, ks)
      use commod, only: nsd
      implicit none

      integer, intent(in) :: eNoN
      real(kind=8), intent(in) :: Nxi(nsd,eNoN), x(nsd,eNoN)
      real(kind=8), intent(out) :: Nx(nsd,eNoN), Jac, ks(nsd,nsd)

      integer :: a
      real(kind=8) :: xXi(nsd,nsd), xiX(nsd,nsd)

      Nx  = 0D0
      xXi = 0D0
      if (nsd .eq. 2) then
         do a=1, eNoN
            xXi(:,1) = xXi(:,1) + x(:,a)*Nxi(1,a)
            xXi(:,2) = xXi(:,2) + x(:,a)*Nxi(2,a)
         end do

         Jac = xXi(1,1)*xXi(2,2) - xXi(1,2)*xXi(2,1)

         xiX(1,1) =  xXi(2,2)/Jac
         xiX(1,2) = -xXi(1,2)/Jac
         xiX(2,1) = -xXi(2,1)/Jac
         xiX(2,2) =  xXi(1,1)/Jac

         ks(1,1) = xiX(1,1)*xiX(1,1) + xiX(2,1)*xiX(2,1)
         ks(1,2) = xiX(1,1)*xiX(1,2) + xiX(2,1)*xiX(2,2)
         ks(2,2) = xiX(1,2)*xiX(1,2) + xiX(2,2)*xiX(2,2)
         ks(2,1) = ks(1,2)

         do a=1, eNoN
            Nx(1,a) = Nx(1,a)+ Nxi(1,a)*xiX(1,1) + Nxi(2,a)*xiX(2,1)
            Nx(2,a) = Nx(2,a)+ Nxi(1,a)*xiX(1,2) + Nxi(2,a)*xiX(2,2)
         end do
      else
         do a=1, eNoN
            xXi(:,1) = xXi(:,1) + x(:,a)*Nxi(1,a)
            xXi(:,2) = xXi(:,2) + x(:,a)*Nxi(2,a)
            xXi(:,3) = xXi(:,3) + x(:,a)*Nxi(3,a)
         end do

         Jac = xXi(1,1)*xXi(2,2)*xXi(3,3) + &
               xXi(1,2)*xXi(2,3)*xXi(3,1) + &
               xXi(1,3)*xXi(2,1)*xXi(3,2) - &
               xXi(1,1)*xXi(2,3)*xXi(3,2) - &
               xXi(1,2)*xXi(2,1)*xXi(3,3) - &
               xXi(1,3)*xXi(2,2)*xXi(3,1)

         xiX(1,1) = (xXi(2,2)*xXi(3,3) - xXi(2,3)*xXi(3,2))/Jac
         xiX(1,2) = (xXi(3,2)*xXi(1,3) - xXi(3,3)*xXi(1,2))/Jac
         xiX(1,3) = (xXi(1,2)*xXi(2,3) - xXi(1,3)*xXi(2,2))/Jac
         xiX(2,1) = (xXi(2,3)*xXi(3,1) - xXi(2,1)*xXi(3,3))/Jac
         xiX(2,2) = (xXi(3,3)*xXi(1,1) - xXi(3,1)*xXi(1,3))/Jac
         xiX(2,3) = (xXi(1,3)*xXi(2,1) - xXi(1,1)*xXi(2,3))/Jac
         xiX(3,1) = (xXi(2,1)*xXi(3,2) - xXi(2,2)*xXi(3,1))/Jac
         xiX(3,2) = (xXi(3,1)*xXi(1,2) - xXi(3,2)*xXi(1,1))/Jac
         xiX(3,3) = (xXi(1,1)*xXi(2,2) - xXi(1,2)*xXi(2,1))/Jac

         ks(1,1) = xiX(1,1)*xiX(1,1)+xiX(2,1)*xiX(2,1)+xiX(3,1)*xiX(3,1)
         ks(1,2) = xiX(1,2)*xiX(1,1)+xiX(2,2)*xiX(2,1)+xiX(3,2)*xiX(3,1)
         ks(1,3) = xiX(1,3)*xiX(1,1)+xiX(2,3)*xiX(2,1)+xiX(3,3)*xiX(3,1)
         ks(2,2) = xiX(1,2)*xiX(1,2)+xiX(2,2)*xiX(2,2)+xiX(3,2)*xiX(3,2)
         ks(2,3) = xiX(1,2)*xiX(1,3)+xiX(2,2)*xiX(2,3)+xiX(3,2)*xiX(3,3)
         ks(3,3) = xiX(1,3)*xiX(1,3)+xiX(2,3)*xiX(2,3)+xiX(3,3)*xiX(3,3)
         ks(2,1) = ks(1,2)
         ks(3,1) = ks(1,3)
         ks(3,2) = ks(2,3)

         do a=1, eNoN
            Nx(1,a) = Nx(1,a) + Nxi(1,a)*xiX(1,1) + &
                                Nxi(2,a)*xiX(2,1) + &
                                Nxi(3,a)*xiX(3,1)

            Nx(2,a) = Nx(2,a) + Nxi(1,a)*xiX(1,2) + &
                                Nxi(2,a)*xiX(2,2) + &
                                Nxi(3,a)*xiX(3,2)

            Nx(3,a) = Nx(3,a) + Nxi(1,a)*xiX(1,3) + &
                                Nxi(2,a)*xiX(2,3) + &
                                Nxi(3,a)*xiX(3,3)
         end do
      end if

      return
      end subroutine GNN

!**************************************************

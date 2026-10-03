!============================================================================
!
      module ang
!- - - - - - - - - -
      use parameters, only : pr

      implicit none

      real(pr), allocatable, private :: FL(:)

      real(pr), private, parameter :: err = 1.e-15_pr
      real(pr), private, parameter :: eps = 1.e-12_pr

      private :: DL, FRTSJL, FRONTL

!
!============================================================================
      contains
!============================================================================
!
      subroutine ka2jl ( ka, j, l )
!- - - - - - - - - -
      implicit none
! .. Input / Output ..
      integer, intent(in) :: ka
      integer, intent(out) :: j, l
!- - - - - - - - - -

        if ( ka /= 0 ) then

          j = 2 * abs ( ka ) - 1

          l = j + abs ( ka ) / ka

        else

          write(*,'(1x,a14)') 'ka2jl : ka = 0'

        endif

       return

      end subroutine ka2jl
!
!============================================================================
!
      real(pr) function CGcoef ( j1, m1, j2, m2, jtot, mtot )
!- - - - - - - - - -
      implicit none
! .. Input / Output ..
      integer, intent(in) :: j1, m1, j2, m2, jtot, mtot
!- - - - - - - - - -

        CGcoef = ipt( ( j1 - j2 + mtot )/2 ) * sqrt( jtot + 1._pr )               &
               * thrj(  j1, j2, jtot, m1, m2, -mtot )

       return

      end function CGcoef
!
!============================================================================
!
      integer function is_trng ( l1, l2, l3 )
!- - - - - - - - - -
      implicit none
! .. Input / Output ..
      integer, intent(in) :: l1, l2, l3
!- - - - - - - - - -

        if ( ( ( l1+l2-l3)>=0 ) .and. ( ( l1-l2+l3)>=0 ) .and. ( (-l1+l2+l3)>=0 ) ) then

          is_trng = 1

        else

          is_trng = 0

        endif

       return

      end function is_trng
!
!============================================================================
!
      logical function if_trng ( l1, l2, l3 )
!- - - - - - - - - -
      implicit none
! .. Input / Output ..
      integer, intent(in) :: l1, l2, l3
!- - - - - - - - - -

        if ( ( ( l1+l2-l3)>=0 ) .and. ( ( l1-l2+l3)>=0 ) .and. ( (-l1+l2+l3)>=0 ) ) then

          if_trng = .true.

        else

          if_trng = .false.

        endif

       return

      end function if_trng
!
!============================================================================
!
      logical function if_even ( J )
!- - - - - - - - - -
      implicit none
! .. Input / Output ..
      integer, intent(in) :: J
!- - - - - - - - - -

        if ( J/2*2 == J ) then

          if_even = .true.

        else

          if_even = .false.

        endif

       return

      end function if_even
!
!============================================================================
!
      integer function ipt ( ll )
!- - - - - - - - - -
      implicit none
! .. Input / Output ..
      integer, intent(in) :: ll
!- - - - - - - - - -

        if ( (ll/2)*2 == ll ) then
          ipt =  1
        else
          ipt = -1
        endif

       return

      end function ipt
!
!============================================================================
!
      integer function inz ( x )
!- - - - - - - - - -
      implicit none
! .. Input / Output ..
      real(pr), intent(in) :: x
!- - - - - - - - - -

        if ( abs(x) > err ) then
          inz = 1
        else
          inz = 0
        endif

       return

      end function inz
!
!============================================================================
!
      integer function ixmax_J ( J )
!- - - - - - - - - -
      implicit none
! .. Input / Output ..
      integer, intent(in) :: J
!- - - - - - - - - -

        ixmax_J = -1
        if ( J == 0 ) ixmax_J = 0

       return

      end function ixmax_J
!
!============================================================================
!
      real(pr) function G_l ( ka1, ka2, ll )
!- - - - - - - - - -
      implicit none
! .. Input / Output ..
      integer, intent(in) :: ka1, ka2, ll
! .. Local ..
      integer :: jj1, ll1, jj2, ll2
      integer :: jj1_jj2_ll1_ll2
!- - - - - - - - - -

        call ka2jl ( ka1, jj1, ll1 )
        call ka2jl ( ka2, jj2, ll2 )

        jj1_jj2_ll1_ll2 = ( jj1 + 1 ) * ( jj2 + 1 )                     &
     &                  * ( ll1 + 1 ) * ( ll2 + 1 )

        G_l = ipt( (jj2+1)/2 )                                          &
     &      * sqrt( real( jj1_jj2_ll1_ll2, pr ) )                       &
     &      * thrj ( ll1, ll, ll2, 0, 0, 0 )                            &
     &      * sixj ( jj1, ll, jj2, ll2, 1, ll1 )

        if ( abs(G_l) < eps ) G_l = 0._pr

       return

      end function G_l
!
!============================================================================
!
      real(pr) function H_xl ( ka1, ka2, lx, ll )
!- - - - - - - - - -
      implicit none
! .. Input / Output ..
      integer, intent(in) :: ka1, ka2, lx, ll
! .. Local ..
      integer :: jj1, ll1, jj2, ll2
      integer :: jj1_jj2_ll1_ll2
      real(pr) :: dninej
!- - - - - - - - - -

        call ka2jl ( ka1, jj1, ll1 )
        call ka2jl ( ka2, jj2, ll2 )

        jj1_jj2_ll1_ll2 = ( jj1 + 1 ) * ( jj2 + 1 )                     &
     &                  * ( ll1 + 1 ) * ( ll2 + 1 )

        dninej = ninej ( jj1, 1, ll1, lx, 2, ll, jj2, 1, ll2 )

        H_xl = ipt( ll1/2 )                                             &
     &       * sqrt( 6._pr * real( jj1_jj2_ll1_ll2, pr ) )              &
     &       * thrj ( ll1, ll, ll2, 0, 0, 0 )                           &
     &       * dninej

        if ( abs(H_xl) < eps ) H_xl = 0._pr

       return

      end function H_xl
!
!============================================================================
!
      real(pr) function C_J ( ka1, ka2, jj )
!- - - - - - - - - - - - - - - - - - - - - - -
!  C_J(ka1,ka2,JJ) = (-1)^J G_L(ka1,ka2,JJ)
!- - - - - - - - - - - - - - - - - - - - - - -
      implicit none
! .. Input / Output ..
      integer, intent(in) :: ka1, ka2, jj
! .. Local ..
      integer :: jj1, ll1, jj2, ll2
      integer :: jj1_jj2
!- - - - - - - - - -

        call ka2jl ( ka1, jj1, ll1 )
        call ka2jl ( ka2, jj2, ll2 )

        C_J = 0._pr

        if ( .not. if_even( (ll1+ll2+jj)/2 ) ) return

        jj1_jj2 = ( jj1 + 1 ) * ( jj2 + 1 )

        C_J = ipt( (jj1+1)/2 )                                          &
     &      * sqrt( real( jj1_jj2, pr ) )                               &
     &      * thrj ( jj2, jj, jj1, 1, 0, -1 )


       return

      end function C_J
!
!============================================================================
!
      real(pr) function S_JL ( ka1, ka2, jj, ll )
!- - - - - - - - - -
      implicit none
! .. Input / Output ..
      integer, intent(in) :: ka1, ka2, jj, ll
! .. Local ..
      integer :: j
!- - - - - - - - - -

        S_JL = 0._pr

        if ( jj == 0 ) then

          if ( ll == 2 ) then

            S_JL = C_J ( -ka2, ka1, 0 )

          endif

        else

          j = jj / 2

          if     ( ll == jj+2 ) then

            S_JL = sqrt( ( j + 1._pr ) / ( jj + 1._pr ) )               &
     &           * ( 1._pr + ( ka1 + ka2 ) / ( j + 1._pr ) )            &
     &           * C_J ( -ka2, ka1, jj )

          elseif ( ll == jj )   then

            S_JL = ( ka1 - ka2 ) / sqrt( j * ( j + 1._pr ) )            &
     &           * C_J ( ka2, ka1, jj )

          elseif ( ll == jj-2 ) then

            S_JL = sqrt( j / ( jj + 1._pr ) )                           &
     &           * ( -1._pr + ( ka1 + ka2 ) / real(j,pr) )              &
     &           * C_J ( -ka2, ka1, jj )

          endif  ! ll

        endif  ! jj

       return

      end function S_JL
!
!============================================================================
!
      real(pr) function ninej ( ja, jb, jc, jd, je, jf, jg, jh, jj )
!- - - - - - - - - - - - - - - -
!  Calculates 9j sybol
!  All ji are doubled
!- - - - - - - - - - - - - - - -
      implicit none
! .. Input / Output ..
      integer, intent(in) :: ja, jb, jc, jd, je, jf, jg, jh, jj
! .. Local ..
      integer :: klo, khi, k
      real(pr) :: phase
!- - - - - - - - - -

        ninej = 0._pr

        klo = max( abs(ja-jj), abs(jf-jb), abs(jh-jd) )
        khi = min(     ja+jj,      jf+jb,      jh+jd  )

        phase = (-1)**klo

        do k = klo, khi, 2
          ninej = ninej + real( k+1, pr )                               &
     &                    * sixj( ja, jb, jc, jf, jj, k  )              &
     &                    * sixj( jd, je, jf, jb, k,  jh )              &
     &                    * sixj( jg, jh, jj, k,  ja, jd )
        enddo  ! k

        ninej = phase * ninej

       return

      end function ninej
!
!============================================================================
!
      real(pr) function SIXJ ( J1D, J2D, J3D, J4D, J5D, J6D )
!- - - - - - - - - - - - - - - -
!  Calculates 6j sybol
!  All JID are doubled
!- - - - - - - - - - - - - - - -
      implicit none
! .. Input / Output ..
      integer, intent(in) :: J1D, J2D, J3D, J4D, J5D, J6D
! .. Local ..
      integer :: JM1, JM2, JM3, JM4, JX1, JX2, JX3
      integer :: JM, JX
      integer :: KM, KX
      integer :: I, I1
      real(pr) :: TERM1, TERM2, TERM
!- - - - - - - - - -

        SIXJ = 0._pr

        if ( (J1D+J2D-J3D) < 0 ) return
        if ( (J2D+J3D-J1D) < 0 ) return
        if ( (J3D+J1D-J2D) < 0 ) return
        if ( (J1D+J2D+J3D)/2*2 /= (J1D+J2D+J3D) ) return

        if ( (J3D+J4D-J5D) < 0 ) return
        if ( (J4D+J5D-J3D) < 0 ) return
        if ( (J5D+J3D-J4D) < 0 ) return
        if ( (J3D+J4D+J5D)/2*2 /= (J3D+J4D+J5D) ) return

        if ( (J1D+J5D-J6D) < 0 ) return
        if ( (J5D+J6D-J1D) < 0 ) return
        if ( (J6D+J1D-J5D) < 0 ) return
        if ( (J1D+J5D+J6D)/2*2 /= (J1D+J5D+J6D) ) return

        if ( (J2D+J4D-J6D) < 0 ) return
        if ( (J4D+J6D-J2D) < 0 ) return
        if ( (J6D+J2D-J4D) < 0 ) return
        if ( (J2D+J4D+J6D)/2*2 /= (J2D+J4D+J6D) ) return

        JM1 = ( J1D + J2D + J3D ) / 2
        JM2 = ( J1D + J5D + J6D ) / 2
        JM3 = ( J4D + J2D + J6D ) / 2
        JM4 = ( J4D + J5D + J3D ) / 2
        JX1 = ( J1D + J2D + J4D + J5D ) / 2
        JX2 = ( J2D + J3D + J5D + J6D ) / 2
        JX3 = ( J3D + J1D + J6D + J4D ) / 2

        JM = JM1
        if ( JM2 > JM ) JM = JM2
        if ( JM3 > JM ) JM = JM3
        if ( JM4 > JM ) JM = JM4

        JX = JX1
        if ( JX2 < JX ) JX = JX2
        if ( JX3 < JX ) JX = JX3

        KM = JM + 1
        KX = JX + 1

        if ( KM > KX ) return

        TERM1 = FRTSJL ( J1D, J2D, J3D, J4D, J5D, J6D )

        do I1 = KM, KX
          I = I1 - 1
          TERM2 = FL(I+2) - FL(I+1-JM1) - FL(I+1-JM2) - FL(I+1-JM3)     &
     &          - FL(I+1-JM4) - FL(JX1-I+1) - FL(JX2-I+1) - FL(JX3-I+1)
          TERM = 0._pr
          TERM = exp( TERM1 + TERM2 ) * real( (-1)**I, pr )
          SIXJ = SIXJ + TERM
        enddo  ! I1

       return

      end function SIXJ
!
!============================================================================
!
      real(pr) function THRJ ( J1D, J2D, J3D, M1D, M2D, M3D )
!- - - - - - - - - - - - - - - -
!  Calculates 3jm sybol
!  All JID and MID are doubled
!- - - - - - - - - - - - - - - -
      implicit none
! .. Input / Output ..
      integer, intent(in) :: J1D, J2D, J3D, M1D, M2D, M3D
! .. Local ..
      real(pr) :: X1, X2, X3, Y1, Y2, Y3
      integer :: KMIN, KMIN1, KMIN2
      integer :: KMAX, KMAX1, KMAX2, KMAX3
      integer :: JMIN, JMAX
      real(pr) :: TERM1, TERM2, TERM
      real(pr) :: SUMM
      integer :: MSIGN
      integer :: I1, I
!- - - - - - - - - -

        THRJ = 0._pr

        X1 = real( J1D, pr ) / 2._pr
        X2 = real( J2D, pr ) / 2._pr
        X3 = real( J3D, pr ) / 2._pr
        Y1 = real( M1D, pr ) / 2._pr
        Y2 = real( M2D, pr ) / 2._pr
        Y3 = real( M3D, pr ) / 2._pr

        if ( abs(M1D) > J1D ) return
        if ( abs(M2D) > J2D ) return
        if ( abs(M3D) > J3D ) return

        if ( (J1D+J2D-J3D) < 0 ) return
        if ( (J2D+J3D-J1D) < 0 ) return
        if ( (J3D+J1D-J2D) < 0 ) return
        if ( -(M1D+M2D+M3D) /= 0 ) return

        KMIN  = ( J3D - J1D - M2D ) / 2
        KMIN1 = KMIN
        KMIN2 = ( J3D - J2D + M1D ) / 2
        if ( KMIN2 < KMIN ) KMIN = KMIN2
        KMIN = -KMIN
        if (KMIN < 0) KMIN = 0

        KMAX = int( X1 + X2 - X3 + 0.1_pr )
        KMAX1 = KMAX
        KMAX2 = int( X1 - Y1 + 0.1_pr )
        KMAX3 = int( X2 + Y2 + 0.1_pr )
        if ( KMAX2 < KMAX ) KMAX = KMAX2
        if ( KMAX3 < KMAX ) KMAX = KMAX3

        if ( KMIN > KMAX ) return

        JMIN = KMIN + 1
        JMAX = KMAX + 1
        TERM1 = FRONTL ( X1, X2, X3, Y1, Y2, Y3 )
        MSIGN = (-1)**( (J1D-J2D-M3D)/2 )
        SUMM = 0._pr

        do I1 = JMIN, JMAX
           I = I1 - 1
           TERM2 = FL(I1) + FL(KMIN1+I1) + FL(KMIN2+I1)                 &
     &           + FL(KMAX1-I+1) + FL(KMAX2-I+1) + FL(KMAX3-I+1)
           TERM = exp( TERM1 - TERM2 )
           TERM = TERM * MSIGN * ( (-1)**I )
           SUMM = SUMM + TERM
        enddo  ! I1

        THRJ = SUMM

       return

      end function THRJ
!
!============================================================================
!
      real(pr) function FRTSJL ( J1D, J2D, J3D, J4D, J5D, J6D )
!- - - - - - - - - -
      implicit none
! .. Input / Output ..
      integer, intent(in) :: J1D, J2D, J3D, J4D, J5D, J6D
!- - - - - - - - - -

        FRTSJL = DL ( J1D, J2D, J3D ) + DL ( J1D, J5D, J6D )            &
     &         + DL ( J4D, J2D, J6D ) + DL ( J4D, J5D, J3D )

       return

      end function FRTSJL
!
!============================================================================
!
      real(pr) function DL ( J1D, J2D, J3D )
!- - - - - - - - - -
      implicit none
! .. Input / Output ..
      integer, intent(in) :: J1D, J2D, J3D
! .. Local ..
      integer :: L1, L2, L3, L4
!- - - - - - - - - -

        L1 = ( J1D + J2D - J3D ) / 2
        L2 = ( J2D + J3D - J1D ) / 2
        L3 = ( J3D + J1D - J2D ) / 2
        L4 = ( J1D + J2D + J3D ) / 2 + 1

        DL = FL(L1+1) + FL(L2+1) + FL(L3+1) - FL(L4+1)
        DL = DL / 2._pr

       return

      end function DL
!
!============================================================================
!
      real(pr) function FRONTL ( X1, X2, X3, Y1, Y2, Y3 )
!- - - - - - - - - -
      implicit none
! .. Input / Output ..
      real(pr), intent(in) :: X1, X2, X3, Y1, Y2, Y3
! .. Local ..
      integer :: L1, L2, L3, L4, L5, L6, L7, L8, L9, L10
!- - - - - - - - - -

        L1  = int( X1 + X2 -X3 + 1.1_pr )
        L2  = int( X2 + X3 -X1 + 1.1_pr )
        L3  = int( X3 + X1 -X2 + 1.1_pr )
        L4  = int( X1 + X2 +X3 + 2.1_pr )
        L5  = int( X1 + Y1 + 1.1_pr )
        L6  = int( X1 - Y1 + 1.1_pr )
        L7  = int( X2 + Y2 + 1.1_pr )
        L8  = int( X2 - Y2 + 1.1_pr )
        L9  = int( X3 + Y3 + 1.1_pr )
        L10 = int( X3 - Y3 + 1.1_pr )

        FRONTL = FL(L1) + FL(L2) + FL(L3) - FL(L4) + FL(L5) + FL(L6)    &
     &         + FL(L7) + FL(L8) + FL(L9) + FL(L10)
        FRONTL = FRONTL / 2._pr

       return

      end function FRONTL
!
!============================================================================
!
      subroutine factorials_setup ( NFL )
!- - - - - - - - - -
      implicit none
! .. Input / Output ..
      integer, intent(in) :: NFL
! .. Local ..
      integer :: i, i1
!- - - - - - - - - -

        if ( allocated(FL) ) deallocate( FL )
        allocate( FL(NFL) )
        FL = 0._pr

        do i = 2, NFL, 1

          i1 = i - 1
          FL(i) = FL(i1) + log( real( i1, pr ) )

        enddo  ! i

       return

      end subroutine factorials_setup


!
!============================================================================
      end module ang
!============================================================================
!

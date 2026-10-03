!============================================================================
!
      module dirac_pnt
!- - - - - - - - - -
      use parameters
      use ang

      implicit none

      integer, parameter, private :: NNmax = 100
      integer, parameter, private :: LLmax = 20

!
!============================================================================
      contains
!============================================================================
!
      subroutine Dirac_pnt_p ( nrx, kax, Zx, p, G, F, beta_o )
!- - - - - - - - -
      implicit none
! .. Input / Output ..
      integer, intent(in) :: nrx, kax
      real(pr), intent(in) :: Zx
      real(pr), intent(in) :: p
      real(pr), intent(out) :: G, F
      real(pr), intent(in), optional :: beta_o
! .. Local ..
      real(pr) :: alphaZ, lmbd, a, eps, Norm
      real(pr) :: C1(0:NNmax), C2(0:NNmax)
      real(pr) :: D(-LLmax:NNmax,0:LLmax)
      integer :: l_p, l_m, lmax
      real(pr) :: sumg, sumf
      integer :: m
!- - - - - - - - -

        call Dirac_pnt_param_prepare ( nrx, kax, Zx, alphaZ, lmbd, a, eps, Norm, l_p, l_m )
        lmax = max( l_p, l_m )

        call C1C2_calc ( nrx, kax, alphaZ, lmbd, a, C1, C2 )

        call Dml_calc ( lmax, nrx, D(-lmax:,:), p, lmbd, a, beta_o=beta_o )

        sumg = 0._pr
        sumf = 0._pr
        do m = 0, nrx, 1
          sumg = sumg + ( C1(m) + C2(m) ) * D(m,l_p)
            !write(*,*) "oekg" , D(m,l_p)
          sumf = sumf + ( C1(m) - C2(m) ) * D(m,l_m)
         ! Write(*,*) "Jjjj" , D(m,l_m) , "l_m = " , l_m  , m
        enddo  ! m
        G = sumg * Norm * sqrt( 1._pr + eps ) * p ! * (2 * a * Zx / acl ) ** (-beta_o)
        F = sumf * Norm * sqrt( 1._pr - eps ) * p !* (2 * a * Zx / acl ) ** (-beta_o)

        if ( kax > 0 ) F = -F

       return

      end subroutine Dirac_pnt_p
!
!============================================================================
!
      subroutine dDirac_pnt_p ( nrx, kax, Zx, p, dG, dF, beta_o )

      implicit none
      integer, intent(in) :: nrx, kax
      real(pr), intent(in) :: Zx
      real(pr), intent(in) :: p
      real(pr), intent(out) :: dG, dF
      real(pr), intent(in), optional :: beta_o

      real(pr) :: alphaZ, lmbd, a, eps, Norm
      real(pr) :: C1(0:NNmax), C2(0:NNmax)
      real(pr) :: D(-LLmax:NNmax,0:LLmax)
      integer :: l_p, l_m, lmax
      real(pr) :: sumg, sumf
      real(pr) :: dG1, dF1
      real(pr) :: dG2, dF2
      real(pr) :: G, F
      integer :: m
      real(pr) :: beta

        call Dirac_pnt_p ( nrx, kax, Zx, p, G, F, beta_o=beta_o )


        call Dirac_pnt_param_prepare ( nrx, kax, Zx, alphaZ, lmbd, a, eps, Norm, l_p, l_m )
        lmax = max( l_p, l_m ) + 1

        dG1 = G * l_p / p
        dF1 = F * l_m / p

        call C1C2_calc ( nrx, kax, alphaZ, lmbd, a, C1, C2 )

        beta = 1._pr
        if ( present(beta_o) ) beta = beta + beta_o
        call Dml_calc ( lmax, nrx, D(-lmax:,:), p, lmbd, a, beta_o=beta )

        sumg = 0._pr
        sumf = 0._pr
        do m = 0, nrx, 1
          sumg = sumg + ( C1(m) + C2(m) ) * D(m,l_p+1)
          sumf = sumf + ( C1(m) - C2(m) ) * D(m,l_m+1)
!          Write(*,*)  "D(m,l_p+1) =", D(m,l_p+1)
!          write(*,*) 'D(m,l_m+1) =' , D(m,l_m+1)
        enddo  ! m
        dG2 = sumg * Norm * sqrt( 1._pr + eps ) * p
        dF2 = sumf * Norm * sqrt( 1._pr - eps ) * p

        dG = dG1 - dG2
        dF =  dF1 - dF2
        !write(*,*) " dG1 df1" , dG1 , dF1
        if ( kax > 0 ) dF = -dF

       return

      end subroutine dDirac_pnt_p
!
!============================================================================
!
subroutine dirac_pnt_XY (nrx, kax, Zx, p, X, Y, beta_o)
!- - - - - - - - -
      implicit none
! .. Input / Output ..
      integer, intent(in) :: nrx, kax
      real(pr), intent(in) :: Zx
      real(pr), intent(in) :: p
      real(pr), intent(out) :: X, Y
      real(pr), intent(in), optional :: beta_o

! .. Local ..
      real(pr) :: alphaZ, lmbd, a, eps, Norm
      real(pr) :: C1(0:NNmax), C2(0:NNmax)
      real(pr) :: D(-LLmax:NNmax,0:LLmax), D_r(-LLmax:NNmax,0:LLmax)
      integer :: l_p, l_m, lmax
      real(pr) :: term_fX, term_fY, term_gX, term_gY, fctrX1, fctrX2, fctrY1, fctrY2, factor1, factor2, beta
      integer :: m

! --- параметры ---
    call Dirac_pnt_param_prepare ( nrx, kax, Zx, alphaZ, lmbd, a, eps, Norm, l_p, l_m )
    lmax = max( l_p, l_m )

    call C1C2_calc ( nrx, kax, alphaZ, lmbd, a, C1, C2 )

    call Dml_calc ( lmax, nrx, D(-lmax:,:), p, lmbd, a, beta_o=beta_o )
    if (present(beta_o)) then
        beta= beta_o + 1.0_pr
        call Dml_calc ( lmax, nrx, D_r(-lmax:,:), p, lmbd, a, beta_o=beta )
    else
        call Dml_calc ( lmax, nrx, D_r(-lmax:,:), p, lmbd, a , 1.0_pr)
    end if

    term_fX = 0.0_pr
    term_gX = 0.0_pr
    term_fY = 0.0_pr
    term_gY = 0.0_pr

    fctrX1 = (real(2*kax, pr) - 1.0_pr + 2 * kax * eps) / 2
    fctrX2 = (1.0_pr - real(2 * kax, pr) ) / 2 * kax
    fctrY1 = (real(2*kax, pr) + 1.0_pr - 2 * kax * eps) / 2
    fctrY2 = (1.0_pr + real(2 * kax, pr) ) / 2 * kax

    do m = 0, nrx, 1
        term_fX = term_fX + ( C1(m) - C2(m) ) * (fctrX1 * D_r(m, l_p) + alphaZ * kax * D(m, l_p))
        term_gX = term_gX + ( C1(m) + C2(m) ) * D(m, l_p) * fctrX2
        term_gY = term_gY + ( C1(m) + C2(m) ) * ( fctrY1 * D_r(m, l_m) - alphaZ * kax * D(m, l_m))
        term_fY = term_fY + ( C1(m) - C2(m) ) * D(m, l_m) * fctrY2
    enddo

    factor1 = Norm * sqrt( 1._pr + eps ) * p
    factor2 = Norm * sqrt( 1._pr - eps ) * p

    term_fX = term_fX * factor2;
    term_gX = term_gX * factor1;
    term_gY = term_gY * factor1 ;
    term_fY = term_fY * factor2 ;

    X = term_fX + term_gX;
    Y = term_fY + term_gY; ! сделать умножение на b_0
    if (kax > 0) Y = -Y
      return
end subroutine dirac_pnt_XY


!============================================================================
!
      subroutine wf_real_fctr_multiply ( G, F, fctr )
!- - - - - - - - -
      implicit none
! .. Input / Output ..
      real(pr), intent(inout) :: G, F
      real(pr), intent(in) :: fctr
!- - - - - - - - -

        G = G * fctr
        F = F * fctr

       return

      end subroutine wf_real_fctr_multiply
!
!============================================================================
!
      subroutine wf_copy ( G, F, G_copy, F_copy )
!- - - - - - - - -
      implicit none
! .. Input / Output ..
      real(pr), intent(in) :: G, F
      real(pr), intent(out) :: G_copy, F_copy
!- - - - - - - - -

        G_copy = G
        F_copy = F

       return

      end subroutine wf_copy
!
!============================================================================
!
      subroutine Dirac_pnt_param_prepare ( nrx, kax, Zx, alphaZ, lmbd, a, eps, Norm, l_p, l_m )
!- - - - - - - - - -
      implicit none
! .. Input / Output ..
      integer, intent(in) :: nrx, kax
      real(pr), intent(in) :: Zx
      real(pr), intent(out) :: alphaZ, lmbd, a, eps, Norm
      integer, intent(out) :: l_p, l_m
! .. Local ..
      integer :: n
      real(pr) :: dkap, dnrx, dn
      real(pr) :: fact
      real(pr) :: gam1, gam2, gam3
      integer :: jx
!- - - - - - - - - -

        alphaZ = Zx / acl

        dkap = real( kax, pr )
        dnrx = real( nrx, pr )

        n = nrx + abs( kax )
        dn = real( n, pr )

        lmbd = sqrt( dkap**2._pr - alphaZ**2._pr )
        a = alphaZ / sqrt( (dnrx+lmbd)**2._pr + alphaZ**2._pr )

        fact = alphaZ / ( dnrx + lmbd )
        eps = 1._pr / sqrt( 1._pr + fact**2._pr )

        call GAMMAFunc ( 2._pr*lmbd+1._pr+dnrx, gam1 )
        call GAMMAFunc ( 2._pr*lmbd+1._pr     , gam2 )
        call GAMMAFunc (            1._pr+dnrx, gam3 )

        Norm = 2._pr * a**5._pr / alphaZ / ( alphaZ - dkap * a ) * gam1 / ( gam2*gam2) / gam3
        Norm = sqrt( Norm )

        call ka2jl (  kax, jx, l_p )
        l_p = l_p / 2
        call ka2jl ( -kax, jx, l_m )
        l_m = l_m / 2

       return

      end subroutine Dirac_pnt_param_prepare
!
!============================================================================
!
      subroutine C1C2_calc ( nrx, kax, alphaZ, lmbd, a, C1, C2 )
!- - - - - - - - - -
      implicit none
! .. Input / Output ..
      integer, intent(in) :: nrx, kax
      real(pr), intent(in) :: alphaZ, lmbd, a
      real(pr), intent(out) :: C1(0:), C2(0:)
! .. Local ..
      real(pr) :: dkap, dnrx
      integer :: m
      real(pr) :: dm, fact
!- - - - - - - - - -

        dkap = real( kax, pr )
        dnrx = real( nrx, pr )

        C1(0) = -dnrx
        C2(0) = alphaZ / a - dkap

        do m = 0, nrx-1, 1
          dm = real( m, pr )
          fact = ( dm + 2._pr * lmbd + 1._pr ) * ( dm + 1._pr )
          C1(m+1) = C1(m) * ( dm + 1._pr - dnrx ) / fact
          C2(m+1) = C2(m) * ( dm - dnrx ) / fact
        enddo  ! m

       return

      end subroutine
!
!============================================================================
!
      subroutine Dml_calc ( lmax, mmax, D, p, lmbd, a, beta_o )
!- - - - - - - - - -
      implicit none
! .. Input / Output ..
      integer, intent(in) :: lmax, mmax
      real(pr), intent(out) :: D(-lmax:,0:)
      real(pr), intent(in) :: p
      real(pr), intent(in) :: lmbd, a
      real(pr), intent(in), optional :: beta_o
! .. Local ..
      real(pr) :: beta
      real(pr) :: a2p2, pwr1, pwr2, gam_m, fctr
      integer :: l, m
!- - - - - - - - - -

        beta = 0._pr
        if ( present(beta_o) ) beta = beta_o

        D = 0._pr

        a2p2 = a**2._pr + p**2._pr
        do m = -lmax, mmax
          pwr1 = real(m,pr) + lmbd
          pwr2 = pwr1 + beta + 1._pr
          call GAMMAFunc ( pwr2, gam_m )
          D(m,0) = ( 2._pr*a )**pwr1 * gam_m * sin( pwr2*atan( p/a ) ) / a2p2**(pwr2/2._pr)

        enddo  ! m

        fctr = 2._pr * pi / a / p
        D(-lmax:mmax,0) = D(-lmax:mmax,0) * fctr
      !  write(*,*) 'D(m,0)' , D(-lmax:mmax,0) , mmax
        do l = 0, lmax-1
          do m = -lmax+1+l, mmax
            D(m,l+1) = ( real(l+m,pr) + lmbd + beta + 1._pr ) * 2._pr * a * D(m-1,l) - a * D(m,l)
            D(m,l+1) = D(m,l+1) / p
          !  Write(*,*) "HHH" , D
          enddo  ! m
        enddo  ! l

       return

      end subroutine Dml_calc
!
!============================================================================
!
      real(pr) function en_Dirac_pnt ( nrx, kax, Zx )
!- - - - - - - - - -
      implicit none
! .. Input / Output ..
      integer, intent(in) :: nrx, kax
      real(pr), intent(in) :: Zx
! .. Local ..
      real(pr) :: alphaZ
      real(pr) :: gamma_p
!- - - - - - - - - -

        alphaZ = Zx / acl

        gamma_p = sqrt( real(kax,pr)**2._pr - alphaZ**2._pr )

        en_Dirac_pnt = 1._pr / sqrt( 1._pr + ( alphaZ / ( gamma_p + real(nrx,pr) ) )**2._pr )

       return

      end function en_Dirac_pnt
!
!============================================================================
      end module dirac_pnt
!============================================================================
!

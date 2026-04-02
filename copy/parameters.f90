!============================================================================
!
      module parameters
!- - - - - - - - - -
      implicit none

     ! integer, parameter :: pr = selected_real_kind(15,307)
      integer, parameter :: pr = selected_real_kind(32)                   ! quadruple precision

      real(pr), parameter :: pi = 4._pr * atan(1._pr)
      real(pr), parameter :: two_pi = 8._pr * atan(1._pr)

      real(pr), parameter :: accuracy_desirable = 1.e-16_pr
      real(pr), parameter :: accuracy_acceptable = 1.5e-15_pr

!      real(pr), parameter :: acl = 137.035999679_pr                       ! CODATA 2006:   137.035 999 679(94)
!      real(pr), parameter :: acl = 137.035999074_pr                       ! CODATA 2010:   137.035 999 074(44)
!      real(pr), parameter :: acl = 137.035999139_pr                       ! CODATA 2014:   137.035 999 139(31)
      real(pr), parameter :: acl = 137.035999084_pr                       ! CODATA 2018:   137.035 999 084(21)
!      real(pr), parameter :: acl = 137.036_pr
      real(pr) , parameter :: alpha = 1._pr / acl
!   Bohr radius for electron:  a_0 = hbar / (alpha*m_e*c)
!      real(pr), parameter :: au_fm = 52917.720859_pr                      ! CODATA 2006:   52 917.720 859(36)
!      real(pr), parameter :: au_fm = 52917.721092_pr                      ! CODATA 2010:   52 917.721 092(17)
!      real(pr), parameter :: au_fm = 52917.721067_pr                      ! CODATA 2014:   52 917.721 067(12)
      real(pr), parameter :: au_fm = 52917.7210903_pr                     ! CODATA 2018:   52 917.721 0903(80)

!      real(pr), parameter :: ru_ev = 0.510998918e6_pr                     ! CODATA 2002:   0.510 998 918(44)    MeV
!      real(pr), parameter :: ru_ev = 0.510998910e6_pr                     ! CODATA 2006:   0.510 998 910(13)    MeV
!      real(pr), parameter :: ru_ev = 0.510998928e6_pr                     ! CODATA 2010:   0.510 998 928(11)    MeV
!      real(pr), parameter :: ru_ev = 0.5109989461e6_pr                    ! CODATA 2014:   0.510 998 9461(31)   MeV
      real(pr), parameter :: ru_ev = 0.51099895000e6_pr                   ! CODATA 2018:   0.510 998 950 00(15) MeV

!      real(pr), parameter :: EHartree = 27.21138386_pr ! eV               ! CODATA 2006:   27.211 383 86(68)      eV
!      real(pr), parameter :: EHartree = 27.21138505_pr ! eV               ! CODATA 2010:   27.211 385 05(60)      eV
!      real(pr), parameter :: EHartree = 27.21138602_pr ! eV               ! CODATA 2014:   27.211 386 02(17)      eV
      real(pr), parameter :: EHartree = 27.211386245988_pr ! eV           ! CODATA 2018:   27.211 386 245 988(53) eV

      real(pr), parameter :: h_si = 6.62607015e-34_pr ! J*s               ! CODATA 2018:   6.626 070 15e-34 (exact)
      real(pr), parameter :: hbar_si = h_si / two_pi
      real(pr), parameter :: SoL_si = 299792458 ! m/s                     ! CODATA 2018:   299 792 458 (exact)
      real(pr), parameter :: mass_si = 9.1093837015e-31_pr ! kg           ! CODATA 2018:   9.109 383 7015(28)e-31
!   Elementary charge times 1 V, i.e., 1 eV in J
      real(pr), parameter :: eV_J = 1.602176634e-19_pr                    ! CODATA 2018:   1.602 176 634e-19 (exact)

!      real(pr), parameter :: fm_ru = acl / au_fm
      real(pr), parameter :: fm_ru = 1.e-15_pr * ( mass_si * SoL_si ) / hbar_si         ! for CODATA 2018

!      real(pr), parameter :: ev_ru = 1._pr / ru_ev
      real(pr), parameter :: ev_ru = eV_J / ( mass_si * SoL_si**2._pr )                 ! for CODATA 2018

      complex(pr), parameter :: ci = (0._pr,1._pr)

!
!============================================================================
      end module parameters
!============================================================================
!


module angINT
use parameters
use ang

contains
!=====================================================================
    real(pr) function Int25 (k_a , k_b , mu)
    integer , intent(in) :: k_a , k_b
    real(pr), intent(in) :: mu
    integer :: j_ax2 , l_ax2 , j_bx2 , l_bx2 , mux2
    real(pr) :: ij3 , ij6 , fctr , j_a , j_b , l
    mux2 = int(mu * 2)
    j_a = abs(k_a) - 0.5_pr
    j_b = abs(k_b) - 0.5_pr
    l = abs(k_a + 0.5_pr) - 0.5_pr
    call ka2jl ( k_a, j_ax2, l_ax2 )
    call ka2jl ( k_b, j_bx2, l_bx2 )
    call factorials_setup(100)
    ij3 = THRJ(j_ax2 , 2 , j_bx2 , - mux2 , 0 , mux2)
    ij6 = SIXJ(j_ax2 , 2 , j_bx2 , 1 , l_ax2 , 1)
    fctr = (-1.0_pr)**(2*j_a - mu - 0.5_pr + l) * sqrt(6.0_pr)
    fctr = fctr * sqrt((2*j_a+1.0_pr)*(2*j_b+1.0_pr))
    int25 = fctr * ij3 * ij6

    end function Int25
!=====================================================================
    complex(pr) function Int26(k_a, k_b , mu)
    integer , intent(in) :: k_a , k_b
    real(pr), intent(in) :: mu
    integer :: j_ax2 , l_ax2 , j_bx2 , l_bx2 , mux2
    real(pr) :: ij3  , j_a , H
    complex(pr) :: i , fctr
    i = (0.0_pr, 1.0_pr)
    mux2 = int(mu * 2)
    j_a = abs(k_a) - 0.5_pr
    call ka2jl ( k_a, j_ax2, l_ax2 )
    call ka2jl ( k_b, j_bx2, l_bx2 )
    call factorials_setup(100)
    H = H_xl(k_a, k_b, 2, 2)
    ij3 = THRJ(j_ax2, 2, j_bx2, -mux2, 0, mux2)
    fctr = (-1.0_pr) ** (j_a - mu) * sqrt(6.0_pr) * i

    Int26 = fctr * ij3 * H
    end function Int26

!=====================================================================
    real(pr) function Y62(k_a , ksi)
    real(pr) , intent(in) :: ksi
    integer , intent(in) :: k_a

    select case (k_a)
    case (-1)
        Y62 = 2._pr * pi
    case (1)
        Y62 = -2._pr * pi / 3._pr * ksi
    case (-2)
        Y62 = 2._pr * pi / 3._pr * ksi
    case (2)
        Y62 = pi * (1._pr - 3._pr * ksi * ksi) / 5._pr
    case (-3)
        Y62 = -pi * (1._pr - 3._pr * ksi * ksi) / 5._pr
    case (3)
        Y62 = pi * (3._pr * ksi - 5._pr * ksi**3) / 7._pr
    case default
        write(*,*) 'Error: k_a out of range'
        stop
    end select

    end function Y62
!=====================================================================
    real(pr) function Y63(k_a , ksi)
    real(pr) , intent(in) :: ksi
    integer , intent(in) :: k_a

    select case (k_a)
    case (-1)
        Y63 = -2._pr * pi / 3._pr * ksi
    case (1)
        Y63 = -2._pr * pi / 3._pr
    case (-2)
        Y63 = pi * (1._pr - 3._pr * ksi * ksi) / 15._pr
    case (2)
        Y63 = -2._pr * pi / 15._pr * ksi
    case (-3)
        Y63 = -pi / 35._pr * (-3._pr * ksi + 5._pr * ksi**3)
    case (3)
        Y63 = -pi / 35._pr * (-1._pr + 3._pr * ksi * ksi)
    case default
        write(*,*) 'Error: k_a out of range'
        stop
    end select

    end function Y63
!=====================================================================
    real(pr) function Y64(k_a , ksi)
    real(pr) , intent(in) :: ksi
    integer , intent(in) :: k_a

    select case (k_a)
    case (-1)
        Y64 = -2._pr * pi / 3._pr
    case (1)
        Y64 = -2._pr * pi / 3._pr * ksi
    case (-2)
        Y64 = -2._pr * pi / 15._pr * ksi
    case (2)
        Y64 = pi * (1._pr - 3._pr * ksi * ksi) / 15._pr
    case (-3)
        Y64 = pi * (1._pr - 3._pr * ksi * ksi) / 35._pr
    case (3)
        Y64 = -pi / 35._pr * (-3._pr * ksi + 5._pr * ksi**3)
    case default
        write(*,*) 'Error: k_a out of range'
        stop
    end select

    end function Y64
!=====================================================================
    real(pr) function Y76(k_a , ksi , p1 , p2)
    real(pr) , intent(in) :: ksi ,  p1 , p2
    integer , intent(in) :: k_a

    select case (k_a)
    case (-1)
        Y76 = - 4 * pi / 3 * (p2 - p1 * ksi)
    case (1)
        Y76 = - 4 * pi / 3 * (p1 - p2 * ksi)
    case (-2)
        Y76 = - 4 * pi / 15 * (p1 + 2 * p2 * ksi - 3 * p1 * ksi**2)
    case (2)
        Y76 = - 4 * pi / 15 * ( p2 + 2 * p1 * ksi - 3 * p2 * ksi**2)
    case (-3)
        Y76 = 6 * pi /  35 * ( p2 - 3*p2 * ksi**2 + p1 * ksi * (-3 + 5*ksi**2))
    case (3)
        Y76 =  6 * pi /  35 * ( p1 - 3*p1 * ksi**2 + p2 * ksi * (-3 + 5*ksi**2))
    case default
        write(*,*) 'Error: k_a out of range'
        stop
    end select

    end function Y76
!=====================================================================
    real(pr) function Y77(k_a , ksi )
    real(pr) , intent(in) :: ksi
    integer , intent(in) :: k_a

    select case(k_a)
    case (-1)
        Y77 = 0.0_pr
    case (1)
        Y77 = 2 * pi / 3 * (1.0_pr - ksi**2)
    case (-2)
        Y77 =  pi / 3 * (1.0_pr - ksi**2)
    case (2)
        Y77 = 3 * pi / 5 * (1.0_pr - ksi**2) * ksi
    case (-3)
        Y77 = 2 * pi / 5 * (1.0_pr - ksi**2) * ksi
    case (3)
        Y77 = - pi / 7 * (1.0_pr - 6 * ksi**2 + 5 * ksi**4 )
    case default
        write(*,*) 'Error: k_a out of range'
        stop

    end select

    end function Y77
!=====================================================================
    real(pr) function Y26(k_a) !! [ sigma x p ]_z
    implicit none
    integer , intent(in) :: k_a

    select case (k_a)
    case (1)
        Y26 = 2.0_pr / 3
    case (-1)
        Y26 = - 2.0_pr / 3
    case (2)
        Y26 = 4.0_pr / 15
    case (-2)
        Y26 = - 4.0_pr / 15
    case (3)
        Y26 = 6.0_pr / 35
    case (-3)
        Y26 = - 6.0_pr / 35
    case default
        write(*,*) 'Error: k_a out of range'
        stop
    end select

    end function Y26
!=====================================================================

    real(pr) function Y34(k_a) !! [ p x nabla ]_z
    implicit none
    integer , intent(in) :: k_a

    select case (k_a)
    case (1)
        Y34 = 2.0_pr / 3
    case (-1)
        Y34 = 0.0_pr
    case (2)
        Y34 = 3.0_pr / 5
    case (-2)
        Y34 = 1.0_pr / 3
    case (3)
        Y34 = 4.0_pr / 7
    case (-3)
        Y34 = -2.0_pr / 5 !??????? знак
    case default
        write(*,*) 'Error: k_a out of range'
        stop
    end select

    end function Y34
!=====================================================================

    real(pr) function Y35(k_a) !! [ sigma x nabla ]_z
    implicit none
    integer , intent(in) :: k_a

    select case (k_a)
    case (1)
        Y35 = 0.0_pr
    case (-1)
        Y35 = - 4.0_pr / 3
     case (2)
        Y35 = -4.0_pr / 15
     case (-2)
        Y35 = -4.0_pr / 5
     case (3)
        Y35 = -12.0_pr / 35
     case (-3)
        Y35 = -24.0_pr / 35
     case default
        write(*,*) 'Error: k_a out of range'
        stop
     end select

    end function Y35
!=====================================================================
    real(pr) function Y25(k_a)
    implicit none
    integer , intent(in) :: k_a
    select case (k_a)
    case (1)
        Y25 = -1.0_pr / 3._pr
    case (-1)
        Y25 =  1.0_pr
     case (2)
        Y25 = -0.2_pr
     case (-2)
        Y25 = 1.0_pr / 3._pr
     case (3)
        Y25 = -1.0_pr / 7._pr
     case (-3)
        Y25 = 0.2_pr
     case default
        write(*,*) 'Error: k_a out of range'
        stop
     end select

    end function Y25

End module angINT

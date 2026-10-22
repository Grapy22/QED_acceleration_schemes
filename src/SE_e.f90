module SE_e
    use parameters
    use prep_to_SE_e
    use quantum_state
    use others
    use Integration_Cuba
    use angINT
    implicit none

contains 

    real(pr) function SE_pot0(p, state)
    real(pr), intent(in) :: p
    type(q_state) , intent(in) :: state
    real(pr) :: A, B, g, f, X, Y, wf1 , wf2, wf3
    associate( nr => state%nr, k => state%k, en => state%en, Znuc => state%Znuc )

    call Dirac_pnt_p ( nr, k, Znuc, p, g, f)
    call dirac_pnt_XY(nr, k, Znuc, p, X, Y)
    call coeff_2024(p, en, A=A, B=B)

    wf1 = g * X - f * Y
    wf2 = g * X + f * Y
    wf3 = g * Y + f * X

    SE_pot0 = A * wf1 + B * (en * wf2 + p * wf3)   
    end associate
    end function SE_pot0

    real(pr) function SE_2pot_via_0pot(p, state)
    real(pr), intent(in) :: p
    type(q_state) , intent(in) :: state
    real(pr) :: ddA, dB, ddB, g, f, X, Y, wf1 , wf2, wf3
    associate( nr => state%nr, k => state%k, en => state%en, Znuc => state%Znuc )
    call Dirac_pnt_p ( nr, k, Znuc, p, g, f, -1._pr)
    call dirac_pnt_XY(nr, k, Znuc, p, X, Y, -1._pr)
    call coeff_2024(p, en, ddA=ddA, dB=dB, ddB = ddB)

    wf1 = g * X - f * Y
    wf2 = g * X + f * Y
    wf3 = g * Y + f * X

    SE_2pot_via_0pot = ddA * wf1 + (en * ddB + 2 * dB) * wf2 + p * ddB * wf3
    end associate
    end function SE_2pot_via_0pot

    
    real(pr) function SE_2pot_via_1pot( y , p1 , p2 , ksi , state)
    real(pr), intent(in) :: y , p1 , p2 , ksi
    type(q_state) , intent(in) :: state
    real(pr) :: g , f , X2, Y2, dF1, dF2, q2, Vq
    associate( nr => state%nr, k => state%k, en => state%en, Znuc => state%Znuc , &
         alphaZ => state%alphaZ , l => state%l , l_ => state%l_)

    call Dirac_pnt_p ( nr, k, Znuc, p1, g, f, -1._pr)
    call dirac_pnt_XY(nr, k, Znuc, p2, X2, Y2, -0._pr)
    call SE_dF_i(y, en, p1, p2, ksi, g, f, X2, Y2, dF1, dF2) !второй вариант диаграммы SE_F_i(y, en, p2, p1, ksi, X2, Y2, g, f, dF1, dF2)
    
    q2 = p1*p1 + p2*p2 - 2 * p1 * p2 * ksi
    Vq = -4 * pi * alphaZ / q2

    SE_2pot_via_1pot = Vq * ( dF1 * P_l(l , ksi) + dF2 * P_l(l_ , ksi) )
        
    end associate
    end  function SE_2pot_via_1pot   

    real(pr) function SE_1pot (y , p1 , p2 , ksi , state)
    real(pr), intent(in) :: y , p1 , p2 , ksi
    type(q_state) , intent(in) :: state
    real(pr) :: g1, f1, g2, f2, Vq, FF1, FF2, q2  
    associate( nr => state%nr, k => state%k, en => state%en, Znuc => state%Znuc , &
         alphaZ => state%alphaZ , l => state%l , l_ => state%l_)

    call Dirac_pnt_p ( nr, k, Znuc, p1, g1, f2)
    call Dirac_pnt_XY ( nr, k, Znuc, p2, g2, f2)
    call SE_F_i(y, en, p1, p2, ksi, g1, f1, g2, f2, FF1, FF2)
    q2 = p1*p1 + p2*p2 - 2 * p1 * p2 * ksi
    Vq = -4 * pi * alphaZ / q2

    SE_1pot =  Vq * ( FF1 * P_l(l , ksi) + FF2 * P_l(l_ , ksi) )
    
    end associate
    end function SE_1pot

    subroutine SE_contrib (state , unit_no , contrib)
        type(q_state) , intent(in) :: state
        integer, intent(in) :: unit_no
        character(len=*), intent(in), optional :: contrib
        real(pr) :: res0 , err0 , fctr0, epsabs0_phys, epsrel, epsabs0  ! 0 pot 
        real(pr) :: res2 , err2 , fctr2, epsabs2_phys, epsabs2 ! 2pot via 0 pot
        real(pr) :: res1 , err1 , fctr1, epsabs1_phys, epsabs1 ! 1pot Sc
        real(pr) :: res3 , err3 , fctr3, epsabs3_phys, epsabs3 ! 1 pot

        associate( j => state%j, k => state%k, en => state%en , alphaZ => state%alphaZ)

        fctr0 = alpha / (4 * pi) / (2 * pi)**3 * k / (j*(j+1._pr)) * 2  ! откуда 2 ???
        fctr2 = alpha / (8 * pi) / (2 * pi)**3 * k / (j*(j+1._pr)) * alphaZ * alphaZ * 2
        fctr1 = alpha / (4 * pi) / (2 * pi)**5 * k / (j*(j+1._pr)) * alphaZ * 2
        fctr3 = alpha / (4 * pi) / (2 * pi)**5
        end associate
        epsrel = 1e-13_pr
        epsabs0_phys = 1e-14_pr
        epsabs2_phys = 1e-14_pr
        epsabs1_phys = 1e-14_pr
        epsabs3_phys = 1e-14_pr

        epsabs0 = epsabs0_phys / fctr0
        epsabs2 = epsabs2_phys / fctr2
        epsabs1 = epsabs1_phys / fctr1
        epsabs3 = epsabs3_phys / fctr3

        if (contrib == "0pot") then
            call Cuba_Int1_Cuhre(SE_pot0, state, res0, err0, epsrel, epsabs0)
            res0 = res0 * fctr0
            err0 = err0 * fctr0

            call print_result("SE: pot0", unit_no, state, res0, err0, epsabs0_phys, 1e6_pr)

            endif
        if (contrib == "2pot_via_0pot") then
            call Cuba_Int1_Cuhre(SE_2pot_via_0pot, state, res2, err2, epsrel, epsabs2)
            res2 = res2 * fctr2
            err2 = err2 * fctr2

            call print_result("SE: 2pot via 0pot", unit_no, state, res2, err2, epsabs2_phys, 1e6_pr)
            endif
        if (contrib == "2pot_via_1pot") then
            call Cuba_Int4_Cuhre(SE_2pot_via_1pot, state, res1, err1, epsrel, epsabs1)
            res1 = res1 * fctr1
            err1 = err1 * fctr1
            call print_result("SE: 2pot via 1pot", unit_no, state, res1, err1, epsabs1_phys, 1e6_pr)
            endif
        if (contrib == "1pot") then
            call Cuba_Int4_Cuhre(SE_1pot, state, res3, err3, epsrel, epsabs3)
            res3 = res3 * fctr3
            err3 = err3 * fctr3
            call print_result("SE: 1pot", unit_no, state, res3, err3, epsabs3_phys, 1e6_pr)
            endif
        end subroutine SE_contrib

end module SE_e
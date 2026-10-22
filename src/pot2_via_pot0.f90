module pot2_via_pot0
    use parameters
    use quantum_state
    use coeff_2004
    use prep_to_1pot
    use angint
    use Integration_Cuba
    use others
    implicit none
    
contains

    real(pr) function ver1_pot2_via_0pot (p, state)
    real(pr) , intent(in) :: p
    type(q_state) , intent(in) :: state
    real(pr) :: X , Y , ddA , dA  , g , f
    associate( nr => state%nr, k => state%k, en => state%en, Znuc => state%Znuc ) 
    call Dirac_pnt_p ( nr, k, Znuc, p, g, f , -1._pr)
    call diff2_Coeff(p , en , ddA = ddA)
    call diff_Coeff(p , en , dA = dA)
    X = g**2 * Y25(k) + f**2 * Y25(-k)
    Y = p * g * f * ((Y25(k) + Y25(-k)) + (Y26(k) + Y26(-k)))

    ver1_pot2_via_0pot = ddA * (en * X + Y) + 2 * dA * X
    end associate
    end function ver1_pot2_via_0pot

    real(pr) function ver2_pot2_via_0pot (p, state)
    real(pr) , intent(in) :: p
    type(q_state) , intent(in) :: state
    real(pr) :: g , f , dg , df , X , Y , Z , W , U
    real(pr) ::  db2 , ddb1 , ddb2 , ddb3
    associate( nr => state%nr, k => state%k, en => state%en, Znuc => state%Znuc ) 
    call Dirac_pnt_p ( nr, k, Znuc, p, G, F , -1.0_pr  )
    call dDirac_pnt_p (nr, k, Znuc, p, dG, dF, -1.0_pr )
    call diff_Coeff(p , en , db2 = db2 )
    call diff2_Coeff(p , en , ddb1 = ddb1 , ddb2 = ddb2 , ddb3 = ddb3 )

    X = (f * dg * Y26(-k) + g * df * Y26(k))
    Y = f * g * (Y35(k) +  Y35(-k))
    W = f * g * p * (Y34(k) + Y34(-k))
    Z = g**2 * Y34(k) + f**2 * Y34(-k)
    U = g**2 * Y34(k) - f**2 * Y34(-k)

    ver2_pot2_via_0pot = ddb1 * ( X + Y / p ) + Z * (2*db2 + en * ddb2) + ddb2 * W + ddb3 * U
    end associate
    end function ver2_pot2_via_0pot
    
    real(pr) function red_pot2_via_0pot (p, state)
    real(pr) , intent(in) :: p
    type(q_state) , intent(in) :: state
    real(pr) :: dda1 , dda2 , dda3 , C1 , C2 , C3 , g , f , da1
    associate( nr => state%nr, k => state%k, en => state%en, Znuc => state%Znuc ) 
    call Dirac_pnt_p ( nr, k, Znuc, p, g, f , -1.0_pr)
    call diff2_Coeff(p , en , dda1 = dda1 , dda2 = dda2 , dda3 = dda3 )

    C1 = g * (en * g + p * f ) + f * (en * f + p * g)
    C2 = g * g + f * f
    C3 = g * g- f * f

    red_pot2_via_0pot = dda1 * C1 + 3 * dda2 * C2 + dda3 * C3
    end associate
    end function red_pot2_via_0pot

        subroutine pot2_via_0pot (state , unit_no ,vr , vert , red)
        type(q_state) , intent(in) :: state
        character(len=*), intent(in), optional :: vr , red , vert
        integer, intent(in) :: unit_no
        real(pr) :: res1 , err1 , fctr1 , epsrel , epsabs1 , epsabs2, res2 , err2 , fctr2 
        real(pr) :: res3 , err3 , epsabs3 , fctr3 , g_D , total_res , total_err
        real(pr) ::  epsabs1_phys , epsabs2_phys ,epsabs3_phys
        logical :: flag_vert , flag_red , flag_vr 
        flag_vert = present(vert) .and. vert == "vert"
        flag_red = present(red) .and. red == "red"
        flag_vr = present(vr) .and. vr == "vr"

        !call Symb(state, unit_no)

        associate( j => state%j, k => state%k, en => state%en , alphaZ => state%alphaZ)
        g_D = k / (j*(j+1._pr)) * (en * k - 0.5_pr)

        fctr1 = 2  * alpha / pi / (2*pi)**3 * alphaZ * alphaZ / 2         ! 1/2 ddE
        fctr2 = 2 * alpha / ((2 * pi)**3 * 4 * pi) * alphaZ * alphaZ / 2
        fctr3 =  - g_D * alpha / (pi ** 4 * 32) * alphaZ * alphaZ / 2
        end associate
        epsrel = 1e-13_pr
        epsabs1_phys = 1e-14_pr
        epsabs2_phys = 1e-14_pr
        epsabs3_phys = 1e-14_pr

        epsabs1 = epsabs1_phys / fctr1
        epsabs2 = epsabs2_phys / fctr2
        epsabs3 = epsabs3_phys / abs(fctr3)


    if (flag_vert) then

        call Cuba_Int1_Cuhre(ver1_pot2_via_0pot, state, res1, err1, epsrel, epsabs1)
        call Cuba_Int1_Cuhre(ver2_pot2_via_0pot, state, res2, err2, epsrel, epsabs2)

        res1 = res1 * fctr1
        err1 = err1 * fctr1

        res2 = res2 * fctr2
        err2 = err2 * fctr2

        total_res = res1 + res2
        total_err = sqrt(err1**2 + err2**2)
        
        call print_result("pot2_via_pot0_vertex", unit_no , state, total_res, total_err, epsabs1_phys, 1e6_pr)

    endif

    if (flag_red) then

        call Cuba_Int1_Cuhre(red_pot2_via_0pot, state, res3, err3, epsrel, epsabs3)

        res3 = res3 * fctr3
        err3 = err3 * fctr3


        call print_result("pot2_via_pot0_red", unit_no , state, res3, err3, epsabs3_phys, 1e6_pr)

    endif

    if (flag_vert .and. flag_red .and. flag_vr) then

        total_res = total_res + res3
        total_err = sqrt(total_err**2 + err3**2)

        call print_result("pot2_via_pot0_vr", unit_no, state, total_res, total_err, epsabs1_phys, 1e6_pr)
         
    endif     
    !write(*,*) repeat("-" , 78)
    end subroutine pot2_via_0pot

    real(pr) function red_pot2_via_1pot (x , y , p1 , p2 , ksi , state)
        real(pr) , intent(in) :: x , y , p1 , p2 , ksi
        type(q_state) , intent(in) :: state
        real(pr) :: N , dN , F1r, F2r , dF1r , dF2r , X1 , Y1 , g2 , f2 , q2 , Vq , ddF1r , ddF2r, d3F1r , d3F2r
        real(pr) :: ddN, fctrdNN, Dfctr, A, B, C, D
        real(pr) :: fctrN , fctr , term1 , term2, term3, term4, term5, term6
        associate( nr => state%nr, k => state%k, en => state%en, Znuc => state%Znuc , &
             alphaZ => state%alphaZ , l => state%l , l_ => state%l_)
    
        q2 = p1*p1 + p2*p2 - 2 * p1 * p2 * ksi
        Vq = -4 * pi * alphaZ / q2
    
        call dirac_pnt_XY ( nr, k, Znuc, p1, X1, Y1 )
        call Dirac_pnt_p ( nr, k, Znuc, p2, g2, f2 , -1.0_pr)
        call prep (x , y, en ,p1 , p2 , ksi , N = N , dN = dN, ddN = ddN)
        call prepare_red_1pot(x , y , en, p1 , p2 , ksi , X1 , Y1 , g2 , f2 , F1r = F1r, &
        F2r = F2r, dF1r = dF1r, dF2r = dF2r , ddF1r = ddF1r , ddF2r = ddF2r, d3F1r = d3F1r, d3F2r = d3F2r)
    
        ! 1. Вспомогательные переменные для дробей и производных N
        fctrN   = 1._pr / N
        fctr    = 2._pr * (x - 1._pr) * (1._pr - 2._pr * en * en * (x - 1._pr) * fctrN)
        Dfctr   = -4._pr * (x - 1._pr)**2 * en * fctrN * (2._pr - en * fctrN * dN)
        fctrdNN = (fctrN * dN)**2 - fctrN * ddN 

        ! 2. Свёрнутые суммы исходных функций и их производных по en
        ! (A = исходная, B = 1-я, C = 2-я, D = 3-я производные)
        A = P_l(l, ksi) * F1r   + P_l(l_, ksi) * F2r
        B = P_l(l, ksi) * dF1r  + P_l(l_, ksi) * dF2r
        C = P_l(l, ksi) * ddF1r + P_l(l_, ksi) * ddF2r
        D = P_l(l, ksi) * d3F1r + P_l(l_, ksi) * d3F2r
        ! Каждая группа слагаемых отвечает за свой порядок производной исходной функции
        red_pot2_via_1pot = Vq * ( D * fctrN  - 3._pr * C * dN * fctrN**2 - 2._pr * B * ddN * fctrN**2 &
            + 5._pr * B * dN**2 * fctrN**3 - B * fctr * fctrN**2 + 2._pr * A * dN * ddN * fctrN**3 &
            - 3._pr * A * dN**3 * fctrN**4 - A * Dfctr * fctrN**2 + 2._pr * A * fctr * dN * fctrN**3 )
        end associate
    end function red_pot2_via_1pot


    ! subroutine pot2_via_1pot(state , unit_no ,vr , vert , red)
    ! type(q_state) , intent(in) :: state
    !     integer, intent(in) :: unit_no
    !     character(len=*), intent(in), optional :: vr , red , vert
    !     real(pr) :: res1 , err1 , fctr1 , epsrel , epsabs1 , epsabs2, res2 , err2 , fctr2 
    !     real(pr) :: res3 , err3 , epsabs3 , fctr3 , g_D , total_res , total_err
    !     real(pr) ::  epsabs1_phys , epsabs2_phys ,epsabs3_phys
    !     logical :: flag_vert , flag_red , flag_vr 
    ! 
    !     flag_vert = present(vert) .and. vert == "vert"
    !     flag_red = present(red) .and. red == "red"
    !     flag_vr = present(vr) .and. vr == "vr"
    ! 
    !     associate( j => state%j, k => state%k, en => state%en , alphaZ => state%alphaZ)
    !     g_D = k / (j*(j+1._pr)) * (en * k - 0.5_pr)
    ! 
    !     fctr1 = (4 * alpha ) / (2 * pi)**6 / pi * alphaZ / 2
    !     fctr2 = 4 * alpha / (2 * pi) ** 7 * alphaZ / 2
    !     fctr3 = g_D * alpha /  (2 * pi)**6 * alphaZ  * alphaZ / 2
    !     end associate
    !     epsrel = 1e-8_pr
    !     epsabs1_phys = 1e-9_pr
    !     epsabs2_phys = 1e-9_pr
    !     epsabs3_phys = 1e-9_pr
    ! 
    !     epsabs1 = epsabs1_phys / fctr1
    !     epsabs2 = epsabs2_phys / fctr2
    !     epsabs3 = epsabs3_phys / fctr3
    ! 
    ! 
    !     if (flag_vert) then
    !         write(*,*) 'Error (pot2_via_1pot): not ready yet'
    !         stop    
    !     endif
    ! 
    !     if (flag_red) then
    !         
    !         call Cuba_Int5_Cuhre(red_pot2_via_1pot , state, res3, err3, epsrel, epsabs3)
    ! 
    !         res3 = res3 * fctr3
    !         err3 = err3 * fctr3
    ! 
    !         call print_result("pot2_red_SC", unit_no ,res3, err3, epsabs3_phys, 1e6_pr)
    !         
    !     endif
    ! 
    !     if (flag_vert .and. flag_red .and. flag_vr) then
    !         write(*,*) 'Error (pot2_via_1pot): not ready yet'
    !         stop 
    !     endif     
    !     
    ! 
    ! end subroutine pot2_via_1pot


    
end module pot2_via_pot0
module pot1_2004
    use parameters
    use prep_to_1pot
    use quantum_state
    use others
    use Integration_Cuba
    use angINT
    implicit none

contains

real(pr) function ver1_pot1 (x , y , p1 , p2 , ksi , state)
    real(pr) , intent(in) :: x , y , p1 , p2 , ksi
    type(q_state) , intent(in) :: state
    real(pr) :: q2 , P_1 , P_2, P_3, P_4 , P_5, P_6 , res , Vq , f1 , g1 , f2 , g2
    associate( nr => state%nr, k => state%k, en => state%en, Znuc => state%Znuc , alphaZ => state%alphaZ)
    q2 = p1**2 + p2**2 - 2 * p1 * p2 * ksi
    Vq = -4 * pi * alphaZ / q2

    call Dirac_pnt_p ( nr, k, Znuc, p1, g1 , f1 )
    call Dirac_pnt_p ( nr, k, Znuc, p2, g2 , f2 )

    call P_i_2004 (x , y , en ,p1 , p2 , ksi , g1 , f1 , g2, f2, P_1 , P_2, P_3, P_4 , P_5, P_6)

    res = P_1 * Y62(k , ksi) + P_2 * Y62(-k , ksi) + p1 * ( P_3 * Y63(k , ksi) &
    + P_4 * Y63(-k , ksi)) + p2 * (P_5 * Y64(k , ksi) + P_6 * Y64(-k , ksi))

    ver1_pot1 = res * Vq
    end associate
    end function ver1_pot1


real(pr) function ver2_pot1 (x , y , p1 , p2 , ksi , state)
    real(pr) , intent(in) :: x , y , p1 , p2 , ksi
    type(q_state) , intent(in) :: state
    real(pr)  :: q2 , gradV_q , f1 , g1 , f2 , g2 , res , R_1 , R_2 , R_3 , R_4 , R_5 , R_6
    associate( nr => state%nr, k => state%k, en => state%en, Znuc => state%Znuc , alphaZ => state%alphaZ)
    q2 =  p1*p1 + p2*p2 - 2 * p1 * p2 * ksi
    gradV_q = 8 * pi * alphaZ / (q2 * q2 )

    call Dirac_pnt_p ( nr, k, Znuc, p1, g1 , f1 )
    call Dirac_pnt_p ( nr, k, Znuc, p2, g2 , f2 )

    call R_i_2004(x , y , en ,p1 , p2 , ksi , g1 , f1 , g2, f2 ,R_1 , R_2 , R_3 , R_4 , R_5 , R_6)

    res = R_1 * Y76(k , ksi , p1 , p2) + R_2 * Y76(-k , ksi , p1 , p2) +   &
    p1 * p2 * ( Y77(k , ksi) * (R_3 + R_4) +  Y77(-k , ksi) * (R_5 + R_6) )

    ver2_pot1 = res * gradV_q
    end associate
    end function ver2_pot1


real(pr) function red_pot1 (x , y , p1 , p2 , ksi , state)
    real(pr) , intent(in) :: x , y , p1 , p2 , ksi
    type(q_state) , intent(in) :: state
    real(pr) :: q2 , Vq , g1 , f1 , g2 , f2 , F1r, F2r, dF1r, dF2r , N , dN

    associate( nr => state%nr, k => state%k, en => state%en, Znuc => state%Znuc , &
        alphaZ => state%alphaZ , l => state%l , l_ => state%l_)
    q2 = p1 * p1 + p2 * p2 - 2 * p1 * p2 * ksi
    Vq = -4 * pi * alphaZ / q2

    call Dirac_pnt_p ( nr, k, Znuc, p1, g1, f1)
    call Dirac_pnt_p ( nr, k, Znuc, p2, g2, f2)
    call prep (x , y , en ,p1 , p2 , ksi , N = N , dN = dN)
    call prepare_red_1pot(x , y , en ,p1 , p2 , ksi , g1 , f1 , g2 , f2 , F1r = F1r, &
    F2r = F2r, dF1r = dF1r, dF2r = dF2r)

    red_pot1 = Vq * (P_l(l , ksi) * (dF1r - F1r * dN / N) + P_l(l_, ksi) * (dF2r - F2r * dN / N)) / N


    end associate
    end function red_pot1

subroutine pot1 (state , unit_no ,vr , vert , red)
    type(q_state) , intent(in) :: state
    integer, intent(in) :: unit_no
    character(len=*), intent(in), optional :: vr , red , vert
    real(pr) :: res1 , err1 , fctr1 , epsrel , epsabs1 , epsabs2, res2 , err2 , fctr2
    real(pr) :: res3 , err3 , epsabs3 , fctr3 , g_D , total_res , total_err
    real(pr) ::  epsabs1_phys , epsabs2_phys ,epsabs3_phys
    logical :: flag_vert , flag_red , flag_vr



    flag_vert = present(vert) .and. vert == "vert"
    flag_red = present(red) .and. red == "red"
    flag_vr = present(vr) .and. vr == "vr"



    associate( j => state%j, k => state%k, en => state%en , alphaZ => state%alphaZ)
    g_D = k / (j*(j+1._pr)) * (en * k - 0.5_pr)

    fctr1 = (4 * alpha ) / (2 * pi)**6 / pi
    fctr2 = 4 * alpha / (2 * pi) ** 7
    fctr3 = g_D * alpha /  (2 * pi)**6
    end associate
    epsrel = 1e-12_pr
    epsabs1_phys = 1e-14_pr
    epsabs2_phys = 1e-14_pr
    epsabs3_phys = 1e-14_pr

    epsabs1 = epsabs1_phys / fctr1
    epsabs2 = epsabs2_phys / fctr2
    epsabs3 = epsabs3_phys / fctr3


    if (flag_vert) then

        call Cuba_Int5_Cuhre(ver1_pot1, state, res1, err1, epsrel, epsabs1)

        call Cuba_Int5_Cuhre(ver2_pot1, state, res2, err2, epsrel, epsabs2)

        res1 = res1 * fctr1
        err1 = err1 * fctr1

        res2 = res2 * fctr2
        err2 = err2 * fctr2

        total_res = res1 + res2
        total_err = sqrt(err1**2 + err2**2)


        call print_result("pot1_vertex", unit_no , total_res, total_err, epsabs1_phys, 1e6_pr)

    endif

    if (flag_red) then

        call Cuba_Int5_Cuhre(red_pot1 , state, res3, err3, epsrel, epsabs3)

        res3 = res3 * fctr3
        err3 = err3 * fctr3

        call print_result("pot1_red", unit_no ,res3, err3, epsabs3_phys, 1e6_pr)

    endif

    if (flag_vert .and. flag_red .and. flag_vr) then

        total_res = total_res + res3
        total_err = sqrt(total_err**2 + err3**2)

        call print_result("pot1_vr", unit_no ,total_res, total_err, epsabs1_phys, 1e6_pr)

    endif


    end subroutine pot1

real(pr) function ver1_pot1_SC (x , y , p1 , p2 , ksi, state)
    real(pr) , intent(in) :: x , y , p1 , p2 , ksi
    type(q_state) , intent(in) :: state
    real(pr) :: Vq , q2, g1 , f1 , g2, f2 , dP_1 , dP_2, dP_3, dP_4 , dP_5, dP_6 , res
    associate( nr => state%nr, k => state%k, en => state%en, Znuc => state%Znuc , alphaZ => state%alphaZ)
    call Dirac_pnt_p ( nr, k, Znuc, p1, g1, f1 , -0.5_pr)
    call Dirac_pnt_p ( nr, k, Znuc, p2, g2, f2 , -0.5_pr)
    call dP_i_2004 (x , y , en , p1 , p2 , ksi , g1 , f1 , g2, f2, dP_1 , dP_2, dP_3, dP_4 , dP_5, dP_6)

    q2 = p1**2 + p2**2 - 2 * p1 * p2 * ksi
    Vq = -4 * pi * alphaZ / q2

    res = dP_1 * Y62(k , ksi)  + dP_2 * Y62(-k , ksi) + p1 * ( dP_3 * Y63(k , ksi) &
    + dP_4 * Y63(-k , ksi)) + p2 * (dP_5 * Y64(k , ksi) + dP_6 * Y64(-k , ksi))

    ver1_pot1_SC = res * Vq
    end associate
    end function ver1_pot1_SC

real(pr) function ver2_pot1_SC (x , y , p1 , p2 , ksi, state)
    real(pr) , intent(in) :: x , y , p1 , p2 , ksi
    type(q_state) , intent(in) :: state
    real(pr) :: res , gradV_q , q2 , g1 , f1 , g2, f2 , dR_1 , dR_2 , dR_3 , dR_4 , dR_5 , dR_6
    associate( nr => state%nr, k => state%k, en => state%en, Znuc => state%Znuc , alphaZ => state%alphaZ)
    call Dirac_pnt_p ( nr, k, Znuc, p1, g1, f1 , -0.5_pr)
    call Dirac_pnt_p ( nr, k, Znuc, p2, g2, f2 , -0.5_pr)
    call dR_i_2004(x , y , en ,p1 , p2 , ksi , g1 , f1 , g2, f2 , dR_1 , dR_2 , dR_3 , dR_4 , dR_5 , dR_6)

    q2 =  p1 * p1 + p2 * p2 - 2 * p1 * p2 * ksi
    gradV_q = 8 * pi * alphaZ / (q2 * q2)

    res = dR_1 * Y76(k , ksi , p1 , p2) + dR_2 * Y76(-k , ksi , p1 , p2) +   &
    p1 * p2 * ( Y77(k , ksi) * (dR_3 + dR_4) +  Y77(-k , ksi) * (dR_5 + dR_6) )

    ver2_pot1_SC = res * gradV_q
    end associate
    end function ver2_pot1_SC

real(pr) function red_pot1_SC(x , y , p1 , p2 , ksi , state)
    real(pr) , intent(in) :: x , y , p1 , p2 , ksi
    type(q_state) , intent(in) :: state
    real(pr) :: N , dN , F1r, F2r , dF1r , dF2r , X1 , Y1 , g2 , f2 , q2 , Vq , ddF1r , ddF2r
    real(pr) :: fctrN , fctr , term1 , term2
    associate( nr => state%nr, k => state%k, en => state%en, Znuc => state%Znuc , &
         alphaZ => state%alphaZ , l => state%l , l_ => state%l_)

    q2 = p1*p1 + p2*p2 - 2 * p1 * p2 * ksi
    Vq = -4 * pi * alphaZ / q2

    call dirac_pnt_XY ( nr, k, Znuc, p1, X1, Y1 , -0.0_pr)
    call Dirac_pnt_p ( nr, k, Znuc, p2, g2, f2 , -1.0_pr)
    call prep (x , y, en ,p1 , p2 , ksi , N = N , dN = dN)
    call prepare_red_1pot(x , y , en, p1 , p2 , ksi , X1 , Y1 , g2 , f2 , F1r = F1r, &
    F2r = F2r, dF1r = dF1r, dF2r = dF2r , ddF1r = ddF1r , ddF2r = ddF2r)

    fctrN = 1._pr / N
    fctr =  2 * (x - 1._pr) * (1._pr - 2 * en * en * (x - 1._pr) * fctrN)

    term1 = P_l(l , ksi) * (dF1r - F1r * dN *fctrN) + P_l(l_, ksi) * (dF2r - F2r * dN * fctrN)
    term1 = -term1 * dN * fctrN
    term2 = P_l(l , ksi) * (ddF1r - (dN * dF1r + F1r * fctr) * fctrN)  &
    + P_l(l_ , ksi) * (ddF2r - (dN * dF2r + F2r * fctr) * fctrN)

    red_pot1_SC = (term1 + term2) * Vq * fctrN
    end associate
    end function red_pot1_SC


subroutine pot1_SC(state , unit_no ,vr , vert , red)
type(q_state) , intent(in) :: state
    integer, intent(in) :: unit_no
    character(len=*), intent(in), optional :: vr , red , vert
    real(pr) :: res1 , err1 , fctr1 , epsrel , epsabs1 , epsabs2, res2 , err2 , fctr2
    real(pr) :: res3 , err3 , epsabs3 , fctr3 , g_D , total_res , total_err
    real(pr) ::  epsabs1_phys , epsabs2_phys ,epsabs3_phys
    logical :: flag_vert , flag_red , flag_vr

    flag_vert = present(vert) .and. vert == "vert"
    flag_red = present(red) .and. red == "red"
    flag_vr = present(vr) .and. vr == "vr"

    associate( j => state%j, k => state%k, en => state%en , alphaZ => state%alphaZ)
    g_D = k / (j*(j+1._pr)) * (en * k - 0.5_pr)

    fctr1 = (4 * alpha ) / (2 * pi)**6 / pi * alphaZ / 2
    fctr2 = 4 * alpha / (2 * pi) ** 7 * alphaZ / 2
    fctr3 = g_D * alpha /  (2 * pi)**6 * alphaZ / 2
    end associate
    epsrel = 1e-9_pr
    epsabs1_phys = 1e-10_pr
    epsabs2_phys = 1e-10_pr
    epsabs3_phys = 1e-10_pr

    epsabs1 = epsabs1_phys / fctr1
    epsabs2 = epsabs2_phys / fctr2
    epsabs3 = epsabs3_phys / fctr3


    if (flag_vert) then

        call Cuba_Int5_Cuhre(ver1_pot1_SC, state, res1, err1, epsrel, epsabs1)

        call Cuba_Int5_Cuhre(ver2_pot1_SC, state, res2, err2, epsrel, epsabs2)

        res1 = res1 * fctr1
        err1 = err1 * fctr1

        res2 = res2 * fctr2
        err2 = err2 * fctr2

        total_res = res1 + res2
        total_err = sqrt(err1**2 + err2**2)


        call print_result("pot1_vertex_SC", unit_no , total_res, total_err, epsabs1_phys, 1e6_pr)

    endif

    if (flag_red) then

        call Cuba_Int5_Cuhre(red_pot1_SC , state, res3, err3, epsrel, epsabs3)

        res3 = res3 * fctr3
        err3 = err3 * fctr3

        call print_result("pot1_red_SC", unit_no ,res3, err3, epsabs3_phys, 1e6_pr)

    endif

    if (flag_vert .and. flag_red .and. flag_vr) then

        total_res = total_res + res3
        total_err = sqrt(total_err**2 + err3**2)

        call print_result("pot1_vr_SC", unit_no ,total_res, total_err, epsabs1_phys, 1e6_pr)

    endif


end subroutine pot1_Sc

end module pot1_2004

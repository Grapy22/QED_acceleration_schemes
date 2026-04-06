module pot0_2004
    use parameters
    use quantum_state
    use coeff_2004
    use angint
    use Integration_Cuba
    use others
    
    implicit none

    
    
contains


    real(pr) function ver1_0pot (p , state)
    real(pr) , intent(in) :: p
    type(q_state) , intent(in) :: state
    real(pr) :: A , g , f 

    associate( nr => state%nr, k => state%k, en => state%en, Znuc => state%Znuc ) 
    call Dirac_pnt_p ( nr,  k, Znuc, p, g, f)
    call Coeff(p , en , A = A )

    ver1_0pot = A * ( en * (g*g * Y25(k) + f*f * Y25(-k) )  &
       + p * g * f * (Y25(k) + Y25(-k) - (Y26(k) + Y26(-k)) ) )
    end associate
    end function ver1_0pot
!=====================================================================

    real(pr) function ver2_0pot (p , state)
    real(pr) , intent(in) :: p
    type(q_state) , intent(in) :: state
    real(pr) :: b1 , b2 , b3 , g , f , X ,  Y , Z , df , dg

    associate( nr => state%nr, k => state%k, en => state%en, Znuc => state%Znuc )
    call Dirac_pnt_p ( nr, k, Znuc, p, g, f)
    call Coeff(p , en , b1 = b1 , b2 = b2 , b3 = b3 )
    call dDirac_pnt_p ( nr, k, Znuc, p, dg , df )

    X = b1 * (f * dg * Y26(-k) + g * df * Y26(k))
    Y = g * f * (b1 / p * (Y35(k) + Y35(-k)) + b2 * p * (Y34(k) + Y34(-k)) )
    Z = en * b2 * (g*g * Y34(k) + f*f * Y34(-k)) + b3 * (g*g * Y34(k) - f**2 * Y34(-k))
    ver2_0pot = X + Y + Z
    end associate
    end function ver2_0pot
!=====================================================================

    real(pr) function red_0pot (p , state)
    real(pr) , intent(in) :: p
    type(q_state) , intent(in) :: state
    real(pr) :: a1 , a2 , a3 , term1 , term2 , term3 , g , f

    associate( nr => state%nr, k => state%k, en => state%en, Znuc => state%Znuc )
    call Dirac_pnt_p ( nr, k, Znuc, p, g, f)
    call Coeff(p , en , a1 = a1 , a2 = a2 , a3 = a3 )

    term1 = a1 * (en * (g*g + f*f) + 2 * p * g * f)
    term2 = a2 * (g*g + f*f)
    term3 = a3 * (g*g - f*f)

    red_0pot = term1 + term2 + term3
    end associate
    end function red_0pot
!=====================================================================

    subroutine pot0 (state , unit_no ,vr , vert , red)
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

        associate( j => state%j, k => state%k, en => state%en)
        g_D = k / (j*(j+1._pr)) * (en * k - 0.5_pr)

        fctr1 = 2._pr * alpha / pi / (2 * pi)**3
        fctr2 = 2._pr * alpha / (2 * pi)**3 / ( 4 * pi)
        fctr3 = - g_D * alpha / (pi ** 4 * 32)
        end associate
        epsrel = 1e-13_pr
        epsabs1_phys = 1e-15_pr
        epsabs2_phys = 1e-15_pr
        epsabs3_phys = 1e-15_pr

        epsabs1 = epsabs1_phys / fctr1
        epsabs2 = epsabs2_phys / fctr2
        epsabs3 = epsabs3_phys / abs(fctr3)


    if (flag_vert) then

        call Cuba_Int1_Cuhre(ver1_0pot, state, res1, err1, epsrel, epsabs1)
        call Cuba_Int1_Cuhre(ver2_0pot, state, res2, err2, epsrel, epsabs2)

        res1 = res1 * fctr1
        err1 = err1 * fctr1

        res2 = res2 * fctr2
        err2 = err2 * fctr2

        total_res = res1 + res2
        total_err = sqrt(err1**2 + err2**2)
        
        call print_result("pot0_vertex", unit_no ,total_res, total_err, epsabs1_phys, 1e6_pr)

    endif

    if (flag_red) then

        call Cuba_Int1_Cuhre(red_0pot, state, res3, err3, epsrel, epsabs3)

        res3 = res3 * fctr3
        err3 = err3 * fctr3


        call print_result("pot0_red", unit_no ,res3, err3, epsabs3_phys, 1e6_pr)

    endif

    if (flag_vert .and. flag_red .and. flag_vr) then

        total_res = total_res + res3
        total_err = sqrt(total_err**2 + err3**2)

        call print_result("pot0_vr", unit_no, total_res, total_err, epsabs1_phys, 1e6_pr)
         
    endif     
    !write(*,*) repeat("-" , 78)
    end subroutine pot0
    
!=====================================================================
    real(pr) function ver1_0pot_SC(p , state)
        real(pr) , intent(in) ::  p
        type(q_state) , intent(in) :: state
        real(pr) :: g , f , A , dA , X , Y
        associate( nr => state%nr, k => state%k, en => state%en, Znuc => state%Znuc ) 
        call Dirac_pnt_p ( nr, k, Znuc, p, g, f , -0.5_pr)
        call Coeff(p , en , A = A)
        call diff_Coeff(p , en , dA = dA)

        X = g**2 * Y25(k) + f**2 * Y25(-k)
        Y = p * g * f * (Y25(k) + Y25(-k) - (Y26(k) + Y26(-k)))

        ver1_0pot_SC = dA * (en * X + Y) + A * X
        end associate
    end function ver1_0pot_SC

    real(pr) function ver2_0pot_SC (p , state)
        real(pr) , intent(in) :: p
        type(q_state) , intent(in) :: state
        real(pr) :: g , f , db1 , db2 , db3 , b2 , X , Y , V, W , U , dg , df
        associate( nr => state%nr, k => state%k, en => state%en, Znuc => state%Znuc )
        call Dirac_pnt_p ( nr, k, Znuc, p, g, f , -0.5_pr)
        call dDirac_pnt_p ( nr, k, Znuc, p, dg , df , -0.5_pr)
        call diff_Coeff(p , en , db1 = db1 , db2 = db2 , db3 = db3)
        call Coeff(p , en , b2 = b2)

        X = (f * dg * Y26(-k) + g * df * Y26(k))
        Y = f * g * (Y35(k) + Y35(-k))
        W = f * g * p * (Y34(k) + Y34(-k))
        V = g**2 * Y34(k) + f**2 * Y34(-k)
        U = g**2 * Y34(k) - f**2 * Y34(-k)

        ver2_0pot_SC = db1 * ( X + Y / p ) + V * (b2 + en * db2) + db2 * W + db3 * U
        end associate
    end function ver2_0pot_SC

!=====================================================================
    real(pr) function red_0pot_SC (p ,state)
        real(pr) , intent(in) :: p
        type(q_state) , intent(in) :: state
        real(pr) :: C1 , C2 , C3 , da1 , da2 , da3 , g , f
        associate( nr => state%nr, k => state%k, en => state%en, Znuc => state%Znuc )
        call Dirac_pnt_p ( nr, k, Znuc, p, g, f , -0.5_pr)
        call diff_Coeff(p , en , da1 = da1 , da2 = da2 , da3 = da3)

        C1 = g * (en * g + p * f ) + f * (en * f + p * g)
        C2 = g*g + f*f
        C3 = g*g - f*f

        red_0pot_SC = da1 * C1 + 2 * da2 * C2 + da3 * C3
        end associate
    end function red_0pot_SC


    subroutine pot0_SC (state , unit_no ,vr , vert , red)
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

        !call Symb(state, unit_no)

        associate( j => state%j, k => state%k, en => state%en , alphaZ => state%alphaZ)
        g_D = k / (j*(j+1._pr)) * (en * k - 0.5_pr)

        fctr1 = 2._pr * alpha / pi / (2 * pi)**3 * alphaZ
        fctr2 = 2._pr * alpha / (2 * pi)**3 / ( 4 * pi) * alphaZ
        fctr3 = - g_D * alpha / (pi ** 4 * 32) * alphaZ
        end associate
        epsrel = 1e-13_pr
        epsabs1_phys = 1e-15_pr
        epsabs2_phys = 1e-15_pr
        epsabs3_phys = 1e-15_pr

        epsabs1 = epsabs1_phys / fctr1
        epsabs2 = epsabs2_phys / fctr2
        epsabs3 = epsabs3_phys / abs(fctr3)


    if (flag_vert) then
        call Cuba_Int1_Cuhre(ver1_0pot_SC, state, res1, err1, epsrel, epsabs1)
        call Cuba_Int1_Cuhre(ver2_0pot_SC, state, res2, err2, epsrel, epsabs2)

        res1 = res1 * fctr1
        err1 = err1 * fctr1

        res2 = res2 * fctr2
        err2 = err2 * fctr2

        total_res = res1 + res2
        total_err = sqrt(err1**2 + err2**2)
        
        call print_result("pot0_vertex_SC", unit_no,total_res, total_err, epsabs1_phys, 1e6_pr)

    endif

    if (flag_red) then

        call Cuba_Int1_Cuhre(red_0pot_SC, state, res3, err3, epsrel, epsabs3)

        res3 = res3 * fctr3
        err3 = err3 * fctr3


        call print_result("pot0_red_SC", unit_no,res3, err3, epsabs3_phys, 1e6_pr)

    endif

    if (flag_vert .and. flag_red .and. flag_vr) then

        total_res = total_res + res3
        total_err = sqrt(total_err**2 + err3**2)

        call print_result("pot0_vr_SC", unit_no ,total_res, total_err, epsabs1_phys, 1e6_pr)
         
    endif     
   ! write(*,*) repeat("-" , 78)
    end subroutine pot0_SC
end module pot0_2004
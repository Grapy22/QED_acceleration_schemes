module prep_to_SE_e
    use parameters
    use coeff_2004

    implicit none
    
contains 

    subroutine coeff_2024 (p , en , A , B, ddA, dB, ddB)
    real(pr) , intent(in) :: p , en
    real(pr) , intent(out), optional :: A, B, ddA, dB, ddB
    real(pr) :: ro , log_ro , x , eps
    real(pr) :: ddAro, dAro, b1, db1, ddb1

    eps = 5.e-2_pr
    x =  p*p - en*en
    ro = 1._pr + x
    log_ro = Log(ro)

    if (present(A)) then
        if (abs(x) > eps) then
            A = 2 * ( 1._pr + 2 * ro * Log_ro / (1._pr - ro) )
        else
            A = -2._pr - 2._pr*x + (2._pr/3._pr)*x**2 - (1._pr/3._pr)*x**3 + (1._pr/5._pr)*x**4 &
                - (2._pr/15._pr)*x**5 + (2._pr/21._pr)*x**6 - (1._pr/14._pr)*x**7 + (1._pr/18._pr)*x**8 &
                - (2._pr/45._pr)*x**9 + (2._pr/55._pr)*x**10 - (1._pr/33._pr)*x**11 + (1._pr/39._pr)*x**12 &
                - (2._pr/91._pr)*x**13 + (2._pr/105._pr)*x**14 - (1._pr/60._pr)*x**15
            endif
        endif
    if (present(B)) then
        call Coeff(p, en, b1=b1)
        B = -b1
        endif 

    if (present(ddA)) then
        if (abs(x) > eps) then
        ddA = - 8._pr / (1._pr - ro) * ( 1._pr + 1._pr / (1._pr - ro) * ( log_ro - 2 * en * en * ( (1._pr + ro) / ro &
        + 2 * log_ro / (1._pr - ro) )) )
        else
        dAro = -2._pr + (4._pr/3._pr)*x - x**2 + (4._pr/5._pr)*x**3 - (2._pr/3._pr)*x**4 &
            + (4._pr/7._pr)*x**5 - (1._pr/2._pr)*x**6 + (4._pr/9._pr)*x**7 - (2._pr/5._pr)*x**8 &
            + (4._pr/11._pr)*x**9 - (1._pr/3._pr)*x**10 + (4._pr/13._pr)*x**11 - (2._pr/7._pr)*x**12 &
            + (4._pr/15._pr)*x**13 - (1._pr/4._pr)*x**14 + (4._pr/17._pr)*x**15
        ddAro =  (4._pr/3._pr) - 2._pr*x + (12._pr/5._pr)*x**2 - (8._pr/3._pr)*x**3 + (20._pr/7._pr)*x**4 &
           - 3._pr*x**5 + (28._pr/9._pr)*x**6 - (16._pr/5._pr)*x**7 + (36._pr/11._pr)*x**8 &
           - (10._pr/3._pr)*x**9 + (44._pr/13._pr)*x**10 - (24._pr/7._pr)*x**11 + (52._pr/15._pr)*x**12 &
           - (7._pr/2._pr)*x**13 + (60._pr/17._pr)*x**14 - (32._pr/9._pr)*x**15
       ddA = ddAro * (2 * en)**2 - 2 * dAro         
        endif         
    endif

    if (present(dB)) then
        call diff_Coeff(p, en, db1 = db1)
        dB = -db1
    endif

    if (present(ddB)) then
        call diff2_Coeff(p, en, ddb1=ddb1)
        ddB = -ddb1
        endif
    end subroutine coeff_2024


    subroutine uvlam_eq(y, en , p1 , p2 , ksi , u , v, lam)
    real(pr) , intent(in) :: y , en , p1 , p2 , ksi
    real(pr) , intent(out) ::  u , v , lam
    real(pr) :: en2 , pp_1 , pp_2

    en2 = en * en
    pp_1 = en2 - p1 * p1
    pp_2 = en2 - p2 * p2

    u = y * y * pp_1 + pp_2 * (1.0_pr - y) ** 2 + 2* y * (1.0_pr - y) * (en2 - p1 * p2 * ksi)

    v = 1.0_pr - y * pp_1 - (1.0_pr - y) * pp_2

    lam = u / v

    end subroutine uvlam_eq

    subroutine dudvdlam_eq(y, en , p1 , p2 , ksi , du , dv, dlam)
    real(pr) , intent(in) :: y, en , p1 , p2 , ksi
    real(pr) , intent(out) :: du , dv, dlam
    real(pr) ::  u , v, lam
    call uvlam_eq(y, en , p1 , p2 , ksi , u , v, lam)
    du = 2 * en * y
    dv = -2 * en * y
    dlam =  (du - lam * dv) / v
    end subroutine dudvdlam_eq

    subroutine C_i_eq(y, en, p1, p2, ksi, C0, C11, C12)
        implicit none
        real(pr), intent(in) :: y, en, p1, p2, ksi
        real(pr), intent(out), optional :: C0, C11, C12
        real(pr) :: u, v, lam, Loga, fctr
    
        if (.not. (present(C0) .or. present(C11) .or. present(C12))) return
    
        call uvlam_eq(y, en, p1, p2, ksi, u, v, lam)
        Loga = log(1._pr + lam)
    
        if (present(C0)) then
            C0 = - Loga / u
        end if
    
        if (present(C11) .or. present(C12)) then
            fctr = (1._pr - Loga / lam) / u
    
            if (present(C11)) C11 = fctr * y
            if (present(C12)) C12 = fctr * (1._pr - y)
        end if
    
    end subroutine C_i_eq

    subroutine C_2i_eq(y, en, p1, p2, ksi, C21, C22, C23, C24)
        implicit none
        real(pr), intent(in) :: y, en, p1, p2, ksi
        real(pr), intent(out), optional :: C21, C22, C23, C24
        real(pr) :: u, v, lam, Loga, fctr
    
        if (.not. (present(C21) .or. present(C22) .or. &
                   present(C23) .or. present(C24))) return
    
        call uvlam_eq(y, en, p1, p2, ksi, u, v, lam)
        Loga = log(1._pr + lam)
    
        if (present(C21) .or. present(C22) .or. present(C23)) then
            fctr = (-0.5_pr + 1._pr / lam - Loga / (lam * lam)) / u
    
            if (present(C21)) C21 = y * y * fctr
            if (present(C22)) C22 = (1._pr - y)**2 * fctr
            if (present(C23)) C23 = y * (1._pr - y) * fctr
        end if
    
        if (present(C24)) C24 = 0._pr
    
    end subroutine C_2i_eq

    subroutine diff_Ci_eq(y, en, p1, p2, ksi, dC0, dC11, dC12)
        implicit none
        real(pr), intent(in) :: y, en, p1, p2, ksi
        real(pr), intent(out), optional :: dC0, dC11, dC12
        real(pr) :: u, v, lam, du, dv, dlam, fctr, Loga
    
        if (.not. (present(dC0) .or. present(dC11) .or. &
                   present(dC12))) return
    
        call uvlam_eq(y, en, p1, p2, ksi, u, v, lam)
        call dudvdlam_eq(y, en, p1, p2, ksi, du, dv, dlam)
    
        Loga = log(1._pr + lam)
    
        if (present(dC0)) then
            dC0 = (dv / v * Loga + dlam / lam * &
                  (Loga - lam / (1._pr + lam))) / u
        end if
    
        if (present(dC11) .or. present(dC12)) then
            fctr = (dv / v * (Loga / lam - 1._pr) + &
                   dlam / lam * (2 * Loga / lam - &
                   (2._pr + lam) / (1._pr + lam))) / u
    
            if (present(dC11)) dC11 = fctr * y
            if (present(dC12)) dC12 = fctr * (1._pr - y)
        end if
    
    end subroutine diff_Ci_eq

    subroutine diff_C2i_eq(y, en, p1, p2, ksi, dC21, dC22, dC23, dC24)
        implicit none
        real(pr), intent(in) :: y, en, p1, p2, ksi
        real(pr), intent(out), optional :: dC21, dC22, dC23, dC24
        real(pr) :: u, v, lam, Loga, fctr, du, dv, dlam, lam2
    
        if (.not. (present(dC21) .or. present(dC22) .or. &
                   present(dC23) .or. present(dC24))) return
    
        call uvlam_eq(y, en, p1, p2, ksi, u, v, lam)
    
        if (present(dC21) .or. present(dC22) .or. &
            present(dC23)) then
    
            call dudvdlam_eq(y, en, p1, p2, ksi, du, dv, dlam)
    
            lam2 = lam * lam
            Loga = log(1._pr + lam)
    
            fctr = (dv / v * (Loga / lam2 - 1._pr / lam + 0.5_pr) + &
                   dlam / lam * (3 * Loga / lam2 + 0.5_pr - &
                   (3._pr + 2 * lam) / (lam * (1._pr + lam)))) / u
    
            if (present(dC21)) dC21 = fctr * y * y
            if (present(dC22)) dC22 = fctr * (1._pr - y)**2
            if (present(dC23)) dC23 = fctr * y * (1._pr - y)
        end if
    
        if (present(dC24)) dC24 = 0._pr
    
    end subroutine diff_C2i_eq

    subroutine ABCD_eq(y, en, p1, p2, ksi, A, B, C, D)
        implicit none
    
        real(pr), intent(in) :: y, en, p1, p2, ksi
        real(pr), intent(out), optional :: A, B, C, D
    
        real(pr) :: C21, C22, C23, C24, C0, C11, C12
        real(pr) :: pp1, pp2
    
        ! Если результаты не запрошены
        if (.not. (present(A) .or. present(B) .or. &
                   present(C) .or. present(D))) return
    
        pp1 = en*en - p1*p1
        pp2 = en*en - p2*p2
    
        ! Коэффициенты C_i нужны для всех A, B, C, D
        call C_i_eq(y, en, p1, p2, ksi, C0=C0, C11=C11, C12=C12)
    
        ! A
        if (present(A)) then
            call C_2i_eq(y, en, p1, p2, ksi, C24=C24)
    
            A = C24 - 2._pr + pp1*C11 + pp2*C12 + &
                4*(en*en - p1*p2*ksi)*(C0 + C11 + C12) + &
                (-2*C0 + C11 + C12)
        end if
    
        ! B
        if (present(B)) then
            call C_2i_eq(y, en, p1, p2, ksi, C21=C21, C23=C23)
    
            B = -4*(C0 + 2*C11 + C12 + C21 + C23)
        end if
    
        ! C
        if (present(C)) then
            call C_2i_eq(y, en, p1, p2, ksi, C22=C22, C23=C23)
    
            C = -4*(C0 + C11 + 2*C12 + C22 + C23)
        end if
    
        ! D
        if (present(D)) then
            D = 2*(C0 + C11 + C12)
        end if
    
    end subroutine ABCD_eq
    
    subroutine diff_ABCD_eq(y, en, p1, p2, ksi, dA, dB, dC, dD)
        implicit none
    
        real(pr), intent(in) :: y, en, p1, p2, ksi
        real(pr), intent(out), optional :: dA, dB, dC, dD
    
        real(pr) :: dC0, dC11, dC12
        real(pr) :: C0, C11, C12
        real(pr) :: dC21, dC22, dC23, dC24
    
        ! Проверяем, запрошен ли хотя бы один результат
        if (.not. (present(dA) .or. present(dB) .or. &
                   present(dC) .or. present(dD))) return
    
        ! Общие коэффициенты: нужны для всех производных
        call diff_Ci_eq(y, en, p1, p2, ksi, &
                        dC0=dC0, dC11=dC11, dC12=dC12)
    
        ! Производная A
        if (present(dA)) then
            call C_i_eq(y, en, p1, p2, ksi, &
                        C0=C0, C11=C11, C12=C12)
    
            call diff_C2i_eq(y, en, p1, p2, ksi, dC24=dC24)
    
            dA = dC24 + 2*en*(C11 + 2*(C0 + C11 + C12)) &
               + (en*en - p1*p1)*dC11 &
               + (en*en - p2*p2)*dC12 &
               + 4*(en*en - p1*p2*ksi)*(dC0 + dC11 + dC12) &
               + (-2*dC0 + dC11 + dC12)
        end if
    
        ! Производная B
        if (present(dB)) then
            call diff_C2i_eq(y, en, p1, p2, ksi, &
                             dC21=dC21, dC23=dC23)
    
            dB = -4*(dC0 + 2*dC11 + dC12 + dC21 + dC23)
        end if
    
        ! Производная C
        if (present(dC)) then
            call diff_C2i_eq(y, en, p1, p2, ksi, &
                             dC22=dC22, dC23=dC23)
    
            dC = -4*(dC0 + dC11 + 2*dC12 + dC22 + dC23)
        end if
    
        ! Производная D
        if (present(dD)) then
            dD = 2*(dC0 + dC11 + dC12)
        end if
    
    end subroutine diff_ABCD_eq

   subroutine SE_dF_i (y, en , p1 , p2 , ksi, g1, f1, g2, f2 ,dF1 , dF2)
   real(pr), intent(in) :: y, en , p1 , p2 , ksi, g1, f1, g2, f2 
   real(pr), intent(out) :: dF1 , dF2
   real(pr) :: dA, dB , dC , dD , C1, B1, B, D, H1, C21, C11, C0, C23, C12

   call C_i_eq(y, en, p1, p2, ksi, C0=C0, C11=C11, C12 = C12)
   call C_2i_eq(y, en, p1, p2, ksi, C21=C21, C23 = C23)
   H1 = 4 * (C0 + 2 * C11)
   B1 = -4 * (C11 + C21)
   C1 = -4 * (C0 + C11 + C12 + C23)
   call ABCD_eq(y, en, p1, p2, ksi, B=B, D=D)
   call diff_ABCD_eq(y, en, p1, p2, ksi, dA=dA, dB=dB, dD=dD, dC = dC)
   dF1 = (dA + en * (B + 4 * dD) + H1) * g1 * g2 + (B1 + en * dB) * (en * g1 + p1 * f1) * g2 &
   + (C1 + en * dC + D) * g1 * (en * g2 + p2 * f2) + dD * (en * g1 + p1 * f1) * (en * g2 + p2 * f2)

   dF2 = (dA + en * (B - 4 * dD) - H1) * f1 * f2 + (B1 + en * dB) * (en * f1 + p1 * g1) * f2 &
   + (C1 + en * dC + D) * f1 * (en * f2 + p2 * g2) + dD * (en * f1 + p1 * g1) * (en * f2 + p2 * g2)
   end subroutine SE_dF_i

   subroutine SE_F_i(y, en, p1, p2, ksi, g1, f1, g2, f2, FF1, FF2)
   
      real(pr), intent(in)  :: y, en, p1, p2, ksi
      real(pr), intent(in)  :: g1, f1, g2, f2
      real(pr), intent(out) :: FF1, FF2
   
      real(pr) :: A, B, D, C   
      ! A is obtained from the existing coefficient routine.
      call ABCD_eq(y, en, p1, p2, ksi, A=A, D=D, B=B, C = C)
   
      ! Eq. (C22)
      FF1 = (A + 4 * en * D)*g1*g2 + en*B*(en*g1 + p1*f1)*g2 + en*C*g1*(en*g2 + p2*f2) + D*(en*g1 + p1*f1)*(en*g2 + p2*f2)  
      ! Eq. (C23)
      FF2 = (A - 4* en* D)*f1*f2 + en*B*(en*f1 + p1*g1)*f2 + en*C*f1*(en*f2 + p2*g2) + D*(en*f1 + p1*g1)*(en*f2 + p2*g2)
   
   end subroutine SE_F_i
end module prep_to_SE_e
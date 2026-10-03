module coeff_2004
    use parameters
    implicit none

    contains

    subroutine Coeff(p , en , A , b1 , b2 , b3 , a1, a2 ,a3)
    real(pr) , intent(in) :: p , en
    real(pr) , intent(out) , optional :: A , b1 , b2 , b3 , a1, a2 ,a3
    real(pr) :: ro , log_ro , x , eps 

    eps = 5.e-2_pr
    x =  p*p - en*en
    ro = 1._pr + x
    log_ro = Log(ro)

    if (present(A)) then
        if (abs(x) > eps) then
            A = 1._pr / (1._pr - ro) * (1._pr + log_ro / (1._pr - ro))
        else
            A = -(1._pr/2._pr) + x**1._pr/3._pr - x**2._pr/4._pr + x**3._pr/5._pr - x**4._pr/6._pr + x**5._pr/7._pr &
             - x**6._pr/8._pr + x**7._pr/9._pr - x**8._pr/10._pr + x**9._pr/11._pr - x**10._pr/12._pr &
              + x**11._pr/13._pr - x**12._pr/14._pr + x**13._pr/15._pr - x**14._pr/16._pr + x**15._pr/17._pr
        end if
    end if

    if (present(b1)) then
        if (abs(x) > eps) then
            b1 = ( 2._pr - ro ) / ( 1._pr - ro ) * ( 1._pr + ro * log_ro / ( 1._pr - ro ) )
        else
            b1 = 1._pr/2._pr - (2._pr * x)/3._pr + x**2._pr/4._pr - (2._pr * x**3._pr)/15._pr + &
            x**4._pr/12._pr - (2._pr * x**5._pr)/35._pr + x**6._pr/24._pr - (2._pr * x**7._pr)/63._pr + &
            x**8._pr/40._pr - (2._pr * x**9._pr)/99._pr + x**10._pr/60._pr - (2._pr * x**11._pr)/143._pr + &
            x**12._pr/84._pr - (2._pr * x**13._pr)/195._pr + x**14._pr/112._pr - (2._pr * x**15._pr)/255._pr
        end if
    end if

    if (present(b2)) then
        if (abs(x) > eps) then
            b2 = -2._pr / ( 1._pr - ro )**2._pr * ( 3._pr - ro + 2._pr * log_ro / ( 1._pr - ro ) )
        else
            b2 = 4._pr/3._pr - x**1._pr + (4._pr * x**2._pr)/5._pr - (2._pr * x**3._pr)/3._pr + &
            (4._pr * x**4._pr)/7._pr - x**5._pr/2._pr + (4._pr * x**6._pr)/9._pr - (2._pr * x**7._pr)/5._pr + &
            (4._pr * x**8._pr)/11._pr - x**9._pr/3._pr + (4._pr * x**10._pr)/13._pr - (2._pr * x**11._pr)/7._pr + &
            (4._pr * x**12._pr)/15._pr - x**13._pr/4._pr + (4._pr * x**14._pr)/17._pr - (2._pr * x**15._pr)/9._pr
        end if
    end if

    if (present(b3)) then
        if (abs(x) > eps) then
            b3 = 8._pr / ( 1._pr - ro ) * ( 1._pr + log_ro / ( 1._pr - ro ) )
        else
            b3 = -4._pr + (8._pr * x)/3._pr - 2._pr * x**2._pr + (8._pr * x**3._pr)/5._pr - &
            (4._pr * x**4._pr)/3._pr + (8._pr * x**5._pr)/7._pr - x**6._pr + (8._pr * x**7._pr)/9._pr - &
            (4._pr * x**8._pr)/5._pr + (8._pr * x**9._pr)/11._pr - (2._pr * x**10._pr)/3._pr + &
            (8._pr * x**11._pr)/13._pr - (4._pr * x**12._pr)/7._pr + (8._pr * x**13._pr)/15._pr - &
            x**14._pr/2._pr + (8._pr * x**15._pr)/17._pr
        end if
    end if

    if (present(a1)) then
        if (abs(x) > eps) then
            a1 = -2._pr * en / ( 1._pr - ro )**2._pr * ( 3._pr - ro + 2._pr * log_ro / ( 1._pr - ro ) )
        else
            a1 = 4._pr/3._pr - x + 4._pr*x**2._pr/5._pr - 2._pr*x**3._pr/3._pr + 4._pr*x**4._pr/7._pr - &
            x**5._pr/2._pr + 4._pr*x**6._pr/9._pr - 2._pr*x**7._pr/5._pr + 4._pr*x**8._pr/11._pr - &
            x**9._pr/3._pr + 4._pr*x**10._pr/13._pr - 2._pr*x**11._pr/7._pr + 4._pr*x**12._pr/15._pr - &
            x**13._pr/4._pr + 4._pr*x**14._pr/17._pr - 2._pr*x**15._pr/9._pr
            a1 = a1 * en
        end if
    end if

    if (present(a2)) then
        if (abs(x) > eps) then
            a2 = 2._pr + ro / ( 1._pr - ro ) * ( 1._pr + ( 2._pr - ro ) * log_ro / ( 1._pr - ro ) )
        else
            a2 = 1._pr/2._pr - 2._pr*x/3._pr + x**2._pr/4._pr - 2._pr*x**3._pr/15._pr + x**4._pr/12._pr - &
            2._pr*x**5._pr/35._pr + x**6._pr/24._pr - 2._pr*x**7._pr/63._pr + x**8._pr/40._pr - &
            2._pr*x**9._pr/99._pr + x**10._pr/60._pr - 2._pr*x**11._pr/143._pr + x**12._pr/84._pr - &
            2._pr*x**13._pr/195._pr + x**14._pr/112._pr - 2._pr*x**15._pr/255._pr
        end if
    end if

    if (present(a3)) then
        if (abs(x) > eps) then
            a3 = 8._pr * en / ( 1._pr - ro ) * ( 1._pr + log_ro / ( 1._pr - ro ) )
        else
            a3 = -4._pr + 8._pr*x/3._pr - 2._pr*x**2._pr + 8._pr*x**3._pr/5._pr - 4._pr*x**4._pr/3._pr + &
            8._pr*x**5._pr/7._pr - x**6._pr + 8._pr*x**7._pr/9._pr - 4._pr*x**8._pr/5._pr + &
            8._pr*x**9._pr/11._pr - 2._pr*x**10._pr/3._pr + 8._pr*x**11._pr/13._pr - 4._pr*x**12._pr/7._pr + &
            8._pr*x**13._pr/15._pr - x**14._pr/2._pr + 8._pr*x**15._pr/17._pr
            a3 = a3 * en
        end if
    end if

end subroutine Coeff

subroutine diff_Coeff(p , en , dA , db1 , db2 , db3 , da1, da2 , da3)
    real(pr) , intent(in) :: p , en
    real(pr) , intent(out) , optional :: dA , db1 , db2 , db3 , da1, da2 , da3
    real(pr) :: ro , log_ro , x , eps , dro , e2

    eps = 5.e-2_pr
    x =  p*p - en*en
    ro = 1._pr + x
    log_ro = Log(ro)
    dro = -2._pr * en
    e2 = en * en

    if (present(dA)) then
        if (abs(x) > eps) then
            dA = dro * ( 1._pr / ro - ro + 2._pr * log_ro  ) / ( 1._pr - ro )**3._pr
        else
            dA = 1._pr/3._pr - x/2._pr + (3._pr * x**2._pr)/5._pr - (2._pr * x**3._pr)/3._pr + &
            (5._pr * x**4._pr)/7._pr - (3._pr * x**5._pr)/4._pr + (7._pr * x**6._pr)/9._pr - &
            (4._pr * x**7._pr)/5._pr + (9._pr * x**8._pr)/11._pr - (5._pr * x**9._pr)/6._pr + &
            (11._pr * x**10._pr)/13._pr - (6._pr * x**11._pr)/7._pr + (13._pr * x**12._pr)/15._pr - &
            (7._pr * x**13._pr)/8._pr + (15._pr * x**14._pr)/17._pr - (8._pr * x**15._pr)/9._pr
            dA = dA * dro
        endif
    endif


    if (present(db1)) then
        if  (abs(x) > eps) then
            db1 = dro *( 3._pr - ro + 2._pr * log_ro / ( 1._pr - ro ) ) / ( 1._pr - ro )**2._pr
        else
            db1 = -(2._pr/3._pr) + x/2._pr - (2._pr * x**2._pr)/5._pr + x**3._pr/3._pr - &
            (2._pr * x**4._pr)/7._pr + x**5._pr/4._pr - (2._pr * x**6._pr)/9._pr + x**7._pr/5._pr - &
            (2._pr * x**8._pr)/11._pr + x**9._pr/6._pr - (2._pr * x**10._pr)/13._pr + x**11._pr/7._pr - &
            (2._pr * x**12._pr)/15._pr + x**13._pr/8._pr - (2._pr * x**14._pr)/17._pr + x**15._pr/9._pr
            db1 = db1 * dro
        endif
    endif


    if (present(db2)) then
        if  (abs(x) > eps) then
            db2 = -2._pr * dro / ( 1._pr - ro )**4._pr * ( 2._pr / ro + 3._pr - 6._pr * ro + ro**2._pr + 6._pr * log_ro )
        else
            db2 = -1._pr + (8._pr * x)/5._pr - 2._pr * x**2._pr + (16._pr * x**3._pr)/7._pr - &
            (5._pr * x**4._pr)/2._pr + (8._pr * x**5._pr)/3._pr - (14._pr * x**6._pr)/5._pr + &
            (32._pr * x**7._pr)/11._pr - 3._pr * x**8._pr + (40._pr * x**9._pr)/13._pr - &
            (22._pr * x**10._pr)/7._pr + (16._pr * x**11._pr)/5._pr - (13._pr * x**12._pr)/4._pr + &
            (56._pr * x**13._pr)/17._pr - (10._pr * x**14._pr)/3._pr + (64._pr * x**15._pr)/19._pr
            db2 = db2 * dro
        endif
    endif


    if (present(db3)) then
        if  (abs(x) > eps) then
            db3 = 8._pr * dro / ( 1._pr - ro )**3._pr * ( 1._pr / ro - ro + 2._pr * log_ro )
        else
            db3 = 8._pr/3._pr - 4._pr * x + (24._pr * x**2._pr)/5._pr - (16._pr * x**3._pr)/3._pr + &
            (40._pr * x**4._pr)/7._pr - 6._pr * x**5._pr + (56._pr * x**6._pr)/9._pr - &
            (32._pr * x**7._pr)/5._pr + (72._pr * x**8._pr)/11._pr - (20._pr * x**9._pr)/3._pr + &
            (8._pr * x**10._pr)/13._pr - (48._pr * x**11._pr)/7._pr + (104._pr * x**12._pr)/15._pr - &
            7._pr * x**13._pr + (120._pr * x**14._pr)/17._pr - (64._pr * x**15._pr)/9._pr
            db3 = db3 * dro
        endif
    endif


    if (present(da1)) then
        if  (abs(x) > eps) then
            da1 =-2._pr / ( 1._pr - ro )**2._pr                                             &
              * ( 3._pr - ro + 2._pr / ( 1._pr - ro )                                     &
                * ( log_ro - e2 * ( 2._pr / ro + 5._pr - ro                              &
                                   + 6._pr * log_ro / ( 1._pr - ro ) ) ) )
        else
            da1 = ( 4._pr/3._pr + 2._pr*e2 )                                                  &
              - ( 1._pr + 16._pr*e2/5._pr )*x + ( 4._pr/5._pr + 4._pr*e2 )*x**2._pr        &
              - ( 2._pr/3._pr + 32._pr*e2/7._pr )*x**3._pr + ( 4._pr/7._pr + 5._pr*e2 )*x**4._pr   &
              - ( 1._pr/2._pr + 16._pr*e2/3._pr )*x**5._pr + ( 4._pr/9._pr + 28._pr*e2/5._pr )*x**6._pr   &
              - ( 2._pr/5._pr + 64._pr*e2/11._pr )*x**7._pr + ( 4._pr/11._pr + 6._pr*e2 )*x**8._pr  &
              - ( 1._pr/3._pr + 80._pr*e2/13._pr )*x**9._pr + ( 4._pr/13._pr + 44._pr*e2/7._pr )*x**10._pr &
              - ( 2._pr/7._pr + 32._pr*e2/5._pr )*x**11._pr
        endif
    endif


    if (present(da2)) then
        if  (abs(x) > eps) then
            da2 =-2._pr * en / ( 1._pr - ro )**2._pr * ( 3._pr - ro + 2._pr * log_ro / ( 1._pr - ro ) )
        else
            da2 = 4._pr/3._pr - x + 4._pr*x**2._pr/5._pr - 2._pr*x**3._pr/3._pr + 4._pr*x**4._pr/7._pr - &
            x**5._pr/2._pr + 4._pr*x**6._pr/9._pr - 2._pr*x**7._pr/5._pr + 4._pr*x**8._pr/11._pr - &
            x**9._pr/3._pr + 4._pr*x**10._pr/13._pr - 2._pr*x**11._pr/7._pr + 4._pr*x**12._pr/15._pr - &
            x**13._pr/4._pr + 4._pr*x**14._pr/17._pr - 2._pr*x**15._pr/9._pr
            da2 = da2 * en
        endif
    endif


    if (present(da3)) then
        if  (abs(x) > eps) then
            da3 = 8._pr / ( 1._pr - ro ) * ( 1._pr + 1._pr / ( 1._pr - ro )                 &
                * ( log_ro - 2._pr * e2 * ( 1._pr / ro                                    &
                  + 1._pr + 2._pr * log_ro / ( 1._pr - ro ) ) ) )
        else
            da3 = - ( 4._pr + 16._pr*e2/3._pr ) + ( 8._pr/3._pr + 8._pr*e2 )*x                   &
              - ( 2._pr + 48._pr*e2/5._pr )*x**2._pr + ( 8._pr/5._pr + 32._pr*e2/3._pr )*x**3._pr                 &
              - ( 4._pr/3._pr + 80._pr*e2/7._pr )*x**4._pr + ( 8._pr/7._pr + 12._pr*e2 )*x**5._pr                 &
              - ( 1._pr + 112._pr*e2/9._pr )*x**6._pr + ( 8._pr/9._pr + 64._pr*e2/5._pr )*x**7._pr                 &
              - ( 4._pr/5._pr + 144._pr*e2/11._pr )*x**8._pr + ( 8._pr/11._pr + 40._pr*e2/3._pr )*x**9._pr             &
              - ( 2._pr/3._pr + 176._pr*e2/13._pr )*x**10._pr + ( 8._pr/13._pr + 96._pr*e2/7._pr )*x**11._pr
        endif
    endif

end subroutine diff_Coeff

subroutine diff2_Coeff(p , en , ddA , ddb1 , ddb2 , ddb3 , dda1, dda2 , dda3)
    real(pr) , intent(in) :: p , en
    real(pr) , intent(out) , optional :: ddA , ddb1 , ddb2 , ddb3 , dda1, dda2 , dda3
    real(pr) :: ro , log_ro , x , eps , dro , e2 , ddAro , dA , ddb1ro , db1 , db2 , db3 , ddb2ro , ddb3ro
    real(pr) :: ddro , a1r , a1rr , a2r , a2rr , a3r , a3rr
    eps = 5.e-2_pr
    x =  p*p - en*en
    ro = 1._pr + x
    log_ro = Log(ro)
    dro = -2._pr * en
    ddro = -2._pr
    e2 = en * en

    if (present(ddA)) then
        if (abs(x) > eps) then
        ddAro = - (3._pr + 1._pr / ro**2 - 6._pr / ro + 2._pr*ro - 6._pr*Log_ro) / (1._pr - ro)**4
        else
        ddAro = -(1._pr/2._pr) + (6._pr * x)/5._pr - 2._pr * x**2  + &
         (20._pr * x**3)/7._pr - (15._pr * x**4)/4._pr + (14._pr * x**5)/3._pr - &
         (28._pr * x**6)/5._pr + (72._pr * x**7)/11._pr - (15._pr * x**8)/2._pr + &
         (110._pr * x**9)/13._pr - (66._pr * x**10)/7._pr + (52._pr * x**11)/5._pr - &
         (91._pr * x**12)/8._pr + (210._pr * x**13)/17._pr - (40._pr * x**14)/3._pr + &
         (272._pr * x**15)/19._pr
        end if
    call diff_Coeff(p , en , dA = dA)
    ddA = ddAro * dro**2 + dA / en    ! ddro / dro = 1 / en
    end if

    if (present(ddb1)) then
        if (abs(x) > eps) then
        ddb1ro = (3._pr + ro * (ro - 6._pr) + 6._pr * log_ro + 2._pr / ro) / (1._pr - ro)**4
        else
        ddb1ro = 1._pr/2._pr - (4._pr * x)/5._pr + x**2 - (8._pr * x**3)/7._pr + &
                 5._pr * x**4/4._pr - (4._pr * x**5)/3._pr + (7._pr * x**6)/5._pr - &
                 (16._pr * x**7)/11._pr + 3._pr * x**8/2._pr - (20._pr * x**9)/13._pr + &
                 (11._pr * x**10)/7._pr - (8._pr * x**11)/5._pr + (13._pr * x**12)/8._pr - &
                 (28._pr * x**13)/17._pr + (5._pr * x**14)/3._pr - (32._pr * x**15)/19._pr
        endif
    call diff_Coeff(p , en , db1 = db1)
    ddb1 = ddb1ro * dro**2 + db1 / en
    end if

    if (present(ddb2)) then
        if (abs(x) > eps) then
        ddb2ro = 4._pr / (1._pr - ro)**5 * (1._pr / ro**2 - 8._pr / ro + 8._pr * ro - ro**2 - 12._pr * log_ro)
        else
        ddb2ro = 8._pr/5._pr - 4._pr * x + (48._pr * x**2)/7._pr - 10._pr * x**3 + &
                 (40._pr * x**4)/3._pr - (84._pr * x**5)/5._pr + (224._pr * x**6)/11._pr - &
                 24._pr * x**7 + (360._pr * x**8)/13._pr - (220._pr * x**9)/7._pr + &
                 (176._pr * x**10)/5._pr - 39._pr * x**11 + (728._pr * x**12)/17._pr - &
                 (140._pr * x**13)/3._pr + (960._pr * x**14)/19._pr - (272._pr * x**15)/5._pr
        endif
    call diff_Coeff(p , en , db2 = db2)
    ddb2 = ddb2ro * dro**2 + db2 / en
    endif

    if (present(ddb3)) then
        if (abs(x) > eps) then
        ddb3ro = 8._pr / (1._pr - ro)**4 * (-3._pr - 1._pr / ro**2 + 6._pr / ro - 2._pr * ro + 6._pr * log_ro)
        else
        ddb3ro = -4._pr + (48._pr * x)/5._pr - 16._pr * x**2 + (160._pr * x**3)/7._pr - &
                 30._pr * x**4 + (112._pr * x**5)/3._pr - (224._pr * x**6)/5._pr + &
                 (576._pr * x**7)/11._pr - 60._pr * x**8 + (880._pr * x**9)/13._pr - &
                 (528._pr * x**10)/7._pr + (416._pr * x**11)/5._pr - 91._pr * x**12 + &
                 (1680._pr * x**13)/17._pr - (320._pr * x**14)/3._pr + (2176._pr * x**15)/19._pr
        endif
    call diff_Coeff(p , en , db3 = db3)
    ddb3 = ddb3ro * dro**2 + db3 / en
    endif

    if (present(dda1)) then
        if (abs(x) > eps) then
        a1r = -2 * (3._pr + ro *(ro-6._pr) + 6 * Log_ro + 2._pr /ro) / (1._pr - ro)**4
        a1rr = -4 * (ro**2 - 8*ro + 12*Log_ro + 8._pr / ro - 1._pr / ro**2) / (1._pr - ro)**5
        else
        a1r = -1._pr + (8._pr/5._pr)*x - 2._pr*x**2 + (16._pr/7._pr)*x**3 &
        - (5._pr/2._pr)*x**4 + (8._pr/3._pr)*x**5 - (14._pr/5._pr)*x**6 &
        + (32._pr/11._pr)*x**7 - 3._pr*x**8 + (40._pr/13._pr)*x**9 &
        - (22._pr/7._pr)*x**10 + (16._pr/5._pr)*x**11 - (13._pr/4._pr)*x**12 &
        + (56._pr/17._pr)*x**13 - (10._pr/3._pr)*x**14 + (64._pr/19._pr)*x**15

        a1rr = (8._pr/5._pr) - 4._pr*x + (48._pr/7._pr)*x**2 - 10._pr*x**3 &
        + (40._pr/3._pr)*x**4 - (84._pr/5._pr)*x**5 + (224._pr/11._pr)*x**6 &
        - 24._pr*x**7 + (360._pr/13._pr)*x**8 - (220._pr/7._pr)*x**9 &
        + (176._pr/5._pr)*x**10 - 39._pr*x**11 + (728._pr/17._pr)*x**12 &
        - (140._pr/3._pr)*x**13 + (960._pr/19._pr)*x**14 - (272._pr/5._pr)*x**15
        endif
    dda1 = (a1rr*dro**2 + a1r*ddro) * en + 2*(a1r)*dro
    endif

    if (present(dda2)) then
        if (abs(x) > eps) then
        a2r =  (3._pr - 4 *ro + ro**2 + 2 * Log_ro) / (1._pr - ro)**3
        a2rr = (3._pr + ro * (ro -6._pr) + 6 * log_ro + 2._pr / ro) / (1._pr - ro)**4
        else
        a2r = -(2._pr/3._pr) + (1._pr/2._pr)*x - (2._pr/5._pr)*x**2 &
        + (1._pr/3._pr)*x**3 - (2._pr/7._pr)*x**4 + (1._pr/4._pr)*x**5 &
        - (2._pr/9._pr)*x**6 + (1._pr/5._pr)*x**7 - (2._pr/11._pr)*x**8 &
        + (1._pr/6._pr)*x**9 - (2._pr/13._pr)*x**10 + (1._pr/7._pr)*x**11 &
        - (2._pr/15._pr)*x**12 + (1._pr/8._pr)*x**13 - (2._pr/17._pr)*x**14 &
        + (1._pr/9._pr)*x**15

        a2rr = 1._pr/2._pr - (4._pr/5._pr)*x + 1._pr*x**2 - (8._pr/7._pr)*x**3 &
        + (5._pr/4._pr)*x**4 - (4._pr/3._pr)*x**5 + (7._pr/5._pr)*x**6 &
        - (16._pr/11._pr)*x**7 + (3._pr/2._pr)*x**8 - (20._pr/13._pr)*x**9 &
        + (11._pr/7._pr)*x**10 - (8._pr/5._pr)*x**11 + (13._pr/8._pr)*x**12 &
        - (28._pr/17._pr)*x**13 + (5._pr/3._pr)*x**14 - (32._pr/19._pr)*x**15
        endif
    dda2 = a2rr * dro**2 + a2r * ddro
    endif

    if (present(dda3)) then
        if  (abs(x) > eps) then
        a3r = 8 * (2 * Log_ro - ro + 1._pr / ro) / (1._pr - ro)**3
        a3rr = -8 * (2 * ro + 3._pr - 6 * log_ro - 6._pr / ro + 1._pr / ro**2) / (1._pr - ro)**4
        else
        a3r = (8._pr/3._pr) - 4._pr*x + (24._pr/5._pr)*x**2 - (16._pr/3._pr)*x**3 &
        + (40._pr/7._pr)*x**4 - 6._pr*x**5 + (56._pr/9._pr)*x**6 &
        - (32._pr/5._pr)*x**7 + (72._pr/11._pr)*x**8 - (20._pr/3._pr)*x**9 &
        + (88._pr/13._pr)*x**10 - (48._pr/7._pr)*x**11 + (104._pr/15._pr)*x**12 &
        - 7._pr*x**13 + (120._pr/17._pr)*x**14 - (64._pr/9._pr)*x**15

        a3rr = -4._pr + (48._pr/5._pr)*x - 16._pr*x**2 + (160._pr/7._pr)*x**3 &
        - 30._pr*x**4 + (112._pr/3._pr)*x**5 - (224._pr/5._pr)*x**6 &
        + (576._pr/11._pr)*x**7 - 60._pr*x**8 + (880._pr/13._pr)*x**9 &
        - (528._pr/7._pr)*x**10 + (416._pr/5._pr)*x**11 - 91._pr*x**12 &
        + (1680._pr/17._pr)*x**13 - (320._pr/3._pr)*x**14 &
        + (2176._pr/19._pr)*x**15
        endif
    dda3 = (a3rr*dro**2 + a3r*ddro) * en + 2*(a3r)*dro
    endif

end subroutine diff2_Coeff

end module  coeff_2004   
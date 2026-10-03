module prep_to_1pot
    use parameters
    use others

    implicit none
    
contains

    subroutine prep (x , y , en , p1 , p2 , ksi , N , a_ , dN, ddN)
    real(pr) , intent(in) :: p1 , p2 , ksi , x ,y , en
    real(pr) , intent(out) , optional :: N , a_ , dN, ddN
    real(pr) :: pp_1 , pp_2 , ro1 , ro2 , b2

    pp_1 = en * en - p1 * p1
    pp_2 = en * en - p2 * p2
    ro1 = 1._pr - pp_1
    ro2 = 1._pr - pp_2

    b2 = y * y * pp_1 + (1._pr - y) ** 2 * pp_2 + 2 * y * (1._pr - y) * (en * en - p1 * p2 * ksi)

    if (present(N)) then
            N = x * b2 + y * ro1 + (1._pr - y) * ro2
    endif

    if (present(a_)) then
        a_ = x * y * (1._pr - x * y) * pp_1 + x * (1._pr - y)*(1._pr - x + x * y) * pp_2 &
        - 2 * (1._pr - x * y) * (1._pr - x + x * y) * (en * en - p1 * p2 * ksi)
    endif

    if (present(dN)) dN = 2 * en * (x - 1._pr)

    if (present(ddN)) ddN = 2 * (x - 1._pr)

    end subroutine prep


    subroutine prepare_1pot(x, y , en, A0 , C1 , C2 , D1, D2, Ff1, Ff2, Gg1, Gg2 , H1 , B , C , D , H  )
    real(pr) , intent(in) :: x , y , en
    real(pr) , intent(out) , optional :: A0 , C1 , C2 , D1, D2, Ff1, Ff2, Gg1, Gg2, H1 , B , C , D , H

    ! Локальные переменные для общих множителей
    real(pr) :: m_1x, m_1xy, m_bracket, m_x1y

    ! --- Группа множителей для экономии ---

    ! 1. Множитель (1 - x). Нужен для A0, Ff1, Ff2, H1, B, C, D, H
    if (present(A0) .or. present(Ff1) .or. present(Ff2) .or. present(H1) .or. &
        present(B)  .or. present(C)   .or. present(D)   .or. present(H)) then
        m_1x = 1.0_pr - x
    end if

    ! 2. Множитель (1 - x*y). Нужен для A0, C1, D1, Ff1, H1, B
    if (present(A0) .or. present(C1) .or. present(D1) .or. present(Ff1) .or. &
        present(H1) .or. present(B)) then
        m_1xy = 1.0_pr - x * y
    end if

    ! 3. множитель (1 - x + x*y). Нужен для D1, D2, C
    if (present(D1) .or. present(D2) .or. present(C)) then
        m_bracket = 1.0_pr - x + x * y
    end if

    ! 4. Множитель x * (1 - y). Нужен для C2, D2, Ff2
    if (present(C2) .or. present(D2) .or. present(Ff2)) then
        m_x1y = x * (1.0_pr - y)
    end if


    if (present(A0))  A0  = 1.0_pr + 2.0_pr * en * en * m_1x * m_1xy   !fake warning
    if (present(C1))  C1  = x * y * m_1xy
    if (present(C2))  C2  = -x * x * y * (1.0_pr - y)!-x**2 *y*(1.0_pr-y) !
    if (present(D1))  D1  = -m_bracket * m_1xy
    if (present(D2))  D2  = m_bracket * m_x1y
    if (present(Ff1)) Ff1 = -en * m_1x * m_1xy
    if (present(Ff2)) Ff2 = -en * m_1x * m_x1y
    if (present(Gg1)) Gg1 = 2.0_pr
    if (present(Gg2)) Gg2 = -2.0_pr
    if (present(H1))  H1  = -2.0_pr * en * m_1x * m_1xy
    if (present(B))   B   = 2.0_pr * m_1xy * m_1x
    if (present(C))   C   = 2.0_pr * m_bracket * m_1x
    if (present(D))   D   = -m_1x
    if (present(H))   H   = -4.0_pr * m_1x

end subroutine prepare_1pot
    

    subroutine P_i_2004 (x , y , en , p1 , p2 , ksi , g1 , f1 , g2, f2, P_1 , P_2, P_3, P_4 , P_5, P_6)
    real(pr) , intent(in) :: x , y , p1 , p2 , ksi , g1 , f1 , g2, f2 , en
    real(pr) , intent(out) :: P_1 , P_2, P_3, P_4 , P_5, P_6
    real(pr) :: A0 , C1 , C2 , D1, D2, Ff1, Ff2, Gg1, Gg2 , H1 , fctr , N

    call prepare_1pot (x, y , en, A0 = A0 , C1 = C1 , C2 = C2 , D1 = D1, D2 = D2, Ff1 = Ff1, Ff2 = Ff2 &
    , Gg1 = Gg1, Gg2 = Gg2 , H1 = H1)
    call prep (x , y , en, p1, p2 , ksi , N = N)
    fctr = (1._pr - y) / (N * N)

    P_1 = fctr * (A0 * g1 * g2 + Ff1 * (en * g1 + p1 * f1)*g2 + Ff2 * g1 * (en * g2 + p2 * f2))

    P_2 = fctr * (A0 * f1 * f2 + Ff1 * (en * f1 + p1 * g1)*f2 + Ff2 * f1 * (en * f2 + p2 * g2))

    P_3 = fctr * (C1 * (en * g1 + p1 * f1)*f2 + C2 * g1 * (en * f2 + p2 * g2) + (H1 - Gg1) * g1 *f2)

    P_4 = fctr * (C1 * (en * f1 + p1 * g1)*g2 + C2 * f1 * (en * g2 + p2 * f2) + (H1 + Gg1) * f1 * g2)

    P_5 = fctr * (D1 * (en * g1 + p1 * f1)*f2 + D2 * g1 * (en * f2 + p2 * g2) - Gg2 * g1 * f2)

    P_6 = fctr * (D1 * (en * f1 + p1 * g1)*g2 + D2 * f1 * (en * g2 + p2 * f2) + Gg2 * f1 * g2)

    end subroutine  P_i_2004


    subroutine R_i_2004(x , y , en , p1 , p2 , ksi , g1 , f1 , g2, f2 ,R_1 , R_2 , R_3 , R_4 , R_5 , R_6)
    real(pr), intent(in) :: x , y , p1 , p2 , ksi, g1 , f1 , g2, f2 , en
    real(pr), intent(out) :: R_1 , R_2 , R_3 , R_4 , R_5 , R_6
    real(pr) :: K1 , K2 , a_ , N , fctr , fctr2 , B, C , D , H

    call prep (x , y , en ,p1 , p2 , ksi , N = N , a_ = a_)
    call prepare_1pot(x, y ,en, B = B, C = C, D = D, H = H)

    K1 = (a_ + 1.0_pr + 2 * x * N + en * H)*g1*g2 + en * B * (en * g1 + p1 * f1)*g2 &
    + en * C * g1 * (en * g2 + p2 * f2) + D * (en * g1 + p1 * f1) * (en * g2 + p2 * f2)

    K2 = (a_ + 1.0_pr + 2 * x * N - en * H)*f1*f2 + en * B * (en * f1 + p1 * g1)*f2 &
    + en * C * f1 * (en * f2 + p2 * g2) + D * (en * f1 +p1 * g1) * (en * f2 + p2 * g2)

    fctr = 1._pr / N
    fctr2 = 2._pr * (1._pr - y) * fctr * fctr

    R_1 = fctr * (en * C * g1 * f2 + D * (en * g1 + p1 * f1)*f2)

    R_2 = fctr * (en * C * f1 * g2 + D * (en * f1 + p1 * g1)*g2)

    R_3 =  -fctr2 * x * y * K1 - 2 * (1.0_pr - x) * g1 * g2 * fctr

    R_4 = fctr2 * (1.0_pr - x + x * y) * K1

    R_5 = -fctr2 * x * y * K2 - 2 * (1.0_pr - x) * fctr * f1 * f2

    R_6 = fctr2 * (1.0_pr - x + x * y) * K2

    end subroutine R_i_2004

    subroutine prepare_red_1pot(x , y , en, p1 , p2 , ksi , g1 , f1 , g2 , f2 , F1r, F2r, dF1r, dF2r, ddF1r, ddF2r, d3F1r, d3F2r)
    real(pr) , intent(in) :: x , y , p1 , p2 , ksi , g1 , f1 , g2 , f2 , en
    real(pr) , intent(out) , optional :: F1r, F2r, dF1r, dF2r, ddF1r, ddF2r, d3F1r, d3F2r
    real(pr) :: A , a_ , N , B , C, D , H , K6 , K7 ,K8 ,K9  , K10, dK10, dN, ddN, fctrN

    if (present(F1r) .or. present(F2r)) then
             call prep (x , y , en , p1 , p2 , ksi , N = N , a_ = a_, dN = dN)
             call prepare_1pot(x, y , en , B = B, C = C, D = D, H = H)
             A = a_ + 1.0_pr - N * ( 0.75_pr + x * log(N))
    endif

    if (present(dF1r) .or. present(dF2r)) then
            K6 = (x - 1._pr) * (2 * x * y - 1._pr) !(B+D)
            K7 = (1._pr - x) * (1._pr + 2 * x * (y - 1._pr)) !(C+D)
            K8 = 0.5_pr * (1._pr - x)* (-8._pr + en * (3._pr + 4 * x) + 4 * en * x * log(N)) !(A' + H + en(B+C))
            K9 = 0.5_pr * (1._pr - x)* ( 8._pr + en * (3._pr + 4 * x) + 4 * en * x * log(N)) !(A' - H + en(B+C))
    endif

    if (present(ddF1r) .or. present(ddF2r)) K10 = 0.5_pr*(1._pr-x) * (3._pr + 4 * x *(1._pr + log(N) &
     + 2 * en**2 *(x - 1._pr) / N)) + 2 * (1._pr - x)**2

    if (present(d3F1r) .or. present(d3F2r)) then
        fctrN = 1._pr / N
        dK10 = 2._pr * x * (1._pr - x) * ( fctrN * dN + 2 * (x - 1._pr) * fctrN * (2 * en - en*en * fctrN * dN) )
    endif    
    if (present(F1r)) F1r = (A + en * H)* g1 * g2 + en * B * (en * g1 + p1 * f1)*g2 + en * C * g1 * (en * g2 + p2 * f2)  &
    + D*(en * g1 + p1 * f1)*(en * g2 + p2 * f2)

    if (present(F2r)) F2r = (A - en * H)* f1 * f2 + en * B * (en * f1 + p1 * g1)*f2 + en * C * f1 * (en * f2 + p2 * g2)  &
    + D*(en * f1 + p1 * g1)*(en * f2 + p2 * g2)

    if (present(dF1r)) dF1r = K8 * g1 * g2 + (en * g1 + p1 * f1)*g2 * K6 + (en * g2 + p2 * f2)*g1 * K7

    if (present(dF2r)) dF2r = K9 * f1 * f2 + (en * f1 + p1 * g1)*f2 * K6 + (en * f2 + p2 * g2)*f1 * K7

    if (present(ddF1r)) ddF1r = K10 *  g1 * g2

    if (present(ddF2r)) ddF2r = K10 * f1 * f2
    
    if (present(d3F1r)) d3F1r = dK10 *  g1 * g2

    if (present(d3F2r)) d3F2r = dK10 * f1 * f2


    end subroutine prepare_red_1pot

    subroutine prepare_dpot1(x, y , en, dA0 , dFf1, dFf2 , dH1, Kp1 , Kp2 , Kp3 , Kp4 , Kp5 , Kp6 , Kp7)
    real(pr) , intent(in) :: x, y , en
    real(pr) , intent(out) , optional ::  dA0 , dFf1, dFf2 , dH1 , Kp1 , Kp2 , Kp3 ,  Kp4 , Kp5 , Kp6 , Kp7
    real(pr) :: term_x, term_xy, common_prod

    ! Вычисляем общие части один раз
    term_x = 1._pr - x
    term_xy = 1._pr - x * y

    ! Общий множитель для dA0, dFf1 и dH1
    common_prod = term_x * term_xy


    if (present(dA0))  dA0  = 4.0_pr * en * common_prod
    if (present(dFf1)) dFf1 = -common_prod
    if (present(dH1))  dH1  = -2.0_pr * common_prod
    if (present(dFf2)) dFf2 = -term_x * x * (1._pr - y)
    if (present(Kp1)) Kp1 = -en * term_x * (-3._pr + x + 2 * x * y )
    if (present(Kp2)) Kp2 = term_x * (3 * x * y - 2._pr)
    if (present(Kp3)) Kp3 = -term_x * (1._pr - x + x * y)
    if (present(Kp4)) Kp4 = -2 * term_x * (2._pr + en * (3 * x - 1))
    if (present(Kp5)) Kp5 = -2 * term_x  * (-2._pr + en * (3 * x - 1))
    if (present(Kp6)) Kp6 = -term_x * (2 * x * y - 1._pr)
    if (present(Kp7)) Kp7 = term_x * (1._pr + 2 * x * (y - 1._pr))

    end subroutine prepare_dpot1

    subroutine dP_i_2004 (x , y , en, p1 , p2 , ksi , g1 , f1 , g2, f2, dP_1 , dP_2, dP_3, dP_4 , dP_5, dP_6)
    real(pr) , intent(in) :: x , y , p1 , p2 , ksi , g1 , f1 , g2, f2 , en
    real(pr) , intent(out) :: dP_1 , dP_2, dP_3, dP_4 , dP_5, dP_6
    real(pr) :: N , dN , fctrNN , fctrdN , fctrN , dFf1, dFf2 , dH1, Kp1 , Kp2 , Kp3
    real(pr) :: P_1 , P_2, P_3, P_4 , P_5, P_6

    call prep (x , y , en, p1, p2 , ksi , N = N , dN = dN)
    call prepare_dpot1(x, y , en ,dFf1 = dFf1 , dFf2 = dFf2 , dH1 = dH1, Kp1 = Kp1 , Kp2 = Kp2 , Kp3 = Kp3)
    call P_i_2004 (x , y , en, p1 , p2 , ksi , g1 , f1 , g2, f2, P_1 , P_2, P_3, P_4 , P_5, P_6)

    fctrN = 1._pr / N
    fctrNN =  (1._pr - y) * fctrN**2
    fctrdN = -2 * fctrN * dN

    dP_1 = fctrdN * P_1 + (g1 * g2 * Kp1 + dFf1 * (en * g1 + p1 * f1) * g2 + dFf2*(en * g2 + p2 * f2)*g1) * fctrNN
    dP_2 = fctrdN * P_2 + (f1 * f2 * Kp1 + dFf1 * (en * f1 + p1 * g1) * f2 + dFf2 * f1 * (en * f2 + p2 * g2)) * fctrNN
    dP_3 = fctrdN * P_3 + g1 * f2 * Kp2 * fctrNN
    dP_4 = fctrdN * P_4 + f1 * g2 * Kp2 *  fctrNN
    dP_5 = fctrdN * P_5 + g1 * f2 * Kp3 * fctrNN
    dP_6 = fctrdN * P_6 + f1 * g2 * Kp3 * fctrNN

    end subroutine dP_i_2004

    subroutine dR_i_2004(x , y , en, p1 , p2 , ksi , g1 , f1 , g2, f2 , dR_1 , dR_2 , dR_3 , dR_4 , dR_5 , dR_6)
    real(pr), intent(in) :: x , y , p1 , p2 , ksi, g1 , f1 , g2, f2 , en
    real(pr), intent(out) :: dR_1 , dR_2 , dR_3 , dR_4 , dR_5 , dR_6
    real(pr) :: N , dN , R_1 , R_2 , R_3 , R_4 , R_5 , R_6 , Kp4 , Kp5 , Kp6 , Kp7 , Kr1 , Kr2
    real(pr) :: fctrN2 , fctrdN , fctrN , fctr1 , fctr2 , fctr3

    call prep (x , y ,en , p1 , p2 , ksi , N = N , dN = dN)
    call prepare_dpot1 (x, y , en , Kp4 = Kp4 , Kp5 = Kp5 , Kp6 = Kp6 , Kp7 = Kp7 )
    call R_i_2004(x , y , en , p1 , p2 , ksi , g1 , f1 , g2, f2 ,R_1 , R_2 , R_3 , R_4 , R_5 , R_6)

    Kr1 = g1 * g2 * Kp4 + p1 * f1 * g2 * Kp6 + p2 * g1 * f2 * Kp7
    kr2 = f1 * f2 * Kp5 + p1 * g1 * f2 * Kp6 + p2 * f1 * g2 * Kp7

    fctrN = 1._pr / N
    fctrN2 = fctrN * fctrN
    fctrdN = -fctrN * dN

    fctr1 = (1._pr - x) * (1._pr + 2 * x * (y - 1._pr))
    fctr2 = 2 * (1._pr - y)*(1._pr - x + x*y)
    fctr3 = 2 * (1._pr - y)*(-x * y)

    dR_1 = fctrdN * R_1 + g1 * f2 * fctr1 * fctrN !take it off
    dR_2 = fctrdN * R_2 + f1 * g2 * fctr1 * fctrN
    dR_3 = 2 * fctrdN * (R_3 + (1._pr - x) * g1 * g2 * fctrN) + fctr3 * Kr1 * fctrN * fctrN
    dR_4 = 2 * fctrdN * R_4 + fctr2 * Kr1 * fctrN * fctrN
    dR_5 = 2 * fctrdN * (R_5 + (1._pr - x) * f1 * f2 * fctrN) + fctr3 * Kr2 * fctrN * fctrN
    dR_6 = 2 * fctrdN * R_6 + fctr2 * Kr2 * fctrN * fctrN

    end subroutine dR_i_2004

end module prep_to_1pot
module Integration_Cuba
    use quantum_state
    use parameters
    implicit none

    abstract interface
        real(pr) function func_interface1(p, state)
            import :: pr, q_state
            real(pr), intent(in) :: p
            type(q_state), intent(in) :: state
        end function func_interface1

        real(pr) function func_interface5(x , y , p1 , p2 , ksi , state)
            import :: pr, q_state
            real(pr), intent(in) :: x , y , p1 , p2 , ksi
            type(q_state), intent(in) :: state
        end function func_interface5

        real(pr) function func_interface4(y , p1 , p2 , ksi , state)
            import :: pr, q_state
            real(pr), intent(in) ::  y , p1 , p2 , ksi
            type(q_state), intent(in) :: state
        end function func_interface4
    end interface

    type :: container
        procedure(func_interface1), pointer, nopass :: func1 => null()
        procedure(func_interface5), pointer, nopass :: func5 => null()
        procedure(func_interface4), pointer, nopass :: func4 => null()
        type(q_state), pointer :: state => null()
    end type


    contains
! ==========================================
! Интеграл от 0 до p0
! ==========================================
integer function integrand1_to_p0(ndim, xx, ncomp, f, userdata)
    implicit none
    integer :: ndim, ncomp
    real(pr) :: xx(ndim), f(ncomp)
    type(container) :: userdata

    real(pr) :: p0, p
    p0 =  2._pr * userdata%state%alphaZ
    p = p0 * xx(1)


    f(1) = userdata%func1( p , userdata%state) * p0

    integrand1_to_p0 = 0
end function integrand1_to_p0

integer function integrand5_to_p0(ndim, xx, ncomp, f, userdata)
    implicit none
    integer:: ndim, ncomp
    real(pr) :: xx(ndim)
    type(container) :: userdata
    real(pr) :: f(ncomp)

    real(pr) :: p0, p, r, abs_r, q, x_var, y_var , jacobian

    p0 = 2.0_pr * userdata%state%alphaZ
    p = p0 * xx(1)
    r = p * xx(2)
    abs_r = abs(r)
    q = (sqrt(2.0_pr) * (abs_r + (p - abs_r) * xx(3)))
    x_var = xx(4)
    y_var = xx(5)
    jacobian = p0 * p * (sqrt(2.0_pr) * (p - abs_r))

    f(1) = (func_prq(x_var, y_var , p , r , q , userdata ) & 
    + func_prq(x_var, y_var , p , -r , q , userdata )) * jacobian

    integrand5_to_p0 = 0
end function integrand5_to_p0

integer function integrand4_to_p0(ndim, xx, ncomp, f, userdata)
    implicit none
    integer:: ndim, ncomp
    real(pr) :: xx(ndim)
    type(container) :: userdata
    real(pr) :: f(ncomp)

    real(pr) :: p0, p, r, abs_r, q, y_var , jacobian

    p0 = 2.0_pr * userdata%state%alphaZ
    p = p0 * xx(1)
    r = p * xx(2)
    abs_r = abs(r)
    q = (sqrt(2.0_pr) * (abs_r + (p - abs_r) * xx(3)))
    y_var = xx(4)
    jacobian = p0 * p * (sqrt(2.0_pr) * (p - abs_r))

    f(1) = (func4_prq( y_var , p , r , q , userdata ) & 
    + func4_prq(y_var , p , -r , q , userdata )) * jacobian

    integrand4_to_p0 = 0
end function integrand4_to_p0


! ==========================================
! Интеграл от p0 до бесконечности
! ==========================================
integer function integrand1_to_inf(ndim, xx, ncomp, f, userdata)
    implicit none
    integer :: ndim, ncomp
    real(pr) :: xx(ndim), f(ncomp)
    type(container) :: userdata
    real(pr) :: p0, p

    p0 =  2._pr * userdata%state%alphaZ
    p = p0 / xx(1)

    f(1) = userdata%func1(p , userdata % state) * p * p / p0

    integrand1_to_inf = 0

  end function integrand1_to_inf


integer function integrand5_to_inf(ndim, xx, ncomp, f, userdata)
    implicit none
    integer:: ndim, ncomp
    real(pr) :: xx(ndim)
    type(container) :: userdata
    real(pr) :: f(ncomp)
    
    real(pr) :: p0, p, r, abs_r, q, x_var, y_var , jacobian
    p0 =  2._pr * userdata%state%alphaZ
    p = p0 / xx(1)
    r = p * xx(2)
    abs_r = abs(r)                              ! r in [-p, p]
    q = sqrt(2.0_pr) * (abs_r + (p - abs_r) * xx(3))
    x_var = xx(4)
    y_var = xx(5)

    jacobian = (p * p / p0) * p * sqrt(2.0_pr) * (p - abs_r)
    
    f(1) = (func_prq(x_var, y_var , p , r , q , userdata ) & 
    + func_prq(x_var, y_var , p , -r , q , userdata )) * jacobian

    integrand5_to_inf = 0

  end function integrand5_to_inf  

integer function integrand4_to_inf(ndim, xx, ncomp, f, userdata)
    implicit none
    integer:: ndim, ncomp
    real(pr) :: xx(ndim)
    type(container) :: userdata
    real(pr) :: f(ncomp)
    
    real(pr) :: p0, p, r, abs_r, q, y_var , jacobian
    p0 =  2._pr * userdata%state%alphaZ
    p = p0 / xx(1)
    r = p * xx(2)
    abs_r = abs(r)                              ! r in [-p, p]
    q = sqrt(2.0_pr) * (abs_r + (p - abs_r) * xx(3))
    y_var = xx(4)

    jacobian = (p * p / p0) * p * sqrt(2.0_pr) * (p - abs_r)
    
    f(1) = (func4_prq (y_var , p , r , q , userdata ) & 
    + func4_prq( y_var , p , -r , q , userdata )) * jacobian

    integrand4_to_inf = 0

  end function integrand4_to_inf 

real(pr) function func_prq(x, y, p, r , q, userdata)
    real(pr), intent(in) :: p, r, q, x, y
    type(container), intent(in) :: userdata
    real(pr) :: p1, p2, ksi

    p1 = (p + r) / sqrt(2.0_pr)
    p2 = (p - r) / sqrt(2.0_pr)
    ksi = (p1*p1 + p2*p2 - q*q) / (2.0_pr * p1 * p2)

    func_prq = userdata%func5(x, y, p1, p2, ksi , userdata%state) * q 

end function func_prq

real(pr) function func4_prq(y, p, r , q, userdata)
    real(pr), intent(in) :: p, r, q, y
    type(container), intent(in) :: userdata
    real(pr) :: p1, p2, ksi

    p1 = (p + r) / sqrt(2.0_pr)
    p2 = (p - r) / sqrt(2.0_pr)
    ksi = (p1*p1 + p2*p2 - q*q) / (2.0_pr * p1 * p2)

    func4_prq = userdata%func4( y, p1, p2, ksi , userdata%state) * q 

end function func4_prq

subroutine Cuba_Int1_Cuhre(func, state, res, err, epsrel, epsabs)
    implicit none
    procedure(func_interface1) :: func  
    type(q_state), target :: state
    real(pr) , intent(in) :: epsrel , epsabs 
    real(pr) , intent(out) :: res , err

    integer, parameter :: ndim = 2, ncomp = 1
    integer :: nregions1, neval1, fail1
    integer :: nregions2, neval2, fail2
    real(pr) :: integral1(ncomp), error1(ncomp), prob1(ncomp)
    real(pr) :: integral2(ncomp), error2(ncomp), prob2(ncomp)
    real(pr) :: res1, res2, err1, err2
    type(container) :: cont
    integer :: flags = 0
    integer :: mineval = 1
    integer :: maxeval = 256000000
    integer :: key = 0 , nvec = 1
    character(len = 64) :: statefile
    integer(pr) :: spin
    spin = -1
    statefile = ""
    !инициализация контейнера функции и структуры
    cont%func1 => func
    cont%state => state

        call cuhre(ndim, ncomp, integrand1_to_p0, cont , nvec,&
                epsrel, epsabs, flags, mineval, maxeval, key, statefile , spin ,  &
                nregions1, neval1, fail1, integral1, error1, prob1)

        res1 = integral1(1)
        err1 = error1(1)


        call cuhre(ndim, ncomp, integrand1_to_inf, cont , nvec, &
                epsrel, epsabs, flags, mineval, maxeval, key,statefile , spin , &
                nregions2, neval2, fail2, integral2, error2, prob2)

        res2 = integral2(1)
        err2 = error2(1)

        res = res1 + res2
        err = sqrt(err1**2 + err2**2)

  end subroutine Cuba_Int1_Cuhre

  subroutine Cuba_Int5_Cuhre(func, state, res, err, epsrel, epsabs)
     implicit none
     procedure(func_interface5) :: func  
     type(q_state), target :: state
     real(pr) , intent(in) :: epsrel , epsabs 
     real(pr) , intent(out) :: res , err

     integer, parameter ::  ndim = 5, ncomp = 1
     integer :: nregions1, neval1, fail1
     integer :: nregions2, neval2, fail2
     real(pr) :: integral1(ncomp), error1(ncomp), prob1(ncomp)
     real(pr) :: integral2(ncomp), error2(ncomp), prob2(ncomp)
     real(pr) :: res1, res2, err1, err2
     type(container) :: cont5
     integer :: flags = 0
     integer :: mineval = 1
     integer :: maxeval = 128000000
     integer :: key = 0 , nvec = 1
     character(len = 64) :: statefile
     integer(pr) :: spin
     spin = -1
     statefile = ""
     !инициализация контейнера функции и структуры
     cont5%func5 => func
     cont5%state => state

         call cuhre(ndim, ncomp, integrand5_to_p0, cont5 , nvec,&
                 epsrel, epsabs, flags, mineval, maxeval, key, statefile , spin ,  &
                 nregions1, neval1, fail1, integral1, error1, prob1)

         res1 = integral1(1)
         err1 = error1(1)


         call cuhre(ndim, ncomp, integrand5_to_inf, cont5 , nvec, &
                 epsrel, epsabs, flags, mineval, maxeval, key,statefile , spin , &
                 nregions2, neval2, fail2, integral2, error2, prob2)

         res2 = integral2(1)
         err2 = error2(1)

         res = res1 + res2
         err = sqrt(err1**2 + err2**2)

   end subroutine Cuba_Int5_Cuhre  


  subroutine Cuba_Int4_Cuhre(func, state, res, err, epsrel, epsabs)
     implicit none
     procedure(func_interface4) :: func  
     type(q_state), target :: state
     real(pr) , intent(in) :: epsrel , epsabs 
     real(pr) , intent(out) :: res , err

     integer, parameter ::  ndim = 4, ncomp = 1
     integer :: nregions1, neval1, fail1
     integer :: nregions2, neval2, fail2
     real(pr) :: integral1(ncomp), error1(ncomp), prob1(ncomp)
     real(pr) :: integral2(ncomp), error2(ncomp), prob2(ncomp)
     real(pr) :: res1, res2, err1, err2
     type(container) :: cont4
     integer :: flags = 0
     integer :: mineval = 1
     integer :: maxeval = 128000000
     integer :: key = 0 , nvec = 1
     character(len = 64) :: statefile
     integer(pr) :: spin
     spin = -1
     statefile = ""
     !инициализация контейнера функции и структуры
     cont4%func4 => func
     cont4%state => state

         call cuhre(ndim, ncomp, integrand4_to_p0, cont4 , nvec,&
                 epsrel, epsabs, flags, mineval, maxeval, key, statefile , spin ,  &
                 nregions1, neval1, fail1, integral1, error1, prob1)

         res1 = integral1(1)
         err1 = error1(1)


         call cuhre(ndim, ncomp, integrand4_to_inf, cont4 , nvec, &
                 epsrel, epsabs, flags, mineval, maxeval, key,statefile , spin , &
                 nregions2, neval2, fail2, integral2, error2, prob2)

         res2 = integral2(1)
         err2 = error2(1)

         res = res1 + res2
         err = sqrt(err1**2 + err2**2)

   end subroutine Cuba_Int4_Cuhre 


end module Integration_Cuba
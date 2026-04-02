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
    end interface

    type :: container
        procedure(func_interface1), pointer, nopass :: func => null()
        type(q_state), pointer :: state => null()
    end type

    contains
! ==========================================
! Интеграл от 0 до p0
! ==========================================
integer function integrand_to_p0(ndim, xx, ncomp, f, userdata)
    implicit none
    integer :: ndim, ncomp
    real(pr) :: xx(ndim), f(ncomp)
    type(container) :: userdata

    real(pr) :: p0, p
    p0 =  2._pr * userdata%state%alphaZ
    p = p0 * xx(1)


    f(1) = userdata%func( p , userdata%state) * p0

    integrand_to_p0 = 0
end function integrand_to_p0

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

    f(1) = userdata%func(p , userdata % state) * p * p / p0

    integrand1_to_inf = 0

  end function integrand1_to_inf

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
    cont%func => func
    cont%state => state

        call cuhre(ndim, ncomp, integrand_to_p0, cont , nvec,&
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


end module Integration_Cuba
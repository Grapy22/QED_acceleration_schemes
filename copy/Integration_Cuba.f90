module Integration_Cuba
    use parameters
    use quantum_state 
    use iso_c_binding
    implicit none

    
    interface

    real(pr) function func_interface1(p , state)
        import :: pr, q_state
        real(pr), intent(in) :: p
        type(q_state), intent(in) :: state
    end function func_interface1

    end interface

    ! ============================
    ! Контекст (передаётся в Cuba)
    ! ============================
    type :: integrand_context
        procedure(func_interface1), pointer , nopass  :: func
        type(q_state) :: state
    end type integrand_context
    
contains

! ==========================================
! Интеграл от 0 до p0
! ==========================================
integer function integrand_to_p0(ndim, xx, ncomp, f, userdata)
    use iso_c_binding
    implicit none

    integer, intent(in) :: ndim, ncomp
    real(pr), intent(in) :: xx(ndim)
    real(pr), intent(out) :: f(ncomp)
    type(c_ptr), value :: userdata

    type(integrand_context), pointer :: ctx

    real(pr) :: p0, p

    call c_f_pointer(userdata, ctx)

    p0 = 2._pr * ctx%state%alphaZ
    p  = p0 * xx(1)

    f(1) = ctx%func(p, ctx%state) * p0

    integrand_to_p0 = 0
end function integrand_to_p0

! ==========================================
! Интеграл от p0 до бесконечности
! ==========================================
integer function integrand_to_inf(ndim, xx, ncomp, f, userdata)
    use iso_c_binding
    implicit none

    integer, intent(in) :: ndim, ncomp
    real(pr), intent(in) :: xx(ndim)
    real(pr), intent(out) :: f(ncomp)
    type(c_ptr), value :: userdata

    type(integrand_context), pointer :: ctx

    real(pr) :: p0, p

    call c_f_pointer(userdata, ctx)

    p0 = 2._pr * ctx%state%alphaZ
    p  = p0 / xx(1)

    f(1) = ctx%func(p, ctx%state) * p * p / p0

    integrand_to_inf = 0
end function integrand_to_inf

! ==========================================
! Основной интерфейс
! ==========================================
subroutine integrate_cuhre_1d(func, state, res, err, epsrel, epsabs)

    use iso_c_binding
    implicit none

    procedure(func_interface1) :: func
    type(q_state), intent(in) :: state

    real(pr), intent(out) :: res, err
    real(pr), intent(in)  :: epsrel, epsabs

    type(integrand_context), pointer :: ctx
    type(c_ptr) :: userdata

    integer, parameter :: ndim = 1, ncomp = 1

    integer :: nregions1, neval1, fail1
    integer :: nregions2, neval2, fail2

    real(pr) :: integral1(ncomp), error1(ncomp), prob1(ncomp)
    real(pr) :: integral2(ncomp), error2(ncomp), prob2(ncomp)

    real(pr) :: res1, res2, err1, err2

    integer :: flags, mineval, maxeval, key, nvec
    integer :: spin
    character(len=64) :: statefile

    ! ============================
    ! Настройки Cuba
    ! ============================
    flags   = 0
    mineval = 1
    maxeval = 256000000
    key     = 0
    nvec    = 1
    spin    = -1
    statefile = ""

    ! ============================
    ! Заполняем контекст
    ! ============================
    allocate(ctx)
    ctx%func  => func
    ctx%state = state

    ! ПОЛУЧАЕМ C-АДРЕС
    userdata = c_loc(ctx)

    ! ============================
    ! Интеграл 0 → p0
    ! ============================
    call cuhre(ndim, ncomp, integrand_to_p0, userdata, nvec, &
               epsrel, epsabs, flags, mineval, maxeval, key, &
               statefile, spin, &
               nregions1, neval1, fail1, integral1, error1, prob1)

    res1 = integral1(1)
    err1 = error1(1)

    ! ============================
    ! Интеграл p0 → ∞
    ! ============================
    call cuhre(ndim, ncomp, integrand_to_inf, userdata, nvec, &
               epsrel, epsabs, flags, mineval, maxeval, key, &
               statefile, spin, &
               nregions2, neval2, fail2, integral2, error2, prob2)

    res2 = integral2(1)
    err2 = error2(1)

    ! ============================
    ! Суммируем
    ! ============================
    res = res1 + res2
    err = sqrt(err1**2 + err2**2)
    deallocate(ctx)

end subroutine integrate_cuhre_1d
    
end module Integration_Cuba
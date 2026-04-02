module others
    use quantum_state
    implicit none
    
contains

subroutine print_result(label, res, err, epsabs, scale)
    implicit none
    character(len=*), intent(in) :: label
    real(pr), intent(in) :: res, err, epsabs, scale

    real(pr) :: eff_eps, actual_err, scaled_res
    integer :: digits, power
    character(len=64) :: fmt
    character(len=16) :: scale_str

    actual_err = abs(err * scale)
    scaled_res = res * scale
    
    ! 1. Определяем точность (худшую из запрошенной и полученной)
    eff_eps = max(actual_err, abs(epsabs * scale))

    if (eff_eps > 1e-20_pr) then
        digits = max(0, int(-log10(eff_eps)) + 1) 
    else
        digits = 15
    end if
    digits = min(digits, 15)

    ! 2. Формируем строку множителя на основе scale
    if (abs(scale - 1.0_pr) < 1e-9_pr) then
        scale_str = "" ! Если scale = 1, ничего не пишем
    else
        power = nint(log10(scale))
        write(scale_str, '(" x10^", I0)') power
    end if

    ! 3. Формируем формат (используем F0.X для адаптивной ширины)
    write(fmt, '("(A, F0.", I0, ", A, F0.", I0, ", A)")') digits, digits

    ! 4. Вывод
    write(*, fmt) trim(label)//' : ', scaled_res, ' ± ', actual_err, trim(scale_str)

end subroutine print_result
!     subroutine print_result(label, res, err, epsabs, scale)
!     implicit none

!     character(len=*), intent(in) :: label
!     real(pr), intent(in) :: res, err, epsabs, scale

!     real(pr) :: eff_eps
!     integer :: digits
!     character(len=64) :: fmt

!     ! ============================
!     ! Эффективная точность
!     ! ============================
!     eff_eps = epsabs * scale

!     if (eff_eps <= 0._pr) then
!         digits = 6
!     else
!         digits = max(0, int(-log10(eff_eps)))
!     end if

!     ! Ограничим (чтобы не было безумия)
!     digits = min(digits, 15)

!     ! ============================
!     ! Формируем формат
!     ! ============================
!     write(fmt, '("(A,F20.",I0,",A,F20.",I0,",A)")') digits, digits

!     ! ============================
!     ! Вывод
!     ! ============================
!     write(*,fmt) trim(label)//' : ', res*scale, ' ± ', err*scale, ' x10^6'

! end subroutine print_result

subroutine Symb(state)
    implicit none
    type(q_state) , intent(in) :: state

    integer :: n, l, j2
    character :: shell
    character(len=10) :: j_str
    character(len=20) :: term
    associate( nr => state%nr, kappa => state%k, Znuc => state%Znuc )

    ! Определяем главное квантовое число n
    n = nr + abs(kappa)

    ! Определяем орбитальное квантовое число l и символ оболочки
    if (kappa > 0) then
        l = kappa
    else
        l = -kappa - 1
    end if

    ! Определяем символ электронной оболочки
    select case(l)
        case(0)
            shell = 's'
        case(1)
            shell = 'p'
        case(2)
            shell = 'd'
        case(3)
            shell = 'f'
        case(4)
            shell = 'g'
        case(5)
            shell = 'h'
        case(6)
            shell = 'i'
        case default
            shell = '?'
    end select

    ! Определяем полный угловой момент j
    if (kappa > 0) then
        j2 = 2 * kappa - 1  ! j = l - 1/2
    else
        j2 = -2 * kappa - 1 ! j = l + 1/2
    end if

    ! Формируем строку для j
    if (mod(j2, 2) == 0) then
        write(j_str, '(I0)') j2/2
    else
        write(j_str, '(I0,"/2")') j2
    end if

    write(term, '(I0,A1,"_",A)') n, shell, trim(j_str)


   write(*, '(A)') trim(term)
   write(*, '("Z = ", I0)') int(Znuc)
   end associate

end subroutine Symb

    
end module others
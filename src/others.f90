module others
    use quantum_state
    implicit none
    integer(8), save :: start_time = 0
    integer(8), save :: rate = 0
    
contains

subroutine print_result(label, unit_no, state, res, err, epsabs, scale)
    implicit none

    character(len=*), intent(in) :: label
    integer, intent(in) :: unit_no
    type(q_state), intent(in) :: state
    real(pr), intent(in) :: res, err, epsabs, scale

    real(pr) :: eff_eps, actual_err, scaled_res
    integer :: digits, power
    character(len=64) :: fmt
    character(len=16) :: scale_str
    character(len=64) :: symb_str

    actual_err = abs(err * scale)
    scaled_res = res * scale

    ! Формируем обозначение терма и зарядовое число
    call Symb(state, symb_str)

    ! Определяем точность
    eff_eps = max(actual_err, abs(epsabs * scale))

    if (eff_eps > 1e-20_pr) then
        digits = max(0, int(-log10(eff_eps)) + 1)
    else
        digits = 15
    end if

    digits = min(digits, 15)

    ! Формируем строку множителя
    if (abs(scale - 1.0_pr) < 1e-9_pr) then
        scale_str = ""
    else
        power = nint(log10(scale))
        write(scale_str, '(" x10^", I0)') power
    end if

    ! Формируем формат результата
    write(fmt, '("(A, F0.", I0, ", A, F0.", I0, ", A)")') &
        digits, digits

    ! Выводим терм, зарядовое число и результат
    write(unit_no, '(A)') trim(symb_str)
    write(unit_no, fmt) trim(label)//' : ', scaled_res, &
        ' ± ', actual_err, trim(scale_str)

    flush(unit_no)

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

subroutine Symb(state, symb_str)
    implicit none
    type(q_state), intent(in) :: state
    character(len=*), intent(out) :: symb_str

    integer :: n, l, j2
    character :: shell
    character(len=10) :: j_str
    character(len=20) :: term
    character(len=16) :: z_str

    associate(nr => state%nr, kappa => state%k, Znuc => state%Znuc)

        ! Главное квантовое число n
        n = nr + abs(kappa)

        ! Орбитальное квантовое число l
        if (kappa > 0) then
            l = kappa
        else
            l = -kappa - 1
        end if

        ! Символ оболочки
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

        ! Удвоенное значение полного углового момента j
        if (kappa > 0) then
            j2 = 2 * kappa - 1
        else
            j2 = -2 * kappa - 1
        end if

        ! Строковое представление j
        if (mod(j2, 2) == 0) then
            write(j_str, '(I0)') j2 / 2
        else
            write(j_str, '(I0,"/2")') j2
        end if

        ! Спектральное обозначение терма
        write(term, '(I0,A1,"_",A)') n, shell, trim(j_str)

        ! Зарядовое число
        write(z_str, '(I0)') int(Znuc)

        ! Формируем общую строку
        write(symb_str, '(A,"  Z = ",A)') trim(term), trim(z_str)

    end associate

end subroutine Symb

    subroutine timer(mode)
    implicit none
    character(len=*), intent(in) :: mode
    integer :: finish, minutes, seconds
    real :: elapsed

    if (mode == "start") then
      call system_clock(start_time, rate)

    else if (mode == "stop") then
      call system_clock(finish)

      elapsed = real(finish - start_time) / real(rate)
      minutes = int(elapsed / 60.0)
      seconds = int(elapsed - minutes * 60)

      write(*,'(A, I0, A, I0, A)') "Execution time: ", minutes, " min ", seconds, " sec"
    end if
    end subroutine timer  
    

    subroutine show_status(step_name, pulse)
        character(len=*), intent(in) :: step_name
        integer, intent(in), optional :: pulse ! Опционально для анимации точек
        character(len=4) :: dots
        integer :: i, num_dots
        
        dots = "    "
        if (present(pulse)) then
            ! Цикл анимации: . -> .. -> ... -> ....
            num_dots = mod(pulse - 1, 4) + 1
            do i = 1, num_dots
                dots(i:i) = "."
            end do
        else
            ! Если pulse нет, просто ставим три точки
            dots = "... "
        end if

        ! Используем \r (char(13)) для возврата каретки и перезаписи строки
        ! Добавляем пробелы в конце для затирки хвоста предыдущей строки
        write(*, '(A, A, A, A, A)', advance='no') &
            char(13), "[WORKING] ", step_name, dots, "    "
        
        call flush(6) 
    end subroutine show_status

    subroutine util_time_print ( Nunit, Nskip_before, Nskip_after )
!- - - - - - - - -
      implicit none
! .. Input / Output ..
      integer, intent(in) :: Nunit
      integer, intent(in) :: Nskip_before, Nskip_after
! .. Local ..
      character(8)  :: ddate
      character(10) :: ttime
      integer :: i
!- - - - - - - - -

        do i = 1, Nskip_before, 1
          write(unit=Nunit,fmt=*)
        enddo  ! i

        call date_and_time(DATE=ddate,TIME=ttime)
        write(unit=Nunit,fmt='(1x,a20)') '===================='
        write(unit=Nunit,fmt='(1x,a4,a1,a2,a1,a2,a2,a2,a1,a2,a1,a2)')                   &
                                     ddate(1:4),'.',ddate(5:6),'.',ddate(7:8),'  ',     &
                                     ttime(1:2),':',ttime(3:4),':',ttime(5:6)
        write(unit=Nunit,fmt='(1x,a20)') '===================='

        do i = 1, Nskip_after, 1
          write(unit=Nunit,fmt=*)
        enddo  ! i

       return

      end subroutine util_time_print

end module others
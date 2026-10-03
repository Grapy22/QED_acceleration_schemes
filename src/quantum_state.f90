module quantum_state
    use parameters
    use dirac_pnt
    implicit none

    type :: q_state
    real(pr) :: Znuc , en , j , alphaZ
    integer :: k , nr , l , l_
    contains
        procedure :: init => init_state
        procedure :: print => print_state
    end type


    contains

   subroutine init_state(this, k , nr , Znuc )
       class(q_state), intent(inout) :: this
       integer, intent(in) :: k , nr
       real(pr), intent(in) :: Znuc
       integer :: jjx, llx

       this%k = k
       this%nr = nr
       this%Znuc = Znuc
       this%alphaZ = alpha * Znuc
       call ka2jl (  k, jjx, llx )
       this%l = llx / 2
       this%j = real(jjx , pr) / 2._pr
       call ka2jl ( -k, jjx, llx )
       this%l_ = llx / 2

       this%en = en_Dirac_pnt ( nr, k, Znuc )

   end subroutine init_state

subroutine print_state(this , unit_no)
    class(q_state) , intent(in) :: this
    integer, intent(in) :: unit_no

    integer :: n, l, j2
    character :: shell
    character(len=10) :: j_str
    character(len=20) :: term
    associate( nr => this%nr, kappa => this%k, Znuc => this%Znuc )

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

   write(unit_no, '(A, A, A, I0, A)') '[', trim(term), ' | Z=', int(Znuc), ']'
   end associate

end subroutine print_state
end module quantum_state



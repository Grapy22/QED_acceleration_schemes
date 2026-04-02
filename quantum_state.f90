module quantum_state
    use parameters
    use dirac_pnt
    implicit none

    type :: q_state
    real(pr) :: Znuc , en , j , alphaZ
    integer :: k , nr , l , l_
    contains
        procedure :: init => init_state
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
end module quantum_state



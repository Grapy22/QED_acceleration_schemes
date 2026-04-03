program main
    use parameters
    use quantum_state
    use pot0_2004
    use pot1_2004
    implicit none

    integer :: k = -1, nr = 0
    real(pr) :: Znuc = 10._pr
    type(q_state) :: state1 , state2

    call init_state(state1, k, nr, Znuc)
    call init_state(state2, 1, 1, 92._pr)

    call pot1(state1 , vr = "vr" , vert = "vert" , red = "red")
    !call pot0_SC (state1 , vr = "vr" , vert = "vert" , red = "red")
    !call pot0_SC (state2 , vr = "vr" , vert = "vert" , red = "red")
    
    
    
    


end program main
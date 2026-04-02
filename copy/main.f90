program main
    use parameters
    use quantum_state
    use pot0_2004
    implicit none

    integer :: k = -1, nr = 0
    real(pr) :: Znuc = 10._pr
    type(q_state) :: state1

    call init_state(state1, k, nr, Znuc)

    call pot0 (state1 , vr = "vr" , vert = "vert" , red = "red")
    
    
    
    


end program main
program main
    use parameters
    use quantum_state
    use pot0_2004
    use pot1_2004
    use omp_lib
    use others
    implicit none
  
    integer :: u1, u2
    type(q_state) :: state1 , state2 , state3 , state4
    write(*,*) "STARTED"
    call util_time_print(6 , 0 , 1)

    open(newunit=u1, file="results_calc/output_state1.txt", status="replace", action="write")
    open(newunit=u2, file="results_calc/output_state2.txt", status="replace", action="write")

    call init_state(state1, -1, 0, 92._pr)
    call init_state(state2, -1, 1, 92._pr)
    
    !call pot0_SC (state1 , u1 , vr = "vr" , vert = "vert" , red = "red")
    !call pot0_SC (state2 , u2 , vr = "vr" , vert = "vert" , red = "red")
    ! call pot1(state1 , vr = "vr" , vert = "vert" , red = "red")
    !call pot0_SC (state1 , vr = "vr" , vert = "vert" , red = "red")
    !call pot0_SC (state2 , vr = "vr" , vert = "vert" , red = "red")
    
    !$omp parallel shared(state1, state2, u1, u2)
        !$omp sections
            !$omp section
                call pot1_SC(state1 , u1 , vert = "vert" , red = "red" , vr = "vr" )
                
            !$omp section
                call pot1_SC(state2 , u2 , vert = "vert" , red = "red" , vr = "vr"  )
                
        !$omp end sections
    !$omp end parallel

    close(u1)
    close(u2)
    
    call util_time_print(6 , 1 , 0)
    write(*,*) "DONE"
end program main
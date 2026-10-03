program main
    use parameters
    use quantum_state
    use pot0_2004
    use pot1_2004
    use omp_lib
    use others
    use dirac_pnt
    implicit none
  
    integer :: u1, u2, u3, u4
    type(q_state) :: state1 , state2 , state3 , state4
    write(*,*) "STARTED"
    call util_time_print(6 , 0 , 1)

    open(newunit=u1, file="results_calc/output_state1.txt", status="replace", action="write")
    open(newunit=u2, file="results_calc/output_state2.txt", status="replace", action="write")
    open(newunit=u3, file="results_calc/output_state3.txt", status="replace", action="write")
    open(newunit=u4, file="results_calc/output_state4.txt", status="replace", action="write")

    call init_state(state1, -1, 0, 10._pr)
    call init_state(state2, -1, 1, 10._pr)
    call init_state(state3, 1, 1, 10._pr)
    call init_state(state4, -2, 0, 10._pr)
    !call pot0(state4 , u4, vr = "vr" , vert = "vert" , red = "red")
    
    !call pot0_SC (state2 , u2 , vr = "vr" , vert = "vert" , red = "red")
    ! call pot1(state1 , vr = "vr" , vert = "vert" , red = "red")
    !call pot1_sc (state1 , vr = "vr" , vert = "vert" , red = "red")
    !call pot0_SC (state2 , vr = "vr" , vert = "vert" , red = "red")
   ! call pot1(state1 , u1, vr = "vr" , vert = "vert" , red = "red")
    call pot1(state2 , u2, vr = "vr" , vert = "vert" , red = "red")
    call pot1(state3 , u3, vr = "vr" , vert = "vert" , red = "red")
    call pot1(state4 , u4, vr = "vr" , vert = "vert" , red = "red")


    close(u1)
    close(u2)
    close(u3)
    close(u4)
    call util_time_print(6 , 1 , 0)
    write(*,*) "DONE"
end program main

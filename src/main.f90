program main
    use parameters
    use quantum_state
    use pot0_2004
    use pot1_2004
    use omp_lib
    use others
    use dirac_pnt
    use pot2_via_pot0
    use SE_e
    implicit none

    integer :: u1, u2, u3, u4
    type(q_state) :: state1 , state2 , state3 , state4
    write(*,*) "STARTED"
    call util_time_print(6 , 0 , 1)
    !call cubacores(0, 0)

    open(newunit=u1, file="results_calc/output_state1.txt", status="replace", action="write")
    open(newunit=u2, file="results_calc/output_state2.txt", status="replace", action="write")
    open(newunit=u3, file="results_calc/output_state3.txt", status="replace", action="write")
    open(newunit=u4, file="results_calc/output_state4.txt", status="replace", action="write")

    call init_state(state1, -1, 0, 10._pr)
    

    call init_state(state2, -1, 1, 10._pr)
    call init_state(state3, 1, 1, 10._pr)
    call init_state(state4, -2, 0, 10._pr)

    call SE_contrib(state1, u1, contrib="1pot")
    call SE_contrib(state2, u2, contrib="1pot")
    call SE_contrib(state3, u3, contrib="1pot")
    call SE_contrib(state4, u4, contrib="1pot")

    ! call SE_contrib(state1, u1, contrib="2pot_via_1pot")
    ! call SE_contrib(state2, u2, contrib="2pot_via_1pot")
    ! call SE_contrib(state3, u3, contrib="2pot_via_1pot")
    ! call SE_contrib(state4, u4, contrib="2pot_via_1pot")


    
! сделать человеское распаралеливание 
   ! call pot1_SC(state1 , u2, red = "red")
    !call pot2_via_1pot(state1, u1, vert = "vert")
   ! call pot1_SC(state1 , u1, red = "red") 
    !call pot1_SC(state2 , u2, vr = "vr" , vert = "vert" , red = "red")
    !call pot1_SC(state3 , u3, vr = "vr" , vert = "vert" , red = "red")
    !call pot2_via_0pot(state4 , u4, vr = "vr" , vert = "vert" , red = "red")
    !call pot0(state4 , u4, vr = "vr" , vert = "vert" , red = "red")

    !call pot0_SC (state2 , u2 , vr = "vr" , vert = "vert" , red = "red")
    ! call pot1(state1 , vr = "vr" , vert = "vert" , red = "red")
    !call pot1_sc (state1 , vr = "vr" , vert = "vert" , red = "red")
    !call pot0_SC (state2 , vr = "vr" , vert = "vert" , red = "red")
   ! call pot1(state1 , u1, vr = "vr" , vert = "vert" , red = "red")
    !call pot1(state2 , u2, red = "red")
   ! call pot1(state3 , u3, vr = "vr" , vert = "vert" , red = "red")
   ! call pot1(state4 , u4, vr = "vr" , vert = "vert" , red = "red")

!    !$omp parallel
!        !$omp sections
!             !$omp section
!             call pot0(state1 , u1, vr = "vr" , vert = "vert" , red = "red")
!             close(u1)
!             !$omp section
!             call pot0(state2 , u2, vr = "vr" , vert = "vert" , red = "red")
!             close(u2)
!             !$omp section
!             call pot0(state3 , u3, vr = "vr" , vert = "vert" , red = "red")
!             close(u3)
!             !$omp section
!             call pot0(state4 , u4, vr = "vr" , vert = "vert" , red = "red")
!             close(u4)
!        !$omp end sections
!    !$omp end parallel

    call util_time_print(6 , 1 , 0)
    write(*,*) "DONE"
end program main

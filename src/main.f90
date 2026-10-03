program main
    use parameters
    use quantum_state
    use pot0_2004
    use pot1_2004
    use omp_lib
    use others
    use dirac_pnt
    use pot2_via_pot0
    implicit none

    integer :: u1, u2, u3, u4, nrx, kax
    type(q_state) :: state1 , state2 , state3 , state4
    real(pr) :: p , X, Y, Zx
    write(*,*) "STARTED"
    call util_time_print(6 , 0 , 1)
    !call cubacores(0, 0)
    ! 
    ! experimenral block
    p = 2.3_pr
    Zx = 10.0_pr
    nrx = 1
    kax = -1
    call dirac_pnt_XY(nrx, kax, Zx, p, X, Y, -1._pr)
    write(*,*) X, Y
    write(*,*) X / Y
    write(*,*) Y / X
    call init_state(state1, -1, 0, 10._pr)
    ! experimental block

!    ! write(*,*) "HERE", red_pot2_via_0pot(2.3_pr, state1) / red_pot2_via_0pot(4.5_pr, state1)
! !    write(*,*) "HERE", red_pot2_via_1pot (0.3_pr , 0.7_pr , 1.8_pr , 2.3_pr , 0.6_pr , state1)
!     write(*,*) "HERE", red_pot1_SC (0.3_pr , 0.7_pr , 1.8_pr , 2.3_pr , 0.6_pr , state1) / & 
!     red_pot1_SC (0.2_pr , 0.4_pr , 2.8_pr , 5.3_pr , 0.2_pr , state1) 
! !    write(*,*) "HERE", red_pot2_via_1pot (0.3_pr , 0.7_pr , 1.8_pr , 2.3_pr , 0.6_pr , state1) /&
! write(*,*) "HERE", red_pot1_SC (0.3_pr , 0.7_pr , 1.8_pr , 2.3_pr , 0.6_pr , state1)
! write(*,*) "HERE", red_pot1_SC (0.2_pr , 0.4_pr , 2.8_pr , 5.3_pr , 0.2_pr , state1) 










    open(newunit=u1, file="results_calc/output_state1.txt", status="replace", action="write")
    open(newunit=u2, file="results_calc/output_state2.txt", status="replace", action="write")
    open(newunit=u3, file="results_calc/output_state3.txt", status="replace", action="write")
    open(newunit=u4, file="results_calc/output_state4.txt", status="replace", action="write")

    call init_state(state1, -1, 0, 10._pr)
    call init_state(state2, -2, 0, 10._pr)
    call init_state(state3, -2, 0, 10._pr)
    call init_state(state4, -3, 0, 10._pr)
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

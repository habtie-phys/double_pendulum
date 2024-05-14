program double_pendulum_test
    use double_pendulum, only: m1, l1, m2, l2
    use double_pendulum, only: rk4_update
    use ogpf
    implicit none
    type(gpf) :: gp
    real, parameter :: d2r = 3.14159265/180.
    ! set number of iterations
    integer, parameter :: n = 2000
    ! set initial angular positions for each pendulum
    real, parameter :: theta1_0 = 15., theta2_0 = 5.
    ! set initial angular velocities for each pendulum
    real, parameter :: w1_0 = 0.,  w2_0 = 0. 
    ! set initial state of the system
    real :: S0(4) = [d2r*theta1_0, d2r*theta2_0, w1_0, w2_0]
    real :: S(4, n), t(n)
    real :: x1(n), x2(n), y1(n), y2(n)
    ! set step size for the rk4 update
    real :: h = 1e-3
    integer :: i
    !mass and length of the upper pendulum
    m1 = 0.5 
    l1 = 0.3
    !mass and length of the lower pendulum
    m2 = 1.
    l2 = 0.6
    ! 
    t(1) = 0
    call rk4_update(h,S0, S(:, 1))
    do i = 1, n-1
        t(i+1) = t(i)+h
        call rk4_update(h,S(:, i), S(:, i+1))
    end do
    x1 = l1*sin(S(1, :))
    x2 = l1*sin(S(1, :))+l2*sin(S(2, :))
    ! 
    y1 = -l1*cos(S(1, :))
    y2 = -l1*cos(S(1, :))-l2*cos(S(2, :))
    !
    call gp%xlabel(' Time (s)')
    call gp%ylabel(' y')
    call gp%title('lower pendulum') 
    call gp%plot(real(a=t, kind=wp), real(a=y2, kind=wp), 'with lines lt 1 lw 2')
end program double_pendulum_test
module double_pendulum
    implicit none
    real, parameter :: g = 9.81
    real :: m1, l1, m2, l2
    contains
    function f(S)
        implicit none
        real, intent(in) :: S(4)
        real :: f(4)
        real :: a1, a2
        real :: b1, b2
        real :: c1, c2, c20, c21, c22
        real :: t1, t2, p1, p2
        ! 
        t1 = S(1)
        t2 = S(2)
        p1 = S(3)
        p2 = S(4)
        ! 
        c1 =p1*p2*sin(t1-t2)/(l1*l2*(m1+m2*sin(t1-t2)**2))
        c20 = l2**2*m2*p1**2+l1**2*(m1+m2)*p2**2
        c21 = l1*l2*m2*p1*p2*cos(t1-t2)
        c22 = 2*l1**2*l2**2*(m1+m2*sin(t1-t2)**2)
        c2 = (c20-c21)/c22*sin(2*t1-2*t2)
        ! 
        a1 = l2*p1-l1*p2*cos(t1-t2)
        a2 = l1**2*l2*(m1+m2*sin(t1-t2)**2)
        ! 
        b1 = l1*(m1+m2)*p2-l2*m2*p1*cos(t1-t2)
        b2 = l1*l2**2*m2*(m1+m2*sin(t1-t2)**2)
        f(1) = a1/a2
        f(2) = b1/b2
        f(3) = -(m1+m2)*g*l1*sin(t1)-c1+c2
        f(4) = -m2*g*l2*sin(t2)+c1-c2
        return
    end function f
    ! 
    subroutine rk4_update(h, S0, S)
        implicit none
        real, intent(in) :: h
        real, intent(in) :: S0(4)
        real, intent(out) :: S(4)
        real :: k1(4), k2(4), k3(4), k4(4)
        k1 = f(S0)
        k2 = f(S0+0.5*h*k1)
        k3 = f(S0+0.5*h*k2)
        k4 = f(S0+h*k3)
        S = S0+h*(k1 + 2*k2 + 2*k3 + k4)/6
    end subroutine rk4_update
end module double_pendulum


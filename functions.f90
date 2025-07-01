module functions

    use constantes

    implicit none 

    contains 

    ! Condition initiale u(0, x), differents cas possibles
    function u_init(x, cas) result(res)

        real(kind=PR), intent(in) :: x
        integer, intent(in) :: cas
        real(kind=PR) :: res 

        SELECT CASE(cas)

        CASE(1)
            res = sin(pi*x)
        CASE(2)
            res = sin(pi*x)
        end SELECT

    end function

    ! Condition de bord (gauche ou droite selon le signe de a), depend du cas choisi pour etre C infini
    function u_bound(t, Lx, Rx, a, cas, C, tau) result(res)

        real(kind=PR), intent(in) :: t, Lx, a, Rx, C, tau
        integer, intent(in) :: cas
        real(kind=PR) :: res, I

        SELECT CASE(cas)

        CASE(1) 

            if (a > 0.0_PR) then 

                res = C * (1 - exp(-t / tau)) + u_init(Lx - a*t, cas) * exp(-t / tau)

            else if (a < 0.0_PR) then

                res = C * (1 - exp(-t / tau)) + u_init(Rx - a*t, cas) * exp(-t / tau)

            end if 

        CASE(2) 

            if (a > 0.0_PR) then 

                I = tau / (1.0_PR + (pi*a*tau)**2) * & 
                    (-exp(-t / tau) * (1.0_PR / tau * cos(pi * (Lx - a*t)) + pi * a * sin(pi * (Lx - a*t))) & 
                    + (1.0_PR / tau * cos(pi * Lx) + pi * a * sin(pi * Lx)))

                res = u_init(Lx - a*t, cas) * exp(-t / tau) + I

            else if (a < 0.0_PR) then 

                I = tau / (1.0_PR + (pi*a*tau)**2) * & 
                    (-exp(-t / tau) * (1.0_PR / tau * cos(pi * (Rx - a*t)) + pi * a * sin(pi * (Rx - a*t))) & 
                    + (1.0_PR / tau * cos(pi * Rx) + pi * a * sin(pi * Rx)))

                res = u_init(Rx - a*t, cas) * exp(-t / tau) + I

            end if 

        end SELECT

    end function

    ! Solution exacte à (t, x), depend du cas, et du signe de a
    function sol_exacte(x, t, Lx, Rx, a, cas, C, tau) result(res)

        real(kind=PR), intent(in) :: x, t, Lx, Rx, a, C, tau
        integer, intent(in) :: cas
        real(kind=PR) :: res, I

        SELECT CASE(cas)

        CASE(1)

        if (a > 0.0_PR) then 

            if (x-a*t >= Lx) then 

                res = u_init(x - a*t, cas) * exp(-t / tau) + C * (1.0_PR - exp(-t / tau))
            
            else if (x-a*t < Lx) then 

                res = u_bound(t - (x - Lx)/a, Lx, Rx, a, cas, C, tau) * exp (-(x - Lx) / (a * tau)) + &
                            C * (1.0_PR - exp (-(x - Lx) / (a * tau)))

            end if 
        
        else if (a < 0.0_PR) then 

            if (x-a*t <= Rx) then 

                res = u_init(x - a*t, cas) * exp(- t / tau) + C * (1.0_PR - exp(-t / tau))
            
            else if (x-a*t > Rx) then 

                res = u_bound(t - (x - Rx)/a, Lx, Rx, a, cas, C, tau) * exp (-(x - Rx) / (a * tau)) + &
                            C * (1.0_PR - exp (-(x - Rx) / (a * tau)))

            end if 

        else if (a == 0.0_PR) then 

            res = u_init(x, cas) * exp(-t / tau) + C * (1.0_PR - exp(-t / tau))

        end if 

        CASE(2)

        if (a > 0.0_PR) then 

             if (x-a*t >= Lx) then 

                I = tau / (1.0_PR + (pi*a*tau)**2) * & 
                    ((1.0_PR / tau * cos(pi * x) + pi * a * sin(pi * x)) - & 
                    exp(-t / tau) * (1.0_PR / tau * cos(pi * (x - a*t)) + pi * a * sin(pi * (x - a*t))))

                res = u_init(x - a*t, cas) * exp(-t / tau) + I 
            
            else if (x-a*t < Lx) then 

                I = tau / (1.0_PR + (pi*a*tau)**2) * & 
                    (1.0_PR / tau * cos(pi * x) + pi * a * sin(pi * x) - &
                    exp(- (x - Lx) / (a * tau)) * (1.0_PR / tau * cos(pi * Lx) + pi * a * sin(pi * Lx)))

                res = u_bound(t - (x - Lx)/a, Lx, Rx, a, cas, C, tau) * exp (-(x - Lx) / (a * tau)) + I

            end if 

        else if (a < 0.0_PR) then 

            if (x-a*t <= Rx) then 

                I = tau / (1.0_PR + (pi*a*tau)**2) * & 
                    ((1.0_PR / tau * cos(pi * x) + pi * a * sin(pi * x)) - & 
                    exp(-t / tau) * (1.0_PR / tau * cos(pi * (x - a*t)) + pi * a * sin(pi * (x - a*t))))

                res = u_init(x - a*t, cas) * exp(-t / tau) + I 
            
            else if (x-a*t > Rx) then 

                I = tau / (1.0_PR + (pi*a*tau)**2) * & 
                    (1.0_PR / tau * cos(pi * x) + pi * a * sin(pi * x) - &
                    exp(- (x - Rx) / (a * tau)) * (1.0_PR / tau * cos(pi * Rx) + pi * a * sin(pi * Rx)))

                res = u_bound(t - (x - Rx)/a, Lx, Rx, a, cas, C, tau) * exp (-(x - Rx) / (a * tau)) + I

            end if 

        end if 
        
        END SELECT
 
    end function

    function fact(k) result(res)

        integer, intent(in) :: k 
        real(kind=PR) :: res
        integer :: i

        res = 1.0_PR

        if (k >= 1) then 

            do i = 1, k 

                res = res * i 

            end do 

        end if 

    end function    

    function Sn_i(x, cas, C) result(res)

        real(kind=PR), intent(in) :: x, C
        integer, intent(in) :: cas 
        real(kind=PR) :: res 

        SELECT CASE(cas)

        CASE(1)  

            res = C 

        CASE(2)

            res = cos(pi*x)
        
        end SELECT

    end function


end module
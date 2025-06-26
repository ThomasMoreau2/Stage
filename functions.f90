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

        end SELECT

    end function

    ! Condition de bord (gauche ou droite selon le signe de a), depend du cas choisi pour etre C infini
    function u_bound(t, Lx, Rx, a, cas, C, tau) result(res)

        real(kind=PR), intent(in) :: t, Lx, a, Rx, C, tau
        integer, intent(in) :: cas
        real(kind=PR) :: res

        SELECT CASE(cas)

        CASE(1) 

            if (a > 0.0_PR) then 

                res = C * (1 - exp(-t / tau)) + u_init(Lx - a*t, cas) * exp(-t / tau)

            else if (a < 0.0_PR) then

                res = C * (1 - exp(-t / tau)) + u_init(Rx - a*t, cas) * exp(-t / tau)

            end if 
        
        end SELECT

    end function

    ! Solution exacte à (t, x), depend du cas, et du signe de a
    function sol_exacte(x, t, Lx, Rx, a, cas, C, tau) result(res)

        real(kind=PR), intent(in) :: x, t, Lx, Rx, a, C, tau
        integer, intent(in) :: cas
        real(kind=PR) :: res 

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

end module
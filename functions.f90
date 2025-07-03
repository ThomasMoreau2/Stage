module functions 

    use constantes
    use Legendre    

    implicit none 

    contains 

    function Maxwell(rho, u, T, v) result(res)

        real(kind=PR), intent(in) :: rho, u, T, v 
        real(kind=PR) :: res 

        res = rho / sqrt(2.0_pr * pi * T) * exp(-(v - u)**2 / (2.0_pr * T))

    end function


    function calcul_rho(nx, nv, p, dv, alpha) result(rho)

        integer, intent(in) :: nx, nv, p
        real(kind=PR), intent(in) :: dv
        real(kind=PR), dimension(p*nx, 0:nv), intent(in) :: alpha
        real(kind=PR), dimension(0:nx, 0:nv) :: f
        real(kind=PR), dimension(0:nx) :: rho
        integer :: i, l

        f = calculate_f(nx, nv, p, alpha)

        rho = 0.0_pr

        do i = 0, nx 

            do l = 0, nv  

                rho(i) = rho(i) + f(i, l) * dv

            end do 

        end do 
     
    end function 


    function calcul_u(nx, nv, p, dv, v, alpha) result(u)

        integer, intent(in) :: nx, nv, p
        real(kind=PR), intent(in) :: dv
        real(kind=PR), dimension(0:nv), intent(in) :: v
        real(kind=PR), dimension(p*nx, 0:nv), intent(in) :: alpha
        real(kind=PR), dimension(0:nx, 0:nv) :: f
        real(kind=pr), dimension(0:nx) :: rho_u, rho, u
        integer :: i, l

        f = calculate_f(nx, nv, p, alpha)
        rho = calcul_rho(nx, nv, p, dv, alpha)

        rho_u = 0.0_pr

        do i = 0, nx 

            do l = 0, nv  

                rho_u(i) = rho_u(i) + f(i, l) * v(l) * dv

            end do 

            u(i) = rho_u(i) / rho(i)

        end do 


    end function 


    function calcul_T(nx, nv, p, dv, v, alpha) result(T)

        integer, intent(in) :: nx, nv, p
        real(kind=PR), intent(in) :: dv
        real(kind=PR), dimension(0:nv), intent(in) :: v
        real(kind=PR), dimension(p*nx, 0:nv), intent(in) :: alpha
        real(kind=PR), dimension(0:nx, 0:nv) :: f
        real(kind=pr), dimension(0:nx) :: T, rho, u, E
        integer :: i, l

        f = calculate_f(nx, nv, p, alpha)
        rho = calcul_rho(nx, nv, p, dv, alpha)
        u = calcul_u(nx, nv, p, dv, v, alpha) 

        E = 0.0_PR

        do i = 0, nx 

            do l = 0, nv  

                E(i) = E(i) + 1.0_PR / 2.0_PR * f(i, l) * v(l)**2 * dv

            end do 

            T(i) = (E(i) - 1.0_pr / 2.0_pr * rho(i) * u(i)**2) / (1.0_pr / 2._pr * rho(i))

        end do 
          
    end function 


    function calcul_tau(nx, nv, p, dv, alpha) result(tau)

        integer, intent(in) :: nx, nv, p
        real(kind=PR), intent(in) :: dv
        real(kind=PR), dimension(p*nx, 0:nv), intent(in) :: alpha
        real(kind=pr), dimension(0:nx) :: rho, tau
        integer :: i

        rho = calcul_rho(nx, nv, p, dv, alpha)

        do i = 0, nx 

            tau(i) = 1.0_PR / rho(i)

        end do 

    end function


    function u_init(x, rho_g, rho_d, T_g, T_d, v) result(res)

        real(kind=PR), intent(in) :: x, rho_g, rho_d, T_g, T_d, v
        real(kind=PR) :: res 

        if (x <= 0.0_PR) then 

            res = Maxwell(rho_g, 0.0_PR, T_g, v) 

        else 

            res = Maxwell(rho_d, 0.0_PR, T_d, v) 

        end if 

    end function

    function u_bound(rho_g, u_g, T_g, rho_d, u_d, T_d, v) result(res)

        real(kind=PR), intent(in) :: rho_g, u_g, T_g, rho_d, u_d, T_d, v
        real(kind=PR) :: res 
        
        if (v > 0.0_PR) then 

            res = Maxwell(rho_g, u_g, T_g, v)

        else if (v < 0.0_PR) then 

            res = Maxwell(rho_d, u_d, T_d, v)

        end if 

    end function


end module
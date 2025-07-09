module sol_exacte 

    use constantes

    implicit none 

    contains
    
    ! Fonction calculant rho
    function rho_ex(t, x, rho_g, rho_d, T_g, T_d) result(res)

        real(PR), intent(in) :: t, x, rho_G, rho_D, T_G, T_D
        real(PR) :: res, xg, xd

        xg = x/(t*sqrt(2.0_PR*T_g))
        xd = x/(t*sqrt(2.0_PR*T_d))

        res = 0.5_pr * (rho_g*erfc(xg) + rho_d*erfc(-xd))

    end function 

    ! Fonction calculant rho*u
    function rho_u_ex(t, x, rho_g, rho_d, T_g, T_d) result(res)

        real(PR), intent(in) :: t, x, rho_G, rho_D, T_G, T_D
        real(PR) :: res, xg, xd

        xg = x/(t*sqrt(2.0_PR*T_g))
        xd = x/(t*sqrt(2.0_PR*T_d))

        res = 0.5_PR * (rho_g*sqrt(2.0_PR*t_g)/sqrt(pi) *exp(-xg**2) - & 
            rho_d*sqrt(2.0_PR*t_d)/sqrt(pi) *exp(-xd**2))
    end function 

    ! Fonction calculant E
    function E_ex(t, x, rho_g, rho_d, T_g, T_d) result(res)

        real(PR), intent(in) :: t, x, rho_G, rho_D, T_G, T_D
        real(PR) :: res, xg, xd 

        xg = x/(t*sqrt(2.0_PR*T_g))
        xd = x/(t*sqrt(2.0_PR*T_d))

        res = 0.5_PR * (rho_g*t_g*( xg/sqrt(pi)*exp(-xg**2) + 0.5_PR*erfc(xg)) &
            + rho_d*t_d*( -xd/sqrt(pi)*exp(-xd**2) + 0.5_PR*erfc(-xd)))

    end function 

    ! Fonction calculant u 
    function u_ex(t, x, rho_g, rho_d, T_g, T_d) result(res)

        real(PR), intent(in) :: t, x, rho_G, rho_D, T_G, T_D
        real(PR) :: res

        res = rho_u_ex(t, x, rho_g, rho_d, T_g, T_d) / rho_ex(t, x, rho_g, rho_d, T_g, T_d)

    end function

    ! Fonction calculant T
    function T_ex(t, x, rho_g, rho_d, T_g, T_d) result(res)

        real(PR), intent(in) :: t, x, rho_G, rho_D, T_G, T_D
        real(PR) :: res

        res = (E_ex(t, x, rho_g, rho_d, T_g, T_d) - 1.0_PR / 2.0_PR * &
            rho_ex(t, x, rho_g, rho_d, T_g, T_d) * u_ex(t, x, rho_g, rho_d, T_g, T_d)**2) / &
            (1.0_PR / 2.0_PR * rho_ex(t, x, rho_g, rho_d, T_g, T_d))

    end function

end module 

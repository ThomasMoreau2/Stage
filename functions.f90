module functions

    use constantes

    implicit none 

    contains 

    ! Condition initiale u(0, x), differents cas possibles
    function u_init(x, cas) result(res)

        real(kind=PR), intent(in) :: x 
        integer, intent(in) :: cas
        real(kind=PR) :: res 

        select case(cas)
        case(1)
            res = cos(5.0_PR*pi*x)
        case(2)
            res = 1.0_PR/2.0_PR + 1.0_PR/2.0_PR*sin(2.0_PR*pi*x)
        case(3)
            res = 1.0_PR + x + x**2
        case(4)
            res = 1.0_PR + x + x**2 + x**3 
        case(5)
            res = 1.0_PR + x + x**2 + x**3 + x**4 
        case(6)
            res = cos(2.0_PR*pi*x)
        end select

    end function

    ! Condition de bord (gauche ou droite selon le signe de a), depend du cas choisi pour etre C infini
    function u_bound(t, Lx, Rx, a, cas) result(res)
        real(kind=PR), intent(in) :: t, Lx, a, Rx
        integer, intent(in) :: cas
        real(kind=PR) :: res

        if (a > 0) then
            res = u_init(Lx - a*t, cas)
        else
            res = u_init(Rx - a*t, cas)
        end if

    end function

    ! Solution exacte à (t, x), depend du cas, et du signe de a
    function sol_exacte(x, t, Lx, Rx, a, cas) result(res)
    real(kind=PR), intent(in) :: x, t, Lx, Rx, a
    integer, intent(in) :: cas
    real(kind=PR) :: res 

    if (a > 0) then 
        if (x - a*t > Lx) then 
            res = u_init(x - a*t, cas)
        else 
            res = u_bound(t - (x - Lx)/a, Lx, Rx, a, cas)
        end if 
    else 
        if (x - a*t < Rx) then 
            res = u_init(x - a*t, cas)
        else 
            res = u_bound(t - (x - Rx)/a, Lx, Rx, a, cas)
        end if 
    end if 
end function
    
end module
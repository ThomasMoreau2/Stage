module functions

    use constantes

    implicit none 

    contains 

    ! Condition initiale u(0, x), différents cas possibles
    function u_init(x, i) result(res)

        real(kind=PR), intent(in) :: x 
        integer, intent(in) :: i
        real(kind=PR) :: res 

        select case(i)
        case(1)
            res = cos(5._PR*pi*x)
        case(2)
            res = 1._PR + x
        case(3)
            res = 1._PR/2._PR + 1._PR/2._PR*sin(2._PR*pi*x)
        end select

    end function

    ! Condition de bord (gauche ou droite selon le signe de a), dépend du cas choisi pour être C infini
    function u_bound(t, Lx, a, i) result(res)
        real(kind=PR), intent(in) :: t, Lx, a
        integer, intent(in) :: i
        real(kind=PR) :: res

        if (a > 0) then
            res = u_init(a*t, i)
        else
            res = u_init(Lx-a*t, i)
        end if
    end function

    ! Solution exacte à (t, x), dépend du cas, et du signe de a
    function sol_exacte(x, t, Lx, a, i) result(res)

        real(kind=PR), intent(in) :: x, t, Lx, a
        integer, intent(in) :: i
        real(kind=PR) :: res 

        if (a>0) then 
            if (x-a*t > 0._PR) then 
                res = u_init(x-a*t, i)
            else 
                res = u_bound(t-x/a, Lx, a, i)
            end if 
        else 
            if (x-a*t < Lx) then 
                res = u_init(x-a*t, i)
            else 
                res = u_bound(t-(x-Lx)/a, Lx, a, i)
            end if 
        end if 

    end function
    
end module
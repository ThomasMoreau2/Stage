module functions

    use constantes

    implicit none 

    contains 

    function u_init(x) result(res)

        real(kind=PR), intent(in) :: x 
        real(kind=PR) :: res 

        res = cos(5._PR*pi*x)
    end function

    function u_L(t, Lx) result(res)

        real(kind=PR), intent(in) :: t, Lx
        real(kind=PR) :: res 

        if (a>0) then 
            res = cos(5._PR*pi*a*t)
        else 
            res = cos(5._PR*pi*(Lx-a*t))
        end if 
    end function

    function sol_exacte(x, t, Lx) result(res)

        real(kind=PR), intent(in) :: x, t, Lx
        real(kind=PR) :: res 

        if (a>0) then 
            if (x-a*t > 0._PR) then 
                res = u_init(x-a*t)
            else 
                res = u_L(t-x/a, Lx)
            end if 
        else 
            if (x-a*t < Lx) then 
                res = u_init(x-a*t)
            else 
                res = u_L(t+(x-Lx)/a, Lx)
            end if 
        end if 

    end function
    
end module
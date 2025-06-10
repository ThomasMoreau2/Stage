module Legendre

    use constantes

    implicit none 

    real(kind=PR), dimension(21), parameter :: weight = [ &
    ! n = 1
    2.0000000000000000_PR, &
    ! n = 2
    1.0000000000000000_PR, 1.0000000000000000_PR, &
    ! n = 3
    0.5555555555555556_PR, 0.8888888888888888_PR, 0.5555555555555556_PR, &
    ! n = 4
    0.3478548451374538_PR, 0.6521451548625461_PR, &
    0.6521451548625461_PR, 0.3478548451374538_PR, &
    ! n = 5
    0.2369268850561891_PR, 0.4786286704993665_PR, 0.5688888888888889_PR, &
    0.4786286704993665_PR, 0.2369268850561891_PR, &
    ! n = 6
    0.1713244923791704_PR, 0.3607615730481386_PR, 0.4679139345726910_PR, &
    0.4679139345726910_PR, 0.3607615730481386_PR, 0.1713244923791704_PR ]


    real(kind=PR), dimension(21), parameter :: points = [ &
    ! n = 1
     0.0000000000000000_PR, &
    ! n = 2
    -0.5773502691896257_PR,  0.5773502691896257_PR, &
    ! n = 3
    -0.7745966692414834_PR,  0.0000000000000000_PR,  0.7745966692414834_PR, &
    ! n = 4
    -0.8611363115940526_PR, -0.3399810435848563_PR, &
     0.3399810435848563_PR,  0.8611363115940526_PR, &
    ! n = 5
    -0.9061798459386640_PR, -0.5384693101056831_PR,  0.0000000000000000_PR, &
     0.5384693101056831_PR,  0.9061798459386640_PR, &
    ! n = 6
    -0.9324695142031521_PR, -0.6612093864662645_PR, -0.2386191860831969_PR, &
     0.2386191860831969_PR,  0.6612093864662645_PR,  0.9324695142031521_PR ]



    contains 

    function Leg(i, x) result(res)

        integer, intent(in) :: i
        real(kind=PR), intent(in) :: x
        real(kind=PR) :: res

        select case(i)
            case default
                res = 1._PR
            case(2)
                res = x
            case(3)
                res = 0.5_PR * (3.0_PR * x**2 - 1.0_PR)
            case(4)
                res = 0.5_PR * (5.0_PR * x**3 - 3.0_PR * x)
            case(5)
                res = (1.0_PR / 8.0_PR) * (35.0_PR * x**4 - 30.0_PR * x**2 + 3.0_PR)
            case(6)
                res = (1.0_PR / 8.0_PR) * (63.0_PR * x**5 - 70.0_PR * x**3 + 15.0_PR * x)
        end select
    end function

    function norme(i) result(res)

        integer, intent(in) :: i 
        real(kind=PR) :: res 
        
        res = 2._PR/(2._PR*(i-1)+1)

    end function

    function calculate_u(alpha, i_max, p, dx) result(u)
         
        real(kind=PR), dimension (:), intent(in) :: alpha
        integer, intent(in) :: i_max, p
        real(kind=PR), intent(in) :: dx
        real(kind=PR), dimension(i_max+1) :: u 
        integer :: i, l
        
        do i=1, i_max
            u(i) = 0._PR
            do l =1, p 
                u(i) = u(i) + alpha((i-1)*p + l)*Leg(l, ((i-1._PR)*dx-(i-0.5_PR)*dx)/(dx/2._PR))
            end do 
        end do 

        u(i_max+1) = 0._PR

        do l = 1, p 
            u(i_max+1) = u(i_max+1) + alpha((i_max-1)*p + l) * Leg(l, (i_max*dx - (i_max-0.5_PR)*dx)/(dx/2._PR))
        end do 
    end function 


end module 
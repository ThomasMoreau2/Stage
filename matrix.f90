module matrix

    use constantes
    use Legendre
    use functions

    implicit none 

    contains 

    function norme(i) result(res)

        integer, intent(in) :: i 
        real(kind=PR) :: res 
        
        res = 2._PR/(2._PR*(i-1)+1)

    end function


    function quad_1(i, j, lambda, p) result(res)
    
        integer, intent(in) :: i, j, p
        real(kind=PR), intent(in) :: lambda
        real(kind=PR) :: res 
        integer :: k, q 

        q = (p-1)/2+2

        res = 0._PR

        do k = q*(q-1)/2+1, q*(q-1)/2+q
            res = res + weight(k)*Leg(i,1._PR-(points(k)+1._PR)/lambda)*Leg(j, points(k))
        end do 

    end function

    function quad_2(i, j, lambda, p) result(res)
    
        integer, intent(in) :: i, j, p
        real(kind=PR), intent(in) :: lambda
        real(kind=PR) :: res 
        integer :: k, q

        q = (p-1)/2+2

        res = 0._PR

        do k = q*(q-1)/2+1, q*(q-1)/2+q
            res = res + weight(k)*Leg(i, -points(k))*Leg(j, (points(k)+1._PR)/lambda-1._PR)
        end do 

    end function

    function quad_3(i, j, lambda, p) result(res)
    
        integer, intent(in) :: i, j, p
        real(kind=PR), intent(in) :: lambda
        real(kind=PR) :: res 
        integer :: k, q

        q = (p-1)/2+2

        res = 0._PR

        do k = q*(q-1)/2+1, q*(q-1)/2+q
            res = res + weight(k)*Leg(i, (points(k)-1._PR)*(1._PR-1._PR/lambda)+1._PR-2._PR/lambda)* &
                    Leg(j, ((points(k)-1._PR)*(1._PR-1._PR/lambda)+1._PR))
        end do 

    end function

    function quad_4(i, j, lambda, p) result(res)
    
        integer, intent(in) :: i, j, p
        real(kind=PR), intent(in) :: lambda
        real(kind=PR) :: res 
        integer :: k, q 

        q = (p-1)/2+2

        res = 0._PR

        do k = q*(q-1)/2+1, q*(q-1)/2+q
            res = res + weight(k)*Leg(i,1._PR-(points(k)+1._PR)/abs(lambda))*Leg(j, -points(k))
        end do 

    end function

    function quad_5(i, j, lambda, p) result(res)
    
        integer, intent(in) :: i, j, p
        real(kind=PR), intent(in) :: lambda
        real(kind=PR) :: res 
        integer :: k, q

        q = (p-1)/2+2

        res = 0._PR

        do k = q*(q-1)/2+1, q*(q-1)/2+q
            res = res + weight(k)*Leg(i, points(k))*Leg(j, (points(k)+1._PR)/abs(lambda)-1._PR)
        end do 

    end function


    function make_L(lambda, p) result(L)

        real(kind=PR), intent(in) :: lambda
        integer, intent(in) :: p
        real(kind=PR), dimension(p, p) :: L
        integer :: i, j

        if (lambda>0) then 
            do i = 1, p 
                do j = 1, p
                    L(i, j) = 1._PR/norme(i)*quad_1(j, i, lambda, p)
                end do 
            end do 
        else 
             do i = 1, p 
                do j = 1, p
                    L(i, j) = 1._PR/norme(i)*quad_4(j, i, lambda, p)
                end do 
            end do 
        end if

    end function

    function make_M(lambda, p) result(M)

        real(kind=PR), intent(in) :: lambda
        integer, intent(in) :: p
        real(kind=PR), dimension(p, p) :: M
        integer :: i, j

        if (lambda>0) then 
            do i = 1, p
                do j = 1, p
                    M(i, j) = 1._PR/norme(i)/lambda*quad_2(j, i, lambda, p)
                end do 
            end do 
        else
            do i = 1, p 
                do j = 1, p
                    M(i, j) = 1._PR/norme(i)*quad_5(j, i, lambda, p)
                end do 
            end do 
        end if 
    end function

    function make_N(lambda, p) result(N)

        real(kind=PR), intent(in) :: lambda
        integer, intent(in) :: p
        real(kind=PR), dimension(p, p) :: N
        integer :: i, j

        do i = 1, p
            do j = 1, p
                N(i, j) = (1._PR-1._PR/lambda)/norme(i)*quad_3(j, i, lambda, p)
            end do 
        end do 

    end function


     function quad_init(i, j, p, dx) result(res)
    
        integer, intent(in) :: i, j, p
        real(kind=PR), intent(in) :: dx
        real(kind=PR) :: res 
        integer :: k, q 

        q = (p-1)/2+2

        res = 0._PR

        do k = q*(q-1)/2+1, q*(q-1)/2+q
            res = res + weight(k)*cos((dx/2._PR*points(k)+(i-0.5_PR)*dx)*5*pi)* &
                    Leg(j, points(k))
        end do 

    end function

    function quad_L(n, j, p, dt, Lx) result(res)
    
        integer, intent(in) :: j, p, n
        real(kind=PR), intent(in) :: dt, Lx
        real(kind=PR) :: res 
        integer :: k, q 

        q = (p-1)/2+2

        res = 0._PR

        do k = q*(q-1)/2+1, q*(q-1)/2+q
            res = res + weight(k)*Leg(j, points(k))*u_L(dt/2._PR*points(k)+(n+0.5_PR)*dt, Lx)
        end do 

    end function



end module 
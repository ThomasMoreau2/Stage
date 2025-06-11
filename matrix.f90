module matrix

    use constantes
    use Legendre
    use functions

    implicit none 

    contains 

    ! Calcule la quadrature utilisée pour la matrice L décrite dans le rapport: dépend de l'indice i et j des polynomes 
    ! de Legendre utilisés, de lambda, et de p (degré du polynome max: p-1)
    function quad_L(i, j, lambda, p) result(res)
    
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

    ! Calcule la quadrature utilisée pour la matrice M
    function quad_M(i, j, lambda, p) result(res)
    
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

    ! Calcule la quadrature utilisée pour la matrice N
    function quad_N(i, j, lambda, p) result(res)
    
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

    ! Crée la matrice L
    function make_L(lambda, p) result(L)

        real(kind=PR), intent(in) :: lambda
        integer, intent(in) :: p
        real(kind=PR), dimension(p, p) :: L
        integer :: i, j

        do i = 1, p 
            do j = 1, p
                L(i, j) = 1._PR/norme(i)*quad_L(j, i, lambda, p)
            end do 
        end do 


    end function

    ! Crée la matrice M
    function make_M(lambda, p) result(M)

        real(kind=PR), intent(in) :: lambda
        integer, intent(in) :: p
        real(kind=PR), dimension(p, p) :: M
        integer :: i, j
   
        do i = 1, p
            do j = 1, p
                M(i, j) = 1._PR/norme(i)/lambda*quad_M(j, i, lambda, p)
            end do 
        end do 
    end function

    ! Crée la matrice N
    function make_N(lambda, p) result(N)

        real(kind=PR), intent(in) :: lambda
        integer, intent(in) :: p
        real(kind=PR), dimension(p, p) :: N
        integer :: i, j

        do i = 1, p
            do j = 1, p
                N(i, j) = (1._PR-1._PR/lambda)/norme(i)*quad_N(j, i, lambda, p)
            end do 
        end do 

    end function

    ! Calcule la quadrature utilisant la condition initiale u(0, x), sert à projeter la condition initiale 
    ! dans la base de Legendre
    function quad_init(i, j, p, dx, case) result(res)
    
        integer, intent(in) :: i, j, p, case
        real(kind=PR), intent(in) :: dx
        real(kind=PR) :: res 
        integer :: k, q 

        q = (p-1)/2+2

        res = 0._PR

        do k = q*(q-1)/2+1, q*(q-1)/2+q
            res = res + weight(k)*u_init((dx/2._PR*points(k)+(i-0.5_PR)*dx), case)* &
                    Leg(j, points(k))
        end do 

    end function

    ! Calcule la quadrature utilisant la condition de bord, sert à projeter la condition de bord
    ! dans la base de Legendre
    function quad_bound(n, j, p, dt, Lx, a, case) result(res)
    
        integer, intent(in) :: j, p, n, case
        real(kind=PR), intent(in) :: dt, Lx, a
        real(kind=PR) :: res 
        integer :: k, q 

        q = (p-1)/2+2

        res = 0._PR
        do k = q*(q-1)/2+1, q*(q-1)/2+q
            res = res + weight(k)*Leg(j, points(k))*u_bound(dt/2._PR*points(k)+(n+0.5_PR)*dt, Lx, a, case)
        end do 

    end function



end module 
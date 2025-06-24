module matrix

    use constantes
    use Legendre
    use functions

    implicit none 

    contains 

    ! Calcule la quadrature utilisée pour la matrice L décrite dans le rapport: depend de l'indice i et j des polynomes 
    ! de Legendre utilises, et de lambda
    function quad_L(i, j, lambda) result(res)
    
        integer, intent(in) :: i, j
        real(kind=PR), intent(in) :: lambda
        real(kind=PR) :: res 
        integer :: k, q 

        q = (i+j-1)/2 + 1 ! Nombre de points a utiliser pour la quadrature 

        res = 0.0_PR

        do k = q*(q-1)/2+1, q*(q-1)/2+q ! Indices permettant de prendre les bons poids/abcisse dans la liste, bases sur la formule n(n+1)/2
            res = res + weight(k)*Leg(j,1.0_PR-(points(k)+1.0_PR)/lambda)*Leg(i, points(k))
        end do 

    end function

    ! Calcule la quadrature utilisee pour la matrice M
    function quad_M(i, j, lambda) result(res)
    
        integer, intent(in) :: i, j
        real(kind=PR), intent(in) :: lambda 
        real(kind=PR) :: res 
        integer :: k, q

        q = (i+j-1)/2 + 1

        res = 0.0_PR

        do k = q*(q-1)/2+1, q*(q-1)/2+q
            res = res + weight(k)*Leg(j, -points(k))*Leg(i, (points(k)+1.0_PR)/lambda-1.0_PR)
        end do 

    end function

    ! Calcule la quadrature utilisee pour la matrice N
    function quad_N(i, j, lambda) result(res)
    
        integer, intent(in) :: i, j
        real(kind=PR), intent(in) :: lambda
        real(kind=PR) :: res 
        integer :: k, q

        q = (i+j-1)/2 + 1

        res = 0.0_PR

        do k = q*(q-1)/2+1, q*(q-1)/2+q
            res = res + weight(k)*Leg(j, points(k)*(1.0_PR-1.0_PR/lambda)-1.0_PR/lambda)* &
                    Leg(i, points(k)*(1.0_PR-1.0_PR/lambda)+1.0_PR/lambda)
        end do 

    end function

    ! Cree la matrice L
    function make_L(lambda, p, a) result(L)

        real(kind=PR), intent(in) :: lambda, a
        integer, intent(in) :: p
        real(kind=PR), dimension(p, p) :: L
        integer :: i, j

        do i = 1, p 
            do j = 1, p
                L(i, j) = 1.0_PR/norme(i)*quad_L(i, j, lambda)

                if (a<0) then 
                    L(i, j) = (-1.0_PR)**(i+1) * L(i, j)
                end if 
            end do 
        end do 


    end function

    ! Cree la matrice M
    function make_M(lambda, p, a) result(M)

        real(kind=PR), intent(in) :: lambda, a 
        integer, intent(in) :: p
        real(kind=PR), dimension(p, p) :: M
        integer :: i, j
   
        do i = 1, p
            do j = 1, p
                M(i, j) = 1.0_PR/norme(i)/lambda*quad_M(i, j, lambda)

                if (a<0) then 
                    M(i, j) = (-1.0_PR)**(j+1) * M(i, j)
                end if 
            end do 
        end do 
    end function

    ! Cree la matrice N
    function make_N(lambda, p) result(N)

        real(kind=PR), intent(in) :: lambda
        integer, intent(in) :: p
        real(kind=PR), dimension(p, p) :: N
        integer :: i, j

        do i = 1, p
            do j = 1, p
                N(i, j) = (1.0_PR-1.0_PR/lambda)/norme(i)*quad_N(i, j, lambda)
            end do 
        end do 

    end function

    ! Cree la matrice O
    function make_O(lambda, p, a) result(O)

        real(kind=PR), intent(in) :: lambda, a 
        integer, intent(in) :: p
        real(kind=PR), dimension(p, p) :: O
        integer :: i, j
   
        if (a>0) then 

            O = make_N(1.0_PR/lambda, p)

        else if (a<0) then 

            do i = 1, p
                do j = 1, p
                    O(i, j) = (-1.0_PR)**(i+j)*(1.0_PR-lambda)/norme(i)*quad_N(i, j, 1.0_PR/lambda)
                end do 
            end do 

        end if 
    end function

    ! Cree la matrice P
    function make_P(lambda, p, a) result(mat_P)

        real(kind=PR), intent(in) :: lambda, a 
        integer, intent(in) :: p
        real(kind=PR), dimension(p, p) :: mat_P
        integer :: i, j
   
        if (a>0) then 

            mat_P = make_M(1.0_PR/lambda, p, a)

        else if (a<0) then 

            do i = 1, p
                do j = 1, p
                    mat_P(i, j) = (-1.0_PR)**(i+1)*lambda/norme(i)*quad_M(i, j, 1.0_PR/lambda)
                end do 
            end do 

        end if 
    end function

    ! Cree la matrice Q
    function make_Q(lambda, p, a) result(Q)

        real(kind=PR), intent(in) :: lambda, a 
        integer, intent(in) :: p
        real(kind=PR), dimension(p, p) :: Q
        integer :: i, j
   
        if (a>0) then 

            Q = make_L(1.0_PR/lambda, p, a)

            else if (a<0) then 

            do i = 1, p
                do j = 1, p
                    Q(i, j) = (-1.0_PR)**(j+1)*1.0_PR/norme(i)*quad_L(i, j, 1.0_PR/lambda)
                end do 
            end do 

        end if 
    end function


    ! Calcule la quadrature utilisant la condition initiale u(0, x), sert a projeter la condition initiale 
    ! dans la base de Legendre
    function quad_init(i, j, p, dx, Lx, cas) result(res)
    
        integer, intent(in) :: i, j, p, cas
        real(kind=PR), intent(in) :: dx, Lx
        real(kind=PR) :: res 
        integer :: k, q 

        q = (p+j-1)/2 + 1

        do k = q*(q-1)/2+1, q*(q-1)/2+q
            res = res + weight(k)*u_init((dx/2.0_PR*points(k) + Lx + (i-0.5_PR)*dx), cas)* &
                    Leg(j, points(k))
        end do 

    end function

    ! Calcule la quadrature utilisant la condition de bord, sert a projeter la condition de bord
    ! dans la base de Legendre
    function quad_bound(n, j, p, dt, Lx, Rx, a, cas) result(res)
    
    integer, intent(in) :: j, p, n, cas
    real(kind=PR), intent(in) :: dt, Lx, Rx, a
    real(kind=PR) :: res 
    integer :: k, q 

    q = (p+j-1)/2 + 1

    res = 0.0_PR
    do k = q*(q-1)/2+1, q*(q-1)/2+q
        res = res + weight(k) * Leg(j, points(k)) * &
              u_bound(dt/2.0_PR*points(k)+(n+0.5_PR)*dt, Lx, Rx, a, cas)
    end do 

end function




end module 
module quad 

    use constantes
    use functions
    use Legendre

    implicit none 

    contains 

    ! Calcule la quadrature utilisée pour la matrice L décrite dans le rapport: depend de l'indice i et j des polynomes 
    ! de Legendre utilises, et de lambda
    function quad_L(i, j, p, lambda, tau, tau_2) result(res)
    
        integer, intent(in) :: i, j, p
        real(kind=PR), intent(in) :: lambda, tau, tau_2
        real(kind=PR) :: res 
        integer :: k, q 

        q = (i+j+p)/2 

        res = 0.0_PR

        do k = q * (q-1)/2 + 1, q * (q-1)/2 + q ! Indices permettant de prendre les bons poids/abcisse dans la liste, bases sur la formule n(n+1)/2

            res = res + weight(k) * Leg(j, 1.0_PR - (points(k) + 1.0_PR) / lambda) * &
                    Leg(i, points(k)) * &
                    exp( -tau_2 / (2.0_PR * tau) * (1.0_PR + points(k)))

        end do 

    end function    

    ! Calcule la quadrature utilisee pour le vecteur V_L
    function quad_V_L(i, p, tau, tau_2) result(res)
    
        integer, intent(in) :: i, p
        real(kind=PR), intent(in) :: tau, tau_2
        real(kind=PR) :: res 
        integer :: k, q 

        q = (i+p-1)/2 + 1 

        res = 0.0_PR

        do k = q * (q-1)/2 + 1, q * (q-1)/2 + q 

            res = res + weight(k) * (1.0_PR - exp( -tau_2 / (2.0_PR * tau) * (1.0_PR + points(k)))) * &
                    Leg(i, (points(k)))

        end do 

    end function    

    ! Calcule la quadrature utilisee pour la matrice M
    function quad_M(i, j, p, lambda, tau, tau_2) result(res)
    
        integer, intent(in) :: i, j, p
        real(kind=PR), intent(in) :: lambda, tau, tau_2
        real(kind=PR) :: res 
        integer :: k, q

        q = (i+j+p)/2 

        res = 0.0_PR

        do k = q * (q-1)/2 + 1, q * (q-1)/2 + q

            res = res + weight(k) * Leg(j, -points(k)) * &
                    Leg(i, (points(k) + 1.0_PR) / lambda - 1.0_PR) * & 
                    exp( -tau_2 / (2.0_PR * tau) * (1.0_PR + points(k)))      
                       
        end do 

    end function

    ! Calcule la quadrature utilisee pour le vecteur V_M
    function quad_V_M(i, p, lambda, tau, tau_2) result(res)
    
        integer, intent(in) :: i, p
        real(kind=PR), intent(in) :: lambda, tau, tau_2
        real(kind=PR) :: res 
        integer :: k, q 

        q = (i+p-1)/2 + 1 

        res = 0.0_PR

        do k = q * (q-1)/2 + 1, q * (q-1)/2 + q 

            res = res + weight(k) * (1.0_PR - exp(-tau_2 / (2.0_PR * tau) * (1.0_PR + points(k)))) * &
                    Leg(i, (points(k) + 1.0_PR) / lambda - 1.0_PR)

        end do 

    end function   

    ! Calcule la quadrature utilisee pour la matrice N
    function quad_N(i, j, lambda, tau, tau_2) result(res)
    
        integer, intent(in) :: i, j
        real(kind=PR), intent(in) :: lambda, tau, tau_2
        real(kind=PR) :: res 
        integer :: k, q

        q = (i+j-1)/2 + 1

        res = 0.0_PR

        do k = q * (q-1)/2 + 1, q * (q-1)/2 + q

            res = res + weight(k) * Leg(j, points(k) * (1.0_PR - 1.0_PR/lambda) -1.0_PR/lambda)* &
                    Leg(i, points(k) * (1.0_PR - 1.0_PR/lambda) + 1.0_PR/lambda) 
        
        end do 

        res = res * exp(-tau_2 / tau)

    end function

    ! Calcule la quadrature utilisee pour le vecteur V_N
    function quad_V_N(i, lambda, tau, tau_2) result(res)
    
        integer, intent(in) :: i
        real(kind=PR), intent(in) :: lambda, tau, tau_2
        real(kind=PR) :: res 
        integer :: k, q 

        q = i/2 + 1

        res = 0.0_PR

        do k = q * (q-1)/2 + 1, q * (q-1)/2 + q 

            res = res + weight(k) * Leg(i, points(k) * (1.0_PR - 1.0_PR/lambda) + 1.0_PR/lambda)

        end do 

        res = res * (1.0_PR - exp(- tau_2 / tau))

    end function  

    function quad_a_0(i, j, tau, tau_2) result(res)

        integer, intent(in) :: i, j
        real(kind=PR), intent(in) :: tau, tau_2
        real(kind=PR) :: res 
        integer :: k, q

        q = (i+j-1)/2 + 1

        res = 0.0_PR

        do k = q * (q-1)/2 + 1, q * (q-1)/2 + q 

            res = res + weight(k) * Leg(i, points(k)) * Leg(j, points(k))

        end do 

        res = res * exp(-tau_2 / tau)
    
        end function

    ! Calcule la quadrature utilisant la condition initiale u(0, x), sert a projeter la condition initiale 
    ! dans la base de Legendre
    function quad_init(i, j, p, dx, Lx, cas) result(res)
    
        integer, intent(in) :: i, j, p, cas
        real(kind=PR), intent(in) :: dx, Lx
        real(kind=PR) :: res 
        integer :: k, q 

        q = (p+j-1)/2 + 1 

        res = 0.0_PR

        do k = q * (q-1)/2 + 1, q * (q-1)/2 + q

            res = res + weight(k) * u_init((dx/2.0_PR * points(k) + Lx + (i-0.5_PR) * dx), cas)* &
                    Leg(j, points(k))

        end do 

    end function

    ! Calcule la quadrature utilisant la condition de bord, sert a projeter la condition de bord
    ! dans la base de Legendre
    function quad_bound(n, j, p, dt, Lx, Rx, a, cas, C, tau) result(res)
    
    integer, intent(in) :: j, p, n, cas
    real(kind=PR), intent(in) :: dt, Lx, Rx, a, C, tau
    real(kind=PR) :: res 
    integer :: k, q 

    q = (j + 2*p)/2 + 2

    res = 0.0_PR

    do k = q * (q-1)/2 + 1, q * (q-1)/2 + q

        res = res + weight(k) * Leg(j, points(k)) * &
              u_bound(dt/2.0_PR * points(k) + (n + 0.5_PR) * dt, Lx, Rx, a, cas, C, tau)
    end do 

    end function

end module
module quad 

    use constantes
    use functions
    use Legendre

    implicit none 

    contains 

    ! Calcule la quadrature utilisée pour la matrice L décrite dans le rapport: depend de l'indice i et j des polynomes 
    ! de Legendre utilises, et de lambda
    function quad_L(i, j, lambda) result(res)
    
        integer, intent(in) :: i, j
        real(kind=PR), intent(in) :: lambda
        real(kind=PR) :: res 
        integer :: k, q 

        q = (i+j-1)/2 + 2 ! Indice que l'on utilisera pour chaque quadrature, permet de calculer l'integrale de polynome de maniere exacte

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

        q = (i+j-1)/2 + 2

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

        q = (i+j-1)/2 + 2

        res = 0.0_PR

        do k = q*(q-1)/2+1, q*(q-1)/2+q
            res = res + weight(k)*Leg(j, points(k)*(1.0_PR-1.0_PR/lambda)-1.0_PR/lambda)* &
                    Leg(i, points(k)*(1.0_PR-1.0_PR/lambda)+1.0_PR/lambda)
        end do 

    end function

end module
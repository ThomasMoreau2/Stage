module Legendre

    use constantes

    implicit none 

    ! Liste des poids pour les quadratures
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

    ! Liste des abcisses pour les quadratures
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

    ! Fonction renvoyant la valeur du polynome de Legendre de degre (i-1) a l'abcisse x 
    function Leg(i, x) result(res)

        integer, intent(in) :: i
        real(kind=PR), intent(in) :: x
        real(kind=PR) :: res

        select case(i)
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
        case(7)
            res = (1.0_PR / 16.0_PR) * (231.0_PR * x**6 - 315.0_PR * x**4 + 105.0_PR * x**2 - 5.0_PR)
        case(8)
            res = (1.0_PR / 16.0_PR) * (429.0_PR * x**7 - 693.0_PR * x**5 + 315.0_PR * x**3 - 35.0_PR * x)
        case(9)
            res = (1.0_PR / 128.0_PR) * (6435.0_PR * x**8 - 12012.0_PR * x**6 + 6930.0_PR * x**4 - 1260.0_PR * x**2 + 35.0_PR)
        case(10)
            res = (1.0_PR / 128.0_PR) * (12155.0_PR * x**9 - 25740.0_PR * x**7 + 18018.0_PR * x**5 - &   
                                            4620.0_PR * x**3 + 315.0_PR * x)
        case default
            res = 1.0_PR
        end select
    end function

    ! Norme d'un polynome de Legendre de degre i-1
    function norme(i) result(res)

        integer, intent(in) :: i 
        real(kind=PR) :: res 
        
        res = 2.0_PR/(2.0_PR*(i-1.0_PR)+1.0_PR)

    end function

    ! Fonction prennant une liste de coeffients alpha (coordonnees dans la base de Legendre locale) et renvoie 
    ! la valeur de u aux points (x_1/2, x_3/2, ..., x_imax+1/2)
    function calculate_u(alpha, i_max, p, a) result(u)
         
    real(kind=PR), dimension(:), intent(in) :: alpha
        integer, intent(in) :: i_max, p
        real(kind=PR), intent(in) :: a
        real(kind=PR), dimension(i_max+1) :: u 
        integer :: i, l

        if (a>0) then ! On distingue a>0 ou a<0

            do i = 1, i_max
                u(i) = 0.0_PR
                do l = 1, p 
                    u(i) = u(i) + alpha((i-1)*p + l) * Leg(l, -1.0_PR) ! L'abcisse (i-1)*dx correpsond à evaluer le polynome en -1
                end do 
            end do 

            u(i_max+1) = 0._PR
            do l = 1, p 
                u(i_max+1) = u(i_max+1) + alpha((i_max-1)*p + l) * Leg(l, 1.0_PR) ! Pour le dernier, il s'agit du polynome en 1
            end do 
        
        else if (a<0) then 

            do i = 1, i_max
                u(i) = 0.0_PR
                do l = 1, p 
                    u(i) = u(i) + alpha((i-1)*p + l) * Leg(l, 1.0_PR) ! Ici c'est l'inverse, car on parcourt de droite a gauche, c'est donc en 1
                end do 
            end do 

            u(i_max+1) = 0._PR
            do l = 1, p 
                u(i_max+1) = u(i_max+1) + alpha((i_max-1)*p + l) * Leg(l, -1.0_PR) ! Puis en -1 pour le dernier
            end do 
        
            
        end if 
    end function



end module 
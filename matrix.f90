module matrix

    use constantes
    use functions
    use Legendre
    use quad

    implicit none 

    contains 

    ! Cree la matrice L
    function make_L(p, lambda, a) result(L)

        real(kind=PR), intent(in) :: lambda, a
        integer, intent(in) :: p
        real(kind=PR), dimension(0:p, p, p) :: L
        integer :: i, j, k

        do k = 0, p
            do i = 1, p 
                do j = 1, p

                    L(k, i, j) = 1.0_PR / norme(i) / fact(k) * quad_L(k, i, j, lambda)

                    if (a < 0.0_PR) then 

                        L(k, i, j) = (-1.0_PR)**(i+1) * L(k, i, j)

                    end if 

                end do 
            end do 
        end do 

    end function

    ! Cree le vecteur V_L
    function make_V_L(p, a) result(V_L)

        real(kind=PR), intent(in) :: a
        integer, intent(in) :: p
        real(kind=PR), dimension(0:p, p) :: V_L
        integer :: i, k

        V_L(0, 1:p) = 0.0_PR

        do k = 1, p
            do i = 1, p

                V_L(k, i) = - 1.0_PR / norme(i) / fact(k) * quad_V_L(k, i)

                if (a < 0.0_PR) then

                    V_L(k, i) = (-1.0_PR)**(i+1) * V_L(k, i)

                end if 

            end do 
        end do 

    end function

    ! Cree la matrice M
    function make_M(p, lambda, a) result(M)

        real(kind=PR), intent(in) :: lambda, a
        integer, intent(in) :: p
        real(kind=PR), dimension(0:p, p, p) :: M
        integer :: i, j, k
        
        do k = 0, p
            do i = 1, p
                do j = 1, p 

                    M(k, i, j) = 1.0_PR / norme(i) / lambda / fact(k) * quad_M(k, i, j, lambda)

                    if (a < 0.0_PR) then 

                        M(k, i, j) = (-1.0_PR)**(j+1) * M(k, i, j)

                    end if 

                end do 
            end do 
        end do 

    end function

    ! Cree le vecteur V_M
    function make_V_M(p, lambda) result(V_M)

        real(kind=PR), intent(in) ::  lambda
        integer, intent(in) :: p
        real(kind=PR), dimension(0:p, p) :: V_M
        integer :: i, k 

        V_M(0, 1:p) = 0.0_PR

        do k = 1, p
            do i = 1, p

                V_M(k, i) = -1.0_PR / norme(i) / lambda / fact(k) * quad_V_M(k, i, lambda)

            end do 
        end do 

    end function

    ! Cree la matrice N
    function make_N(p, lambda) result(N)

        real(kind=PR), intent(in) :: lambda
        integer, intent(in) :: p
        real(kind=PR), dimension(p, p) :: N
        integer :: i, j

        do i = 1, p
            do j = 1, p 

                N(i, j) = (1.0_PR - 1.0_PR/lambda) / norme(i) * quad_N(i, j, lambda) 

            end do 
        end do 

    end function

    ! Cree le vecteur V_N
    function make_V_N(p, lambda) result(V_N)

        real(kind=PR), intent(in) :: lambda
        integer, intent(in) :: p
        real(kind=PR), dimension(p) :: V_N
        integer :: i

        do i = 1, p

            V_N(i) = (1.0_PR - 1.0_PR/lambda) / norme(i) * quad_V_N(i, lambda)

        end do 

    end function

    ! Cree la matrice O
    function make_O(p, lambda, a) result(O)

        real(kind=PR), intent(in) :: lambda, a 
        integer, intent(in) :: p
        real(kind=PR), dimension(p, p) :: O
        integer :: i, j

        O = make_N(p, 1.0_PR/lambda)

        if (a < 0.0_PR) then 

            do i = 1, p
                do j = 1, p

                    O(i, j) = (-1.0_PR)**(i+j) * O(i, j)

                end do 
            end do 

        end if 

    end function

    ! Cree le vecteur V_O
    function make_V_O(p, lambda, a) result(V_O)

        real(kind=PR), intent(in) :: lambda, a 
        integer, intent(in) :: p
        real(kind=PR), dimension(p) :: V_O
        integer :: i

        V_O = make_V_N(p, 1/lambda)

        if (a < 0.0_PR) then 

            do i=1, p 

                V_O(i) = (-1.0_PR)**(i+1) * V_O(i) 

            end do 

        end if  

    end function

    ! Cree la matrice P
    function make_P(p, lambda, a) result(mat_P)

        real(kind=PR), intent(in) :: lambda, a
        integer, intent(in) :: p
        real(kind=PR), dimension(0:p, p, p) :: mat_P
        integer :: i, j, k

        mat_P = make_M(p, 1.0_PR/lambda, abs(a))

        if (a < 0.0_PR) then 

            do k = 0, p
                do i = 1, p
                    do j = 1, p

                        mat_P(k, i, j) = (-1.0_PR)**(i+1) * mat_P(k, i, j)
                        
                    end do 
                end do 
            end do 

        end if 
    end function

    ! Cree le vecteur V_P
    function make_V_P(p, lambda, a) result(V_P)

        real(kind=PR), intent(in) :: lambda, a 
        integer, intent(in) :: p
        real(kind=PR), dimension(0:p, p) :: V_P
        integer :: i, k

        V_P = make_V_M(p, 1/lambda)

        if (a < 0.0_PR) then 

            do k = 0, p
                do i=1, p 

                    V_P(k, i) = (-1.0_PR)**(i+1) * V_P(k, i) 

                end do 
            end do 

        end if  

    end function

    ! Cree la matrice Q
    function make_Q(p, lambda, a) result(Q)

        real(kind=PR), intent(in) :: lambda, a
        integer, intent(in) :: p
        real(kind=PR), dimension(0:p, p, p) :: Q
        integer :: i, j, k

        Q = make_L(p, 1.0_PR/lambda, abs(a))

        if (a < 0.0_PR) then 

            do k = 0, p
                do i = 1, p
                    do j = 1, p

                        Q(k, i, j) = (-1.0_PR)**(j+1) * Q(k, i, j)

                    end do 
                end do 
            end do

        end if 
    end function

     ! Cree le vecteur V_Q
    function make_V_Q(p, a) result(V_Q)

        real(kind=PR), intent(in) :: a
        integer, intent(in) :: p
        real(kind=PR), dimension(0:p, p) :: V_Q

        V_Q = make_V_L(p, abs(a))

    end function

    ! Cree la matrice lorsque a=0
    function make_a_0(p) result(mat_a_0)

        integer, intent(in) :: p
        real(kind=PR), dimension(p, p) :: mat_a_0
        integer :: i, j

        do i = 1, p
            do j = 1, p

                mat_a_0(i, j) = 1.0_PR / norme(i) * quad_a_0(i, j)

            end do 
        end do 

    end function

    ! Cree le vecteur lorsque a=0
    function make_V_a_0(p) result(V_a_0)

        integer, intent(in) :: p
        real(kind=PR), dimension(p) :: V_a_0
        integer :: i

        do i = 1, p

            V_a_0(i) = 1.0_PR / norme(i) * quad_V_a_0(i)

        end do 

    end function



end module 
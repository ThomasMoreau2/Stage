module matrix

    use constantes
    use functions
    use Legendre
    use quad

    implicit none 

    contains 

    ! Cree la matrice L
    function make_L(p, lambda, tau, tau_2, a) result(L)

        real(kind=PR), intent(in) :: lambda, a, tau, tau_2
        integer, intent(in) :: p
        real(kind=PR), dimension(p, p) :: L
        integer :: i, j

        do i = 1, p 
            do j = 1, p

                L(i, j) = 1.0_PR / norme(i) * quad_L(i, j, p, lambda, tau, tau_2)

                if (a < 0.0_PR) then 
                    L(i, j) = (-1.0_PR)**(i+1) * L(i, j)
                end if 

            end do 
        end do 

    end function

    ! Cree le vecteur V_L
    function make_V_L(p, tau, tau_2, a) result(V_L)

        real(kind=PR), intent(in) :: tau, tau_2, a
        integer, intent(in) :: p
        real(kind=PR), dimension(p) :: V_L
        integer :: i

        do i = 1, p

            V_L(i) = 1.0_PR / norme(i) * quad_V_L(i, p, tau, tau_2)

            if (a < 0.0_PR) then 
                V_L(i) = (-1.0_PR)**(i+1) * V_L(i)
            end if 

        end do 

    end function

    ! Cree la matrice M
    function make_M(p, lambda, tau, tau_2, a) result(M)

        real(kind=PR), intent(in) :: lambda, a, tau, tau_2
        integer, intent(in) :: p
        real(kind=PR), dimension(p, p) :: M
        integer :: i, j
   
        do i = 1, p
            do j = 1, p 

                M(i, j) = 1.0_PR / norme(i) / lambda * quad_M(i, j, p, lambda, tau, tau_2)

                if (a < 0.0_PR) then 
                    M(i, j) = (-1.0_PR)**(j+1) * M(i, j)
                end if 

            end do 
        end do 

    end function

    ! Cree le vecteur V_M
    function make_V_M(p, lambda, tau, tau_2) result(V_M)

        real(kind=PR), intent(in) :: tau, tau_2, lambda
        integer, intent(in) :: p
        real(kind=PR), dimension(p) :: V_M
        integer :: i

        do i = 1, p

            V_M(i) = 1.0_PR / norme(i) / lambda * quad_V_M(i, p, lambda, tau, tau_2)

        end do 

    end function

    ! Cree la matrice N
    function make_N(p, lambda, tau, tau_2) result(N)

        real(kind=PR), intent(in) :: lambda, tau, tau_2
        integer, intent(in) :: p
        real(kind=PR), dimension(p, p) :: N
        integer :: i, j

        do i = 1, p
            do j = 1, p 

                N(i, j) = (1.0_PR - 1.0_PR/lambda) / norme(i) * quad_N(i, j, lambda, tau, tau_2) 

            end do 
        end do 

    end function

    ! Cree le vecteur V_N
    function make_V_N(p, lambda, tau, tau_2) result(V_N)

        real(kind=PR), intent(in) :: tau, tau_2, lambda
        integer, intent(in) :: p
        real(kind=PR), dimension(p) :: V_N
        integer :: i

        do i = 1, p

            V_N(i) = (1.0_PR - 1.0_PR/lambda) / norme(i) * quad_V_N(i, lambda, tau, tau_2)

        end do 

    end function

    ! Cree la matrice O
    function make_O(p, lambda, tau, tau_2, a) result(O)

        real(kind=PR), intent(in) :: lambda, a , tau, tau_2
        integer, intent(in) :: p
        real(kind=PR), dimension(p, p) :: O
        integer :: i, j

        O = make_N(p, 1.0_PR/lambda, tau, tau_2)

        if (a < 0.0_PR) then 

            do i = 1, p
                do j = 1, p
                    O(i, j) = (-1.0_PR)**(i+j) * O(i, j)
                end do 
            end do 

        end if 

    end function

    ! Cree le vecteur V_O
    function make_V_O(p, lambda, tau, tau_2, a) result(V_O)

        real(kind=PR), intent(in) :: tau, tau_2, lambda, a 
        integer, intent(in) :: p
        real(kind=PR), dimension(p) :: V_O
        integer :: i

        V_O = make_V_N(p, 1/lambda, tau, tau_2)

        if (a < 0.0_PR) then 

            do i=1, p 
                V_O(i) = (-1.0_PR)**(i+1) * V_O(i) 
            end do 

        end if  

    end function

    ! Cree la matrice P
    function make_P(p, lambda, tau, tau_2, a) result(mat_P)

        real(kind=PR), intent(in) :: lambda, a, tau, tau_2
        integer, intent(in) :: p
        real(kind=PR), dimension(p, p) :: mat_P
        integer :: i, j

        mat_P = make_M(p, 1.0_PR/lambda, tau, tau_2, abs(a))

        if (a < 0.0_PR) then 

            do i = 1, p
                do j = 1, p
                    mat_P(i, j) = (-1.0_PR)**(i+1) * mat_P(i, j)
                end do 
            end do 

        end if 
    end function

    ! Cree le vecteur V_P
    function make_V_P(p, lambda, tau, tau_2, a) result(V_P)

        real(kind=PR), intent(in) :: tau, tau_2, lambda, a 
        integer, intent(in) :: p
        real(kind=PR), dimension(p) :: V_P
        integer :: i

        V_P = make_V_M(p, 1/lambda, tau, tau_2)

        if (a < 0.0_PR) then 

            do i=1, p 
                V_P(i) = (-1.0_PR)**(i+1) * V_P(i) 
            end do 

        end if  

    end function

    ! Cree la matrice Q
    function make_Q(p, lambda, tau, tau_2, a) result(Q)

        real(kind=PR), intent(in) :: lambda, a, tau, tau_2
        integer, intent(in) :: p
        real(kind=PR), dimension(p, p) :: Q
        integer :: i, j

        Q = make_L(p, 1.0_PR/lambda, tau, tau_2, abs(a))

        if (a < 0.0_PR) then 

            do i = 1, p
                do j = 1, p
                    Q(i, j) = (-1.0_PR)**(j+1) * Q(i, j)
                end do 
            end do 

        end if 
    end function

     ! Cree le vecteur V_P
    function make_V_Q(p, tau, tau_2, a) result(V_Q)

        real(kind=PR), intent(in) :: tau, tau_2, a
        integer, intent(in) :: p
        real(kind=PR), dimension(p) :: V_Q

        V_Q = make_V_L(p, tau, tau_2, abs(a))

    end function

    ! Cree la matrice lorsque a=0
    function make_a_0(p, tau, tau_2) result(mat_a_0)

        real(kind=PR), intent(in) :: tau, tau_2
        integer, intent(in) :: p
        real(kind=PR), dimension(p, p) :: mat_a_0
        integer :: i, j

        do i = 1, p
            do j = 1, p
                mat_a_0(i, j) = 1.0_PR / norme(i) * quad_a_0(i, j, tau, tau_2)
            end do 
        end do 

    end function


end module 
program main 

    use constantes
    use functions
    use Legendre
    use quad
    use matrix

    implicit none 

    integer :: i_max, n_max, p, i, j, r, k, cas
    real(kind=PR) :: dx, dt, lambda, t_final, Lx, Rx, a, err_N2, Mn_i, tau, t_exacte
    real(kind=PR), dimension(:), allocatable :: alpha, alpha_np1, beta, u
    real(kind=PR), dimension(:), allocatable :: V_L, V_M, V_N, V_O, V_P, V_Q
    real(kind=PR), dimension(:, :), allocatable :: mat_L, mat_M, mat_N, mat_O, mat_P, mat_Q, mat_a_0


    ! Definition des parametres
    p = 3
    cas = 1
    Lx = -1.0_PR
    Rx = 1.0_PR
    t_final = 1.0_PR
    a = 0.0_PR 
    lambda = 2.0
    tau = 1.0_PR
    Mn_i = 0.0_PR

    ! Allocation des matrices L, M, N et vecteurs V_L, V_M, V_N 
    if (a /= 0.0_PR) then 

        if (abs(lambda) >= 1.0_PR) then 

            allocate(mat_L(p, p), mat_M(p, p), mat_N(p, p))
            allocate(V_L(p), V_M(p), V_N(p))

        else if ((abs(lambda) < 1.0_PR)) then 

            allocate(mat_O(p, p), mat_P(p, p), mat_Q(p, p))
            allocate(V_O(p), V_P(p), V_Q(p))

        end if 

    else 

        allocate(mat_a_0(p, p))

    end if 

    ! Ouverture du fichier pour l'ordre 
    open(unit=30, file="erreur.dat", action='write')
  
    
    ! Boucle pour le calcul d'ordre
    do k = 0, 4

        ! Allocation des parametres pour le calcul d'ordre
        dx = 0.5_PR / (2.0_PR**k)
        i_max = int((Rx-Lx)/dx)

        if (a /= 0.0_PR) then 
            dt = lambda * dx / abs(a)
        else 
            dt = lambda * dx 
        end if 

        n_max = int(t_final/dt)
        t_exacte = n_max * dt

        ! Creation des matrices/vecteurs
        if (a /= 0.0_PR) then 

            if (abs(lambda) >= 1.0_PR) then 

                mat_L = make_L(p, lambda, tau, abs(dx / a), a)
                mat_M = make_M(p, lambda, tau, abs(dx / a), a)
                mat_N = make_N(p, lambda, tau, abs(dx / a))

                V_L = make_V_L(p, tau, abs(dx / a), a)
                V_M = make_V_M(p, lambda, tau, abs(dx / a))
                V_N = make_V_N(p, lambda, tau, abs(dx / a))

            else if ((abs(lambda) < 1.0_PR)) then 

                mat_O = make_O(p, lambda, tau, dt, a)
                mat_P = make_P(p, lambda, tau, dt, a)
                mat_Q = make_Q(p, lambda, tau, dt, a)

                V_O = make_V_O(p, lambda, tau, dt, a)
                V_P = make_V_P(p, lambda, tau, dt, a)
                V_Q = make_V_Q(p, tau, dt, a)
                
            end if 
        
        else

            mat_a_0 = make_a_0(p, tau, dt)

        end if 
    
        ! Allocation des tableaux
        allocate(alpha(p*i_max), alpha_np1(p*i_max), beta(p*(i_max+1)))
        allocate(u(i_max+1))

        ! Projection de la condition initiale
        do i = 1, i_max
            do j = 1, p 
                alpha((i-1)*p + j) = quad_init(i, j, p, dx, Lx, cas)/norme(j)
            end do 
        end do 

        ! Projection de la condition de bord, depend du signe de a
        if (a > 0.0_PR) then 
            do j = 1, p 
                beta(j) = quad_bound(0, j, p, dt, Lx, Rx, a, cas, Mn_i, tau)/norme(j)
            end do 
        else if (a < 0.0_PR) then 
            do j = 1, p 
                beta(p*i_max+j) = quad_bound(0, j, p, dt, Lx, Rx, a, cas, Mn_i, tau)/norme(j)
            end do 
        end if 


        ! Implementation du schema 
        do r=1, n_max ! Boucle en temps 

            do i = 1, i_max  ! Boucle en espace 

                if (a > 0.0_PR) then ! Cas a positif: on parcourt le vecteur beta de gauche a droite

                    if (abs(lambda) >= 1.0_PR) then 

                        alpha_np1((i-1)*p + 1 : i*p) = matmul(mat_L, beta((i-1)*p + 1 : i*p)) + & 
                                                        Mn_i * V_L

                        beta(i*p + 1 : (i+1)*p) = matmul(mat_M, alpha((i-1)*p + 1 : i*p)) + & 
                                                    matmul(mat_N, beta((i-1)*p + 1 : i*p)) + &
                                                    Mn_i * (V_M + V_N)
                    else if (abs(lambda) < 1.0_PR ) then  

                        alpha_np1((i-1)*p + 1 : i*p) = matmul(mat_O, alpha((i-1)*p + 1 : i*p)) + &
                                                        matmul(mat_P, beta((i-1)*p + 1 : i*p)) + & 
                                                        Mn_i * (V_O + V_P)

                        beta(i*p + 1 : (i+1)*p) = matmul(mat_Q, alpha((i-1)*p + 1 : i*p)) + & 
                                                    Mn_i * V_Q

                    end if 

                else if (a < 0.0_PR) then ! Cas a negatif: on parcourt le vecteur beta de droite a gauche 

                    j = i_max - i + 1 ! Changement d'indice pour parcourir de droite à gauche

                    if (abs(lambda) >= 1.0_PR) then 

                        alpha_np1((j-1)*p + 1 : j*p) = matmul(mat_L, beta(j*p + 1 : (j+1)*p)) + &  
                                                        Mn_i * V_L

                        beta((j-1)*p + 1 : j*p) = matmul(mat_M, alpha((j-1)*p + 1 : j*p)) + &
                                                    matmul(mat_N, beta(j*p + 1 : (j+1)*p)) + & 
                                                    Mn_i * (V_M + V_N)
                    else if (abs(lambda) < 1.0_PR) then  

                        alpha_np1((j-1)*p + 1 : j*p) = matmul(mat_O, alpha((j-1)*p + 1 : j*p)) + &
                                                        matmul(mat_P, beta(j*p + 1 : (j+1)*p)) + & 
                                                        Mn_i * (V_O + V_P)

                        beta((j-1)*p + 1 : j*p) = matmul(mat_Q, alpha((j-1)*p + 1 : j*p)) + &
                                                    Mn_i * V_Q

                    end if 

                else if (a == 0.0_PR) then

                    alpha_np1((i-1)*p + 1 : i*p) = matmul(mat_a_0, alpha((i-1)*p + 1 : i*p)) + &
                                                        Mn_i * 2.0_PR * (1 - exp(-dt / tau))

                end if 
            end do 

            ! Re-projection de la condition de bord à chaque pas de temps 
            if (a > 0.0_PR) then 
                do j = 1, p 
                    beta(j) = quad_bound(r, j, p, dt, Lx, Rx, a, cas, Mn_i, tau)/norme(j)
                end do 
            else if (a < 0.0_PR) then 
                do j = 1, p 
                    beta(p*i_max+j) = quad_bound(r, j, p, dt, Lx, Rx, a, cas, Mn_i, tau)/norme(j)
                end do 
            end if 

            alpha = alpha_np1 ! Mise a jour du vecteur alpha
        end do 


        ! Calcul d'erreur en norme 2 
        u = calculate_u(alpha, i_max, p)
        err_N2 = 0._PR
        do i = 1, i_max+1
            err_N2 = err_N2 + (u(i) - sol_exacte(Lx+dx*(i-1.0_PR), t_exacte, Lx, Rx, a, cas, Mn_i, tau))**2
        end do 
        err_N2 = sqrt((err_N2*dx))

        write(30, *) dx, err_N2


        ! Calcul de la solution et ecriture dans un fichier
        open(unit=20, file="resultats.dat", status='replace', action='write')
        ! u = calculate_u(alpha, i_max, p, dx, Lx)
        do i = 1, i_max+1
            write (20, *) Lx+dx*(i-1.0_PR), u(i), sol_exacte(Lx+dx*(i-1.0_PR), t_exacte, Lx, Rx, a, cas, Mn_i, tau)
        end do 
        close(20)

        ! Deallocation des tableaux
        deallocate(u)
        deallocate(alpha, alpha_np1, beta)

        ! Deallocation des matrices
        if (a /= 0.0_PR) then 
            if (abs(lambda) > 1.0_PR) then
                deallocate(mat_L, mat_M, mat_N) 
            else if ((abs(lambda) < 1.0_PR)) then
                deallocate(mat_O, mat_P, mat_Q)
            end if 
        else 
            deallocate(mat_a_0)
        end if 


    end do 


    ! Fermeture du fichier pour l'ordre 
    close(30)

end program
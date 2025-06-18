program main 

    use constantes
    use functions
    use Legendre
    use quad
    use matrix

    implicit none 

    integer :: i_max, n_max, p, i, j, r, k, cas
    real(kind=PR) :: dx, dt, lambda, t_final, Lx, Rx, a, err_N2
    real(kind=PR), dimension(:), allocatable :: alpha, alpha_np1, beta, u
    real(kind=PR), dimension(:, :), allocatable :: mat_L, mat_M, mat_N, mat_O, mat_P, mat_Q

    ! Lecture et definition des parametres 
    ! open(unit=10, file="parametres.dat", action='read')
    ! read (10, *) p, i_max, n_max, t_final, Lx, Rx, a, cas
    ! close(10) 

    ! dx = Lx/i_max
    ! dt = t_final/n_max

    ! lambda = a*dt/dx
    ! print *, lambda


    ! Parametres pour le calcul d'ordre
    p = 3
    cas = 3
    Lx = 0.0_PR
    Rx = 1.0_PR
    t_final = 1.0_PR
    a = -1.0_PR 
    lambda = 0.25_PR

    ! Creation des matrices L, M, N
    if (abs(lambda)>1) then 

        allocate(mat_L(p, p), mat_M(p, p), mat_N(p, p))

        mat_L = make_L(abs(lambda), p, a)
        mat_M = make_M(abs(lambda), p, a)
        mat_N = make_N(abs(lambda), p)

    else if (abs(lambda)<1) then 

        allocate(mat_O(p, p), mat_P(p, p), mat_Q(p, p))

        mat_O = make_O(abs(lambda), p, a)
        mat_P = make_P(abs(lambda), p, a)
        mat_Q = make_Q(abs(lambda), p, a)

    end if 

    ! Ouverture du fichier pour l'ordre 
    open(unit=30, file="erreur.dat", action='write')
  
    
    ! Boucle pour le calcul d'ordre
    do k = 0, 0

        ! Allocation des parametres pour le calcul d'ordre
        dx = 0.0625_PR / (2**k)
        i_max = int((Rx-Lx)/dx)
        dt = lambda * dx / abs(a)
        n_max = int(t_final/dt)

        print *, dx, i_max, dt, n_max
    
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
        if (a>0) then 
            do j = 1, p 
                beta(j) = quad_bound(0, j, p, dt, Lx, Rx, a, cas)/norme(j)
            end do 
        else if (a<0) then 
            do j = 1, p 
                beta(p*i_max+j) = quad_bound(0, j, p, dt, Lx, Rx, a, cas)/norme(j)
            end do 
        end if 


        ! Implementation du schema 
        do r=1, n_max ! Boucle en temps 

            do i = 1, i_max  ! Boucle en espace 

                if (a>0) then ! Cas a positif: on parcourt le vecteur beta de gauche a droite
                    if (abs(lambda)>1) then 
                        alpha_np1((i-1)*p + 1 : i*p) = matmul(mat_L, beta((i-1)*p + 1 : i*p))
                        beta(i*p + 1 : (i+1)*p) = matmul(mat_M, alpha((i-1)*p + 1 : i*p)) + & 
                                                    matmul(mat_N, beta((i-1)*p + 1 : i*p))
                    else if (abs(lambda)<1) then     
                        alpha_np1((i-1)*p + 1 : i*p) = matmul(mat_O, alpha((i-1)*p + 1 : i*p)) + &
                                                        matmul(mat_P, beta((i-1)*p + 1 : i*p))
                        beta(i*p + 1 : (i+1)*p) = matmul(mat_Q, alpha((i-1)*p + 1 : i*p))
                    end if 

                else if (a<0) then ! Cas a negatif: on parcourt le vecteur beta de droite a gauche 

                    j = i_max - i + 1 ! Changement d'indice pour parcourir de droite à gauche

                    if (abs(lambda)>1) then 
                        alpha_np1((j-1)*p + 1 : j*p) = matmul(mat_L, beta(j*p + 1 : (j+1)*p))
                        beta((j-1)*p + 1 : j*p) = matmul(mat_M, alpha((j-1)*p + 1 : j*p)) + &
                                                    matmul(mat_N, beta(j*p + 1 : (j+1)*p))
                    else if (abs(lambda)<1) then     
                        alpha_np1((j-1)*p + 1 : j*p) = matmul(mat_O, alpha((j-1)*p + 1 : j*p)) + &
                                                        matmul(mat_P, beta(j*p + 1 : (j+1)*p))
                        beta((j-1)*p + 1 : j*p) = matmul(mat_Q, alpha((j-1)*p + 1 : j*p))
                    end if 
                end if 
            end do 

            ! Re-projection de la condition de bord à chaque pas de temps 
            if (a>0) then 
                do j = 1, p 
                    beta(j) = quad_bound(r, j, p, dt, Lx, Rx, a, cas)/norme(j)
                end do 
            else if (a<0) then 
                do j = 1, p 
                    beta(p*i_max+j) = quad_bound(r, j, p, dt, Lx, Rx, a, cas)/norme(j)
                end do 
            end if 

            alpha = alpha_np1 ! Mise a jour du vecteur alpha
        end do 


        ! Calcul d'erreur en norme 2 
        u = calculate_u(alpha, i_max, p)
        err_N2 = 0._PR
        do i = 1, i_max+1
            err_N2 = err_N2 + (u(i) - sol_exacte(Lx+dx*(i-1.0_PR), t_final, Lx, Rx, a, cas))**2
        end do 
        err_N2 = sqrt((err_N2*dx))

        write(30, *) dx, err_N2


        ! Calcul de la solution et ecriture dans un fichier
        open(unit=20, file="resultats.dat", status='replace', action='write')
        ! u = calculate_u(alpha, i_max, p, dx, Lx)
        do i = 1, i_max+1
            write (20, *) Lx+dx*(i-1.0_PR), u(i), sol_exacte(Lx+dx*(i-1.0_PR), t_final, Lx, Rx, a, cas)
        end do 
        close(20)

        ! Deallocation des tableaux
        deallocate(u)
        deallocate(alpha, alpha_np1, beta)

    end do 


    ! Fermeture du fichier pour l'ordre 
    close(30)


    ! Deallocation des matrices
    if (abs(lambda)>1) then
        deallocate(mat_L, mat_M, mat_N) 
    else if (abs(lambda)<1) then
        deallocate(mat_O, mat_P, mat_Q)
    end if 

end program
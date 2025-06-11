program main 

    use constantes
    use Legendre
    use functions
    use matrix

    implicit none 

    integer :: i_max, n_max, p, i, j, r, k, case
    real(kind=PR) :: dx, dt, lambda, t_final, Lx, a, err_N2
    real(kind=PR), dimension(:), allocatable :: alpha, alpha_np1, beta, u
    real(kind=PR), dimension(:, :), allocatable :: mat_L, mat_M, mat_N, mat_O, mat_P, mat_Q

    ! Lecture et definition des parametres 
    ! open(unit=10, file="parametres.dat", action='read')
    ! read (10, *) p, i_max, n_max, t_final, Lx, a, case
    ! close(10)

    ! dx = Lx/i_max
    ! dt = t_final/n_max

    ! lambda = a*dt/dx
    ! print *, lambda


    ! Paramètres pour le calcul d'ordre
    p = 2
    case = 3
    Lx = 1._PR
    t_final = 1._PR
    a = -1._PR 
    lambda = 4._PR

    ! Creation des matrices L, M, N
    if (abs(lambda)>1) then 
        allocate(mat_L(p, p), mat_M(p, p), mat_N(p, p))
        mat_L = make_L(abs(lambda), p)
        mat_M = make_M(abs(lambda), p)
        mat_N = make_N(abs(lambda), p)
    else if (abs(lambda)<1) then 
        allocate(mat_O(p, p), mat_P(p, p), mat_Q(p, p))
        mat_O = make_N(1._PR/abs(lambda), p)
        mat_P = make_M(1._PR/abs(lambda), p)
        mat_Q = make_L(1._PR/abs(lambda), p)
    end if 

    ! Ouverture du fichier pour l'ordre 
    open(unit=30, file="erreur.dat", action='write')
  
    
    ! Boucle pour le calcul d'erreur
    do k = 0, 4

        ! Ordre en espace
        dx = 0.0625_PR / (2**k)
        i_max = int(Lx/dx)
        dt = lambda * dx / abs(a)
        n_max = int(t_final/dt)

        ! Ordre en temps 
        ! dt = 0.1_PR / (2**k)
        ! n_max = int(t_final/dt)
        ! dx = abs(a) * dt / lambda
        ! i_max = int(Lx/dx)

        print *, dx, i_max, dt, n_max
    
        ! Allocation des tableaux
        allocate(alpha(p*i_max), alpha_np1(p*i_max), beta(p*(i_max+1)))
        allocate(u(i_max+1))

        ! Projection des conditions initiales
        do i = 1, i_max
            do j = 1, p 
                alpha((i-1)*p + j) = quad_init(i, j, p, dx, case)/norme(j)
            end do 
        end do 

        do j = 1, p 
            beta(j) = quad_bound(0, j, p, dt, Lx, a, case)/norme(j)
        end do 


        ! Implementation du schéma 
        do r=1, n_max
            do i = 1, i_max  
                if (abs(lambda)>1) then 
                    alpha_np1((i-1)*p + 1 : i*p) = matmul(mat_L, beta((i-1)*p + 1 : i*p))
                    beta(i*p + 1 : (i+1)*p) = matmul(mat_M, alpha((i-1)*p + 1 : i*p)) + matmul(mat_N, beta((i-1)*p + 1 : i*p))
                else if (abs(lambda)<1) then     
                    alpha_np1((i-1)*p + 1 : i*p) = matmul(mat_O, alpha((i-1)*p + 1 : i*p)) + matmul(mat_P, beta((i-1)*p + 1 : i*p))
                    beta(i*p + 1 : (i+1)*p) = matmul(mat_Q, alpha((i-1)*p + 1 : i*p))
                end if 
            end do 

            do j = 1, p 
                beta(j) = quad_bound(r, j, p, dt, Lx, a, case)/norme(j)
            end do 

            alpha = alpha_np1
        end do 


        ! Calcul d'erreur en norme 2 
        u = calculate_u(alpha, i_max, p, dx)
        err_N2 = 0._PR
        do i = 1, i_max+1
            err_N2 = err_N2 + (u(i) - sol_exacte(dx*(i-1), t_final, Lx, a, case))**2
        end do 
        err_N2 = sqrt((err_N2*dx))

        write(30, *) dx, err_N2


        ! Calcul de la solution et ecriture dans un fichier
        open(unit=20, file="resultats.dat", status='replace', action='write')
        ! u = calculate_u(alpha, i_max, p, dx)
        do i = 1, i_max+1
            if (a>0) then 
                write (20, *) dx*(i-1), u(i), sol_exacte(dx*(i-1), t_final, Lx, a, case)
            else 
                write (20, *) dx*(i_max-i+1), u(i), sol_exacte(dx*(i_max-i+1), t_final, Lx, a, case)
            end if 
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
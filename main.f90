program main 

    use constantes
    use Legendre
    use functions
    use matrix

    implicit none 

    integer :: i_max, n_max, p, i, k, j, r
    real(kind=PR) :: dx, dt, lambda, t_final, Lx
    real(kind=PR), dimension(:), allocatable :: alpha, alpha_np1, beta, u, u_np1
    real(kind=PR), dimension(:, :), allocatable :: mat_L, mat_M, mat_N, mat_O, mat_P, mat_Q

    ! Lecture et definition des parametres 
    open(unit=10, file="parametres.dat", action='read')
    read (10, *) p, i_max, n_max, t_final, Lx
    close(10)
    k = 1

    dx = Lx/i_max
    dt = t_final/n_max

    lambda = a*dt/dx
    print *, lambda
    
    ! Allocation des tableaux
    allocate(alpha(p*i_max), alpha_np1(p*i_max), beta(p*(i_max+1)))
    allocate(u(i_max*k+1), u_np1(i_max*k+1))

    ! Projection des conditions initiales
    do i = 1, i_max
        do j = 1, p 
            alpha((i-1)*p + j) = quad_init(i, j, p, dx)/norme(j)
        end do 
    end do 

    if (a>0) then 
        do j = 1, p 
            beta(j) = quad_L(0, j, p, dt, Lx)/norme(j)
        end do 
    else 
        do j = 1, p 
            beta(p*i_max+j) = quad_L(0, j, p, dt, Lx)/norme(j)
        end do 
    end if 
        
    ! Creation des matrices
    if (a>0) then 
        if (abs(lambda)>1) then 
            allocate(mat_L(p, p), mat_M(p, p), mat_N(p, p))
            mat_L = make_L(lambda, p)
            mat_M = make_M(lambda, p)
            mat_N = make_N(lambda, p)
        else if (abs(lambda)<1) then 
            allocate(mat_O(p, p), mat_P(p, p), mat_Q(p, p))
            mat_O = make_N(1._PR/lambda, p)
            mat_P = make_M(1._PR/lambda, p)
            mat_Q = make_L(1._PR/lambda, p)
        end if 
    else 
        if (abs(lambda)>1) then 
            allocate(mat_L(p, p), mat_M(p, p), mat_N(p, p))
            mat_L = make_L(lambda, p)
            mat_M = make_M(lambda, p)
            mat_N = make_N(abs(lambda), p)
        else if (abs(lambda)<1) then 
            allocate(mat_O(p, p), mat_P(p, p), mat_Q(p, p))
            mat_O = make_N(1._PR/lambda, p)
            mat_P = make_M(1._PR/lambda, p)
            mat_Q = make_L(1._PR/lambda, p)
        end if 
    end if 


    ! Implementation du schéma 
    do r=1, n_max
        do i = 1, i_max  
            if (a>0) then 
                if (abs(lambda)>1) then 
                    alpha_np1((i-1)*p + 1 : i*p) = matmul(mat_L, beta((i-1)*p + 1 : i*p))
                    beta(i*p + 1 : (i+1)*p) = matmul(mat_M, alpha((i-1)*p + 1 : i*p)) + matmul(mat_N, beta((i-1)*p + 1 : i*p))
                else if (abs(lambda)<1) then     
                    alpha_np1((i-1)*p + 1 : i*p) = matmul(mat_O, alpha((i-1)*p + 1 : i*p)) + matmul(mat_P, beta((i-1)*p + 1 : i*p))
                    beta(i*p + 1 : (i+1)*p) = matmul(mat_Q, alpha((i-1)*p + 1 : i*p))
                end if 
            else 
                j = i_max - i + 1
                if (abs(lambda)>1) then 
                    alpha_np1((j-1)*p + 1 : j*p) = matmul(mat_L, beta(j*p + 1 : (j+1)*p))
                    beta((j-1)*p + 1 : j*p) = matmul(mat_M, alpha((j-1)*p + 1 : j*p)) + matmul(mat_N, beta(j*p + 1 : (j+1)*p))
                else if (abs(lambda)<1) then     
                    alpha_np1((j-1)*p + 1 : j*p) = matmul(mat_O, alpha((j-1)*p + 1 : j*p)) + matmul(mat_P, beta((j-1)*p + 1 : j*p))
                    beta(j*p + 1 : (j+1)*p) = matmul(mat_Q, alpha((j-1)*p + 1 : j*p))
                end if 
            end if 
        end do 

        if (a>0) then 
            do j = 1, p 
                beta(j) = quad_L(r, j, p, dt, Lx)/norme(j)
            end do 
        else 
            do j = 1, p 
                beta(p*i_max+j) = quad_L(r, j, p, dt, Lx)/norme(j)
            end do 
        end if 

        alpha = alpha_np1
    end do 


    ! Calcul de la solution et ecriture dans un fichier
    open(unit=20, file="resultats.dat", status='replace', action='write')
    u = calculate_u(alpha, i_max, p, k, dx)
    u_np1 = calculate_u(alpha_np1, i_max, p, k, dx)
    do i = 1, i_max*k+1
        write (20, *) dx/k*(i-1), u_np1(i), sol_exacte(dx/k*(i-1), t_final, Lx)
    end do 
    close(20)

    ! Deallocation des tableaux
    if (abs(lambda)>1) then
        deallocate(mat_L, mat_M, mat_N) 
    else if (abs(lambda)<1) then
        deallocate(mat_O, mat_P, mat_Q)
    end if 
    deallocate(alpha, alpha_np1, beta)
    deallocate(u, u_np1)

end program
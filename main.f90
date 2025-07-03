program main 

    use constantes
    use functions
    use Legendre
    use quad
    use matrix

    implicit none 

    integer :: nx, p, i, j, k, n, nv, l, nt
    real(kind=PR) :: dx, dt, lambda, t_final, Lx, Rx, Mn_i, a, tau 
    real(kind=PR) :: rho_g, rho_d, u_g, u_d, T_g, T_d, eps, temps
    real(kind=PR) :: v_min, v_max, dv
    real(kind=PR), dimension(:, :), allocatable :: alpha, alpha_np1, beta, f 
    real(kind=PR), dimension(:), allocatable :: rho, u, T, vec_tau, v
    real(kind=PR), dimension(:), allocatable :: V_N, V_O, V_a_0
    real(kind=PR), dimension(:, :), allocatable :: V_L, V_M, V_P, V_Q
    real(kind=PR), dimension(:, :), allocatable ::  mat_N, mat_O, mat_a_0
    real(kind=PR), dimension(:, :, :), allocatable :: mat_L, mat_M, mat_P, mat_Q

    ! Choix de l'ordre
    p = 5

    ! Definition des donnes
    rho_g = 1.0_PR ; u_g = 0.0_PR ; T_g = 1.0_PR 
    rho_d = 0.125_PR ; u_d = 0.0_PR ; T_d = 0.8_PR
    eps = 0.1_PR

    ! Maillage en espace et en temps 
    Lx = -1.0_PR
    Rx = 1.0_PR
    nx = 100
    dx = (Rx - Lx)/nx 

    t_final = 0.2_PR
    temps = 0.0_PR

    ! Allocation des tableaux des quantites macro 
    allocate(rho(0:nx), u(0:nx), T(0:nx), vec_tau(0:nx))

    ! Maillage en vitesse
    v_min = u_g - 5 * sqrt(T_g)
    v_max = u_d + 5 * sqrt(T_d)
    nv = 50
    dv = (v_max - v_min)/nv

    ! Allocation du vecteur vitesse et initialisation
    allocate(v(0:nv))
    
    do l = 0, nv 

        v(l) = v_min + dv * l

    end do 

    dt = dx / MAXVAL(abs(v))
    nt = int(t_final/dt)

    ! Allocation des matrices L, M, N et vecteurs V_L, V_M, V_N 
    allocate(mat_L(0:p, p, p), mat_M(0:p, p, p), mat_N(p, p))
    allocate(V_L(0:p, p), V_M(0:p, p), V_N(p))
    allocate(mat_O(p, p), mat_P(0:p, p, p), mat_Q(0:p, p, p))
    allocate(V_O(p), V_P(0:p, p), V_Q(0:p, p))
    allocate(mat_a_0(p, p), V_a_0(p))

    ! Allocation des tableaux de coefficients
    allocate(alpha(p*nx, 0:nv), alpha_np1(p*nx, 0:nv), beta(p*(nx+1), 0:nv))

    ! Projection de la condition initiale
    do l = 0, nv
        do i = 1, nx
            do j = 1, p 
 
                alpha((i-1)*p + j, l) = quad_init(i, j, dx, Lx, rho_g, rho_d, T_g, T_d, v(l))/norme(j)

            end do 
        end do 
    end do 

    ! Initialisation des quantites macro
    rho = calcul_rho(nx, nv, p, dv, alpha)
    vec_tau = calcul_tau(nx, nv, p, dv, alpha)
    u = calcul_u(nx, nv, p, dv, v, alpha)
    T = calcul_T(nx, nv, p, dv, v, alpha)

    ! Boucle en temps 
    do n = 0, nt - 1 

        alpha_np1 = 0.0_PR  

        do l = 0, nv ! Boucle en vitesse 

            beta = 0.0_PR

            ! Définition de lambda et a 
            lambda = abs(v(l)) * dt / dx 
            a = v(l)

            ! initialisation des matrices
            if (a /= 0.0_PR) then 

                if (abs(lambda) >= 1.0_PR) then 

                    mat_L = make_L(p, lambda, abs(dx / a), a)
                    mat_M = make_M(p, lambda, abs(dx / a), a)
                    mat_N = make_N(p, lambda)

                    V_L = make_V_L(p, abs(dx / a), a)
                    V_M = make_V_M(p, lambda, abs(dx / a))
                    V_N = make_V_N(p, lambda)

                else if ((abs(lambda) < 1.0_PR)) then 

                    mat_O = make_O(p, lambda, a)
                    mat_P = make_P(p, lambda, dt, a)
                    mat_Q = make_Q(p, lambda, dt, a)

                    V_O = make_V_O(p, lambda, a)
                    V_P = make_V_P(p, lambda, dt, a)
                    V_Q = make_V_Q(p, dt, a)

                    
                end if 
        
            else

                mat_a_0 = make_a_0(p)
                V_a_0 = make_V_a_0(p)

            end if  

            ! Projection des conditions de bords 
            if (a > 0.0_PR) then 
                do j = 1, p 
                    beta(j, l) = quad_bound(j,rho_g, u_g, T_g, rho_d, u_d, T_d, a)/norme(j)
                end do 
            else if (a < 0.0_PR) then 
                do j = 1, p 
                    beta(p*nx + j, l) = quad_bound(j, rho_g, u_g, T_g, rho_d, u_d, T_d, a)/norme(j)
                end do 
            end if

            ! Schema GRP espace-temps

            do i = 1, nx  ! Boucle en espace 

                tau = eps * vec_tau(i)
                Mn_i = (Maxwell(rho(i-1), u(i-1), abs(T(i-1)), a) + Maxwell(rho(i), u(i), abs(T(i)), a)) / 2.0_PR  
        
                if (a > 0.0_PR) then ! Cas a positif: on parcourt le vecteur beta de gauche a droite

                    beta(i*p + 1 : (i+1)*p, l) = 0.0_PR

                    if (abs(lambda) >= 1.0_PR) then ! Cas lambda > 1 

                        do k = 0, p  ! Boucle pour la serie entiere de exp

                            alpha_np1((i-1)*p + 1 : i*p, l) = alpha_np1((i-1)*p + 1 : i*p, l) + 1.0_PR / tau**k * & 
                                                        matmul(mat_L(k, 1:p, 1:p), beta((i-1)*p + 1 : i*p, l)) + & 
                                                        Mn_i / tau**k * V_L(k, 1:p)

                            beta(i*p + 1 : (i+1)*p, l) = beta(i*p + 1 : (i+1)*p, l) + 1.0_PR / tau**k * & 
                                                    matmul(mat_M(k, 1:p, 1:p), alpha((i-1)*p + 1 : i*p, l)) + & 
                                                    Mn_i / tau**k * V_M(k, 1:p)

                        end do 

                        beta(i*p + 1 : (i+1)*p, l) = beta(i*p + 1 : (i+1)*p, l) + &
                                                exp(- dx / a / tau) * matmul(mat_N, beta((i-1)*p + 1 : i*p, l)) + &
                                                Mn_i * (1 - exp(- dx / a / tau)) * V_N

                    else if (abs(lambda) < 1.0_PR ) then ! Cas lambda > 1 

                        do k = 0, p 

                        alpha_np1((i-1)*p + 1 : i*p, l) = alpha_np1((i-1)*p + 1 : i*p, l) + 1.0_PR / tau**k * &
                                                    matmul(mat_P(k, 1:p, 1:p), beta((i-1)*p + 1 : i*p, l)) + & 
                                                    Mn_i / tau**k * V_P(k, 1:p)

                        beta(i*p + 1 : (i+1)*p, l) = beta(i*p + 1 : (i+1)*p, l) + 1.0_PR / tau**k * &
                                                matmul(mat_Q(k, 1:p, 1:p), alpha((i-1)*p + 1 : i*p, l)) + & 
                                                Mn_i / tau**k * V_Q(k, 1:p)

                        end do 

                        alpha_np1((i-1)*p + 1 : i*p, l) = alpha_np1((i-1)*p + 1 : i*p, l) + & 
                                                exp(-dt / tau) *matmul(mat_O, alpha((i-1)*p + 1 : i*p, l)) + & 
                                                Mn_i * (1 - exp(-dt / tau)) * V_O

                    end if 

                else if (a < 0.0_PR) then ! Cas a negatif: on parcourt le vecteur beta de droite a gauche 

                    j = nx - i + 1 ! Changement d'indice pour parcourir de droite à gauche

                    beta((j-1)*p + 1 : j*p, l) = 0.0_PR

                    if (abs(lambda) >= 1.0_PR) then 

                        do k = 0, p 

                            alpha_np1((j-1)*p + 1 : j*p, l) = alpha_np1((j-1)*p + 1 : j*p, l) + 1.0_PR / tau**k * & 
                                                            matmul(mat_L(k, 1:p, 1:p), beta(j*p + 1 : (j+1)*p, l)) + & 
                                                            Mn_i / tau**k * V_L(k, 1:p)

                            beta((j-1)*p + 1 : j*p, l) = beta((j-1)*p + 1 : j*p, l) + 1.0_PR / tau**k * & 
                                                    matmul(mat_M(k, 1:p, 1:p), alpha((j-1)*p + 1 : j*p, l)) + & 
                                                    Mn_i / tau**k * V_M(k, 1:p)

                        end do 

                        beta((j-1)*p + 1 : j*p, l) = beta((j-1)*p + 1 : j*p, l) + &
                                                exp(- dx / abs(a) / tau) * matmul(mat_N, beta(j*p + 1 : (j+1)*p, l)) + &
                                                Mn_i * (1 - exp(- dx / abs(a) / tau)) * V_N


                    else if (abs(lambda) < 1.0_PR) then  

                        do k = 0, p 

                            alpha_np1((j-1)*p + 1 : j*p, l) = alpha_np1((j-1)*p + 1 : j*p, l) + 1.0_PR / tau**k * &
                                                        matmul(mat_P(k, 1:p, 1:p), beta(j*p + 1 : (j+1)*p, l)) + & 
                                                        Mn_i / tau**k * V_P(k, 1:p)

                            beta((j-1)*p + 1 : j*p, l) = beta((j-1)*p + 1 : j*p, l) + 1.0_PR / tau**k * &
                                                    matmul(mat_Q(k, 1:p, 1:p), alpha((j-1)*p + 1 : j*p, l)) + & 
                                                    Mn_i / tau**k * V_Q(k, 1:p)

                        end do 

                        alpha_np1((j-1)*p + 1 : j*p, l) = alpha_np1((j-1)*p + 1 : j*p, l) + & 
                                                exp(-dt / tau) *matmul(mat_O, alpha((j-1)*p + 1 : j*p, l)) + & 
                                                Mn_i * (1 - exp(-dt / tau)) * V_O
                    end if 

                else if (a == 0.0_PR) then

                    alpha_np1((i-1)*p + 1 : i*p, l) = matmul(mat_a_0, alpha((i-1)*p + 1 : i*p, l)) * exp(-dt / tau) + &
                                                        Mn_i * (1 - exp(-dt / tau)) * V_a_0

                end if 

            end do 

        end do 

        ! Mise a jour de alpha 
        alpha = alpha_np1

        if (any(alpha /= alpha)) then
            print *, "❌ NaN détecté à n =", n
            stop 
        end if

        ! Mise a jour des quantites macro
        rho = calcul_rho(nx, nv, p, dv, alpha)
        vec_tau = calcul_tau(nx, nv, p, dv, alpha)
        u = calcul_u(nx, nv, p, dv, v, alpha)
        T = calcul_T(nx, nv, p, dv, v, alpha)

        ! Mise a jour du temps 
        temps = temps + dt 

    end do 

    ! Ouverture des fichiers pour tracer les quantites macro a tfinal
    open(unit=10, file='rho.dat', status='replace')
    open(unit=11, file='u.dat',   status='replace')
    open(unit=12, file='T.dat',   status='replace')

    ! Écriture ligne par ligne : x, valeur
    do i = 0, nx

        write(10, *) Lx + i*dx, rho(i)
        write(11, *) Lx + i*dx, u(i)
        write(12, *) Lx + i*dx, T(i)
        
    end do

    ! Fermeture des fichiers
    close(10)
    close(11)
    close(12)

end program
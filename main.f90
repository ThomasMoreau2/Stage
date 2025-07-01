program main 

    use constantes
    use functions
    use Legendre
    use quad
    use matrix

    implicit none 

    integer :: nx, p, i, j, r, k, n_max, nv, l
    real(kind=PR) :: dx, dt, lambda, t_final, Lx, Rx, a, err_N2, Mn_i, t_exacte, C
    real(kind=PR) :: rho_g, rho_d, u_g, u_d, T_g, T_d, eps, temps, n_iter
    real(kind=PR) :: v_min, v_max, dv, x
    real(kind=PR), dimension(:, :), allocatable :: alpha, alpha_np1, beta
    real(kind=PR), dimension(:), allocatable :: rho, u, T, tau, v
    real(kind=PR), dimension(:), allocatable :: V_N, V_O, V_a_0
    real(kind=PR), dimension(:, :), allocatable :: V_L, V_M, V_P, V_Q
    real(kind=PR), dimension(:, :), allocatable ::  mat_N, mat_O, mat_a_0
    real(kind=PR), dimension(:, :, :), allocatable :: mat_L, mat_M, mat_P, mat_Q

    ! Choix du degré
    p = 2

    ! Definition des donnes
    rho_g = 1.0_pr ; u_g = 0.0_pr ; T_g = 1.0_pr 
    rho_d = 0.125_pr ; u_d = 0.0_pr ; T_d = 0.8_pr 
    eps = 0.1_pr 
    temps = 0.0_pr
    n_iter = 0

    ! Maillage en espace et temps final
    Lx = -1.0_PR
    Rx = 1.0_PR
    nx = 10
    dx = (Rx - Lx)/nx 
    t_final = 0.2_pr

    ! Allocation des tableaux des quantites macro 
    allocate(rho(0:nx), u(0:nx), T(0:nx), tau(0:nx))

    ! Maillage en vitesse
    v_min = u_g - 5 * sqrt(T_g)
    v_max = u_d + 5 * sqrt(T_d)
    nv = 5
    dv = (v_max - v_min)/nv

    ! Allocation des vitesses, et des tableaux de coefficients
    allocate(v(0:nv))
    allocate(alpha(p*nx, 0:nv), alpha_np1(p*nx, 0:nv), beta(p*(nx+1), 0:nv))

    ! Initialisation du vecteur des vitesses
    do l = 0, nv 

        v(l) = v_min + dv * l

    end do 

    ! Projection de la condition initiale
    do l = 0, nv
        do i = 1, nx
            do j = 1, p 
 
                alpha((i-1)*p + j, l) = quad_init(i, j, p, dx, Lx, rho_g, rho_d, T_g, T_d, v(l))/norme(j)

            end do 
        end do 
    end do 

    ! Initialisation des quantites macro
    rho = calcul_rho(nx, nv, p, dv, alpha)
    tau = calcul_tau(nx, nv, p, dv, alpha)
    u = calcul_u(nx, nv, p, dv, v, alpha)
    T = calcul_T(nx, nv, p, dv, v, alpha)

    ! Allocation des matrices L, M, N et vecteurs V_L, V_M, V_N 
    if (a /= 0.0_PR) then 

        if (abs(lambda) >= 1.0_PR) then 

            allocate(mat_L(0:p, p, p), mat_M(0:p, p, p), mat_N(p, p))
            allocate(V_L(0:p, p), V_M(0:p, p), V_N(p))

        else if ((abs(lambda) < 1.0_PR)) then 

            allocate(mat_O(p, p), mat_P(0:p, p, p), mat_Q(0:p, p, p))
            allocate(V_O(p), V_P(0:p, p), V_Q(0:p, p))

        end if 

    else 

        allocate(mat_a_0(p, p), V_a_0(p))

    end if 

    ! Boucle en temps 

  
end program
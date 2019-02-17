
! Code by: Jostein Brandshoi

! ============================================================================================
! Compile/run on Windows:
! > gfortran -o nonlinear_swe.exe nonlinear_swe.f90
! > nonlinear_swe.exe <output_filename>
!
! Compile/run on Linux:
! > gfortran -o nonlinear_swe nonlinear_swe.f90
! > ./nonlinear_swe <output_filename>
!
! Command line arguments:
!     output_filename : Give a filename (*.dat) to store numerical solutions
!
! Purpose of program:
!   Solve nonlinear rotating shallow water equations
!
! Equations:
!   dU/dt + u * du/dx - fV = -gH_0 * dh/dx
!   dV/dt + u * dv/dx fU = 0
!   dh/dt + h * du/dx + u * dh/dx = 0
!
! Scheme (method of characteristics):
!   u_star = (u_P + 2*c_P + u_Q - 2*c_Q + ((f*dt)/2)*(2*v_R -&
!       f*dt*u_R + v_Q + v_P)) / (2 + ((f*f*dt*dt)/2))
!   v_star = v_R - ((f*dt)/2)*(u_star + u_R)
!   c_star = 0.5*(u_star - u_Q + 2*c_Q - ((f*dt)/2)*(v_star + v_Q))
!
!   u(j, n_next) = (1 - alpha_j) * u_star + alpha_j * 0.0
!   v(j, n_next) = (1 - alpha_j) * v_star + alpha_j * 0.0
!   c(j, n_next) = (1 - alpha_j) * c_star + alpha_j * dsqrt(g * H)
! ============================================================================================

PROGRAM advection_solver
    IMPLICIT NONE
    ! ------------------------------- VARIABLE DECLARATIONS -------------------------------
    double precision, parameter :: g = 10.0              ! Accelereation of gravity [m/s^2]
    double precision, parameter :: f = 1E-4              ! Coriolis parameter [s^(-1)]
    double precision, parameter :: H = 100.0             ! Equillibrium height [m]
    double precision, parameter :: dH = H / 2.0          ! Coefficient in initial condition [m]
    double precision, parameter :: L = 2000E+3           ! Length of solution domain [m]
    double precision, parameter :: kappa = 10.0 / L      ! Coefficient in initial condition [m^(-1)]

    double precision, parameter :: dx = L / 100.0        ! Distance between spatial gridpoints [m]
    double precision, parameter :: dt = 60.0             ! Time between temporal gridpoints [s]

    integer :: j, n                                      ! Counter variables (j: space, n: time)
    integer :: n_temp, n_curr = 1, n_next = 2                          ! Indices for curr and next time
    integer, parameter :: FRS_length = 10                              ! Number of points in FRS-zone
    integer, parameter :: j_max = idnint(L / dx) + 1 + 2 * FRS_length  ! Total number of space points
    integer, parameter :: n_max = 50000                  ! Some arbitrary upper bound for simulation time

    double precision :: C_Pjn, C_Qjn, C_Rjn              ! To hold local "Courant" numbers
    double precision :: u_P, v_P, c_P                    ! To hold variable values at point P
    double precision :: u_Q, v_Q, c_Q                    ! To hold variable values at point Q
    double precision :: u_R, v_R                         ! To hold variable values at point R
    double precision :: alpha_j                          ! "Relaxation function"
    double precision :: u_star, v_star, c_star           ! To hold predictor values

    double precision, dimension(j_max, 2) :: u = 0.0     ! To store u at times n and n + 1 (also sets IC)
    double precision, dimension(j_max, 2) :: v = 0.0     ! To store v at times n and n + 1 (also sets IC)
    double precision, dimension(j_max, 2) :: c = 0.0     ! To store c at times n and n + 1

    character(len = 32) :: output_filename               ! Variable to hold output filename for results
    integer, parameter :: unit_number = 10               ! Unit number for identifying output file
    ! -------------------------------------------------------------------------------------

    ! ---------- GETTING ARGS FROM CMD-LINE AND OPEN .DAT FILE FOR RESULT OUTPUT ----------
    CALL GET_COMMAND_ARGUMENT(1, output_filename)                ! Get cmd-line arg for output file
    open(unit = unit_number, file = output_filename, form = "formatted")  ! Open file that will store results
    ! -------------------------------------------------------------------------------------

    print*, j_max, n_max

    ! ------------------------------ INITIAL CONDITION FOR c ------------------------------
    do j = 1, j_max
        c(j, n_curr) = dsqrt(g * (H + dH*dtanh(kappa*(-FRS_length*dx - L/2 + (j - 1)*dx))))
    end do
    ! -------------------------------------------------------------------------------------

    ! ----------------------------- METHOD OF CHARACTERISTICS -----------------------------
    do n = 0, n_max - 1         ! Loop, in time to compute all three variables u, v, c
        do j = 1, j_max
            C_Pjn = (dt / dx) * dabs(u(j, n_curr) + c(j, n_curr))         ! Coefficient for point P
            C_Qjn = (dt / dx) * dabs((u(j, n_curr) - c(j, n_curr)))   ! Coefficient for point Q
            C_Rjn = (dt / dx) * dabs(u(j, n_curr))                    ! Coefficient for point R

            if (C_Pjn > 1 .or. C_Qjn > 1 .or. C_Rjn > 1) then    ! In case we cannot interpolate closest
                print*, "x_P, x_Q or x_R outside! Exiting..."    ! Give user some info what's happening
                call exit(1)                                     ! Terminate execution of program
            end if

            u_P = (1 - C_Pjn) * u(j, n_curr) + C_Pjn * u(j - 1, n_curr)   ! Value of u at point P
            v_P = (1 - C_Pjn) * v(j, n_curr) + C_Pjn * v(j - 1, n_curr)   ! Value of v at point P
            c_P = (1 - C_Pjn) * c(j, n_curr) + C_Pjn * c(j - 1, n_curr)   ! Value of c at point P

            u_Q = (1 - C_Qjn) * u(j, n_curr) + C_Qjn * u(j + 1, n_curr)   ! Value of u at point Q
            v_Q = (1 - C_Qjn) * v(j, n_curr) + C_Qjn * v(j + 1, n_curr)   ! Value of v at point Q
            c_Q = (1 - C_Qjn) * c(j, n_curr) + C_Qjn * c(j + 1, n_curr)   ! Value of c at point Q

            if (u(j, n_curr) >= 0) then     ! If characteristic has positive slope
                u_R = (1 - C_Rjn) * u(j, n_curr) + C_Rjn * u(j - 1, n_curr)   ! Value of u at point R
                v_R = (1 - C_Rjn) * v(j, n_curr) + C_Rjn * v(j - 1, n_curr)   ! Value of v at point R
            else                            ! If characteristic has negative slope
                u_R = (1 - C_Rjn) * u(j, n_curr) + C_Rjn * u(j + 1, n_curr)   ! Value of u at point R
                v_R = (1 - C_Rjn) * v(j, n_curr) + C_Rjn * v(j + 1, n_curr)   ! Value of v at point R
            end if

            if (j >= 2 .and. j <= j_max - 1) then   ! Compute predictor for all points except outer bounds.
                u_star = (u_P + 2*c_P + u_Q - 2*c_Q + ((f*dt)/2)*(2*v_R -&
                    f*dt*u_R + v_Q + v_P)) / (2 + ((f*f*dt*dt)/2))                  ! Predictor for u
                v_star = v_R - ((f*dt)/2)*(u_star + u_R)                            ! Predictor for u
                c_star = 0.5*(u_star - u_Q + 2*c_Q - ((f*dt)/2)*(v_star + v_Q))     ! Predictor for u
            end if

            if (j < FRS_length + 1) then            ! For left FRS-zone
                alpha_j = 1 - dtanh((dble(j) - 1.0) / 2.0)      ! From 1 at outer bound. towards zero approchnig OB
            else if (j > j_max - FRS_length) then   ! For right FRS-zone
                alpha_j = 1 - dtanh((dble(j_max - j)) / 2.0)    ! From 1 at outer bound. towards zero approchnig OB
            else                                    ! For "real" points
                alpha_j = 0.0   ! Zero for all "real" points including open boundaries
            end if

            u(j, n_next) = (1 - alpha_j) * u_star + alpha_j * 0.0            ! Corrector step for u
            v(j, n_next) = (1 - alpha_j) * v_star + alpha_j * 0.0            ! Corrector step for v

            if (j < FRS_length + 1) then
                c(j, n_next) = (1 - alpha_j) * c_star + alpha_j * dsqrt(g * (H-dH))  ! Corrector c in left FRS-zone
            else if (j > j_max - FRS_length) then
                c(j, n_next) = (1 - alpha_j) * c_star + alpha_j * dsqrt(g * (H+dH))  ! Corrector c in right FRS-zone
            else
                c(j, n_next) = c_star   ! Corrector c in left "real" domain
            end if
        end do

        if (modulo(n, 10) == 0) then    ! Write line-wise to file. Write only every 10'th time level
            write(unit = 10, fmt = "(303f20.14)") u(FRS_length + 1:j_max -&
                FRS_length, n_curr),&
                v(FRS_length + 1:j_max - FRS_length, n_curr),&
                c(FRS_length + 1:j_max - FRS_length, n_curr)
        end if

        n_temp = n_next     ! Set help index to next time step column column
        n_next = n_curr     ! Set next time step index to current ime step column
        n_curr = n_temp     ! Set current time step index to next time step column
    end do
    ! -------------------------------------------------------------------------------------

    close(unit = unit_number)    ! Close file after writing

END PROGRAM advection_solver


! Code by: Jostein Brandshoi

! This program uses the upwind scheme together with cyclic boundary conditions and MPDATA method to solve the
! advection equation for the region x in [0, j_max] and t > 0. MPDATA is only used when specified by the user
! through the command line. Below is info on how to compile (with gfortran) and run the program (for Windows).
!
! > gfortran -o advection_solver.exe advection_solver.f90
! > advection_solver.exe <method_type> <differentation_type> <output_filename>
!
! Command line args (in order):
!     method_type         : Choose either "mpdata" or an arbitrary non-empy string. Determines if MPDATA is
!                           used or not. something different than "mpdata" is given, ordinary upwind is used
!     differentation_type : Choose between "forward", "centered" or "backward". Determines what FDA is
!                           used for computing the anti-diffusive velocity u_star.
!     output_filename     : Give a filename (*.dat) to store the theta values at releveant time levels
!                           for plotting.

PROGRAM advection_solver
    IMPLICIT NONE
    ! ------------------------------- VARIABLE DECLARATIONS -------------------------------
    double precision, parameter :: L = 50000.0                      ! Location (L/2) where initial bell is max [m]
    double precision, parameter :: sigma = L / 10.0                 ! Measure of with of initial Gaussian bell [m]
    double precision, parameter :: theta_0 = 10.0                   ! Max temperature value of initial bell [deg C]
    double precision, parameter :: u_0 = 1.0                        ! Advection speed (constant here) [m/s]
    double precision, parameter :: S_c = 1.1                        ! Scaling factor for use in MPDATA
    double precision, parameter :: epsilon = 1E-15                  ! Smalln to avoid zero division in MPDATA

    double precision, parameter :: C = 0.50                         ! Courant number [u_0 * delta_t / delta_x] [m/s]
    double precision, parameter :: delta_x = sigma / 10.0           ! Distance between spatial gridpoints [m]
    double precision, parameter :: delta_t = C * delta_x / u_0      ! Time between temporal gridpoints [s]

    integer :: i, j, n                                              ! Counter variables (i: iters, j: space, n: time)
    integer :: n_temp, n_curr = 1, n_next = 2                       ! Indices to represent curr and next time stetp
    integer, parameter :: num_iterations = 1
    integer, parameter :: j_max = idint(L / delta_x + 1)            ! Max number of space points in [0: j = 1, L: j = j_max]
    integer, parameter :: n_max = idint((25 * L) / (u_0 * delta_t)) ! Max number of time points needed for 25 cycles
    integer, parameter :: n_05c = idint((5 * L) / (u_0 * delta_t))  ! Number of time points needed for 5 cycles
    integer, parameter :: n_10c = idint((10 * L) / (u_0 * delta_t)) ! Number of time points needed for 10 cycles
    integer, parameter :: n_15c = idint((15 * L) / (u_0 * delta_t)) ! Number of time points needed for 15 cycles
    integer, parameter :: n_20c = idint((20 * L) / (u_0 * delta_t)) ! Number of time points needed for 20 cycles

    double precision, dimension(j_max) :: theta_star = 0.0, u_star = 0.0, F_star = 0.0
    double precision, dimension(j_max, 2) :: theta = 0.0            ! 2D-array to store solution at time points n and n + 1

    character(len = 32) :: method_type                              ! String to specify which algorithm to use
    character(len = 32) :: differentation_type                      ! String to specify which FDA to use in MPDATA
    character(len = 32) :: output_filename                          ! Variable to hold output filename for results

    logical :: use_scaling = .TRUE.                                ! Boolean that specifies if scaling is used
    ! -------------------------------------------------------------------------------------

    ! -------------------------------- INITIAL CONDITION ----------------------------------
    ! Compute initial condition for theta and where x_j = (j - 1) * delta_x.
    do j = 1, j_max
        theta(j, n_curr) = theta_0 * dexp(-((2 * ((j - 1) * delta_x) - L) / sigma) ** 2)
    end do
    ! -------------------------------------------------------------------------------------

    ! ---------- GETTING ARGS FROM CMD-LINE AND OPEN .DAT FILE FOR RESULT OUTPUT ----------
    CALL GET_COMMAND_ARGUMENT(1, method_type)                    ! Get cmd-line arg for method type
    CALL GET_COMMAND_ARGUMENT(2, differentation_type)            ! Get cmd-line arg for diff type file
    CALL GET_COMMAND_ARGUMENT(3, output_filename)                ! Get cmd-line arg for output file
    open(unit = 10, file = output_filename, form = "formatted")  ! Open file that will store results
    ! -------------------------------------------------------------------------------------

    ! ---------------- USING UPWIND SCHEME, BC'S AND MPDATA IF REQUESTED ------------------
    do n = 0, n_max     ! Loop, in time as long as necessary to complete requested number of cycles
        theta(1, n_curr) = theta(j_max, n_curr)     ! BC theta at the left boaundary x = 0

        ! Use upwind scheme to compute a predictor to theta at time level n + 1
        do j = 2, j_max
            theta_star(j) = theta(j, n_curr) - C * (theta(j, n_curr) - theta(j - 1, n_curr))
        end do

        if (method_type == "mpdata") then   ! Employ MPDATA method when specified
            do i = 1, num_iterations        ! Use MPDATA for num_iterations iterations
                theta_star(1) = theta_star(j_max)   ! BC for theta_star for each iteration
                u_star(1) = u_star(j_max)           ! BC for u_star for each iteration
                F_star(1) = F_star(j_max)           ! BC for F_star for each iteration

                ! Compute u_star (at all interior points) if theta_star > 0 with requested FDA.
                do j = 2, j_max - 1
                    if (theta_star(j) > 0) then
                        if (differentation_type == "centered") then
                            u_star(j) = 0.25 * (1 - C) * u_0 * ((theta_star(j + 1) - &
                                theta_star(j - 1)) / (theta_star(j) + epsilon))
                        else if (differentation_type == "forward") then
                            u_star(j) = 0.50 * (1 - C) * u_0 * ((theta_star(j + 1) - &
                                theta_star(j)) / (theta_star(j) + epsilon))
                        else if (differentation_type == "backward") then
                            u_star(j) = 0.50 * (1 - C) * u_0 * ((theta_star(j) - &
                                theta_star(j - 1)) / (theta_star(j) + epsilon))
                        end if
                    else
                        u_star(j) = 0.0     ! u_star set to zero when theta_star <= 0
                    end if

                    ! Scale anti-diffusive velocity when scaled approach is used.
                    if (use_scaling .eqv. .TRUE.) then
                        u_star(j) = u_star(j) * S_c
                    end if
                end do

                ! Compute u_star (at the j_max using BC) if theta_star > 0.
                if (theta_star(j_max) > 0) then
                    if (differentation_type == "centered") then
                        u_star(j_max) = 0.25 * (1 - C) * u_0 * ((theta_star(1) - &
                            theta_star(j_max - 1)) / (theta_star(j_max) + epsilon))
                    else if(differentation_type == "forward") then
                        u_star(j_max) = 0.50 * (1 - C) * u_0 * ((theta_star(1) - &
                        theta_star(1)) / (theta_star(1) + epsilon))
                    else if (differentation_type == "backward") then
                        u_star(j_max) = 0.50 * (1 - C) * u_0 * ((theta_star(1) - &
                            theta_star(j_max - 1)) / (theta_star(1) + epsilon))
                    end if
                else
                    u_star(j_max) = 0.0     ! u_star set to zero when theta_star <= 0
                end if

                ! Scale anti-diffusive velocity when scaled approach is used.
                if (use_scaling .eqv. .TRUE.) then
                    u_star(j) = u_star(j) * S_c
                end if

                ! Compute F_star at all interior points.
                do j = 2, j_max - 1
                    F_star(j) = 0.5 * ((u_star(j) + dabs(u_star(j))) * &
                        theta_star(j) + (u_star(j + 1) - dabs(u_star(j + 1))) * &
                        theta_star(j + 1)) * (delta_t / delta_x)
                end do

                ! Compute F_star at the j_max using BC.
                F_star(j_max) = 0.5 * ((u_star(j_max) + dabs(u_star(j_max))) * &
                    theta_star(j_max) + (u_star(1) - dabs(u_star(1))) * &
                    theta_star(1)) * (delta_t / delta_x)

                ! Update theta_star to be the corrected value
                do j = 2, j_max
                    theta_star(j) = theta_star(j) - (F_star(j) - F_star(j - 1))
                end do

                if (use_scaling .eqv. .TRUE.) then
                    exit    ! Exit iteration loop after one iteration if scaling approach is used
                end if
            end do
        end if

        do j = 2, j_max
            theta(j, n_next) = theta_star(j)    ! Set theta at time level n + 1 equal theta_star
        end do

        ! Write theta to file for time steps corresponding to the requested cycles.
        if (n == 0 .or. n == n_05c .or. n == n_10c .or. n == n_15c .or. n == n_20c .or. n == n_max) then
            write(unit = 10, fmt = "(102f20.14)") float(n), theta(:, n_curr) / theta_0  ! Write line-wise
        end if

        n_temp = n_next     ! Set help index to second theta column (2)
        n_next = n_curr     ! Set next time step index to first theta column (1)
        n_curr = n_temp     ! Set current time step index to second column (2)
    end do
    ! -------------------------------------------------------------------------------------

    close(unit = 10)    ! Close file after writing

END PROGRAM advection_solver

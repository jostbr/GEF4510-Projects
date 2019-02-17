
! Code by: Jostein Brandshoi

PROGRAM advection_solver
    IMPLICIT NONE
    ! ------------------------------- VARIABLE DECLARATIONS -------------------------------
    double precision :: upwind_scheme                               ! Function declaration for upwind scheme

    double precision, parameter :: L = 50000.0                      ! Location (L/2) where initial bell is max [m]
    double precision, parameter :: sigma = L / 10.0                 ! Measure of with of initial Gaussian bell [m]
    double precision, parameter :: theta_0 = 10.0                   ! Max temperature value of initial bell [deg C]
    double precision, parameter :: u_0 = 1.0                        ! Advection speed (constant here) [m/s]
    double precision, parameter :: epsilon = 1E-15

    double precision, parameter :: C = 0.50                         ! Courant number [u_0 * delta_t / delta_x] [m/s]
    double precision, parameter :: delta_x = sigma / 10.0           ! Distance between spatial gridpoints [m]
    double precision, parameter :: delta_t = (C * delta_x) / u_0    ! Time between temporal gridpoints [s]
    double precision :: x_j                                         ! Helper variable to increase readability in IC-loop

    integer :: j, n                                                 ! Counter variables (j: space, n: time)
    integer :: n_temp, n_curr = 1, n_next = 2                       ! Indices to represent curr and next time stetp
    integer, parameter :: j_max = idint((L / delta_x) + 1)          ! Max number of space points in [0: j = 1, L: j = j_max]
    integer, parameter :: n_max = idint((25 * L) / (u_0 * delta_t)) ! Max number of time points needed for 25 cycles
    integer, parameter :: n_05c = idint((5 * L) / (u_0 * delta_t))  ! Number of time points needed for 5 cycles
    integer, parameter :: n_10c = idint((10 * L) / (u_0 * delta_t)) ! Number of time points needed for 10 cycles
    integer, parameter :: n_15c = idint((15 * L) / (u_0 * delta_t)) ! Number of time points needed for 15 cycles
    integer, parameter :: n_20c = idint((20 * L) / (u_0 * delta_t)) ! Number of time points needed for 20 cycles

    double precision :: theta_jp2_star, theta_jp1_star, theta_j_star                ! Helper variables for MPDATA method
    double precision :: theta_jm1_star, theta_jm2_star              ! Helper variables for MPDATA method
    double precision :: u_jp1_star, u_j_star, u_jm1_star, F_j_star, F_jm1_star  ! Helper variables for MPDATA method
    double precision, dimension(j_max, 2) :: theta = 0.0            ! 2D-array to store solution at time points n and n + 1

    character(len = 32) :: output_filename                          ! Variable to hold output filename for results
    character(len = 32) :: method_type                              ! String to specify which algorithm to use
    ! -------------------------------------------------------------------------------------

    ! -------- GETTING FILENAME FROM CMD-LINE AND OPEN .DAT FILE FOR RESULT OUTPUT --------
    CALL GET_COMMAND_ARGUMENT(1, output_filename)                ! Get cmd-line arg for output file
    CALL GET_COMMAND_ARGUMENT(2, method_type)                    ! Get cmd-line arg for output file
    open(unit = 10, file = output_filename, form = "formatted")  ! Open file that will store results
    ! -------------------------------------------------------------------------------------
    
    ! -------------------------------- INITIAL CONDITION ----------------------------------
    do j = 1, j_max                 ! Loop over all space points in the defined region [0, L]
        x_j = (j - 1) * delta_x     ! x-value corresponding to the j'th gridpoint in space
        theta(j, n_curr) = theta_0 * dexp(-((2 * x_j - L) / sigma) ** 2)    ! Initial condition
    end do
    ! -------------------------------------------------------------------------------------

    ! --------------------------- USING UPWIND SCHEME AND BC'S ----------------------------
    do n = 0, n_max     ! Loop, in time as long as necessary to complete requested number of cycles
        theta(1, n_curr) = theta(j_max, n_curr)     ! Cyclic boundary condition (theta(0, t) = theta(L, t))

        do j = 2, j_max ! Loop over all space points in the defined region [0, L]
            theta_j_star = upwind_scheme(theta, j, j - 1, n_curr, j_max, C)  ! Upwind scheme
            F_j_star = 0.0; F_jm1_star = 0.0;   ! Set to zero in case MPDATA is not used

            ! Employ MPDATA method when specified through cmd-line argument
            if (method_type == "mpdata" .and. theta(j, n_curr) > 0) then
                if (j == 2) then
                    theta_jp2_star = upwind_scheme(theta, 4, 3, n_curr, j_max, C)
                    theta_jp1_star = upwind_scheme(theta, 3, 2, n_curr, j_max, C)
                    theta_jm1_star = upwind_scheme(theta, 1, j_max, n_curr, j_max, C)
                    theta_jm2_star = upwind_scheme(theta, j_max, j_max - 1, n_curr, j_max, C)
                else if (j == 3) then
                    theta_jp2_star = upwind_scheme(theta, 5, 4, n_curr, j_max, C)
                    theta_jp1_star = upwind_scheme(theta, 4, 3, n_curr, j_max, C)
                    theta_jm1_star = upwind_scheme(theta, 2, 1, n_curr, j_max, C)
                    theta_jm2_star = upwind_scheme(theta, 1, j_max, n_curr, j_max, C)
                else if (j == j_max) then
                    theta_jp2_star = upwind_scheme(theta, 2, 1, n_curr, j_max, C)
                    theta_jp1_star = upwind_scheme(theta, 1, j_max, n_curr, j_max, C)
                    theta_jm1_star = upwind_scheme(theta, j_max - 1, j_max - 2, n_curr, j_max, C)
                    theta_jm2_star = upwind_scheme(theta, j_max - 2, j_max - 3, n_curr, j_max, C)
                else
                    theta_jp2_star = upwind_scheme(theta, j + 2, j + 1, n_curr, j_max, C)
                    theta_jp1_star = upwind_scheme(theta, j + 1, j, n_curr, j_max, C)
                    theta_jm1_star = upwind_scheme(theta, j - 1, j - 2, n_curr, j_max, C)
                    theta_jm2_star = upwind_scheme(theta, j - 2, j - 3, n_curr, j_max, C)
                end if

                u_jp1_star = 0.25 * (1 - C) * u_0 * ((theta_jp2_star - theta_j_star) / (theta_jp1_star + epsilon))
                u_j_star = 0.25 * (1 - C) * u_0 * ((theta_jp1_star - theta_jm1_star) / (theta_j_star + epsilon))
                u_jm1_star = 0.25 * (1 - C) * u_0 * ((theta_j_star - theta_jm2_star) / (theta_jm1_star + epsilon))

                F_j_star = 0.5 * ((u_j_star + dabs(u_j_star)) * theta_j_star + (u_jp1_star - &
                    dabs(u_jp1_star)) * theta_jp1_star) * (delta_t / delta_x)
                F_jm1_star = 0.5 * ((u_jm1_star + dabs(u_jm1_star)) * theta_jm1_star + &
                    (u_j_star - dabs(u_j_star)) * theta_j_star) * (delta_t / delta_x)
            end if

            theta(j, n_next) = theta_j_star - (F_j_star - F_jm1_star)
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

double precision function upwind_scheme(theta, j, jm1, n_curr, j_max, C)
    double precision, dimension (j_max, 2) :: theta  ! 2D theta array
    double precision C                               ! Courant number
    integer j, jm1                                   ! j'th and (j - 1)'th point

    upwind_scheme = theta(j, n_curr) - C * (theta(j, n_curr) - theta(jm1, n_curr)) ! Upwind scheme
    return
end



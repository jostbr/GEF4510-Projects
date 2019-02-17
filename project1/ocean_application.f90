
! Program by: Jostein Brandshoi

PROGRAM exercise_2d
    ! This program solves the diffusion equation numerically for the oceanic boundary layer with the initial- and
    ! boundary conditions given in the exercise. After execution a file (with filename specified as a command line
    ! argument) will be outputed with the results for specific time-values. Thi program can be compiled and executed
    ! by itself just fine, but the script visualize_results.py automates the process and provides plot of the results.
    IMPLICIT NONE

    ! ------------------------------- VARIABLE DECLARATIONS -------------------------------
    integer :: j, n                                                     ! Counter variables used in loops below (j: space, n: time)
    integer, parameter :: j_max = 27, n_max = 1201                      ! Maximum number of gridpoints in space and time

    double precision, parameter :: D = 30.0                                         ! Hight of the oceanic boundary layer [m]
    double precision, parameter :: psi_0 = 10.0                                     ! Steady max temperature at top of ocean [deg C]
    double precision, parameter :: K = 0.45                                         ! K = (kappa * delta_t) / delta_z [dim-less]
    double precision, parameter :: kappa_O = 30.0 * 1E-4                            ! Diffusion coefficient for the ocean [m^2/s]
    double precision, parameter :: delta_z = D / (j_max - 1.0)                      ! Increment in space (j_max points implies (j_max - 1) steps) [m]
    double precision, parameter :: delta_t = (K * delta_z ** 2) / kappa_O           ! Increment in time [s]
    double precision, parameter :: gmma = 1.5, t_c = 6.0 * 24 * 3600                ! Gamma parameter [dim-less], Reference time of 6 days [s]

    double precision, dimension(1 : j_max) :: psi_n, psi_np1, psi_temp
    double precision :: t_n                                                         ! t-value in [0, H] to represent time. Used in boundary con. [m]

    character(len = 32) :: output_filename                              ! Variable to hold output filename specified by cmd-line arg
    ! -------------------------------------------------------------------------------------

    ! -------------------------- INITIAL AND BOUNDARY CONDITIONS --------------------------
    psi_n(:) = 0.0         ! Initially the water column n has 0 deg C for all the z-values
    ! -------------------------------------------------------------------------------------

    ! ------- GETTING FILENAME FROM CMD-LINE AND WRITE RESULTS TO .DAT-FILE. EACH LINE CORRESPONDS TO ONE N-VALUES -------
    CALL GET_COMMAND_ARGUMENT(1, output_filename)                  ! Get command line argument that specifies output file
    open(unit = 10, file = output_filename, form = "formatted")    ! Open file that will store results based on cmd-line arg
    ! -------------------------------------------------------------------------------------

    ! ------------------------ SOLVING THE EQUATION WITH FTCS-SCHEME ----------------------
    do n = 0, n_max - 2
        do j = 2, j_max - 1
            t_n = n * delta_t                                     ! Compute t-value at n'th time point
            psi_np1(1) = 0.0                                      ! The bottom (z = -D) is held at 0 deg C for all times
            psi_np1(j_max) = psi_0 * tanh((gmma * t_n) / t_c)     ! Compute psi at top (z = 0) for each t_n according to BC.
            psi_np1(j) = psi_n(j) + K * (psi_n(j + 1) - 2 * psi_n(j) + psi_n(j - 1))
        end do

        ! Write psi for the requested n-values (plus some more) to a .dat file where each line is on the form "n psi(1) psi(2) ..."
        if (n == 0 .or. n == 100 .or. n == 200 .or. n == 400 .or. n == 600 .or. n == 800 .or. n == 1000 .or. n == 1199) then
            write(unit = 10, fmt = "(28f20.14)") float(n), psi_n(:) / psi_0
        end if

        psi_n = psi_np1     ! Update psi_n to be the value of the newly computed timestep
    end do
    ! -------------------------------------------------------------------------------------

    close(unit = 10)! Close file after writing

END PROGRAM exercise_2d

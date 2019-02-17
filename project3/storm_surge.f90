
! Code by: Jostein Brandshoi

! ---------------------------------------------------------------------------------------------------------------
! Compile/run on Windows:
! > gfortran -o storm_surge.exe storm_surge.f90
! > storm_surge.exe <output_filename>
!
! Compile/run on Linux:
! > gfortran -o storm_surge storm_surge.f90
! > ./storm_surge <output_filename>
!
! Command line args:
!     output_filename : Give a filename (*.dat) to store numerical solutions
!
! Purpose of program:
!   Solve exercise k) (linear rotating (with wind stress and friction) shallow water
!   equations applied to a storm surge problem) on Problem Set 10 in GEF4510
!
! Equations:
!   dU/dt - fV = -gH_0 * dh/dx + (tau_sx - tau_bx) / rho_0
!   dV/dt + fU = (tau_sy - tau_by) / rho_0
!   dh/dt \ dU/dx = 0
!
! Scheme (FBTCS):
!   U(j, n+1) = U(j, n) + dt * f * V(j, n) - ((g * H_0 * dt) / dx) *&
!       (h(j + 1, n) - h(j, n)) + (dt / rho_0) * (tau_sx - tau_bx_coeff * U(j, n))
!
!   V(j, n+1) = V(j, n) - dt * f * U(j, n+1) + (dt / rho_0) * &
!       (tau_sy - tau_by_coeff * V(j, n))
!
!   h(j, n+1) = h(j, n) - (dt / dx) * (U(j, n+1) - U(j - 1, n+1))
! ---------------------------------------------------------------------------------------------------------------

PROGRAM advection_solver
    IMPLICIT NONE
    ! ------------------------------- VARIABLE DECLARATIONS -------------------------------
    double precision, parameter :: g = 10.0              ! Accelereation of gravity [m/s^2]
    double precision, parameter :: rho_0 = 1E+3          ! Uniform constant density [kg/m^3]
    double precision, parameter :: f = 1E-4              ! Coriolis parameter [s^(-1)]
    double precision, parameter :: R = 2.4E-3            ! Factor used in bottom stress [m/s]
    double precision, parameter :: H_0 = 300.0           ! Equillibrium height [m]
    double precision, parameter :: c_0 = dsqrt(g * H_0)  ! Phase speed [m/s]
    double precision, parameter :: L_R = c_0 / f         ! Rossby's deformation radius [m]
    double precision, parameter :: L = L_R * 10.0        ! Max distance from coast [m]
    double precision, parameter :: tau_sx = 0.0          ! Surface wind stress (x-dir) [Pa]
    double precision, parameter :: tau_sy = 0.1          ! Surface wind stress (y-dir) [Pa]
    double precision, parameter :: tau_bx_coeff = (rho_0 * R) / H_0  ! Bottom stress coeff (x-dir) [kg/m^3s]
    double precision, parameter :: tau_by_coeff = (rho_0 * R) / H_0  ! Bottom stress coeff (y-dir) [kg/m^3s]

    double precision, parameter :: C = 0.80              ! Courant number [c_0 * delta_t / delta_x] [m/s]
    double precision, parameter :: dx = L_R / 10.0       ! Distance between spatial gridpoints [m]
    double precision, parameter :: dt = C * dx / c_0     ! Time between temporal gridpoints [s]

    integer :: j, n                                      ! Counter variables (j: space, n: time)
    integer :: n_temp, n_curr = 1, n_next = 2            ! Indices to represent curr and next time
    integer, parameter :: j_max = idnint(L / dx) + 1     ! Total number of space points in domain
    integer, parameter :: n_max = 1000                   ! Some arbitrary upper bound for simulation time

    double precision, dimension(j_max, 2) :: U = 0.0     ! To store U at times n and n + 1 (also sets IC)
    double precision, dimension(j_max, 2) :: V = 0.0     ! To store V at times n and n + 1 (also sets IC)
    double precision, dimension(j_max, 2) :: h = H_0     ! To store h at times n and n + 1 (also sets IC)

    character(len = 32) :: output_filename               ! Variable to hold output filename for results
    integer, parameter :: unit_number = 10               ! Unit number for identifying output file
    ! -------------------------------------------------------------------------------------

    ! ---------- GETTING ARGS FROM CMD-LINE AND OPEN .DAT FILE FOR RESULT OUTPUT ----------
    CALL GET_COMMAND_ARGUMENT(1, output_filename)                ! Get cmd-line arg for output file
    open(unit = unit_number, file = output_filename, form = "formatted")  ! Open file that will store results
    ! -------------------------------------------------------------------------------------

    ! ---------------------------- USING FBTCS SCHEME AND BC'S ----------------------------
    do n = 0, n_max - 1         ! Loop, in time to compute all three variables U, V, h
        U(j_max, n_next) = 0.0  ! BC: No flow through coast
        h(1, n_next) = H_0      ! BC: Ekman solution at position far away from coast

        do j = 1, j_max - 1     ! Loop in space to compute U_j^(n+1)
            U(j, n_next) = U(j, n_curr) + dt * f * V(j, n_curr) - ((g * H_0 * dt) / dx) *&
                (h(j + 1, n_curr) - h(j, n_curr)) + (dt / rho_0) * (tau_sx - tau_bx_coeff * U(j, n_curr))
        end do

        do j = 1, j_max         ! Loop in space to compute V_j^(n+1)
            V(j, n_next) = V(j, n_curr) - dt * f * U(j, n_next) + (dt / rho_0) * &
                (tau_sy - tau_by_coeff * V(j, n_curr))
        end do

        do j = 2, j_max         ! Loop in space to compute h_j^(n+1)
            h(j, n_next) = h(j, n_curr) - (dt / dx) * (U(j, n_next) - U(j - 1, n_next))
        end do

        ! Write line-wise to file. Write only solution for half the domain (closest to coast) to file.
        write(unit = 10, fmt = "(157f20.14)") float(n), U(50:, n_curr), V(50:, n_curr), h(50:, n_curr)

        n_temp = n_next     ! Set help index to next time step column column
        n_next = n_curr     ! Set next time step index to current ime step column
        n_curr = n_temp     ! Set current time step index to next time step column
    end do
    ! -------------------------------------------------------------------------------------

    close(unit = unit_number)    ! Close file after writing

END PROGRAM advection_solver

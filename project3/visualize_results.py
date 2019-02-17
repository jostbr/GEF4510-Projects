
# Code by: Jostein Brandshoi

"""
This program is written while using Windows, but should work in other OS'es as well, perhaps with a minor change in how
to run the executable resulting from compiling the fortran source. Another note: Here Python 2.x is used, however it should
run without problems in a Python 3.x environment as well.

This script that tells the operating system to compile the requested (by first command line argument) .f90 file with with the desired
command line arguments (see below). Then, if compilation successful, it extracts the result data from the specified file and plots the
results using matplotlib. Extraction of the data from the file is based on knowing the format (line-wise) outputed by the .f90-file.
Below is a table, with description, over the command line arguments for this script:

sys.argv[0] : Name of the program itself (visualize_results.py)
sys.argv[1] : Name of the .f90-file without the extension
sys.argv[2] : Name of the result output file, preferrably *.dat

Below is an overview on how to run the script with appropriate command line arguments.
> python visualize_results.py storm_surge <output_filename>
"""

import os, sys, time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Function that plots the contours from the given 2D-array psi in xt-space.
def plot_contours(x, t, psi, x_scale, t_scale, num_contours, color_map, title):
    plt.figure()

    domain = (x[0] / x_scale, x[-1] / x_scale, t[0] / t_scale, t[-1] / t_scale)

    background_image = plt.imshow(psi, interpolation = "bilinear", origin = "lower",
        cmap = color_map, extent = domain, aspect = "auto")

    contour_plot = plt.contour(x / x_scale, t / t_scale, psi, num_contours,
        origin = "lower", linewidths = 2, extent = domain)

    plt.clabel(contour_plot, inline = 1, fontsize = 12)
    plt.xlabel("Space coordinate: $x$ [km]", fontsize = 16)
    plt.ylabel("Time: $t$ [h]", fontsize = 16)
    plt.title(title, fontsize = 18)

# Function that computes the analtyic solutions based on exercise d)
def get_analytical_solutions(x, t):
    x, t = np.meshgrid(x, t)            # Generate mesh in x,t-space
    U = np.zeros((n_max, j_max))
    V = np.zeros((n_max, j_max))
    h = np.zeros((n_max, j_max))

    U_E = (tau_sy / (rho_0 * f))        # Ekman transport

    U = U_E * (1 - np.exp(x / L_R))
    V = f * t * U_E * np.exp(x / L_R)
    h = H_0 * (1 + ((t * U_E) / (L_R * H_0)) * np.exp(x / L_R))

    return U, V, h

if (__name__ == "__main__"):
    try:
        fortran_filename = sys.argv[1]
        output_filename = sys.argv[2]

    except Exception:
        sys.exit("\nInvalid use of command line arguments. See docstring in source.\n")

    if (not os.system("gfortran -o {0}.exe {0}.f90".format(fortran_filename))):
        time_start = time.clock()
        os.system("{0}.exe {1}".format(fortran_filename, output_filename))
        print("Execution time: {0} seconds".format(time.clock() - time_start))

        g = 10.0; rho_0 = 1000.0; H_0 = 300.0; c_0 = np.sqrt(g * H_0); f = 1E-4; L_R = c_0 / f;
        tau_sy = 0.1; L = L_R * 10.0; C = 0.80; dx = L_R / 10.0; dt = C * dx / c_0;
        j_max = 52; n_max = 1000;

        x = np.linspace(-L / 2.0, 0, j_max)     # Only plot for half the domain
        t = np.linspace(0, n_max * dt, n_max)   # Plot for all times
        U = np.zeros((n_max, j_max))
        V = np.zeros((n_max, j_max))
        h = np.zeros((n_max, j_max))

        with open(output_filename, "r") as results_file:
            for line_num, current_line in enumerate(results_file):
                list_of_nUVh = current_line.split()                   # List of n, U, V and h strings
                current_n = int(float(list_of_nUVh[0]))               # Read current n-value
                U[line_num, :] = np.array([float(U_j) for U_j in list_of_nUVh[1:j_max + 1]])                # U at time level n
                V[line_num, :] = np.array([float(V_j) for V_j in list_of_nUVh[j_max + 1 : 2 * j_max + 1]])  # V at time level n
                h[line_num, :] = np.array([float(h_n) for h_n in list_of_nUVh[2 * j_max + 1:]])             # h at time level n

        plot_contours(x, t, U, 1000.0, 3600.0, 15, cm.hot, "Numerical solution for $U(x,t)$")    # Plot numerical U
        plot_contours(x, t, V, 1000.0, 3600.0, 15, cm.hot, "Numerical solution for $V(x,t)$")    # Plot numerical V
        plot_contours(x, t, h, 1000.0, 3600.0, 15, cm.BuPu, "Numerical solution for $h(x,t)$")   # Plot numerical h

        x_a = np.linspace(-L / 4.0, 0, j_max)
        U_a, V_a, h_a = get_analytical_solutions(x_a, t)

        plot_contours(x_a, t, U_a, 1000.0, 3600.0, 15, cm.hot, "Analytic solution for $U(x)$")     # Plot analytic U
        plot_contours(x_a, t, V_a, 1000.0, 3600.0, 15, cm.hot, "Analytic solution for $V(x,t)$")   # Plot analytic V
        plot_contours(x_a, t, h_a, 1000.0, 3600.0, 15, cm.BuPu, "Analytic solution for $h(x,t)$")  # Plot analytic h

        plt.show()

    else:
        print("\n---------- Compilation of '{0}.f90' failed ----------\n".format(sys.argv[1]))

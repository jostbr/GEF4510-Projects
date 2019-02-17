
# Code by: Jostein Brandshoi

"""
This program is written while using Windows, but should work in other OS'es as well, perhaps with a
minor change in how to run the executable resulting from compiling the fortran source. Another note:
Here Python 2.x is used, however it should run without problems in a Python 3.x environment as well.

This script that tells the operating system to compile the requested (by first command line argument)
.f90 file with with the desired command line arguments (see below). Then, if compilation successful,
it extracts the result data from the specified file and plots the results using matplotlib. Extraction
of the data from the file is based on knowing the format (line-wise) outputed by the .f90-file. Below
is a table, with description, over the command line arguments for this script:

sys.argv[0] : Name of the program itself (visualize_results.py)
sys.argv[1] : Name of the .f90-file without the extension
sys.argv[2] : Name of the result output file, preferrably *.dat

Below is an overview on how to run the script with appropriate command line arguments.
> python visualize_results.py nonlinear_swe <output_filename>
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

        g = 10.0
        L = 2000E+3;
        dx = L / 100.0
        dt = 60.0
        j_max = int(L / dx) + 1
        n_max = 50000
        num_lines = 5000

        x = np.linspace(- L / 2.0, L / 2.0, j_max)
        t = np.linspace(0, n_max * dt, num_lines)
        u = np.zeros((num_lines, j_max))
        v = np.zeros((num_lines, j_max))
        c = np.zeros((num_lines, j_max))


        plt.ion()

        with open(output_filename, "r") as results_file:
            for line_num, current_line in enumerate(results_file):
                list_of_nuvh = current_line.split()
                u[line_num, :] = np.array([float(u_j) for u_j in list_of_nuvh[0 : j_max]])
                v[line_num, :] = np.array([float(v_j) for v_j in list_of_nuvh[j_max: 2 * j_max]])
                c[line_num, :] = np.array([float(c_j) for c_j in list_of_nuvh[2 * j_max:]])
                plt.axis([-L/2, L/2, -12, 10])
                plt.plot(x, u[line_num, :], linewidth = 2.0, label = "Time: {0:.1f} hours".format(line_num * 10 * dt / 3600))
                plt.legend()
                plt.pause(0.005)
                plt.clf()

        #plot_contours(x, t, u, 1000.0, 3600.0, 8, cm.hot, "Numerical solution for $u(x,t)$")
        #plot_contours(x, t, v, 1000.0, 3600.0, 12, cm.hot, "Numerical solution for $v(x,t)$")
        #plot_contours(x, t, c**2/g, 1000.0, 3600.0, 17, cm.BuPu, "Numerical solution for $h(x,t)$")
        #plt.show()

    else:
        print("\n---------- Compilation of '{0}.f90' failed ----------\n".format(sys.argv[1]))

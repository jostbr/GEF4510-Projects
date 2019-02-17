
# Code by: Jostein Brandshoi

"""This program is written while using Windows, but should work in other OS'es as well, perhaps with a minor change in how
to run the executable resulting from compiling the fortran source. Another note: Here Python 2.x is used, however it should
run without problems in a Python 3.x environment as well.

This script that tells the operating system to compile the requested (by first command line argument) .f90 file with with the desired
command line arguments (see below). Then, if compilation successful, it extracts the result data from the specified file and plots the
results using matplotlib. Extraction of the data from the file is based on knowing the format (line-wise) outputed by the .f90-file.
Below is a table, with description, over the command line arguments for this script:

sys.argv[0] : Name of the program itself (visualize_results.py)
sys.argv[1] : Name of the .f90-file without the extension
sys.argv[2] : Name of method/algortihm used to solve the advection equation
sys.argv[3] : Type of FDA used for computing anti-diffusive velocity in MPDATA
sys.argv[4] : Name of the result output file, preferrably *.dat

Below is an overview on how to run the script with appropriate command line arguments.

> python visualize_results.py advection_solver <method_type> <diff_type> <output_filename>

where method_type can be set to 'mpdata' if one wishes to use this method. Otherwise the argument can be whatever, however non-empty,
and the regular upwind scheme without any flux correction will be used. <diff_type> may be either "forward", "centered" or "backward".
"""

import os, sys, time
import numpy as np
import matplotlib.pyplot as plt

if (__name__ == "__main__"):
    try:
        fortran_filename = sys.argv[1]
        method_type = sys.argv[2]
        diff_type = sys.argv[3]
        output_filename = sys.argv[4]
    except Exception:
        sys.exit("\nInvalid use of command line arguments. See docstring in source.")

    if (not os.system("gfortran -o {0}.exe {0}.f90".format(fortran_filename))):
        time_start = time.clock()
        os.system("{0}.exe {1} {2} {3}".format(fortran_filename, method_type, diff_type, output_filename))
        print("Execution time: {0} seconds".format(time.clock() - time_start))

        L = 50000.0; sigma = L / 10; u_0 = 1.0; C = 0.50    # Defining parameters
        delta_x = sigma / 10; delta_t = C * delta_x / u_0;  # Defining parameters
        j_max = int(round((L / delta_x) + 1))               # Defining parameters
        x = np.linspace(0, int(L / L), j_max)               # Dim-less space-array

        with open(output_filename, "r") as theta_values_file:
            for line_num, current_line in enumerate(theta_values_file):
                list_theta_j = current_line.split()                          # List of theta strings
                current_n = int(float(list_theta_j[0]))                      # Read current n-value
                current_cycle = int(round((current_n * u_0 * delta_t) / L))  # Compute current cycle
                theta_current_n = np.array([float(theta_j) for theta_j in list_theta_j[1:]])
                plt.plot(x, theta_current_n, linewidth = 1.5, label = "n = {0}, cycle: {1}".format(current_n, current_cycle))
                plt.hold("on")

        plt.xlabel("Dimension-less space coordinate: $x / L$", fontsize = 16)                   # Set label on the x-axis of the figure
        plt.ylabel("Dimension-less temperature: $\\theta / \\theta_0$", fontsize = 16)          # Set label on the y-axis of the figure
        plt.title("Numerical solution of advection equation, C = {}".format(C), fontsize = 17)  # Set title of the figure
        plt.legend(loc = 2)                                                                     # Show labels on all plotted curves
        plt.show()

    else:
        print("\n---------- Compilation of '{0}.f90' failed ----------\n".format(sys.argv[1]))

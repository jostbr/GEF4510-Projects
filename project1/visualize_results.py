
# Program by: Jostein Brandshoi

"""This program is written and tested while using Windows, but should work in other OS'es as well, perhaps with a minor change in
how to run the executable resulting from compiling the fortran source. Another note: Here Python 2.x is used, however it should
run without problems in a Python 3.x environment as well.

This script that tells the operating system to compile the requested (by first command line argument) .f90 file with the requested
name for the result output file (specified by the second command line argument). Then, if compilation successful, it extracts the
result data from the specified file and plots the results using matplotlib. Extraction of the data from the file is based on knowing
the format outputed by the .f90-file. Below is a table, with description, over the command line arguments for this script:

sys.argv[0] : Name of the program itself (visualize_results.py)
sys.argv[1] : Name of the .f90-file without the extension
sys.argv[2] : Name of the result output file, preferrably *.dat

This script can be run for both the atmosphere and ocean case, meaning one only needs to run the script with different command line
arguments corresponding to whether it's exercise 2d) or 2e). Below is an example on how to run the script for the two exercises.

> python visualize_results.py atmosphere_application atmosphere_results.dat
> python visualize_results.py ocean_application ocean_results.dat
"""

import os
import sys
import numpy as np
from matplotlib import pyplot

def visualize_analytical(n_value):
	analytical_solution = atmosphere_analytical

	H = 270.0; kappa = 30.0; j_max = 27; K = 0.45; psi_0 = 10.0 	# Defining parameters
	delta_z = H / (j_max - 1.0)										# Space step distance
	delta_t = (K * (delta_z ** 2)) / kappa							# Time step distance
	t = n_value * delta_t											# t-value corresponding to n

	z = np.linspace(0, H, j_max)										# Array of z-values
	psi = analytical_solution(z, t, psi_0, kappa, H)					# Array of psi-values
	pyplot.plot(psi / psi_0, z / H, label = "n = {}".format(n_value))	# Plot dim-less psi

def atmosphere_analytical(z, t, psi_0, kappa, H):
	return psi_0 * np.exp(-(((np.pi ** 2) * kappa) / (H ** 2)) * t) * np.sin((np.pi * z) / H)	# Analytical solution

def ocean_analytical(z, t, psi_0, kappa, D):
	return psi_0 * (z / D + 1)		# Steday-state solution

if (__name__ == "__main__"):
	# If the operating system succesfully compiles the Fortran file.
	if (not os.system("gfortran -o {0}.exe {0}.f90".format(sys.argv[1]))):
		os.system("{0}.exe {1}".format(sys.argv[1], sys.argv[2]))	# Execute the program to generate .dat-file with results

		if (sys.argv[1] == "atmosphere_application"):
			length_scale = "H"					# Length scale for labeling x-axis
			y_start = 0; y_stop = 1;			# Start and stop values for dim-less height

		elif (sys.argv[1] == "ocean_application"):
			length_scale = "(-D)"				# Length scale for labeling x-axis
			y_start = -1; y_stop = 0;			# Start and stop values for dim-less height

		n_values = []							# List to store n-values for analytical solution
		z = np.linspace(y_start, y_stop, 27)	# Create dimensionless height array

		# Loop over all lines in file where each line is a psi array for a given n-value.
		with open(sys.argv[2], "r") as data_file:
			for line_num, current_line in enumerate(data_file):
				list_of_current_values = current_line.split()									# Put each value on current_line in a list
				n_value = str(int(float(list_of_current_values[0])))							# Value of n for this line of psi-values
				n_values.append(int(n_value))													# Update list with current n_value
				psi_given_n = [float(psi_string) for psi_string in list_of_current_values[1:]]	# Convert psi-values from string to float
				pyplot.plot(np.array(psi_given_n), z, label = "n = {}".format(n_value))			# Plot current psi-array with label

				# Plot steady-state solution in the same figure as the numerical solution.
				if (sys.argv[1] == "ocean_application" and line_num == 7):
					pyplot.plot(ocean_analytical(np.linspace(-30, 0, 27), None, 10.0, None, 30) / 10.0, z, label = "SS-solution")

				pyplot.hold("on")

		pyplot.xlabel("Dimension-less temperature $\psi/\psi_0$")							# Set label on the x-axis of the figure
		pyplot.ylabel("Dimension-less height $z/{0}$".format(length_scale))					# Set label on the y-axis of the figure
		pyplot.title("Numerical solution for the {0}".format(sys.argv[1].split("_")[0]))	# Set title of the figure
		pyplot.legend()																		# Show labels on all plotted curves

		if (sys.argv[1] == "atmosphere_application"):
			pyplot.figure()																	# Generate new figure to plot in
			for current_n in n_values:
				visualize_analytical(current_n)		# Plot analytical solution for same n's
				pyplot.hold("on")					# Hold plot for more curves coming later

			pyplot.xlabel("Dimension-less temperature $\psi/\psi_0$")							# Set label on the x-axis of the figure
			pyplot.ylabel("Dimension-less height $z/{0}$".format(length_scale))					# Set label on the y-axis of the figure
			pyplot.title("Analytical solution for the {0}".format(sys.argv[1].split("_")[0]))	# Set title of the figure
			pyplot.legend()																		# Show labels on all plotted curves

		pyplot.show()

	else:
		print("\n---------- Compilation of '{0}' failed ----------\n".format(sys.argv[1]))	# Give error message if compilations fails

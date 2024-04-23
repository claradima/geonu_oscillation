import numpy as np
import matplotlib.pyplot as plt
import sys
import psutil
import os
print("imports successful")


# alternative calculation suggested by chatgpt
grid_count_crust = 640  # would need 637.1 for 20 km spacing
# 1275 for 10 km spacing
coordinates = np.linspace(-6371, 6371, grid_count_crust)
grid_1d_size_crust = coordinates[1] - coordinates[0]

# Generate the grid coordinates using meshgrid
x_coordinates, y_coordinates, z_coordinates = np.meshgrid(coordinates, coordinates, coordinates)

# Calculate the distance of each point from the origin using vectorized operations
distances_squared = x_coordinates ** 2 + y_coordinates ** 2 + z_coordinates ** 2

# Find indices where distance is less than or equal to the radius squared
crust_indices = np.logical_and(6369 ** 2 < distances_squared, distances_squared <= 6371 ** 2)

# Extract valid coordinates using boolean indexing
crust_grid = np.stack((x_coordinates[crust_indices], y_coordinates[crust_indices], z_coordinates[crust_indices]),
                      axis=-1)

print("crust grid created successfully")

grid_count_mantle = 250  # would need 637.1 for 20 km spacing
# 320 for 40 km spacing
# 300 for 42.5 km spacing
coordinates = np.linspace(-6371, 6371, grid_count_mantle)
grid_1d_size_mantle = coordinates[1] - coordinates[0]

# Generate the grid coordinates using meshgrid
x_coordinates, y_coordinates, z_coordinates = np.meshgrid(coordinates, coordinates, coordinates)

# Calculate the distance of each point from the origin using vectorized operations
distances_squared = x_coordinates ** 2 + y_coordinates ** 2 + z_coordinates ** 2
mantle_indices = np.logical_and(3470 ** 2 < distances_squared, distances_squared <= 6369 ** 2)
mantle_grid = np.stack((x_coordinates[mantle_indices], y_coordinates[mantle_indices], z_coordinates[mantle_indices]),
                       axis=-1)

print("mantle grid created successfully")
# apparently i had these in my code already

print("setting densities and abundances")

A_Th_c = 5 * (10 ** (-6))
A_Th_m = 8 * (10 ** (-8))
A_U_c = 1 * (10 ** (-6))
A_U_m = 2 * (10 ** (-8))

rho_c = 2.7
rho_m = 4.5  # g/cm^3 probably
print("densities and abundances set successfully")

# we have pretty much everything to compute things separately!
# compute for crust
# set position of SNO+

print("setting detector position")

SNO_r = np.array([0, 0, 6369])


# function to calculate relative distance to SNO+ from points in earth grid

def calc_relative_dist(points_array):
    # Calculate the Euclidean distance using vectorized operations
    relative_distances = np.linalg.norm(points_array - SNO_r, axis=1)
    return relative_distances


# define Delta function for an array of points and energies
# inputs: appropriate deltam _ij^2, energy array [MeV], points array
# relative distance calculated in km, hence Delta scaled up by a factor of 1000 to account for m-km conversion

def Delta_ij(energy_array, points_array, delta_m_ij_squared):
    # Calculate relative distances
    relative_distance_array = calc_relative_dist(points_array)

    # Reshape energy_array to perform element-wise division
    energy_array_reshaped = energy_array.reshape(-1, 1)

    # Calculate Delta using vectorized operations
    Delta = (1.27 * delta_m_ij_squared * relative_distance_array * 1000) / energy_array_reshaped

    return Delta


def P_ee_full_optimized(energy_array, points_array, theta_12, theta_13, delta_m_21_squared):
    # P_ee = np.empty((len(energy_array), len(points_array))) # np.empty more efficient
    # than np.zeros or whatever
    # might not need this at all though
    relative_distance_array = calc_relative_dist(points_array)

    Delta_31 = Delta_ij(energy_array, points_array, delta_m_31_squared)
    print("Delta_31 computed successfully")
    Delta_12 = Delta_ij(energy_array, points_array, delta_m_21_squared)
    print("Delta_12 computed successfully")

    A = (np.cos(theta_13)) ** 4 * (np.sin(2 * theta_12)) ** 2
    B = np.sin(2 * theta_13) ** 2

    sin_squared_Delta_31 = np.sin(Delta_31) ** 2
    sin_squared_Delta_12 = np.sin(Delta_12) ** 2

    print("values computed successfully")

    P_ee = 1 - (A * sin_squared_Delta_12 + B * sin_squared_Delta_31)
    print("P_ee computed successfully")

    return P_ee


# set standard oscillation parameters
print("setting standard oscillation parameters")

theta_12 = 0.5903  # rad
theta_23 = 0.8430  # rad
theta_13 = 0.1503  # rad

delta_m_21_squared = 7.39 * 10 ** (-5)  # eV^2
delta_m_32_squared = 2.449 * 10 ** (-3)  # eV^2
delta_m_31_squared = delta_m_32_squared
print("defining energy array")
energy_array = np.linspace(1.8, 3.3, 100)


## optimized definition of integral, optimized P_ee
def integral_over_positions_Th_U(points_array, energy_array, grid_1d_size, theta_12, delta_m_21_squared, A_Th, A_U,
                                 rho):
    dV = grid_1d_size ** 3

    relative_distance_array = calc_relative_dist(points_array)
    print("Relative distance array computed successfully")
    P_ee_array = P_ee_full_optimized(energy_array, points_array, theta_12, theta_13, delta_m_21_squared)
    print("P_ee_array computed successfully")

    # Compute sum_Th
    sum_Th = np.sum(P_ee_array * ((A_Th * rho) / (4 * np.pi * (relative_distance_array ** 2)))[np.newaxis, :] * dV,
                    axis=1)

    print("sum_Th computed successfully")

    # Compute sum_U
    sum_U = np.sum(P_ee_array * ((A_U * rho) / (4 * np.pi * (relative_distance_array ** 2)))[np.newaxis, :] * dV,
                   axis=1)
    print("sum_U computed successfully")

    return sum_Th, sum_U


print("computing integral values for standard oscillation paramaters (crust contribution only)")
Th_integral_values_c, U_integral_values_c = integral_over_positions_Th_U(crust_grid, energy_array, grid_1d_size_crust,
                                                                         theta_12, delta_m_21_squared, A_Th_c, A_U_c,
                                                                         rho_c)
print("done! moving on")


# #make plots
# print("making plots")
# #Thorium
# plt.plot(energy_array, Th_integral_values_c, marker='o', linestyle='-', color='b')
# plt.xlabel('Energy [MeV]')
# plt.ylabel('Integral value')
# plt.title('Integral over emission position, Thorium')
# #plt.savefig('Integral over emission position, Thorium, standard oscillation parameters.pdf', format='pdf')
# plt.show()
#
# #Uranium
# plt.plot(energy_array, U_integral_values_c, marker='o', linestyle='-', color='r')
# plt.xlabel('Energy [MeV]')
# plt.ylabel('Integral value')
# plt.title('Integral over emission position, Uranium')
# #plt.savefig('Integral over emission position, Uranium, standard oscillation parameters.pdf', format='pdf')
# plt.show()
#
# #Both together
#
# plt.plot(energy_array, Th_integral_values_c, marker='o', linestyle='-', color='b', label='Thorium')
# plt.plot(energy_array, U_integral_values_c, marker='o', linestyle='-', color='r', label='Uranium')
# plt.xlabel('Energy [MeV]')
# plt.ylabel('Integral value')
# plt.title('Integral over emission position')
# plt.legend()
# #plt.savefig('Integral over emission position, standard oscillation parameters.pdf', format='pdf')
# plt.show()
def integral_over_positions_Th_U_constant_P_ee(points_array, energy_array, grid_1d_size, A_Th, A_U, rho):
    dV = grid_1d_size ** 3

    relative_distance_array = calc_relative_dist(points_array)
    P_ee = 0.553
    P_ee_array = np.full((len(energy_array), len(points_array)), P_ee)
    # Compute sum_Th
    sum_Th = np.sum(P_ee_array * ((A_Th * rho) / (4 * np.pi * (relative_distance_array ** 2)))[np.newaxis, :] * dV,
                    axis=1)

    # Compute sum_U
    sum_U = np.sum(P_ee_array * ((A_U * rho) / (4 * np.pi * (relative_distance_array ** 2)))[np.newaxis, :] * dV,
                   axis=1)

    return sum_Th, sum_U


print("computing integral values for standard P_ee (crust contribution only")
Th_integral_values_constant_P_ee_c, U_integral_values_constant_P_ee_c = integral_over_positions_Th_U_constant_P_ee(
    crust_grid, energy_array, grid_1d_size_crust, A_Th_c, A_U_c, rho_c)
print("done! moving on")


# #make plots, just for integrals with constant P_ee
# print("making plots")
# #Thorium
# plt.plot(energy_array, Th_integral_values_constant_P_ee_c, marker='o', linestyle='-', color='b')
# plt.xlabel('Energy [MeV]')
# plt.ylabel('Integral value')
# plt.title('Integral over emission position, Thorium')
# plt.show()
#
# #Uranium
# plt.plot(energy_array, U_integral_values_constant_P_ee_c, marker='o', linestyle='-', color='r')
# plt.xlabel('Energy [MeV]')
# plt.ylabel('Integral value')
# plt.title('Integral over emission position, Uranium')
# plt.show()
#
# #Both together
#
# plt.plot(energy_array, Th_integral_values_constant_P_ee_c, marker='o', linestyle='-', color='b', label='Thorium')
# plt.plot(energy_array, U_integral_values_constant_P_ee_c, marker='o', linestyle='-', color='r', label='Uranium')
# plt.xlabel('Energy [MeV]')
# plt.ylabel('Integral value')
# plt.title('Integral over emission position')
# plt.legend()
# plt.show()
# #make plots to compare integral for constant P_ee and non-constant(dependent on energy and oscillation parameters)
#
# #Thorium
# plt.plot(energy_array, Th_integral_values_constant_P_ee_c, marker='o', linestyle='-', color='b', label='constant P_ee')
# plt.plot(energy_array, Th_integral_values_c, marker='*', linestyle='-', color='b', label='non-constant P_ee, standard osc params')
# plt.xlabel('Energy [MeV]')
# plt.ylabel('Integral value')
# plt.title('Integral over emission position, Thorium')
# plt.legend()
# #plt.savefig('Thorium Integral, standard osc params vs constant P_ee.pdf', format='pdf')
# plt.show()
#
# #Uranium
# plt.plot(energy_array, U_integral_values_constant_P_ee_c, marker='o', linestyle='-', color='r')
# plt.plot(energy_array, U_integral_values_c, marker='*', linestyle='-', color='r', label = 'non-constant P_ee, standard osc params')
# plt.xlabel('Energy [MeV]')
# plt.ylabel('Integral value')
# plt.title('Integral over emission position, Uranium')
# plt.legend()
# #plt.savefig('Uranium Integral, standard osc params vs constant P_ee.pdf', format='pdf')
# plt.show()

def sigma_IBD(energy_array):
    # optimized with numpy

    m_e = 0.511  # MeV
    m_p = 938  # MeV
    m_n = 941  # MeV
    E_th = 1.8  # MeV

    sigma = ((energy_array - E_th + m_e) ** 2) * ((1 - (m_e) ** 2 / ((energy_array - E_th + m_e) ** 2)) ** (1 / 2))

    return sigma


# all these + sources in Earth model more bins -3 (for example)

print("setting lambda and mu values")

lambda_U = 4.916
lambda_Th = 1.563
mu_U = 235
mu_Th = 232
# get fluxes
print("getting emission fluxes")
energy_array_U = []
dn_dE_U = []

with open('U238_spectrum.txt', 'r') as file:
    for line in file:
        # Split each line into columns
        columns = line.split()

        # Convert the elements to float and append to arrays
        energy_array_U.append(float(columns[0]))
        dn_dE_U.append(float(columns[1]))

# Scale down all energies by a factor of 1000
energy_array_U = np.array(energy_array_U) / 1000
dn_dE_U = np.array(dn_dE_U)

# plt.plot(energy_array_U, dn_dE_U)
# plt.xlabel('E_nu [MeV]')
# plt.yscale('log')
# plt.ylabel('Intensity (some units, doesn\'t matter)')
# plt.title('Energy spectrum of geonus from U238 decay')
#
# # Add shaded region between 1.8 MeV and 3.2 MeV
# plt.axvspan(1.8, 3.2, alpha=0.3, color='gray')
# # Enable minor ticks on x-axis
# plt.minorticks_on()
# plt.show()

print("done for Uranium, moving on to Thorium")

energy_array_Th = []
dn_dE_Th = []
with open('Th232_spectrum.txt', 'r') as file:
    for line in file:
        # Split each line into columns
        columns = line.split()

        # Convert the elements to float and append to arrays
        energy_array_Th.append(float(columns[0]))
        dn_dE_Th.append(float(columns[1]))

# Scale down all energies by a factor of 1000
energy_array_Th = np.array(energy_array_Th) / 1000
dn_dE_Th = np.array(dn_dE_Th)


# plt.plot(energy_array_Th, dn_dE_Th)
# plt.xlabel('E_nu [MeV]')
# plt.yscale('log')
# plt.ylabel('Intensity (some units, doesn\'t matter)')
# plt.title('Energy spectrum of geonus from Th232 decay')
# # Add shaded region between 1.8 MeV and 3.2 MeV
# plt.axvspan(1.8, 3.2, alpha=0.3, color='gray')
# # Enable minor ticks on x-axis
# plt.minorticks_on()
# plt.show()
# Plot U238 decay data (blue line)
# plt.plot(energy_array_U, dn_dE_U, label='U238 decays', color='blue')
#
# # Plot Th232 decay data (red line)
# plt.plot(energy_array_Th, dn_dE_Th, label='Th232 decays', color='red')
#
# plt.xlabel('E_nu [MeV]')
# plt.yscale('log')
# plt.ylabel('Intensity (some units, doesn\'t matter)')
# plt.title('Energy spectrum of geonus')
#
# # Add shaded region between 1.8 MeV and 3.3 MeV
# plt.axvspan(1.8, 3.3, alpha=0.3, color='gray')
#
# # Enable minor ticks on x-axis
# plt.minorticks_on()
#
# plt.legend(loc='upper right')
#
# plt.show()
def rebin_counts(initial_bins, counts_in_initial_bins, final_bin_midpoints):
    """
    Rebins counts data.

    Parameters:
    - initial_bins: array, bin edges of the initial data
    - counts_in_initial_bins: array, counts in each initial bin
    - final_bin_midpoints: array, midpoints of the final desired bins

    Returns:
    - counts_in_final_bins: array, counts in each final bin
    """

    # Calculate bin midpoints of the initial bins
    bin_midpoints = (initial_bins[:-1] + initial_bins[1:]) / 2

    # Use np.histogram to calculate counts in final bins
    counts_in_final_bins, _ = np.histogram(initial_bins, bins=np.concatenate(
        [initial_bins, [2 * initial_bins[-1] - initial_bins[-2]]]), weights=counts_in_initial_bins)

    # Interpolate the counts to the final bin midpoints
    counts_in_final_bins = np.interp(final_bin_midpoints, bin_midpoints, counts_in_final_bins[:-1])

    return counts_in_final_bins


print("rebin to match energy array")

dn_dE_rebinned_U = rebin_counts(energy_array_U, dn_dE_U, energy_array)
dn_dE_rebinned_Th = rebin_counts(energy_array_Th, dn_dE_Th, energy_array)
# compute total fluxes, coming from crust only

print("computing IBD cross section")
sigma = sigma_IBD(energy_array)
print("computing fluxes at detector (standard oscillation parameters, crust contribution only")
N_Th_c = ((lambda_Th) / (mu_Th)) * sigma * dn_dE_rebinned_Th * Th_integral_values_c
N_U_c = ((lambda_U) / (mu_U)) * sigma * dn_dE_rebinned_U * U_integral_values_c

# # Plot histogram for both
# plt.step(energy_array, N_U_c + N_Th_c, where='mid', label='Total', color='green')
# plt.step(energy_array, N_U_c, where='mid', label='U238 decays', color='blue')
# plt.step(energy_array, N_Th_c, where='mid', label='Th232 decays', color='red')
#
# plt.xlabel('E_nu [MeV]')
# plt.yscale('log')
# plt.ylabel('Expected counts (some units, doesn\'t matter)')
# plt.title('Expected geonus')
#
# # Set lower limit for y-axis to 10^-7
# plt.ylim(bottom=1e-12)
#
# # Enable minor ticks on x-axis
# plt.minorticks_on()
#
# # Add legend
# plt.legend(loc='upper right')
#
# plt.show()

print("computing fluxes at detector (constant P_ee, crust contribution only")
N_Th_constant_P_ee_c = ((lambda_Th) / (mu_Th)) * sigma * dn_dE_rebinned_Th * Th_integral_values_constant_P_ee_c
N_U_constant_P_ee_c = ((lambda_U) / (mu_U)) * sigma * dn_dE_rebinned_U * U_integral_values_constant_P_ee_c

# # Plot histogram for both, constant P_ee only
# plt.step(energy_array, N_U_constant_P_ee_c + N_Th_constant_P_ee_c, where='mid', label='Total, constant P_ee', color='green')
# plt.step(energy_array, N_U_constant_P_ee_c, where='mid', label='U238 decays', color='blue')
# plt.step(energy_array, N_Th_constant_P_ee_c, where='mid', label='Th232 decays', color='red')
#
# plt.xlabel('E_nu [MeV]')
# plt.yscale('log')
# plt.ylabel('Expected counts (some units, doesn\'t matter)')
# plt.title('Expected geonus')
#
# # Set lower limit for y-axis to 10^-7
# plt.ylim(bottom=1e-12)
#
# # Enable minor ticks on x-axis
# plt.minorticks_on()
#
# # Add legend
# plt.legend(loc='upper right')
#
# plt.show()
# plot to compare total for constant P_ee vs non-constant P_ee with standard osc parameters
#
# plt.step(energy_array, N_U_constant_P_ee_c + N_Th_constant_P_ee_c, where='mid', label='Total, constant P_ee', color='green')
# plt.step(energy_array, N_U_c + N_Th_c, where='mid', label='Total, non-constant P_ee', color='lime')
#
# plt.xlabel('E_nu [MeV]')
# plt.yscale('log')
# plt.ylabel('Expected counts (some units, doesn\'t matter)')
# plt.title('Expected geonus')
#
# # Set lower limit for y-axis to 10^-7
# plt.ylim(bottom=1e-12)
#
# # Enable minor ticks on x-axis
# plt.minorticks_on()
#
# # Add legend
# plt.legend(loc='upper right')
# #plt.savefig('Expected counts, standard params vs constant P_ee.pdf', format='pdf')
#
# plt.show()
#
#
# plt.plot(energy_array, (N_U_constant_P_ee_c + N_Th_constant_P_ee_c) / (N_U_c + N_Th_c), label='ratio', color='green')
# plt.xlabel('E_nu [MeV]')
# plt.ylabel('ratio')
# plt.title('Ratio of expected geonus, const P_ee / non-const')
#
# # Set lower limit for y-axis to 10^-7
# #plt.ylim(bottom=1e-12)
#
# # Enable minor ticks on x-axis
# plt.minorticks_on()
#
# # Add legend
# plt.legend(loc='upper right')
# #plt.savefig('Ratio of expected geonus, standard params vs constant P_ee.pdf', format='pdf')
#
# plt.show()
# we did this for the crust, now want to do this for the mantle as well
# mantle grid too big, restrict distances
# adjust size as required

print("Now getting into the complicated stuff")
print("Starting mantle contribution calculations ... ")
print("Separating mantle points")

mantle_grid_distances = calc_relative_dist(mantle_grid)
max_mantle_distance = 8000
mask1 = mantle_grid_distances <= max_mantle_distance
new_mantle_grid1 = np.array(mantle_grid)[mask1]

mask2 = mantle_grid_distances > max_mantle_distance
new_mantle_grid2 = np.array(mantle_grid)[mask2]

print("len(mantle_grid) : " + str(len(mantle_grid)))
print("len(new_mantle_grid1) : " + str(len(new_mantle_grid1)))
# Print the size in memory in gigabytes (GB)
# for max distance = 2R_Earth, we should have the same number of points
size_gb = new_mantle_grid1.nbytes / (1024 ** 3)
print("Size of new_mantle_grid1 in memory:", size_gb, "GB")

print("ratio is " + str(len(mantle_grid) / (len(new_mantle_grid1))))

print("computing integrals for standard oscillation parameters, mantle contribution only")

Th_integral_values_m1, U_integral_values_m1 = integral_over_positions_Th_U(new_mantle_grid1, energy_array,
                                                                           grid_1d_size_mantle, theta_12,
                                                                           delta_m_21_squared, A_Th_m, A_U_m, rho_m)

print("firs component done")
Th_integral_values_m2, U_integral_values_m2 = integral_over_positions_Th_U(new_mantle_grid2, energy_array,
                                                                           grid_1d_size_mantle, theta_12,
                                                                           delta_m_21_squared, A_Th_m, A_U_m, rho_m)
print("second component done")
Th_integral_values_m = Th_integral_values_m1 + Th_integral_values_m2
U_integral_values_m = U_integral_values_m1 + U_integral_values_m2
print("total contribution done")

print("computing integrals for constant P_ee, mantle contribution only")
Th_integral_values_constant_P_ee_m1, U_integral_values_constant_P_ee_m1 = integral_over_positions_Th_U_constant_P_ee(
    new_mantle_grid1, energy_array, grid_1d_size_mantle, A_Th_m, A_U_m, rho_m)

print("first component done")
Th_integral_values_constant_P_ee_m2, U_integral_values_constant_P_ee_m2 = integral_over_positions_Th_U_constant_P_ee(
    new_mantle_grid2, energy_array, grid_1d_size_mantle, A_Th_m, A_U_m, rho_m)

print("second component done")
Th_integral_values_constant_P_ee_m = Th_integral_values_constant_P_ee_m1 + Th_integral_values_constant_P_ee_m2
U_integral_values_constant_P_ee_m = U_integral_values_constant_P_ee_m1 + U_integral_values_constant_P_ee_m2

print("total contribution done")
# #make plots to compare integral for constant P_ee and non-constant(dependent on energy and oscillation parameters)
#
# #Thorium
# plt.plot(energy_array, Th_integral_values_constant_P_ee_m, marker='o', linestyle='-', color='b', label='constant P_ee')
# plt.plot(energy_array, Th_integral_values_m, marker='*', linestyle='-', color='b', label='non-constant P_ee, standard osc params')
# plt.xlabel('Energy [MeV]')
# plt.ylabel('Integral value')
# plt.title('Integral over emission position, Thorium')
# plt.legend()
# #plt.savefig('Thorium Integral, standard osc params vs constant P_ee.pdf', format='pdf')
# plt.show()
#
# #Uranium
# plt.plot(energy_array, U_integral_values_constant_P_ee_m, marker='o', linestyle='-', color='r')
# plt.plot(energy_array, U_integral_values_m, marker='*', linestyle='-', color='r', label = 'non-constant P_ee, standard osc params')
# plt.xlabel('Energy [MeV]')
# plt.ylabel('Integral value')
# plt.title('Integral over emission position, Uranium')
# plt.legend()
# #plt.savefig('Uranium Integral, standard osc params vs constant P_ee.pdf', format='pdf')
# plt.show()

print("computing fluxes at detector location")

N_Th_m = ((lambda_Th) / (mu_Th)) * sigma * dn_dE_rebinned_Th * Th_integral_values_m
N_U_m = ((lambda_U) / (mu_U)) * sigma * dn_dE_rebinned_U * U_integral_values_m

N_Th_constant_P_ee_m = ((lambda_Th) / (mu_Th)) * sigma * dn_dE_rebinned_Th * Th_integral_values_constant_P_ee_m
N_U_constant_P_ee_m = ((lambda_U) / (mu_U)) * sigma * dn_dE_rebinned_U * U_integral_values_constant_P_ee_m

# plot to compare total for constant P_ee vs non-constant P_ee with standard osc parameters

print("plotting comparison between constant P_ee and standard oscillation parameters")

plt.figure()

plt.step(energy_array, N_U_constant_P_ee_c + N_Th_constant_P_ee_c, where='mid',
         label='Total, constant P_ee, crust only', color='pink')
plt.step(energy_array, N_U_c + N_Th_c, where='mid', label='Total, non-constant P_ee, crust only', color='red')
plt.step(energy_array, N_U_constant_P_ee_c + N_Th_constant_P_ee_c + N_U_constant_P_ee_m + N_Th_constant_P_ee_m,
         where='mid', label='Total, constant P_ee, c+m', color='green')
plt.step(energy_array, N_U_c + N_Th_c + N_U_m + N_Th_m, where='mid', label='Total, non-constant P_ee, c+m',
         color='lime')

plt.xlabel('E_nu [MeV]')
plt.yscale('log')
plt.ylabel('Expected counts (some units, doesn\'t matter)')
plt.title('Expected geonus')

# Set lower limit for y-axis to 10^-7
plt.ylim(bottom=1e-12)

# Enable minor ticks on x-axis
plt.minorticks_on()

# Add legend
plt.legend(loc='upper right')
plt.savefig('Expected counts, standard params vs constant P_ee, max mantle distance ' + str(
    max_mantle_distance) + " Different Grids.pdf", format='pdf')

# plt.show()

plt.figure()
plt.plot(energy_array, (N_U_constant_P_ee_c + N_Th_constant_P_ee_c) / (N_U_c + N_Th_c), label='ratio, crust only',
         color='red')
plt.plot(energy_array, (N_U_constant_P_ee_c + N_Th_constant_P_ee_c + N_U_constant_P_ee_m + N_Th_constant_P_ee_m) / (
            N_U_c + N_Th_c + N_U_m + N_Th_m), label='ratio, crust + mantle', color='green')
plt.xlabel('E_nu [MeV]')
plt.ylabel('ratio')
plt.title('Ratio of expected geonus, const P_ee / non-const')

# Set lower limit for y-axis to 10^-7
# plt.ylim(bottom=1e-12)

# Enable minor ticks on x-axis
plt.minorticks_on()

# Add legend
plt.legend(loc='upper right')
plt.savefig('Ratio of expected geonus, standard params vs constant P_ee, max mantle distance ' + str(
    max_mantle_distance) + " Different Grids.pdf", format='pdf')

# plt.show()
plt.figure()
plt.plot(energy_array, (N_U_constant_P_ee_c + N_Th_constant_P_ee_c + N_U_constant_P_ee_m + N_Th_constant_P_ee_m) / (
            N_U_c + N_Th_c + N_U_m + N_Th_m), label='ratio, crust + mantle', color='green')
plt.xlabel('E_nu [MeV]')
plt.ylabel('ratio')
plt.title('Ratio of expected geonus, const P_ee / non-const')

# Set lower limit for y-axis to 10^-7
# plt.ylim(bottom=1e-12)

# Enable minor ticks on x-axis
plt.minorticks_on()

# Add legend
plt.legend(loc='upper right')
plt.savefig('Ratio of expected geonus, standard params vs constant P_ee, Different Grids cm.pdf', format='pdf')

# plt.show()

print("plots saved")

print("setting alternative oscillation parameters")

theta_12_low = theta_12  - 0.0132  # rad; do +- 1sigma now
print("theta_12_low: " + str(theta_12_low))
theta_12_high = theta_12  + 0.0136
print("theta_12_high: " + str(theta_12_high))

delta_m_21_squared_low = 5 * 10 ** (-5)  # eV^2
print("delta_m_21_squared_low: " + str(delta_m_21_squared_low))
delta_m_21_squared_high = 1.2 * 10 ** (-4)
print("delta_m_21_squared_high: " + str(delta_m_21_squared_high))

print("computing contribution from crust for low theta")
Th_integral_values_theta_low_c, U_integral_values_theta_low_c = integral_over_positions_Th_U(crust_grid, energy_array,
                                                                                             grid_1d_size_crust,
                                                                                             theta_12_low,
                                                                                             delta_m_21_squared, A_Th_c,
                                                                                             A_U_c, rho_c)
print("done!! moving on ...")

print("computing contribution from mantle for low theta")
Th_integral_values_theta_low_m1, U_integral_values_theta_low_m1 = integral_over_positions_Th_U(new_mantle_grid1,
                                                                                               energy_array,
                                                                                               grid_1d_size_mantle,
                                                                                               theta_12_low,
                                                                                               delta_m_21_squared,
                                                                                               A_Th_m, A_U_m, rho_m)
print("first component done")
Th_integral_values_theta_low_m2, U_integral_values_theta_low_m2 = integral_over_positions_Th_U(new_mantle_grid2,
                                                                                               energy_array,
                                                                                               grid_1d_size_mantle,
                                                                                               theta_12_low,
                                                                                               delta_m_21_squared,
                                                                                               A_Th_m, A_U_m, rho_m)

print("second component done")
Th_integral_values_theta_low_m = Th_integral_values_theta_low_m1 + Th_integral_values_theta_low_m2
U_integral_values_theta_low_m = U_integral_values_theta_low_m1 + U_integral_values_theta_low_m2
print("total contribution from mantle done!! moving on ...")

print("computing contribution from crust for high theta")
Th_integral_values_theta_high_c, U_integral_values_theta_high_c = integral_over_positions_Th_U(crust_grid, energy_array,
                                                                                               grid_1d_size_crust,
                                                                                               theta_12_high,
                                                                                               delta_m_21_squared,
                                                                                               A_Th_c, A_U_c, rho_c)
print("done!! moving on ...")

print("computing contribution from mantle for high theta")
Th_integral_values_theta_high_m1, U_integral_values_theta_high_m1 = integral_over_positions_Th_U(new_mantle_grid1,
                                                                                                 energy_array,
                                                                                                 grid_1d_size_mantle,
                                                                                                 theta_12_high,
                                                                                                 delta_m_21_squared,
                                                                                                 A_Th_m, A_U_m, rho_m)
print("first component done")
Th_integral_values_theta_high_m2, U_integral_values_theta_high_m2 = integral_over_positions_Th_U(new_mantle_grid2,
                                                                                                 energy_array,
                                                                                                 grid_1d_size_mantle,
                                                                                                 theta_12_high,
                                                                                                 delta_m_21_squared,
                                                                                                 A_Th_m, A_U_m, rho_m)
print("second component done")
Th_integral_values_theta_high_m = Th_integral_values_theta_high_m1 + Th_integral_values_theta_high_m2
U_integral_values_theta_high_m = U_integral_values_theta_high_m1 + U_integral_values_theta_high_m2
print("total contribution from mantle done!! moving on ...")

print("computing fluxes at detector location for low theta")

N_Th_theta_low_c = ((lambda_Th) / (mu_Th)) * sigma * dn_dE_rebinned_Th * Th_integral_values_theta_low_c
N_U_theta_low_c = ((lambda_U) / (mu_U)) * sigma * dn_dE_rebinned_U * U_integral_values_theta_low_c

N_Th_theta_low_m = ((lambda_Th) / (mu_Th)) * sigma * dn_dE_rebinned_Th * Th_integral_values_theta_low_m
N_U_theta_low_m = ((lambda_U) / (mu_U)) * sigma * dn_dE_rebinned_U * U_integral_values_theta_low_m

print("done! moving on ...")
print("computing fluxes at detector location for high theta")
N_Th_theta_high_c = ((lambda_Th) / (mu_Th)) * sigma * dn_dE_rebinned_Th * Th_integral_values_theta_high_c
N_U_theta_high_c = ((lambda_U) / (mu_U)) * sigma * dn_dE_rebinned_U * U_integral_values_theta_high_c

N_Th_theta_high_m = ((lambda_Th) / (mu_Th)) * sigma * dn_dE_rebinned_Th * Th_integral_values_theta_high_m
N_U_theta_high_m = ((lambda_U) / (mu_U)) * sigma * dn_dE_rebinned_U * U_integral_values_theta_high_m

print("done! moving on ...")
# plot comparisons
print("plot comparisons for different theta")
# plot to compare total for constant P_ee vs non-constant P_ee with standard osc parameters
plt.figure()
plt.step(energy_array, N_U_theta_low_c + N_Th_theta_low_c + N_U_theta_low_m + N_Th_theta_low_m, where='mid',
         label=('theta = ' + str(theta_12_low) + '[rad]'), color='green')
plt.step(energy_array, N_U_c + N_Th_c + N_U_m + N_Th_m, where='mid', label='theta = 0.5903[rad]', color='lime')
plt.step(energy_array, N_U_theta_high_c + N_Th_theta_high_c + N_U_theta_high_m + N_Th_theta_high_m, where='mid',
         label=('theta = ' + str(theta_12_high) + '[rad]'), color='olive')

plt.xlabel('E_nu [MeV]')
plt.yscale('log')
plt.ylabel('Expected counts (some units, doesn\'t matter)')
plt.title('Expected geonus, delta_m_21^2 = 7.39e-5, different theta_12 values, different grids')

# Set lower limit for y-axis to 10^-7
plt.ylim(bottom=1e-11)

# Enable minor ticks on x-axis
plt.minorticks_on()

# Add legend
plt.legend(loc='upper right')
plt.savefig(
    'Expected geonu, standard delta m, max mantle distance ' + str(max_mantle_distance) + " Different Grids.pdf",
    format='pdf')

# plt.show()

plt.figure()
plt.plot(energy_array,
         (N_U_theta_low_c + N_Th_theta_low_c + N_U_theta_low_m + N_Th_theta_low_m) / (N_U_c + N_Th_c + N_U_m + N_Th_m),
         label='ratio', color='green')
plt.xlabel('E_nu [MeV]')
plt.ylabel('ratio')
plt.title('Ratio of expected geonus, low theta / standard theta')

# Set lower limit for y-axis to 10^-7
# plt.ylim(bottom=1e-7)

# Enable minor ticks on x-axis
plt.minorticks_on()

# Add legend
plt.legend(loc='upper right')
plt.savefig('Ratio of expected geonu, standard delta m, low theta, max mantle distance ' + str(
    max_mantle_distance) + " Different Grids.pdf", format='pdf')

# plt.show()
plt.figure()
plt.plot(energy_array, (N_U_theta_high_c + N_Th_theta_high_c + N_U_theta_high_m + N_Th_theta_high_m) / (
            N_U_c + N_Th_c + N_U_m + N_Th_m), label='ratio', color='olive')
plt.xlabel('E_nu [MeV]')
plt.ylabel('ratio')
plt.title('Ratio of expected geonus, high theta / standard theta')

# Set lower limit for y-axis to 10^-7
# plt.ylim(bottom=1e-7)

# Enable minor ticks on x-axis
plt.minorticks_on()

# Add legend
plt.legend(loc='upper right')
plt.savefig('Ratio of expected geonu, standard delta m, high theta, max mantle distance ' + str(
    max_mantle_distance) + " Different Grids.pdf", format='pdf')

# plt.show()

print("plots saved! yay!!")

print("moving on to different delta_m_21^2 for standard theta ... ")

print("computing contribution for low delta_m^2 (crust only)")
Th_integral_values_delta_m_low_c, U_integral_values_delta_m_low_c = integral_over_positions_Th_U(crust_grid,
                                                                                                 energy_array,
                                                                                                 grid_1d_size_crust,
                                                                                                 theta_12,
                                                                                                 delta_m_21_squared_low,
                                                                                                 A_Th_c, A_U_c, rho_c)

print("done!! moving on ... ")

print("computing contribution for high delta_m^2 (mantle only)")
Th_integral_values_delta_m_low_m1, U_integral_values_delta_m_low_m1 = integral_over_positions_Th_U(new_mantle_grid1,
                                                                                                   energy_array,
                                                                                                   grid_1d_size_mantle,
                                                                                                   theta_12,
                                                                                                   delta_m_21_squared_low,
                                                                                                   A_Th_m, A_U_m, rho_m)
print("first component done")
Th_integral_values_delta_m_low_m2, U_integral_values_delta_m_low_m2 = integral_over_positions_Th_U(new_mantle_grid2,
                                                                                                   energy_array,
                                                                                                   grid_1d_size_mantle,
                                                                                                   theta_12,
                                                                                                   delta_m_21_squared_low,
                                                                                                   A_Th_m, A_U_m, rho_m)
print("second component done")
Th_integral_values_delta_m_low_m = Th_integral_values_delta_m_low_m1 + Th_integral_values_delta_m_low_m2
U_integral_values_delta_m_low_m = U_integral_values_delta_m_low_m1 + U_integral_values_delta_m_low_m2
print("total contribution done")

print("computing contribution for high delta_m^2 (crust only)")
Th_integral_values_delta_m_high_c, U_integral_values_delta_m_high_c = integral_over_positions_Th_U(crust_grid,
                                                                                                   energy_array,
                                                                                                   grid_1d_size_crust,
                                                                                                   theta_12,
                                                                                                   delta_m_21_squared_high,
                                                                                                   A_Th_c, A_U_c, rho_c)
print("done!! moving on ... ")

print("computing contribution for high delta_m^2 (mantle only)")
Th_integral_values_delta_m_high_m1, U_integral_values_delta_m_high_m1 = integral_over_positions_Th_U(new_mantle_grid1,
                                                                                                     energy_array,
                                                                                                     grid_1d_size_mantle,
                                                                                                     theta_12,
                                                                                                     delta_m_21_squared_high,
                                                                                                     A_Th_m, A_U_m,
                                                                                                     rho_m)
print("first component done")
Th_integral_values_delta_m_high_m2, U_integral_values_delta_m_high_m2 = integral_over_positions_Th_U(new_mantle_grid2,
                                                                                                     energy_array,
                                                                                                     grid_1d_size_mantle,
                                                                                                     theta_12,
                                                                                                     delta_m_21_squared_high,
                                                                                                     A_Th_m, A_U_m,
                                                                                                     rho_m)
print("second component done")
Th_integral_values_delta_m_high_m = Th_integral_values_delta_m_high_m1 + Th_integral_values_delta_m_high_m2
U_integral_values_delta_m_high_m = U_integral_values_delta_m_high_m1 + U_integral_values_delta_m_high_m2
print("total contribution done")

print("computing fluxes at detector location, low delta m^2")
N_Th_delta_m_low_c = ((lambda_Th) / (mu_Th)) * sigma * dn_dE_rebinned_Th * Th_integral_values_delta_m_low_c
N_U_delta_m_low_c = ((lambda_U) / (mu_U)) * sigma * dn_dE_rebinned_U * U_integral_values_delta_m_low_c

N_Th_delta_m_low_m = ((lambda_Th) / (mu_Th)) * sigma * dn_dE_rebinned_Th * Th_integral_values_delta_m_low_m
N_U_delta_m_low_m = ((lambda_U) / (mu_U)) * sigma * dn_dE_rebinned_U * U_integral_values_delta_m_low_m

print("computing fluxes at detector location, high delta m^2")
N_Th_delta_m_high_c = ((lambda_Th) / (mu_Th)) * sigma * dn_dE_rebinned_Th * Th_integral_values_delta_m_high_c
N_U_delta_m_high_c = ((lambda_U) / (mu_U)) * sigma * dn_dE_rebinned_U * U_integral_values_delta_m_high_c

N_Th_delta_m_high_m = ((lambda_Th) / (mu_Th)) * sigma * dn_dE_rebinned_Th * Th_integral_values_delta_m_high_m
N_U_delta_m_high_m = ((lambda_U) / (mu_U)) * sigma * dn_dE_rebinned_U * U_integral_values_delta_m_high_m

print("almost done! now making comparison plots ... ")
# plot comparisons

# plot to compare total for constant P_ee vs non-constant P_ee with standard osc parameters
plt.figure()
plt.step(energy_array, N_U_delta_m_low_c + N_Th_delta_m_low_c + N_U_delta_m_low_m + N_Th_delta_m_low_m, where='mid',
         label=('delta_m^2 = ' + str(delta_m_21_squared_low) + '[rad]'), color='green')
plt.step(energy_array, N_U_c + N_Th_c + N_U_m + N_Th_m, where='mid', label='delta_m_21^2 = 7.39e-5 eV^2', color='lime')
plt.step(energy_array, N_U_delta_m_high_c + N_Th_delta_m_high_c + N_U_delta_m_high_m + N_Th_delta_m_high_m, where='mid',
         label=('delta_m^2 = ' + str(delta_m_21_squared_high) + '[rad]'), color='olive')

plt.xlabel('E_nu [MeV]')
plt.yscale('log')
plt.ylabel('Expected counts (some units, doesn\'t matter)')
plt.title('Expected geonus, theta = 0.5903, different delta_m_12^2 values, different grids')

# Set lower limit for y-axis to 10^-7
plt.ylim(bottom=1e-11)

# Enable minor ticks on x-axis
plt.minorticks_on()

# Add legend
plt.legend(loc='upper right')
plt.savefig('Expected geonu, standard theta, max mantle distance ' + str(max_mantle_distance) + " Different Grids.pdf",
            format='pdf')

# plt.show()

plt.figure()
plt.plot(energy_array, (N_U_delta_m_low_c + N_Th_delta_m_low_c + N_U_delta_m_low_m + N_Th_delta_m_low_m) / (
            N_U_c + N_Th_c + N_U_m + N_Th_m), label='ratio', color='green')
plt.xlabel('E_nu [MeV]')
plt.ylabel('ratio')
plt.title('Ratio of expected geonus, low delta_m^2 / standard delta_m^2')

# Set lower limit for y-axis to 10^-7
# plt.ylim(bottom=1e-7)

# Enable minor ticks on x-axis
plt.minorticks_on()

# Add legend
plt.legend(loc='upper right')
plt.savefig('Ratio of expected geonu, standard theta, low delta m, max mantle distance ' + str(
    max_mantle_distance) + " Different Grids.pdf", format='pdf')

# plt.show()
plt.figure()
plt.plot(energy_array, (N_U_delta_m_high_c + N_Th_delta_m_high_c + N_U_delta_m_high_m + N_Th_delta_m_high_m) / (
            N_U_c + N_Th_c + N_U_m + N_Th_m), label='ratio', color='olive')
plt.xlabel('E_nu [MeV]')
plt.ylabel('ratio')
plt.title('Ratio of expected geonus, high delta m / standard delta m')

# Set lower limit for y-axis to 10^-7
# plt.ylim(bottom=1e-7)

# Enable minor ticks on x-axis
plt.minorticks_on()

# Add legend
plt.legend(loc='upper right')
plt.savefig('Ratio of expected geonu, standard theta, high delta m,, max mantle distance ' + str(
    max_mantle_distance) + " Different Grids.pdf", format='pdf')

# plt.show()

print("plots saved")

print("   ")

print("All done! Congratulations!")

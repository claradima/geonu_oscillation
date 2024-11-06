import numpy as np
import matplotlib.pyplot as plt
import sys
import gc
import psutil
import os
import argparse


### Prints current memory usage ###

def get_memory_usage():

    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    # Convert bytes to megabytes
    mem_usage_mb = mem_info.rss / (1024 ** 2)
    print(f"Current memory usage: {mem_usage_mb:.2f} MB")
    return mem_usage_mb
    
### Generates cube grid
### inputs : grid_counts = the number of points along one dimension
###          halfside = half side length of cube [km] ;
###          default radius of Earth 6371
### returns : x_coords, y_coords, z_coords , dist_sq

def create_3d_grid(grid_counts, halfside = 6371):

    coords = np.linspace(- halfside, halfside, grid_counts)
    grid_1d_size = coords[1] - coords[0]

    x_coords, y_coords, z_coords = np.meshgrid(coords, coords, coords)
    dist_sq = x_coords**2 + y_coords**2 + z_coords**2
    grid_specs = [x_coords, y_coords, z_coords, dist_sq]

    return grid_specs

### Creates spherical shell from 3d Cartesian grid
### inputs : inner_rad; outer_rad, sublayers (number of subls)
###          equal height; grid_specs = x, y, z from meshgrid
###          dist^2 = x^2 + y^2 + z^2
### returns : x_coords, y_coords, z_coords , dist_sq

def cut_shell(inner_rad, outer_rad, sublayers, grid_specs):
    x_coords, y_coords, z_coords, dist_sq = grid_specs
    b_values = np.linspace(inner_rad, outer_rad, sublayers + 1)

    shell_grids = np.empty(len(b_values) - 1, dtype = object)
    shell_indices_array = np.empty(len(b_values) - 1, dtype = object)

    for i in range(len(b_values) - 1):
        # print(i)
        shell_indices = np.logical_and(b_values[i] ** 2 < dist_sq, dist_sq <= b_values[i + 1] ** 2)
        # Attempt to create the DM_grid
        shell_grid = np.stack((x_coords[shell_indices], y_coords[shell_indices], z_coords[shell_indices]),
                                   axis=-1)

        # Store the DM_grid and DM_indices in the object arrays
        shell_grids[i] = shell_grid
        shell_indices_array[i] = shell_indices
                # good_DM_distances_array[i] = distances_squared[DM_indices]**(0.5)

        print(f"Subshell {i + 1} grid done (out of {sublayers}); limits: {b_values[i]} , {b_values[i + 1]}")

    return shell_grids

### Set abundances and densities for crust, mantle, CLM, DM, EM
###
### inputs : set_name : pick out of low, mid, high
### returns :  A_Th_c, A_U_c, A_Th_CLM, A_U_CLM, A_Th_DM,
###            A_U_DM, A_Th_EM, A_U_EM, rho_c, rho_CLM, rho_DM,
###            rho_EM

### TO DO :
###
### Should put these sets in an external file and have them read
### off from there to make things more flexible in terms of
### adding new Earth models;

def set_abund (set_name):
    if set_name not in ['low', 'mid', 'high']:
        print('invalid abundance set name; select from mid, low or high')
        print('can add new sets in function def in functions.py')
        print('you can even make this nicer and move the data in an external file, have abundances read off from there together with other details of the earth model :)')
        print('up to you tho no pressure')
        sys.exit()

    # Placeholder definitions for IDE recognition
    A_Th_c = 0
    A_U_c = 0
    A_Th_CLM = 0
    A_U_CLM = 0
    A_Th_DM = 0
    A_U_DM = 0
    A_Th_EM = 0
    A_U_EM = 0
    rho_c = 0
    rho_CLM = 0
    rho_DM = 0
    rho_EM = 0

    low = {
        "A_Th_c": 5 * (10 ** (-6)),
        "A_U_c": 1 * (10 ** (-6)),
        "A_Th_CLM": 53 * (10 ** (-9)),
        "A_U_CLM": 13 * (10 ** (-9)),
        "A_Th_DM": 17.52 * (10 ** (-9)),
        "A_U_DM": 6.4 * (10 ** (-9)),
        "A_Th_EM": 90 * (10 ** (-9)),
        "A_U_EM": 15 * (10 ** (-9)),
        "rho_c": 2.7,
        "rho_CLM": 3.3,
        "rho_DM": 3.3,
        "rho_EM": 3.3
    }

    mid = {
        "A_Th_c": 5 * (10 ** (-6)),
        "A_U_c": 1 * (10 ** (-6)),
        "A_Th_CLM": 147 * (10 ** (-9)),
        "A_U_CLM": 33 * (10 ** (-9)),
        "A_Th_DM": 21.9 * (10 ** (-9)),
        "A_U_DM": 8 * (10 ** (-9)),
        "A_Th_EM": 147 * (10 ** (-9)),
        "A_U_EM": 33 * (10 ** (-9)),
        "rho_c": 2.7,
        "rho_CLM": 3.3,
        "rho_DM": 3.3,
        "rho_EM": 3.3
    }

    high = {
        "A_Th_c": 5 * (10 ** (-6)),
        "A_U_c": 1 * (10 ** (-6)),
        "A_Th_CLM": 427 * (10 ** (-9)),
        "A_U_CLM": 83 * (10 ** (-9)),
        "A_Th_DM": 26.28 * (10 ** (-9)),
        "A_U_DM": 9.6 * (10 ** (-9)),
        "A_Th_EM": 221 * (10 ** (-9)),
        "A_U_EM": 57 * (10 ** (-9)),
        "rho_c": 2.7,
        "rho_CLM": 3.3,
        "rho_DM": 3.3,
        "rho_EM": 3.3
    }

    abundances = set_name
    globals().update(abundances)

    return A_Th_c, A_U_c, A_Th_CLM, A_U_CLM, A_Th_DM, A_U_DM, A_Th_EM, A_U_EM, rho_c, rho_CLM, rho_DM, rho_EM


### Set fixed parameters : this code is intended to run with
### fixed values for all oscillation parameters except theta_12
### and delta_m_21^2; run this function at the start of your code
###
### NOTE : all standard params from PDG 2020 :
###        theta_12 = 0.5903 #rad
###        theta_23 = 0.8430 #rad
###        theta_13 = 0.1503 #rad
###
###        delta_m_21_squared = 7.39 * 10**(-5) #eV^2
###        delta_m_32_squared = 2.449 * 10**(-3) #eV^2
###        delta_m_31_squared = delta_m_32_squared

def set_fixed_params():
    global theta_23, theta_13, delta_m_32_squared, delta_m_31_squared
    print("setting standard oscillation parameters")
    theta_23 = 0.8430  # rad
    theta_13 = 0.1503  # rad

    delta_m_32_squared = 2.449 * 10 ** (-3)  # eV^2
    delta_m_31_squared = delta_m_32_squared

### Set particle masses and threshold energy for IBD : there is
### no reason for these to change really; just run this at the start
### of your code
###
### NOTE : don't really need these globally, just for sigma_IBD
###        they used to be defined inside that function initially
###        but I liked them better here, I don't think it matters
###        much

def set_masses():
    global m_e, m_p, m_n, E_th
    print('setting masses and E_th for IBD')
    m_e = 0.511  # MeV
    m_p = 938  # MeV
    m_n = 941  # MeV
    E_th = 1.806  # MeV

def set_lambda_mu():
    global lambda_U, lambda_Th, mu_U, mu_Th
    print('setting lambda and mu values for U and Th')
    lambda_U = 4.916
    lambda_Th = 1.563
    mu_U = 235
    mu_Th = 232

### Setting energy array
###
### NOTE : The limits of the range should be unchanged
###        I already defined E_Th globally so I guess that's
###        there already, I might as well use it;
###        SO : only use this after running set_masses()
###
### NOTE : Upper limit fixed 3.3 [MeV]; can change below but don't

def set_energy_array(no_E_bins):
    global energy_array
    print('setting energy array')
    energy_array = np.linspace(E_th, 3.3, no_E_bins)

### Calculate distance between a set of points and SNO+
### inputs : points_array (must be a 2D array where each element
###          is (x, y, z) with the coords of each point
###          SNO_r = position of SNO+ in Cartesian coords;
###          default : np.array([0, 0, 6369])
### returns: relative_distances : array of the same length as
###          points_array, same index - corresponding point

def calc_relative_dist(points_array, SNO_r = np.array([0, 0, 6369])):
    # Calculate the Euclidean distance using vectorized operations
    print("   ")
    if SNO_r == np.array([0, 0, 6369]):
        print('Position of SNO+ set to (0, 0, 6369) by default')
    relative_distances = np.linalg.norm(points_array - SNO_r, axis=1)

    print("Computed relative distances from Earth grid points to SNO+")
    print("   ")

    return relative_distances

### Calculate sigma_IBD for each energy in given energy_array
###
### NOTE : Not sure using the global var is the best way to go
###        about this; will have to do for now
### TO DO (MAYBE) : change things slightly

def calc_sigma_IBD(energy_array):
    global sigma
    sigma = ((energy_array - E_th + m_e) ** 2) * ((1 - (m_e) ** 2 / ((energy_array - E_th + m_e) ** 2)) ** (1 / 2))

### Calculate Delta_ij for given points and energies
### inputs : points_array (must be a 2D array where each element
###          is (x, y, z) with the coords of each point
###          energy_array
###          delta_m_ij_squared : default values given (best fit
###          values from PDG (I think 2020)
###          SNO_r : SNO+ position in Cartesian coordinates
###          default np.array([0, 0, 6369])
### returns: Delta : 2D array, one element corresponds to one
###          energy value and one fix point;
###          size(Delta) = (len(energy_array), len(points_array)

### TO DO :
###
### I think it's very inefficient to have Delta take points_array
### as an argument and calculate their relative position from
### SNO+ each time, especially given that we're not actually working
### with that many different sets of points
###
### Will add new function to calculate Delta given just the
### relative distance

### TO DO :
###
### Write default values for standard oscillation parameters

def Delta_ij(energy_array, points_array, delta_m_ij_squared, SNO_r= np.array([0, 0, 6369])):
    # Calculate relative distances
    relative_distance_array = calc_relative_dist(points_array, SNO_r)

    # Reshape energy_array to perform element-wise division
    energy_array_reshaped = energy_array.reshape(-1, 1)

    # Calculate Delta using vectorized operations
    Delta = (1.27 * delta_m_ij_squared * relative_distance_array * 1000) / (energy_array_reshaped)

    print('computed; deleting useless stuff')
    get_memory_usage()
    del relative_distance_array
    del energy_array_reshaped
    print('deleted')
    get_memory_usage()

    return Delta


### Calculate average distance per layer
### inputs : points_array (must be a 2D array where each element
###          is (x, y, z) with the coords of each point
###          energy_array
###          SNO+ position in Cartesian coordinates
###          default np.array([0, 0, 6369])

def layer_avg_dist (points_array, SNO_r = np.array([0, 0, 6369])):

    print("   ")
    if SNO_r == np.array([0, 0, 6369]):
        print('Position of SNO+ set to (0, 0, 6369) by default')
    relative_distance_array = calc_relative_dist(points_array, SNO_r)

    avg_dist = np.mean(relative_distance_array)

    return avg_dist

### Calculate overall average distance
###
### NOTE: This is super non-flexible, it relies on the Earth model
###       A LOT (and how I split into layers and sublayers]
###       This is because the DM and EM layers being split
###       resulted in them being packed in objects that are arrays of
###       arrays, and the latter have different sizes
###       I'm sure it's not too hard to change, but I couldn't
###       be bothered, I can't imagine having to use this in
###       the future tbh
###
### NOTE : because of the previous note, you must've already
###        made the objects with names matching the input names
###
### NOTE : It is assumed that the whole mantle has equally spaced
###        points
###
### inputs : crust_grid, CLM_grid, DM_grids, EM_grids,
##           grid_1d_size_crust, grid_1d_size_mantle
### returns : average_distance

def avg_dist(crust_grid, CLM_grid, DM_grids, EM_grids, grid_1d_size_crust, grid_1d_size_mantle):

    crust_relative_distances = calc_relative_dist(crust_grid)
    CLM_relative_distances= calc_relative_dist(CLM_grid)

    DM_relative_distances_array = np.empty(len(DM_grids), dtype=object)

    for i in range(len(DM_grids)):
        DM_relative_distances_array[i] = calc_relative_dist(DM_grids[i])

    EM_relative_distances_array = np.empty(len(EM_grids), dtype=object)
    for i in range(len(EM_grids)):
        EM_relative_distances_array[i] = calc_relative_dist(EM_grids[i])

    EM_relative_distances = np.concatenate(EM_relative_distances_array)
    DM_relative_distances = np.concatenate(DM_relative_distances_array)

    mantle_relative_distances = np.concatenate((CLM_relative_distances, DM_relative_distances, EM_relative_distances))

    crust_avg_dist = np.mean(crust_relative_distances)
    mantle_avg_dist = np.mean(mantle_relative_distances)

    print(f'Crust average distance : {crust_avg_dist}')
    print(f'Mantle average distance (grid 1d size mantle : {grid_1d_size_mantle}) : {mantle_avg_dist}')

    # To get the overall average, need to weight by volume element

    grid_volume_el_crust = grid_1d_size_crust ** 3
    grid_volume_el_mantle = grid_1d_size_mantle ** 3

    average_distance = (grid_volume_el_mantle * np.sum(mantle_relative_distances) + grid_volume_el_crust * np.sum(crust_relative_distances)) / (grid_volume_el_mantle * len(mantle_relative_distances) + grid_volume_el_crust * len(crust_relative_distances))

    print(f'Average distance (grid 1d size mantle : {grid_1d_size_mantle}) : {average_distance}')

    return average_distance

### Calculate survival prob for given points and energies
### inputs : points_array (must be a 2D array where each element
###          is (x, y, z) with the coords of each point
###          energy_array
###          theta_12, delta_m_21_squared : default values given (best fit
###          values from PDG (I think 2020)
### NOTE : other oscillation parameters should stay fixed;
###        define them at the start of running any code by
###        running set_fixed_params();
###
### returns: P_ee : 2D array, one element corresponds to one
###          energy value and one fix point;
###          size(P_ee) = (len(energy_array), len(points_array)

### TO DO
###
### When you look at the Delta TO DO, think about whether any of
### the things calculated here can also be taken out and calculated
### separately to avoid repeating operations


### TO DO
### Add a version of this function that takes as input the
### relative distances directly rather than the points array
### Can do this because we only care about distance really

def calc_P_ee(energy_array, points_array, theta_12, theta_13, delta_m_21_squared):

    relative_distance_array = calc_relative_dist(points_array)

    Delta_31 = Delta_ij(energy_array, points_array, delta_m_31_squared)
    print("Delta_31 computed successfully")
    Delta_12 = Delta_ij(energy_array, points_array, delta_m_21_squared)
    print("Delta_12 computed successfully")

    A = (np.cos(theta_13)) ** 4 * (np.sin(2 * theta_12)) ** 2
    B = np.sin(2 * theta_13) ** 2

    sin_squared_Delta_31 = np.sin(Delta_31) ** 2
    sin_squared_Delta_12 = np.sin(Delta_12) ** 2

    print("terms computed successfully")

    P_ee = 1 - (A * sin_squared_Delta_12 + B * sin_squared_Delta_31)
    print("P_ee computed successfully")

    return P_ee

### Calculate volume integral from geonu flux formula, custom
### oscillation parameters
###
### NOTE : this has to be calculated separately for the two
### emmitting species because the other factors depend on
### emitter as well;
###
### inputs : points_array (must be a 2D array where each element
###          is (x, y, z) with the coords of each point
###          energy_array
###          grid_1d_size : distance betweem two consecutive
###          grid points in 1d; used to compute volume element
###          theta_12, delta_m_21_squared : default values given (best fit
###          values from PDG (I think 2020)
###          A_Th, A_U, rho : abundances and density in layer
###
### NOTE : other oscillation parameters should stay fixed;
###        define them at the start of running any code by
###        running set_fixed_params();
###
### NOTE : this is intended to be used to compute the integral in
###        one sublayer, where abundances and density are assumed
###        to be constant
###
### returns: sum_Th, sum_U = values of integral for each emitter
###          and each energy value :
##           len(sum_Th) = len(sum_U) = len(energy_array)

### TO DO
###
### When you look at the Delta TO DO, think about whether any of
### the things calculated here can also be taken out and calculated
### separately to avoid repeating operations

### TO DO
### Add a version of this function that takes as input the
### relative distances directly rather than the points array
### Can do this because we only care about distance really

def vol_integral(points_array, energy_array, grid_1d_size, theta_12, delta_m_21_squared, A_Th, A_U,
                                 rho):
    dV = grid_1d_size ** 3

    relative_distance_array = calc_relative_dist(points_array)
    print("Relative distance array computed successfully")
    P_ee_array = calc_P_ee(energy_array, points_array, theta_12, delta_m_21_squared)
    print("P_ee_array computed successfully")

    # Compute sum_Th
    sum_Th = np.sum(P_ee_array * ((A_Th * rho) / (4 * np.pi * (relative_distance_array ** 2)))[np.newaxis, :] * dV,
                    axis=1)

    print("sum_Th computed successfully")

    # Compute sum_U
    sum_U = np.sum(P_ee_array * ((A_U * rho) / (4 * np.pi * (relative_distance_array ** 2)))[np.newaxis, :] * dV,
                   axis=1)
    print("sum_U computed successfully")

    print('computed; deleting useless stuff')
    get_memory_usage()
    del P_ee_array
    del relative_distance_array
    print('deleted')
    get_memory_usage()

    return sum_Th, sum_U

### Calculate volume integral from geonu flux formula, use
### P_ee = constant as an approximation
###
### NOTE : this has to be calculated separately for the two
### emmitting species because the other factors depend on
### emitter as well;
###
### inputs : points_array (must be a 2D array where each element
###          is (x, y, z) with the coords of each point
###          energy_array
###          grid_1d_size : distance betweem two consecutive
###          grid points in 1d; used to compute volume element
###          A_Th, A_U, rho : abundances and density in layer
###
### NOTE : don't need oscillation parameters here because P_ee
###        is just assumed to be constant ; The value of this should
###        be computed separately and just passed as an argument
###
### NOTE : this is intended to be used to compute the integral in
###        one sublayer, where abundances and density are assumed
###        to be constant
###
### returns: sum_Th, sum_U  - three different versions of these:
###          mid, low, high

### TO DO
###
### I mean ... this is just ugly as hell, can always make it better

### TO DO
###
### Version of this that does error propagation correctly
### instead of just calculating values at midpoint value +- stdev

### TO DO
### Add a version of this function that takes as input the
### relative distances directly rather than the points array
### Can do this because we only care about distance really

def vol_integral_const_P_ee(points_array, energy_array, grid_1d_size, A_Th, A_U, rho, P_ee_mid, P_ee_stdev):

    dV = grid_1d_size**3
    relative_distance_array = calc_relative_dist(points_array)

    P_ee_low = P_ee_mid - P_ee_stdev
    P_ee_high = P_ee_mid + P_ee_stdev
    P_ee_array_mid = np.full((len(energy_array), len(points_array)), P_ee_mid)
    P_ee_array_low = np.full((len(energy_array), len(points_array)), P_ee_low)
    P_ee_array_high = np.full((len(energy_array), len(points_array)), P_ee_high)

    # Compute sum_Th
    sum_mid_Th = np.sum(
        P_ee_array_mid * ((A_Th * rho) / (4 * np.pi * (relative_distance_array ** 2)))[np.newaxis, :] * dV, axis=1)
    sum_low_Th = np.sum(
        P_ee_array_low * ((A_Th * rho) / (4 * np.pi * (relative_distance_array ** 2)))[np.newaxis, :] * dV, axis=1)
    sum_high_Th = np.sum(
        P_ee_array_high * ((A_Th * rho) / (4 * np.pi * (relative_distance_array ** 2)))[np.newaxis, :] * dV, axis=1)

    # Compute sum_U
    sum_mid_U = np.sum(
        P_ee_array_mid * ((A_U * rho) / (4 * np.pi * (relative_distance_array ** 2)))[np.newaxis, :] * dV, axis=1)
    sum_low_U = np.sum(
        P_ee_array_low * ((A_U * rho) / (4 * np.pi * (relative_distance_array ** 2)))[np.newaxis, :] * dV, axis=1)
    sum_high_U = np.sum(
        P_ee_array_high * ((A_U * rho) / (4 * np.pi * (relative_distance_array ** 2)))[np.newaxis, :] * dV, axis=1)

    return sum_mid_Th, sum_mid_U, sum_low_Th, sum_low_U, sum_high_Th, sum_high_U


### Rebins counts data.
###
### inputs : initial_bins: array, bin edges of the initial data
###         counts_in_initial_bins: array, counts in each
###         initial bin
###         final_bin_midpoints: array, midpoints of the final
###         desired bins
###
### returns : counts_in_final_bins: array, counts in each final
### bin
###
### NOTE : used to get emission fluxes

def rebin_counts(initial_bins, counts_in_initial_bins, final_bin_midpoints):

    # Calculate bin midpoints of the initial bins
    bin_midpoints = (initial_bins[:-1] + initial_bins[1:]) / 2

    # Use np.histogram to calculate counts in final bins
    counts_in_final_bins, _ = np.histogram(initial_bins, bins=np.concatenate([initial_bins, [2 * initial_bins[-1] - initial_bins[-2]]]), weights=counts_in_initial_bins)

    # Interpolate the counts to the final bin midpoints
    counts_in_final_bins = np.interp(final_bin_midpoints, bin_midpoints, counts_in_final_bins[:-1])

    return counts_in_final_bins

### Get emission fluxes and rebin to match desired energy bins
### (with centers on the points from energy_arrray)
###
### Note: this function required the files 'U238_spectrum.txt'
### and 'Th232_spectrum.txt'; check explanation doc for source
###
### inputs : energy_array
###          plot_spectrum : bool; say if you want to plot
###          the 'raw' spectrum (not rebinned)
###
### returns : dn_dE_rebinned_U, dn_dE_rebinned_Th - arrays, same
###           len as energy_array

### TO DO : for this and all other plotting bits, figure out a
###         nice flexible way to say if you want to save them or not

### TO DO : also for this and all other plotting bits, should
### probably save results in csv files or something like that
### as well

def get_emission_fluxes(energy_array, plot_spectrum = False):
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

    if plot_spectrum:
        print('plotting emission spectrum')
        # Plot U238 decay data (blue line)
        plt.plot(energy_array_U, dn_dE_U, label='U238 decays', color='blue')
        # Plot Th232 decay data (red line)
        plt.plot(energy_array_Th, dn_dE_Th, label='Th232 decays', color='red')
        plt.xlabel('E_nu [MeV]')
        plt.yscale('log')
        plt.ylabel('Intensity (arbitrary units)')
        plt.title('Geonu emission spectrum')
        # Add shaded region between 1.8 MeV and 3.3 MeV
        plt.axvspan(1.8, 3.3, alpha=0.3, color='gray')
        plt.minorticks_on()
        plt.legend(loc='upper right')
        plt.show()
        # plt.savefig("Emission.pdf", format='pdf')

    print("rebin to match energy array")

    dn_dE_rebinned_U = rebin_counts(energy_array_U, dn_dE_U, energy_array)
    dn_dE_rebinned_Th = rebin_counts(energy_array_Th, dn_dE_Th, energy_array)

    return dn_dE_rebinned_U, dn_dE_rebinned_Th
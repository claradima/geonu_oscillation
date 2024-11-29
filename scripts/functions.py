import numpy as np
import matplotlib.pyplot as plt
import sys
import gc
import psutil
import os
import argparse
import csv
import pandas as pd


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
    # would need 637.1 for 20 km spacing
    # 1275 for 5 km spacing
    # 320 for 40 km spacing
    # 285 for 45 km spacing
    # 215 for 60km spacing
    # 160 for 80km spacing

    coords = np.linspace(- halfside, halfside, grid_counts)
    grid_1d_size = coords[1] - coords[0]
    print('1d coords created')
    get_memory_usage()
    print('meshgriding now')
    x_coords, y_coords, z_coords = np.meshgrid(coords, coords, coords)
    print('done')
    get_memory_usage()
    print('don\'t need 1d coords anymore; delete')
    del coords
    print('deleted')
    get_memory_usage()
    print('computing distances and putting grid specs together')
    dist_sq = x_coords**2 + y_coords**2 + z_coords**2
    grid_specs = [x_coords, y_coords, z_coords, dist_sq]
    print('grid specs put together')
    get_memory_usage()
    print('checking if deleting separate x_coords does anything')
    del x_coords
    del y_coords
    del z_coords
    del dist_sq
    print('deleted')
    get_memory_usage()

    return grid_specs, grid_1d_size

### Creates spherical shell from 3d Cartesian grid
### inputs : inner_rad; outer_rad, sublayers (number of subls)
###          equal height; grid_specs = x, y, z from meshgrid
###          dist^2 = x^2 + y^2 + z^2
### returns : x_coords, y_coords, z_coords , dist_sq

### TO DO : I think defining cut_shells in this way avoids my issue of
###         having different data types if you have one or multiple sublayers

### Layered Earth Model : crust 6350 - 6371 km
###                       CLM 6196 - 6350 km
###                       DM 4202 - 6196 km
###                       EM 3480 - 4202 km

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

def trim_shell(inner_rad, outer_rad, points_array):
    print(f'len of initial array: {len(points_array)}')
    # Compute the squared radii once to avoid redundant calculations
    inner_rad_sq = inner_rad**2
    outer_rad_sq = outer_rad**2

    # Use a boolean mask to filter points based on the condition
    mask = np.sum(points_array**2, axis=1) >= inner_rad_sq
    mask &= np.sum(points_array**2, axis=1) <= outer_rad_sq

    # Return only the points that satisfy the condition, using the mask

    trimmed_arr = points_array[mask]
    print(f'len of trimmed array: {len(trimmed_arr)}')
    del points_array
    del mask
    del inner_rad_sq
    del outer_rad_sq

    return trimmed_arr

### Compute the total volume of a shell

def get_shell_volume(shell_grids, grid_1d_size):
    dV = grid_1d_size**3
    print(f'distance between points : {grid_1d_size}')
    print(f'volume element[km^3] is : {dV}')

    print('computing number of points in shell')

    no_points = 0
    for i in range(len(shell_grids)):
        subshell_points = len(shell_grids[i])
        print(f'there are {subshell_points} points in subshell {i} out of {len(shell_grids)}')
        no_points += subshell_points
        print('points from subshell added to sum')
        print(f'total points (intermediate): {no_points}')

    print(' ')
    print(f'total number of points in shell : {no_points}')

    volume = dV * no_points
    print(f'total shell volume : {volume}')

    return no_points, volume



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

def set_abund(set_name):
    # Check if the set_name is valid
    if set_name not in ['low', 'mid', 'high']:
        print('Invalid abundance set name; select from low, mid, or high.')
        print('You can add new sets in the function definition in functions.py.')
        print('Consider storing the data in an external file for flexibility.')
        sys.exit()

    # Abundance data for each set
    abundances = {
        "low": {
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
        },
        "mid": {
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
        },
        "high": {
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
    }

    # Get the correct abundance set
    if set_name not in abundances:
        print("Error: Invalid abundance set.")
        sys.exit()

    # Retrieve the abundances from the selected set
    selected_abundance = abundances[set_name]

    # Print selected abundance set for debugging
    print(f"Selected {set_name} abundance: {selected_abundance}")

    # Return all the values from the selected set
    return (selected_abundance["A_Th_c"], selected_abundance["A_U_c"], selected_abundance["A_Th_CLM"],
            selected_abundance["A_U_CLM"], selected_abundance["A_Th_DM"], selected_abundance["A_U_DM"],
            selected_abundance["A_Th_EM"], selected_abundance["A_U_EM"], selected_abundance["rho_c"],
            selected_abundance["rho_CLM"], selected_abundance["rho_DM"], selected_abundance["rho_EM"])


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
###
###        params from James' fit (constrained)
###
###        theta_12 = 33.7 +- 0.8 [deg]
###        theta_12 = 0.5882 +- 0.0140 [rad]
###        delta_m_21_squared = 7.58 (+0.18; -0.17) [10**(-5) eV^2]

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
    lambda_U = 4.916 #[10^(-18) s^(-1)]
    lambda_Th = 1.563 #[10^(-18) s^(-1)]
    mu_U = 235 #[g/mol]
    mu_Th = 232 #[g/mol]

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
    energy_array = np.linspace(round(E_th, 2), 3.3, no_E_bins)

    # this would be an issue if E_th would be rounded down

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
    if np.array_equal(SNO_r, np.array([0, 0, 6369])):
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

def calc_sigma_IBD():
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

def calc_P_ee(points_array, theta_12, delta_m_21_squared):

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

    get_memory_usage()
    print('deleting intermediary quantities after calculating P_ee')
    del Delta_12
    del Delta_31
    del sin_squared_Delta_12
    del sin_squared_Delta_31
    del A
    del B
    print('deleted')
    get_memory_usage()

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
###          energy_array defined globally! no input
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

def vol_integral(points_array, grid_1d_size, theta_12,
                 delta_m_21_squared, A_Th, A_U, rho):
    print('test test') #TEMP
    dV = grid_1d_size ** 3
    relative_distance_array = calc_relative_dist(points_array)
    print(f'relative_distance_array : {relative_distance_array}') # TEMP
    print("Relative distance array computed successfully")
    P_ee_array = calc_P_ee(points_array, theta_12, delta_m_21_squared)

    print("P_ee_array computed successfully")
    print(f'P_ee_array : {P_ee_array}')#TEMP
    print(f'shape P_ee_array : {P_ee_array.shape}')#TEMP

    print(f'P_ee_array shape: {P_ee_array.shape}')
    print(f'rho: {rho}')
    print(f'dV: {dV}')
    # Compute sum_Th
    sum_Th = np.sum(P_ee_array * ((A_Th * rho) / (4 * np.pi * (relative_distance_array ** 2)))[np.newaxis, :] * dV,
                    axis=1)
    print(f'sum_Th : {sum_Th}')#TEMP

    print("sum_Th computed successfully")

    # Compute sum_U
    sum_U = np.sum(P_ee_array * ((A_U * rho) / (4 * np.pi * (relative_distance_array ** 2)))[np.newaxis, :] * dV,
                   axis=1)
    print("sum_U computed successfully")
    print(f'sum_U : {sum_U}')
    print('computed; deleting useless stuff')
    get_memory_usage()
    del P_ee_array
    del relative_distance_array
    print('deleted')
    get_memory_usage()

    return sum_Th, sum_U


### Adds up volume integrals for sublayers in a layer
###
### inputs : same as vol_integrals, except points_array;
###          instead, array of sets of points (so
###          grids_array[i] has the same format as points_array
###
### returns : total_Th, total_U : values of integral for Th,
###           U respectively for whole layer made up of
###           sublayers
###
### NOTE : Layers should be quite easy too add, there shouldn't
###        be that many of them so not really worth adding a
###        function to add them up? (perhaps TO DO)
###
### TO DO : Test this when you get the final script!! compare
###         plots with plots made with Optimized_Memory_v4.ipynb

def add_vol_integrals(grids_array, grid_1d_size, theta_12, delta_m_21_squared, A_Th, A_U, rho):

    print('adding contributions from sublayers')
    total_Th = np.zeros(energy_array.shape)
    total_U = np.zeros(energy_array.shape)


    for i in range(len(grids_array)):
        integral_Th, integral_U = vol_integral(grids_array[i], grid_1d_size, theta_12, delta_m_21_squared,A_Th, A_U, rho)
        print(f'integral_Th : {integral_Th}')
        print(f'integral_U : {integral_U}')


        total_Th = total_Th + integral_Th
        total_U = total_U + integral_U

        print(f'added contribution from sublayer {i}')

    print(' ')
    print('full layer contribution computed')

    return total_Th, total_U


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

def vol_integral_const_P_ee(points_array, grid_1d_size, A_Th, A_U, rho, P_ee_mid, P_ee_stdev):

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

### Adds up volume integrals for sublayers, constant P_ee
###
### inputs : same as vol_integral_conts_P_ee, except points_array;
###          instead, array of sets of points (so
###          grids_array[i] has the same format as points_array
###
### returns : total_mid_Th, total_mid_U, total_low_Th,
###           total_low_U, total_high_U, total_high_Th :
###           values of integral for Th, U respectively for
###           whole layer made up of sublayers
###
### NOTE : Layers should be quite easy too add, there shouldn't
###        be that many of them so not really worth adding a
###        function to add them up? (perhaps TO DO)
###
### TO DO : Test this when you get the final script!! compare
###         plots with plots made with Optimized_Memory_v4.ipynb

def add_vol_integrals_const_P_ee(grids_array, grid_1d_size, A_Th, A_U, rho, P_ee_mid, P_ee_stdev):

    print('adding contributions from sublayers')
    total_mid_Th = np.zeros(energy_array.shape)
    total_mid_U = np.zeros(energy_array.shape)
    total_low_Th = np.zeros(energy_array.shape)
    total_low_U = np.zeros(energy_array.shape)
    total_high_Th = np.zeros(energy_array.shape)
    total_high_U = np.zeros(energy_array.shape)

    for i in range(len(grids_array)):
        integral_mid_Th, integral_mid_U, integral_low_Th, integral_low_U, integral_high_Th, integral_high_U = vol_integral_const_P_ee(grids_array[i], grid_1d_size, A_Th, A_U, rho, P_ee_mid, P_ee_stdev)
        total_mid_Th = total_mid_Th + integral_mid_Th
        total_mid_U = total_mid_U + integral_mid_U
        total_low_Th = total_low_Th + integral_low_Th
        total_low_U = total_low_U + integral_low_U
        total_high_Th = total_high_Th + integral_high_Th
        total_high_U = total_high_U + integral_high_U

        print(f'added contribution from sublayer {i}')

    print(' ')
    print('full layer contribution computed')

    return total_mid_Th, total_mid_U, total_low_Th, total_low_U, total_high_Th, total_high_U

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
### creates globals : dn_dE_rebinned_U, dn_dE_rebinned_Th -
### arrays, same len as energy_array

### TO DO : for this and all other plotting bits, figure out a
###         nice flexible way to say if you want to save them or not

### TO DO : also for this and all other plotting bits, should
### probably save results in csv files or something like that
### as well

def get_emission_fluxes(plot_spectrum = False):

    global dn_dE_rebinned_U, dn_dE_rebinned_Th

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

### Calculate expected fluxes (arbitrary units)
### To use this, need lambda values, mu values to be set and
### sigma, volume integral and rebinned emission spectrum to
### be computed beforehand
###
### NOTE : sigma is a predefined global var, so don't need
###        to add as input; lambdas, mus and emissons also;
###        volume ints are assigned to variable names
###
### TO DO : (maybe) set energy array as global variable?
###         OR set the others as not global; I think this
###         would be better, so should change emission and
###         sigma

def calc_exp_spec(U_vol_int, Th_vol_int):
    N_Th = ((lambda_Th)/(mu_Th)) * sigma * dn_dE_rebinned_Th * Th_vol_int
    N_U = ((lambda_U)/(mu_U)) * sigma * dn_dE_rebinned_U * U_vol_int

    return N_Th, N_U

### Calculate expected geonu counts
### To get the number of geonu per second expected in each
### energy bin, we have to multiply by 3.268 * 10^(-6) * dE
### where dE is the width of the energy bin that we are using
###
### To get the total number of expected geonu for some
### total livetime, simply multiply by the livetime in s
###
### Further details about how the constant above was obtained
### will be added to explanation document (if not there already)
###
### inputs : N_U, N_Tot - spectra in arbitrary units,
###          obtained with function above
###          energy_array : to calc dE; defined globally
###          livetime : total observation time [s]
### returns : N_U_scaled, N_Th_scaled : scaled spectra from
###           each emitter
###           geonus_tot : total number of geonus expected
###           in the given livetime across all energies
###
### Note : can use the plotting functions for scaled or
###        'unscaled' spectra, but they have to be consistent
###        if we compute flux ratios

def calc_exp_spec_scale(N_Th, N_U, livetime):
    scaling_ct = 0.3268 * 1000
    # TO DO : need to figure out scaling factor !!!!
    dE = energy_array[1] - energy_array[0]

    N_Th_scaled = N_Th * scaling_ct * dE * livetime
    N_U_scaled = N_U * scaling_ct * dE * livetime

    print(f'scaling_ct : {scaling_ct}') #TEMP
    print(f'dE : {dE}') #TEMP
    print(f'livetime : {livetime}') #TEMP
    print(f'Overall scaling : {scaling_ct * dE * livetime}') #TEMP

    print(f'Computed scaled spectra for a livetime of {livetime} s = {livetime/86400} days')

    geonus_tot_Th = np.sum(N_Th_scaled)
    geonus_tot_U = np.sum(N_U_scaled)
    geonus_tot = geonus_tot_Th + geonus_tot_U

    print(f'Total expected geonus from Th decays : {geonus_tot_Th}')
    print(f'Total expected geonus from U decays : {geonus_tot_U}')
    print(f'Total expected geonus from any decays: {geonus_tot}')

    return N_Th_scaled, N_U_scaled, geonus_tot

##########################################################
### PLOTTING FUNCTIONS ###################################
##########################################################

### NOTE : spec_save can be a bool set with argsparse
###        if set to True, save all plots + data as csv
### ADD THIS TO ALL PLOTTING FUNCTIONS

### Plot spectrum
###
### inputs : energy_array defined globally
###          N_Th, N_U : arays of same len
###          spec_save : bool (set true if you want plot saved)
###          grid_1d_size_crust, grid_1d_size_mantle : floats
###          abd_set : str 'mid', 'low' or 'high'
###          title_prefix : str, to identify plot (e.g. std_osc_params
###          high_theta, const_P_ee, etc)
###
###          makes dir to save plot as pdf and data as csv
###          4 columns for csv : energy_arr, N_U, N_Th, N_U + N_Th
###
### NAMING SCHEME : [title_prefix]_spec_mid_100E20C50M
###                 spec + abundance set name + # E bins + E +
###                 crust spacing in km + C + mantle spacing
###                 in km + M

def plot_spec(N_Th, N_U, spec_save, plot_show, grid_1d_size_crust, grid_1d_size_mantle, abd_set, title_prefix=""):
    # Calculate total number of geonus
    # WARNING : if you don't use scaled spectra, this doesn't mean anything

    geonus_tot_Th = np.sum(N_Th)
    geonus_tot_U = np.sum(N_U)
    geonus_tot = geonus_tot_Th + geonus_tot_U

    # Plotting the data
    plt.step(energy_array, N_Th, where='mid', label='Th232', color='blue')
    plt.step(energy_array, N_U, where='mid', label='U238', color='red')
    plt.step(energy_array, N_U + N_Th, where='mid', label='total', color='green')

    # Add grid
    plt.grid(True)  # Show grid
    plt.grid(which = 'both', color='lightgray', linestyle='--', linewidth=0.5)

    # Adding geonus_tot label to the plot
    label_text = f'geonus_tot = {geonus_tot:.2f}\ngeonus_tot_Th = {geonus_tot_Th:.2f}\ngeonus_tot_U = {geonus_tot_U:.2f}'
    plt.text(0.95, 0.95, label_text, transform=plt.gca().transAxes, ha='right', va='top', fontsize=10,
             bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', boxstyle='round,pad=0.5'))

    plt.xlabel('E_nu [MeV]')
    plt.yscale('log')
    plt.ylabel('Expected geonu count')
    plt.title(f'Expected geonus \n {title_prefix}')

    plt.ylim(bottom=np.max(N_U + N_Th) / 20)
    plt.minorticks_on()
    plt.legend(loc='lower left')

    # Save plot if spec_save is True
    if spec_save:
        # Construct the directory path based on the naming scheme
        dir_name = f"{abd_set}_{len(energy_array)}E{int(np.floor(grid_1d_size_crust))}C{int(np.floor(grid_1d_size_mantle))}M"
        save_dir = os.path.join("..", "plots", dir_name)

        # Check if directory exists, create it if not
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        # Define the full file path with the title prefix
        file_path = os.path.join(save_dir, f"{title_prefix}_spec_{abd_set}_{len(energy_array)}E{int(np.floor(grid_1d_size_crust))}C{int(np.floor(grid_1d_size_mantle))}M.pdf")

        # Print the save location
        print(f"Saving plot in {file_path}")

        # Save the plot
        plt.savefig(file_path, format='pdf')

        # Save data to CSV
        csv_path = os.path.join(save_dir,
                                f"{title_prefix}_spec_data_{abd_set}_{len(energy_array)}E{int(np.floor(grid_1d_size_crust))}C{int(np.floor(grid_1d_size_mantle))}M.csv")

        with open(csv_path, mode='w', newline='') as csv_file:
            writer = csv.writer(csv_file)
            # Write header
            writer.writerow(['Energy [MeV]', 'N_Th', 'N_U', 'Total'])
            # Write data rows
            for e, th, u in zip(energy_array, N_Th, N_U):
                writer.writerow([e, th, u, th + u])

        # Print confirmation of CSV save location
        print(f"Data saved in {csv_path}")

    if plot_show:
        plt.show()

### Plot spectrums and ratio
###
### inputs : energy_array, N_Th, N_U : arays of same len
###          two copies of each N
###          spec_save : bool (set true if you want plot saved)
###          grid_1d_size_crust, grid_1d_size_mantle : floats
###          abd_set : str 'mid', 'low' or 'high'
###          two copies, one for each spectrum
###          title_prefix : str, to identify plot (e.g. std_osc_params
###          high_theta, const_P_ee, etc)
###          two copies, one for each spectrum
###
###          makes dir to save ratio plot as pdf and data as csv
###          2 columns for spectrum csvs: energy_arr, N_U + N_Th
###          4 columns for ratio csv : energy_arr, N_U_1 + N_Th_1,
###          N_U_2 + N_Th_2, total_2 / total_2
###
###
### NAMING SCHEME : spectrum plots
###                 [title_prefix]_spec_mid_100E20C50M
###                 spec + abundance set name + # E bins + E +
###                 crust spacing in km + C + mantle spacing
###                 in km + M
###
### TO DO : might need to change to accommodate spectra of
###         different 1d grid sizes (for mantle mostly)

def plot_rat(N_Th_1, N_U_1, N_Th_2, N_U_2, spec_save, plot_show, grid_1d_size_crust, grid_1d_size_mantle, abd_set_1, abd_set_2, title_prefix_1 ="", title_prefix_2 =""):
    # Calculate total number of geonus
    # WARNING : if you don't use scaled spectra, this doesn't mean anything

    geonus_tot_Th_1 = np.sum(N_Th_1)
    geonus_tot_U_1 = np.sum(N_U_1)
    geonus_tot_1 = geonus_tot_Th_1 + geonus_tot_U_1

    geonus_tot_Th_2 = np.sum(N_Th_2)
    geonus_tot_U_2 = np.sum(N_U_2)
    geonus_tot_2 = geonus_tot_Th_2 + geonus_tot_U_2

    # Plotting the data for spectra
    plt.step(energy_array, N_U_1 + N_Th_1, where='mid', label='total', color='green')

    # Add grid
    plt.grid(True)  # Show grid
    plt.grid(which = 'both', color='lightgray', linestyle='--', linewidth=0.5)

    # Adding geonus_tot label to the plot
    label_text = f'geonus_tot = {geonus_tot_1:.2f}\ngeonus_tot_Th = {geonus_tot_Th_1:.2f}\ngeonus_tot_U = {geonus_tot_U_1:.2f}'
    plt.text(0.95, 0.95, label_text, transform=plt.gca().transAxes, ha='right', va='top', fontsize=10,
             bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', boxstyle='round,pad=0.5'))

    plt.xlabel('E_nu [MeV]')
    plt.yscale('log')
    plt.ylabel('Expected geonu count')
    plt.title(f'Expected geonus \n {title_prefix_1}')

    #plt.ylim(bottom=6e-11)
    plt.ylim(bottom=np.max(N_U_1 + N_Th_1) / 20)
    plt.minorticks_on()
    plt.legend(loc='lower left')

    # Save plot if spec_save is True
    if spec_save:
        # Construct the directory path based on the naming scheme
        # the dir is named after the abd set, energy bins and grid
        # spacing of the _2 fluxes = reference
        #
        # in the main script, this will be the one with standard
        # oscillation params from James/Tony fit

        dir_name = f'{abd_set_2}_{len(energy_array)}E{int(np.floor(grid_1d_size_crust))}C{int(np.floor(grid_1d_size_mantle))}M'
        save_dir = os.path.join("..", "plots", dir_name)

        # Check if spectrum plot directory exists, create it if not
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        # Define the full file path with the title prefix
        file_path = os.path.join(save_dir,
                                 f"{title_prefix_1}_spec_{abd_set_1}_{len(energy_array)}E{int(np.floor(grid_1d_size_crust))}C{int(np.floor(grid_1d_size_mantle))}M.pdf")

        # Print the save location
        print(f"Saving plot in {file_path}")

        # Save the plot
        plt.savefig(file_path, format='pdf')

        # Save data to CSV
        csv_path = os.path.join(save_dir,
                                f"{title_prefix_1}_spec_data_{abd_set_1}_{len(energy_array)}E{int(np.floor(grid_1d_size_crust))}C{int(np.floor(grid_1d_size_mantle))}M.csv")

        with open(csv_path, mode='w', newline='') as csv_file:
            writer = csv.writer(csv_file)
            # Write header
            writer.writerow(['Energy [MeV]', 'Total'])
            # Write data rows
            for e, th, u in zip(energy_array, N_Th_1, N_U_1):
                writer.writerow([e, th + u])

        # Print confirmation of CSV save location
        print(f"Data saved in {csv_path}")

    if plot_show:
        plt.show()

    # second spectrum

    plt.step(energy_array, N_U_2 + N_Th_2, where='mid', label='total', color='green')

    # Add grid
    plt.grid(True)  # Show grid
    plt.grid(which = 'both', color='lightgray', linestyle='--', linewidth=0.5)  # Customize grid appearance

    # Adding geonus_tot label to the plot
    label_text = f'geonus_tot = {geonus_tot_2:.2f}\ngeonus_tot_Th = {geonus_tot_Th_2:.2f}\ngeonus_tot_U = {geonus_tot_U_2:.2f}'
    plt.text(0.95, 0.95, label_text, transform=plt.gca().transAxes, ha='right', va='top', fontsize=10,
             bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', boxstyle='round,pad=0.5'))

    plt.xlabel('E_nu [MeV]')
    plt.yscale('log')
    plt.ylabel('Expected geonu count')
    plt.title(f'Expected geonus \n {title_prefix_2}')

    #plt.ylim(bottom=6e-11)
    plt.ylim(bottom=np.max(N_U_2 + N_Th_2) / 20)
    plt.minorticks_on()
    plt.legend(loc='lower left')

    # Save plot if spec_save is True
    if spec_save:
        # Construct the directory path based on the naming scheme
        dir_name = f'{abd_set_2}_{len(energy_array)}E{int(np.floor(grid_1d_size_crust))}C{int(np.floor(grid_1d_size_mantle))}M'
        save_dir = os.path.join("..", "plots", dir_name)

        # Check if spectrum plot directory exists, create it if not
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        # Define the full file path with the title prefix
        file_path = os.path.join(save_dir,
                                 f"{title_prefix_2}_spec_{abd_set_2}_{len(energy_array)}E{int(np.floor(grid_1d_size_crust))}C{int(np.floor(grid_1d_size_mantle))}M.pdf")

        # Print the save location
        print(f"Saving plot in {file_path}")

        # Save the plot
        plt.savefig(file_path, format='pdf')

        # Save data to CSV
        csv_path = os.path.join(save_dir,
                                f"{title_prefix_2}_spec_data_{abd_set_2}_{len(energy_array)}E{int(np.floor(grid_1d_size_crust))}C{int(np.floor(grid_1d_size_mantle))}M.csv")

        with open(csv_path, mode='w', newline='') as csv_file:
            writer = csv.writer(csv_file)
            # Write header
            writer.writerow(['Energy [MeV]', 'Total'])
            # Write data rows
            for e, th, u in zip(energy_array, N_Th_2, N_U_2):
                writer.writerow([e, th + u])

        # Print confirmation of CSV save location
        print(f"Data saved in {csv_path}")

    if plot_show:
        plt.show()

    # Ratio plot


    plt.step(energy_array, (N_U_1 + N_Th_1) / (N_U_2 + N_Th_2), where='mid', label='total', color='green')

    # Add grid
    plt.grid(True)  # Show grid
    plt.grid(which = 'both', color='lightgray', linestyle='--', linewidth=0.5)  # Customize grid appearance

    plt.xlabel('E_nu [MeV]')
    plt.ylabel('Expected geonu count ratio')
    plt.title(f'Ratio of expected geonus, \n {title_prefix_1}_{abd_set_1} / {title_prefix_2}_{abd_set_2}')

    #plt.ylim(bottom=0.92, top=1.15)
    plt.minorticks_on()
    plt.legend(loc='lower left')

    # Save plot if spec_save is True
    if spec_save:
        # Construct the directory path based on the naming scheme
        dir_name = f'{abd_set_2}_{len(energy_array)}E{int(np.floor(grid_1d_size_crust))}C{int(np.floor(grid_1d_size_mantle))}M'
        save_dir = os.path.join("..", "plots", dir_name)

        # Check if spectrum plot directory exists, create it if not
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        # Define the full file path with the title prefix
        file_path = os.path.join(save_dir,
                                 f"ratio_{title_prefix_1}_{abd_set_1}_{title_prefix_2}_{abd_set_2}_{len(energy_array)}E{int(np.floor(grid_1d_size_crust))}C{int(np.floor(grid_1d_size_mantle))}M.pdf")

        # Print the save location
        print(f"Saving plot in {file_path}")

        # Save the plot
        plt.savefig(file_path, format='pdf')

        # Save data to CSV
        csv_path = os.path.join(save_dir,
                                f"{title_prefix_1}_{abd_set_1}_{title_prefix_2}_{abd_set_2}_{len(energy_array)}E{int(np.floor(grid_1d_size_crust))}C{int(np.floor(grid_1d_size_mantle))}M.csv")

        with open(csv_path, mode='w', newline='') as csv_file:
            writer = csv.writer(csv_file)
            # Write header
            writer.writerow(['Energy [MeV]', 'Total'])
            # Write data rows
            for e, th1, u1, th2, u2 in zip(energy_array, N_Th_1, N_U_1, N_Th_2, N_U_2):
                writer.writerow([e, th1 + u1, th2 + u2, (th1 + u1) / (th2 + u2)])

        # Print confirmation of CSV save location
        print(f"Data saved in {csv_path}")

    if plot_show:
        plt.show()


    ### Plotting comparison with constant P_ee with error
    ### bands (very hand wavy)
    ###
    ### vol_integral_const_P_ee() returns sum_mid_Th,
    ### sum_mid_U, sum_low_Th, sum_low_U, sum_high_Th,
    ### sum_high_U; the N_U and N_Th for each are computed
    ### with calc_exp_spec
    ###
    ### TO DO (not super useful and quite complicated)

##########################################################
### Calculate average P_ee
###
### This section includes a bunch of additional functions
### estimating the average P_ee to be used as a constant
### for calculating the expected flux

### Set arrays of energy weights for each emitter
### See explanation doc for details - these are the f(E)'s
### (Sec. 7.5)
### Each element of an array corresponds to one energy
###
### No inputs, but requires some globally set variables

def set_energy_weights():
    f_E_Th = ((lambda_Th) / (mu_Th)) * sigma * dn_dE_rebinned_Th
    f_E_U = ((lambda_U) / (mu_U)) * sigma * dn_dE_rebinned_U

    print("Shape of f(E)_Th:", f_E_Th.shape)
    print("Shape of f(E)_U:", f_E_U.shape)

    return f_E_Th, f_E_U


### Calculate the position weights; these are the w(r)'s
### from the explanation doc (Sec 7.5)


def get_position_weights(points_array, grid_1d_size,  A_Th, A_U, rho):
    dV = grid_1d_size ** 3

    relative_distance_array = calc_relative_dist(points_array)
    print("Relative distance array computed successfully")

    w_pos_Th = (A_Th * rho * dV) / (4 * np.pi * relative_distance_array ** 2)
    w_pos_U = (A_U * rho * dV) / (4 * np.pi * relative_distance_array ** 2)

    return w_pos_U, w_pos_Th


### TO DO : go through all functions and delete local vars
### when testing code, make sure you check with memory print
### statements if it actually did anything

### To calculate average prob, must calculate terms for
### numerator and denominator separately for each emitter
### and for each layer;
###
### The dE is always the same so it simplifies, but we must
### keep volume elements in explicitly because they are not
### the same in all layers

### One term for the numerator is :
###   f(E) * w(r) * P(E, r) * dV
### One term for the denimunator is :
###   f(E) * w(r) * dV

### Calculate numerator for layer made of sublayers
### Calculate denominator for layer made of sublayers
###
### Keeping both of these in the same function allows us
### to only calculate the position weights once!
###
### inputs : grids (each grids[i] is a 2D array where each elem
###          is (x, y, z) with the coords of each point
###          energy_array
###          grid_1d_size : distance betweem two consecutive
###          grid points in 1d; used to compute volume element
###          theta_12, delta_m_21_squared : default values given (best fit
###          values from PDG (I think 2020)
###          A_Th, A_U, rho : abundances and density in layer

def calc_avg_P_ee_num_den(grids, theta_12, delta_m_21_squared, grid_1d_size, A_Th, A_U, rho, f_E_Th, f_E_U):

    print(' ')
    print('computing numerators')

    numerator_Th_total = 0
    numerator_U_total = 0

    for i in range(len(grids)):
        get_memory_usage()
        P_ee = calc_P_ee(grids[i], theta_12, delta_m_21_squared)
        print(f'average P_ee for layer = {np.mean(P_ee)}')
        w_pos_Th, w_pos_U = get_position_weights(grids[i], grid_1d_size, A_Th, A_U, rho)
        w_pos_Th = np.nan_to_num(w_pos_Th, nan=0)
        w_pos_U = np.nan_to_num(w_pos_U, nan=0)
        print(f'computing numerators for layer {i}')
        numerator_Th = np.nansum(w_pos_Th[:, np.newaxis] * f_E_Th * P_ee.T)
        numerator_U = np.nansum(w_pos_U[:, np.newaxis] * f_E_U * P_ee.T)
        print(f'numerator values computed for layer {i} out of {len(grids)}')
        get_memory_usage()
        print('deleting P_ee and pos weights for layer')
        del P_ee
        del w_pos_Th
        del w_pos_U
        print('deleted')
        get_memory_usage()

        print(f'numerator_Th_total : {numerator_Th_total}')
        print(f'numerator_U_total : {numerator_U_total}')
        numerator_Th_total += numerator_Th
        numerator_U_total += numerator_U
        print(f'numerator_Th_total : {numerator_Th_total}')
        print(f'numerator_U_total : {numerator_U_total}')


    print('numerator values for the entire layer computed')
    print('computing denominator values for each sublayer')

    denominator_Th_total = 0
    denominator_U_total = 0

    for i in range(len(grids)):
        w_pos_Th, w_pos_U = get_position_weights(grids[i], grid_1d_size, A_Th, A_U, rho)
        w_pos_Th = np.nan_to_num(w_pos_Th, nan=0)
        w_pos_U = np.nan_to_num(w_pos_U, nan=0)

        denominator_Th = np.nansum(w_pos_Th[:, np.newaxis] * f_E_Th)
        denominator_U = np.nansum(w_pos_U[:, np.newaxis] * f_E_U)
        print(f'denominator values computed for layer {i} out of {len(grids)}')
        print(' ')

        get_memory_usage()
        print('deleting pos weights for layer')
        del w_pos_Th
        del w_pos_U
        print('deleted')
        get_memory_usage()

        print(f'denominator_Th_total : {denominator_Th_total}')
        print(f'denominator_U_total : {denominator_U_total}')
        denominator_Th_total += denominator_Th
        denominator_U_total += denominator_U
        print(f'denominator_Th_total : {denominator_Th_total}')
        print(f'denominator_U_total : {denominator_U_total}')


    print('denominator values for the entire layer computed')


    print('deleting more stuff')
    del numerator_U
    del numerator_Th
    print('deleted')
    get_memory_usage()


    return numerator_Th_total, numerator_U_total, denominator_Th_total, denominator_U_total


### Calculate numerator and denominator terms for simple layer

def calc_avg_P_ee_num_den_simple(energy_array, points_array, theta_12, theta_13, delta_m_21_squared, grid_1d_size, A_Th, A_U, rho, f_E_Th, f_E_U):
    print('calculating survival probabilities for layer')

    P_ee = calc_P_ee(points_array, theta_12, delta_m_21_squared)
    print('survival probability computed')
    print('calculating position weights')

    w_pos_Th, w_pos_U = get_position_weights(points_array, grid_1d_size, A_Th, A_U, rho)

    print('position weights computed')
    print(' ')
    print('computing numerators')

    numerator_Th = np.nansum(w_pos_Th[:, np.newaxis] * f_E_Th * P_ee.T)
    numerator_U = np.nansum(w_pos_U[:, np.newaxis] * f_E_U * P_ee.T)
    print(f'numerator values computed for layer')
    print(f'Shape of numerator_Th_c : ', numerator_Th.shape)
    print(f'Shape of numerator_U_c : ', numerator_U.shape)
    print()

    print('computing denominator values')
    denominator_Th = np.nansum(w_pos_Th[:, np.newaxis] * f_E_Th)
    denominator_U = np.nansum(w_pos_U[:, np.newaxis] * f_E_U)
    print(f'denominator values computed for layer')
    print(f'Shape of denominator_Th_c : ', denominator_Th.shape)
    print(f'Shape of denominator_U_c : ', denominator_U.shape)

    return numerator_Th, numerator_U, denominator_Th, denominator_U

### TO DO : keeping layers with sublayers in an array of arrays
###         is not great, should use list instead given that
###         we just iterate over each sublayer and don't use
###         vectorized operations
### TO DO : after fixing that, should merge the two functions
###         above into one - check if it's a list, then if it
###         is do the loop, otherwise idk wrap array into a list
###         with one elem and keep the loop anyway


### Calculate average P_ee with stdev for Layered Earth Model
### Important : this function assumes that we are working with
### an Earth model with 4 layers : C (crust), CLM, DM, EM
### and that the numerator and denominator terms associated
### with each of these and each emitter has been computed

def calc_P_ee_avg(numerator_Th_c, numerator_U_c, numerator_Th_CLM, numerator_U_CLM, numerator_Th_DM, numerator_U_DM, numerator_Th_EM, numerator_U_EM, denominator_Th_c, denominator_U_c, denominator_Th_CLM, denominator_U_CLM, denominator_Th_DM, denominator_U_DM, denominator_Th_EM, denominator_U_EM):
    print('calculating average P_ee for given Earth model')
    P_ee_average = ((numerator_Th_c + numerator_U_c +
                    numerator_Th_CLM + numerator_U_CLM +
                    numerator_Th_DM + numerator_U_DM +
                    numerator_Th_EM + numerator_U_EM) /
                    (denominator_Th_c + denominator_U_c +
                     denominator_Th_CLM + denominator_U_CLM +
                     denominator_Th_DM + denominator_U_DM +
                     denominator_Th_EM + denominator_U_EM))

    print(f'Average P_ee computed: P_ee = {P_ee_average}')
    return P_ee_average

### To calculate the variance, the numerators are different
### and need the computed average value

def calc_var_P_ee_num(grids, theta_12, delta_m_21_squared, grid_1d_size, A_Th, A_U, rho, f_E_Th, f_E_U, P_ee_average):
    print('calculating numerator term for current layer')

    numerator_Th_total = 0
    numerator_U_total = 0

    # Process each sublayer individually
    for i, grid in enumerate(grids):
        print(f'Processing layer {i + 1} of {len(grids)}')

        # Calculate survival probabilities
        P_ee = calc_P_ee(grid, theta_12, delta_m_21_squared)
        print(f'average P_ee for layer = {np.mean(P_ee)}')

        # Compute position weights
        w_pos_Th, w_pos_U = get_position_weights(grid, grid_1d_size, A_Th, A_U, rho)
        w_pos_Th = np.nan_to_num(w_pos_Th, nan=0)
        w_pos_U = np.nan_to_num(w_pos_U, nan=0)

        # Compute numerator terms for the current layer
        numerator_Th = np.nansum(w_pos_Th[:, np.newaxis] * f_E_Th * (P_ee.T - P_ee_average) ** 2)
        numerator_U = np.nansum(w_pos_U[:, np.newaxis] * f_E_U * (P_ee.T - P_ee_average) ** 2)
        get_memory_usage()
        print('computed numerators for layer; deleting P_ee and position weights')
        del P_ee
        del w_pos_Th
        del w_pos_U
        print('deleted')
        get_memory_usage()
        # Accumulate totals
        numerator_Th_total += numerator_Th
        numerator_U_total += numerator_U

        print(f'Layer {i + 1}: Th Total = {numerator_Th_total}, U Total = {numerator_U_total}')

    print('Completed numerator calculations for all layers')
    return numerator_Th_total, numerator_U_total

def calc_var_P_ee_num_v2(grids, theta_12, delta_m_21_squared, grid_1d_size, A_Th, A_U, rho, f_E_Th, f_E_U, P_ee_average):
    print('calculating numerator term for current layer')

    numerator_Th_total = 0
    numerator_U_total = 0

    # Process each sublayer individually
    for i in range(len(grids)):
        print(f'Processing layer {i + 1} of {len(grids)}')

        # Calculate survival probabilities
        P_ee = calc_P_ee(grids[i], theta_12, delta_m_21_squared)
        print(f'average P_ee for layer = {np.mean(P_ee)}')

        # Compute position weights
        w_pos_Th, w_pos_U = get_position_weights(grids[i], grid_1d_size, A_Th, A_U, rho)
        w_pos_Th = np.nan_to_num(w_pos_Th, nan=0)
        w_pos_U = np.nan_to_num(w_pos_U, nan=0)

        # Compute numerator terms for the current layer
        numerator_Th = np.nansum(w_pos_Th[:, np.newaxis] * f_E_Th * (P_ee.T - P_ee_average) ** 2)
        numerator_U = np.nansum(w_pos_U[:, np.newaxis] * f_E_U * (P_ee.T - P_ee_average) ** 2)
        get_memory_usage()
        print('computed numerators for layer; deleting P_ee and position weights')
        del P_ee
        del w_pos_Th
        del w_pos_U
        print('deleted')
        get_memory_usage()
        # Accumulate totals
        numerator_Th_total = numerator_Th_total + numerator_Th
        numerator_U_total += numerator_U

        print(f'Layer {i + 1}: Th Total = {numerator_Th_total}, U Total = {numerator_U_total}')

    print('Completed numerator calculations for all layers')
    return numerator_Th_total, numerator_U_total

### Calc numerators for var for single layer

def calc_avg_P_ee_num_den_simple(points_array, theta_12, delta_m_21_squared, grid_1d_size, A_Th, A_U, rho, f_E_Th, f_E_U, P_ee_average):
    print('calculating survival probabilities for layer')

    P_ee = calc_P_ee(points_array, theta_12, delta_m_21_squared)
    print('survival probability computed')
    print('calculating position weights')

    w_pos_Th, w_pos_U = get_position_weights(points_array, grid_1d_size, A_Th, A_U, rho)

    print('position weights computed')
    print(' ')
    print('computing numerators')

    numerator_var_Th = np.nansum(w_pos_Th[:, np.newaxis] * f_E_Th * (P_ee.T - P_ee_average))
    numerator_var_U = np.nansum(w_pos_U[:, np.newaxis] * f_E_U * (P_ee.T - P_ee_average))

    print(f'numerator values computed for layer')
    print(f'Shape of numerator_var_Th_c : ', numerator_var_Th.shape)
    print(f'Shape of numerator_var_U_c : ', numerator_var_U.shape)

    return numerator_var_Th, numerator_var_U

### Calculate stdev for aerage P_ee given Earth model

def calc_P_ee_stdev(numerator_var_Th_c, numerator_var_U_c, numerator_var_Th_CLM, numerator_var_U_CLM, numerator_var_Th_DM, numerator_var_U_DM,
                  numerator_var_Th_EM, numerator_var_U_EM, denominator_Th_c, denominator_U_c, denominator_Th_CLM,
                  denominator_U_CLM, denominator_Th_DM, denominator_U_DM, denominator_Th_EM, denominator_U_EM):
    print('calculating stdev P_ee for given Earth model')
    P_ee_var = (numerator_var_Th_c + numerator_var_U_c + numerator_var_Th_CLM + numerator_var_U_CLM + numerator_var_Th_DM + numerator_var_U_DM + numerator_var_Th_EM + numerator_var_U_EM) / (
                               denominator_Th_c + denominator_U_c + denominator_Th_CLM + denominator_U_CLM + denominator_Th_DM + denominator_U_DM + denominator_Th_EM + denominator_U_EM)

    print(f'P_ee variance computed = {P_ee_var}')
    P_ee_stdev = np.sqrt(P_ee_var)
    print(f'P_ee_stdev computed = {P_ee_stdev}')

    return P_ee_stdev

def days_to_seconds(days):
    seconds_per_day = 86400  # 24 * 60 * 60
    return days * seconds_per_day

### Given 2 csv files with spectra, extract spectra info, plot and save ratio

def extract_columns(csv_title, path=None):
    # Extract the folder title from the file title
    try:
        folder_name = '_'.join(csv_title.split('_')[1:]).split('.')[0]
        if path is None:
            print(f'Looking for file in folder: ../plots/{folder_name}')
    except IndexError:
        raise ValueError(
            "Invalid file title format. "
            "Ensure it's in the format 'blahblah_abd_eEcCmM.csv', "
            "where abd is the abundance set name (e.g., mid, low, high), and e, c, m "
            "represent the number of energy bins, crust spacing, and mantle spacing."
        )

    # Construct the file path if not provided
    if path is None:
        file_path = os.path.join('..', 'plots', folder_name, csv_title)
    else:
        file_path = os.path.join(path, csv_title)

    # Check if the file exists
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")

    # Read the CSV file
    try:
        data = pd.read_csv(file_path)
    except Exception as e:
        raise ValueError(f"Error reading the CSV file: {e}")

    # Extract the first and last columns
    try:
        energy_array_ext = np.array(data.iloc[:, 0])
        N_tot_ext = np.array(data.iloc[:, -1])
    except IndexError:
        raise ValueError("The CSV file does not have the expected structure.")

    return energy_array_ext, N_tot_ext
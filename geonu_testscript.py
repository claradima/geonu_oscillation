############################################################################
# This file mostly copies what's in the latest notebook so it can be ran 
# from a terminal; I mostly want to use it to run this code on hpc and see 
# I can make it work;
#
# I want to submit a batch job
#
# I want to add saving statements for the plots; I will put them in a folder
# that I'll maybe add to gitignore
#############################################################################

import numpy as np
import matplotlib.pyplot as plt
import math
import os #need this to save the plots

# make earth grid
# adjust grid size here -- bigger grid_count, finer grid


grid_count = 250
coordinates = np.linspace(-6371, 6371, grid_count)
x_coordinates, y_coordinates, z_coordinates = np.meshgrid(coordinates, coordinates, coordinates)

earth_grid = np.zeros((grid_count**3, 3))
index = 0

for i in range(grid_count):
    for j in range(grid_count):
        for k in range(grid_count):
            if x_coordinates[i, j, k]**2 + y_coordinates[i, j, k]**2 + z_coordinates[i, j, k]**2 <= 6371**2:
                earth_grid[index] = np.array([x_coordinates[i, j, k], y_coordinates[i, j, k], z_coordinates[i, j, k]])
                index += 1

earth_grid = earth_grid[:index]

#separate mantle from core

core_radius = 3486
mantle_crust_mask = (earth_grid[:, 0]**2 + earth_grid[:, 1]**2 + earth_grid[:, 2]**2 >= core_radius**2)

mantle_crust_grid = earth_grid[mantle_crust_mask]
core_grid = earth_grid[~mantle_crust_mask]

# reshape the arrays to have three columns
mantle_crust_grid = mantle_crust_grid.reshape(-1, 3)
core_grid = core_grid.reshape(-1, 3)

#set abundances and densities; weight average density by shell volume

rho_core = 0
A_Th_core = 0
A_U_core = 0

rho_mantle = 4.5 #g/cm^3 -- I just looked this up
#I'm not sure what the units are supposed to be but this will just be a scaling factor, we mostly care about the shape of 
#the spectrum, not the absolute value 

A_Th_c = 5 * (10**(-6))
A_Th_m = 8 * (10**(-8))
A_U_c = 1 * (10**(-6))
A_U_m = 2 * (10**(-8))

r_outer_crust = 6371
r_inner_crust = 6371 - 20
r_outer_mantle = r_inner_crust
r_inner_mantle = r_outer_mantle - 2900

def shell_volume(r_inner, r_outer):
    return (4/3)*math.pi*(r_outer**2 + r_outer*r_inner + r_inner**2)*(r_outer - r_inner)

crust_volume = shell_volume(r_inner_crust, r_outer_crust)
mantle_volume = shell_volume(r_inner_mantle, r_outer_mantle)

#print("Crust volume is "+str(crust_volume))
#print("Mantle volume is "+str(mantle_volume))

A_Th_mantle = (A_Th_m*mantle_volume + A_Th_c*crust_volume)/(crust_volume + mantle_volume)
A_U_mantle = (A_U_m*mantle_volume + A_U_c*crust_volume)/(crust_volume + mantle_volume)

#set position of SNO+

SNO_r = np.array([0, 0, 6369])

# function to calculate relative distance to SNO+ from points in earth grid

def calc_relative_dist(points_array):
    relative_distances = np.zeros(len(points_array))
    for i in range(len(points_array)):
        # Index each point individually
        point = points_array[i]
        relative_distances[i] = np.sqrt((point[0] - SNO_r[0])**2 + (point[1] - SNO_r[1])**2 + (point[2] - SNO_r[2])**2)

    return relative_distances
#distances calculated in km

#make abundances and density into arrays; not super useful now but might be easier if you have different regions

A_Th_array = np.full(len(mantle_crust_grid), A_Th_mantle)
A_U_array = np.full(len(mantle_crust_grid), A_U_mantle)

rho_array = np.full(len(mantle_crust_grid), rho_mantle*1000) # x1000 to convert from g/cm^3 

#define Delta function for an array of points and energies
#inputs: appropriate deltam _ij^2, energy array [MeV], points array
#relative distance calculated in km, hence Delta scaled up by a factor of 1000 to account for m-km conversion


def Delta_ij(energy_array, points_array, delta_m_ij_squared):
    
    Delta = np.zeros((len(energy_array), len(points_array)))
    #want positions relative to sno+
    relative_distance_array = calc_relative_dist(points_array)
    
    for i in range(len(energy_array)):
        for j in range(len(points_array)):
            Delta[i][j] = (1.27 * delta_m_ij_squared * relative_distance_array[j] * 1000) / (energy_array[i])
            
    return Delta
    
#define FULL survival probability formula for arrays of points and energies

def P_ee_full(energy_array, points_array, theta_12, theta_13, delta_m_21_squared):
    
    P_ee = np.zeros((len(energy_array), len(points_array)))
    #want positions relative to sno+
    relative_distance_array = calc_relative_dist(points_array)
    
    Delta_31 = Delta_ij(energy_array, points_array, delta_m_31_squared)
    #Delta_32 = Delta_ij(energy_array, points_array, delta_m_32_squared)
    #the two things above are equal but keep separate for clarity/flexibility
    #uncomment line if you need
    Delta_32 = Delta_31
    Delta_12 = Delta_ij(energy_array, points_array, delta_m_21_squared)
    
    for i in range(len(energy_array)):
        for j in range(len(points_array)):
            P_ee[i][j] = 1 - ( (np.cos(theta_13))**4 * (np.sin(2 * theta_12))**2 * (np.sin(Delta_12[i][j]))**2 + (np.cos(2*theta_12) **2 * np.sin(Delta_31[i][j])**2  + np.sin(2*theta_12) **2 * np.sin(Delta_32[i][j])**2)* (np.sin(2 * theta_13))**2)

    #because the two deltas with 3 are equal here, all this does is replace
    #the factor of 1/2 in front of sin^2(2theta_13) with a factor of
    #sin^2(delta_31)

    
    
    return P_ee
    
#define simplified survival probability for arrays of points and energies

def P_ee_simplified(energy_array, points_array, theta_12, theta_13, delta_m_21_squared):
    
    P_ee = np.zeros((len(energy_array), len(points_array)))
    #want positions relative to sno+
    relative_distance_array = calc_relative_dist(points_array)
    
        
    Delta_12 = Delta_ij(energy_array, points_array, delta_m_21_squared)
    
    for i in range(len(energy_array)):
        for j in range(len(points_array)):
            P_ee[i][j] = 1 - ( (np.cos(theta_13))**4 * (np.sin(2 * theta_12))**2 * (np.sin(Delta_12[i][j]))**2 + (1/2) * (np.sin(2 * theta_13))**2)
            
    return P_ee
    
#set standard oscillation parameters

theta_12 = 0.5903 #rad
theta_23 = 0.8430 #rad
theta_13 = 0.1503 #rad

delta_m_21_squared = 7.39 * 10**(-5) #eV^2
delta_m_32_squared = 2.449 * 10**(-3) #eV^2
delta_m_31_squared = delta_m_32_squared

#define energy array

energy_array = np.linspace(1.8, 3.3, 100)

#define Thorium integral using approximate formula in terms of points and energies

def integral_over_positions_Th(points_array, energy_array, grid_1d_size, theta_12, delta_m_21_squared, A_Th_array, rho_array):
    sum = np.zeros(len(energy_array))
    dV = grid_1d_size**3
    
    relative_distance_array = calc_relative_dist(points_array)
    P_ee_array = P_ee_full(energy_array, points_array, theta_12, theta_13, delta_m_21_squared)
    
    for j in range(len(energy_array)):
        for i in range(len(points_array)):
            sum[j] += P_ee_array[j][i] * ((A_Th_array[i] * rho_array[i]) / (4 * np.pi * (relative_distance_array[i]**2))) * dV
            #print("computed for energy "+str(j)+" out of "+str(len(energy_array)))
            
    return sum
    
#define Uranium integral using approximate formula in terms of points and energies

def integral_over_positions_U(points_array, energy_array, grid_1d_size, theta_12, delta_m_21_squared, A_U_array, rho_array):
    sum = np.zeros(len(energy_array))
    dV = grid_1d_size**3
    
    relative_distance_array = calc_relative_dist(points_array)
    P_ee_array = P_ee_full(energy_array, points_array, theta_12, theta_13, delta_m_21_squared)
    
    for j in range(len(energy_array)):
        for i in range(len(points_array)):
            sum[j] += P_ee_array[j][i] * ((A_U_array[i] * rho_array[i]) / (4 * np.pi * (relative_distance_array[i]**2))) * dV
            
    return sum #this is an array of length = len(energy_array)
    #can identify integral for each energy by index
    
def integral_over_positions_Th_U(points_array, energy_array, grid_1d_size, theta_12, delta_m_21_squared, A_Th_array, A_U_array, rho_array):
    sum_Th = np.zeros(len(energy_array))
    sum_U = np.zeros(len(energy_array))
    dV = grid_1d_size**3
    
    relative_distance_array = calc_relative_dist(points_array)
    P_ee_array = P_ee_full(energy_array, points_array, theta_12, theta_13, delta_m_21_squared)

    for j in range(len(energy_array)):
        for i in range(len(points_array)):
            sum_Th[j] += P_ee_array[j][i] * ((A_Th_array[i] * rho_array[i]) / (4 * np.pi * (relative_distance_array[i]**2))) * dV
            #print("computed for energy "+str(j)+" out of "+str(len(energy_array)))
       
    for j in range(len(energy_array)):
        for i in range(len(points_array)):
            sum_U[j] += P_ee_array[j][i] * ((A_U_array[i] * rho_array[i]) / (4 * np.pi * (relative_distance_array[i]**2))) * dV
          

    return sum_Th, sum_U

#integral naming schemes

#Th_integral_values for full formula, standard oscillation parameters
#Th_integral_values_low_theta for some lower theta, full formula; same for mid and high
#Th_integral_values_constant_P_ee for integrating using constant survival probability

coordinates = np.linspace(-6371, 6371, grid_count)
grid_1d_size = coordinates[1] - coordinates[0]
print(grid_1d_size)  

#calculate integral arrays

Th_integral_values, U_integral_values = integral_over_positions_Th_U(mantle_crust_grid, energy_array, grid_1d_size, theta_12, delta_m_21_squared, A_Th_array, A_U_array, rho_array)

#make plots

#Thorium
plt.plot(energy_array, Th_integral_values, marker='o', linestyle='-', color='b')
plt.xlabel('Energy [MeV]')
plt.ylabel('Integral value')
plt.title('Integral over emission position, Thorium')
# Create the folder if it doesn't exist
if not os.path.exists('plots'):
    os.makedirs('plots')

# Save the plot with the title as the filename
plt.savefig('plots/Th_Int_standard_grid_' + str(grid_count) + '.png')

#Uranium
plt.plot(energy_array, U_integral_values, marker='o', linestyle='-', color='r')
plt.xlabel('Energy [MeV]')
plt.ylabel('Integral value')
plt.title('Integral over emission position, Uranium')
plt.savefig('plots/U_Int_standard_grid' + str(grid_count) + '.png')

#Both together

plt.plot(energy_array, Th_integral_values, marker='o', linestyle='-', color='b', label='Thorium')
plt.plot(energy_array, U_integral_values, marker='o', linestyle='-', color='r', label='Uranium')
plt.xlabel('Energy [MeV]')
plt.ylabel('Integral value')
plt.title('Integral over emission position')
plt.legend()
plt.savefig('plots/U_Th_Int_standard_grid' + str(grid_count) + '.png')


### IMPORTANT !!! ###
### complete this notebook after testing to see if it runs on apollo 2

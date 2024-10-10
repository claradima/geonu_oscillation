#!/usr/bin/env python
# coding: utf-8

# Here, I'm trying to figure out why I don't get the same result when splitting into sublayers and summing over that. I will just make one layer, and sum over its elements then split into only 2 sublayers. Let's see
# 
# Based on crust but increase a bit more

# In[84]:


import numpy as np
import matplotlib.pyplot as plt
import sys
print("imports successful")

grid_count_crust = 640 # would need 637.1 for 20 km spacing
#1275 for 5 km spacing
coordinates = np.linspace(-6371, 6371, grid_count_crust)
grid_1d_size_crust = coordinates[1]-coordinates[0]

# Generate the grid coordinates using meshgrid
x_coordinates, y_coordinates, z_coordinates = np.meshgrid(coordinates, coordinates, coordinates)

# Calculate the distance of each point from the origin using vectorized operations
distances_squared = x_coordinates**2 + y_coordinates**2 + z_coordinates**2

# Find indices where distance is less than or equal to the radius squared
crust_indices = np.logical_and(6300**2 < distances_squared, distances_squared <= 6371**2)

# Extract valid coordinates using boolean indexing
crust_grid = np.stack((x_coordinates[crust_indices], y_coordinates[crust_indices], z_coordinates[crust_indices]), axis=-1)

print("crust grid created successfully")


# In[85]:


# now make smaller grids
# Find indices where distance is less than or equal to the radius squared
crust_indices_1 = np.logical_and(6300**2 < distances_squared, distances_squared <= 6335**2)

# Extract valid coordinates using boolean indexing
crust_grid_1 = np.stack((x_coordinates[crust_indices_1], y_coordinates[crust_indices_1], z_coordinates[crust_indices_1]), axis=-1)

print("crust grid 1 created successfully")

# Find indices where distance is less than or equal to the radius squared
crust_indices_2 = np.logical_and(6335**2 < distances_squared, distances_squared <= 6371**2)

# Extract valid coordinates using boolean indexing
crust_grid_2 = np.stack((x_coordinates[crust_indices_2], y_coordinates[crust_indices_2], z_coordinates[crust_indices_2]), axis=-1)

print("crust grid 2 created successfully")


# In[86]:


print(len(crust_grid))
print(len(crust_grid_1) + len(crust_grid_2))

print(len(crust_grid_1))
print(len(crust_grid_2))


# In[ ]:





# In[87]:


# Assuming crust_grid_1 and crust_grid_2 are 2D arrays, where each row is a 3D coordinate
combined_grid = np.concatenate((crust_grid_1, crust_grid_2), axis=0)

# Find the unique rows (3D coordinates)
union_array = np.unique(combined_grid, axis=0)

# Print the length of the union array
print(len(union_array))


# In[88]:


# Ensure both arrays have unique rows
unique_union = np.unique(union_array, axis=0)
unique_crust_grid = np.unique(crust_grid, axis=0)

# Sort both arrays to ignore order
sorted_union = np.sort(unique_union, axis=0)
sorted_crust_grid = np.sort(unique_crust_grid, axis=0)

# Check if they are equal
are_equal = np.array_equal(sorted_union, sorted_crust_grid)

if are_equal:
    print("The arrays contain the same elements.")
else:
    print("The arrays do NOT contain the same elements.")
    


# In[89]:


#So then what the fuck is the issue???? Must be a later step .....


# In[90]:


print("setting densities and abundances")
# Define the parameter sets

# Placeholder definitions for IDE recognition
A_Th_c = 5 * (10**(-6)) 
A_U_c = 1 * (10**(-6))

rho_c = 2.7
print("densities and abundances set successfully")

#we have pretty much everything to compute things separately!
# compute for crust
#set position of SNO+

print("setting detector position")

SNO_r = np.array([0, 0, 6369])


# In[91]:


def calc_relative_dist(points_array):
    # Calculate the Euclidean distance using vectorized operations
    relative_distances = np.linalg.norm(points_array - SNO_r, axis=1)
    
    print("   ")
    print("Computed relative distances from Earth grid points to SNO+")
    print("   ")
    
    return relative_distances


# In[92]:


def Delta_ij(energy_array, points_array, delta_m_ij_squared):
    # Calculate relative distances
    relative_distance_array = calc_relative_dist(points_array)
    
    # Reshape energy_array to perform element-wise division
    energy_array_reshaped = energy_array.reshape(-1, 1)
    
    # Calculate Delta using vectorized operations
    Delta = (1.27 * delta_m_ij_squared * relative_distance_array * 1000) / (energy_array_reshaped)
    
    return Delta


# In[93]:


def P_ee_full_optimized(energy_array, points_array, theta_12, theta_13, delta_m_21_squared):
    
    #P_ee = np.empty((len(energy_array), len(points_array))) # np.empty more efficient
    #than np.zeros or whatever
    #might not need this at all though
    relative_distance_array = calc_relative_dist(points_array)
    
    Delta_31 = Delta_ij(energy_array, points_array, delta_m_31_squared)
    print("Delta_31 computed successfully")
    Delta_12 = Delta_ij(energy_array, points_array, delta_m_21_squared)
    print("Delta_12 computed successfully")

    A = (np.cos(theta_13))**4 * (np.sin(2 * theta_12))**2
    B = np.sin(2 * theta_13)**2

    sin_squared_Delta_31 = np.sin(Delta_31) ** 2
    sin_squared_Delta_12 = np.sin(Delta_12) ** 2
    
    print("values computed successfully")

    P_ee = 1 - (A * sin_squared_Delta_12 + B * sin_squared_Delta_31)
    print("P_ee computed successfully")
    
    return P_ee


# In[94]:


#set standard oscillation parameters
print("setting standard oscillation parameters")

theta_12 = 0.5903 #rad
theta_23 = 0.8430 #rad
theta_13 = 0.1503 #rad

delta_m_21_squared = 7.39 * 10**(-5) #eV^2
delta_m_32_squared = 2.449 * 10**(-3) #eV^2
delta_m_31_squared = delta_m_32_squared
print("defining energy array")
energy_array = np.linspace(1.8, 3.3, 20)


# In[95]:


## optimized definition of integral, optimized P_ee
def integral_over_positions_Th_U(points_array, energy_array, grid_1d_size, theta_12, delta_m_21_squared, A_Th, A_U, rho):
    dV = grid_1d_size**3
    
    relative_distance_array = calc_relative_dist(points_array)
    print("Relative distance array computed successfully")
    print(f'size : {len(relative_distance_array)}')
    print("   ")
    
    P_ee_array = P_ee_full_optimized(energy_array, points_array, theta_12, theta_13, delta_m_21_squared)
    print("P_ee_array computed successfully")
    print(f'size : {len(P_ee_array)}')
    print("   ")

    # Compute sum_Th
    sum_Th = np.sum(P_ee_array * ((A_Th * rho) / (4 * np.pi * (relative_distance_array**2)))[np.newaxis, :] * dV, axis=1)
    
    print("sum_Th computed successfully")
    print(f'size : {len(sum_Th)}')
    print("   ")
    
    # Compute sum_U
    sum_U = np.sum(P_ee_array * ((A_U * rho) / (4 * np.pi * (relative_distance_array**2)))[np.newaxis, :] * dV, axis=1)
    print("sum_U computed successfully")
    print(f'size : {len(sum_Th)}')
    print("   ")

    return sum_Th, sum_U, relative_distance_array, P_ee_array


# In[96]:


print("computing integral values for standard oscillation paramaters, whole crust")    
Th_integral_values_c, U_integral_values_c, relative_distance_array_c, P_ee_array_c = integral_over_positions_Th_U(crust_grid, energy_array, grid_1d_size_crust, theta_12, delta_m_21_squared, A_Th_c, A_U_c, rho_c)
print("done! moving on")
print('   ')


# In[97]:


print("computing integral values for standard oscillation paramaters, crust layer 1")    
Th_integral_values_c1, U_integral_values_c1, relative_distance_array_c1, P_ee_array_c1 = integral_over_positions_Th_U(crust_grid_1, energy_array, grid_1d_size_crust, theta_12, delta_m_21_squared, A_Th_c, A_U_c, rho_c)
print("done! moving on")
print('   ')


# In[98]:


print("computing integral values for standard oscillation paramaters, crust layer 2")    
Th_integral_values_c2, U_integral_values_c2, relative_distance_array_c2, P_ee_array_c2 = integral_over_positions_Th_U(crust_grid_2, energy_array, grid_1d_size_crust, theta_12, delta_m_21_squared, A_Th_c, A_U_c, rho_c)
print("done! moving on")
print('   ')


# In[99]:


print(Th_integral_values_c - Th_integral_values_c1 - Th_integral_values_c2)
print(U_integral_values_c - U_integral_values_c1 - U_integral_values_c2)


# In[100]:


print(Th_integral_values_c.dtype)


# In[101]:


# Step 1: Check dimensions
if P_ee_array_c.shape[0] != 20 or P_ee_array_c1.shape[0] != 20 or P_ee_array_c2.shape[0] != 20:
    raise ValueError("All arrays must have 20 rows.")

# Step 2: Concatenate P_ee_array_c1 and P_ee_array_c2
P_ee_combined = np.concatenate((P_ee_array_c1, P_ee_array_c2), axis=1)

# Step 3: Compare the arrays
if np.array_equal(P_ee_array_c, P_ee_combined):
    print("P_ee_array_c is a valid concatenation of P_ee_array_c1 and P_ee_array_c2.")
else:
    print("P_ee_array_c is NOT a valid concatenation of P_ee_array_c1 and P_ee_array_c2.")

# Optional: To see the shapes for debugging
print(f"P_ee_array_c shape: {P_ee_array_c.shape}")
print(f"P_ee_array_c1 shape: {P_ee_array_c1.shape}")
print(f"P_ee_array_c2 shape: {P_ee_array_c2.shape}")
print(f"Combined shape: {P_ee_combined.shape}")


# In[102]:


#Thing above cares about element order; if we don't:


# In[103]:


# Step 1: Check dimensions
if P_ee_array_c.shape[0] != 20 or P_ee_array_c1.shape[0] != 20 or P_ee_array_c2.shape[0] != 20:
    raise ValueError("All arrays must have 20 rows.")

# Step 2: Combine arrays into sets
# Flatten the arrays and convert to sets of tuples (since they are 2D arrays)
set_c = set(map(tuple, P_ee_array_c.T))   # Transpose to consider columns as individual points
set_c1 = set(map(tuple, P_ee_array_c1.T))
set_c2 = set(map(tuple, P_ee_array_c2.T))

# Step 3: Check if P_ee_array_c contains all elements from P_ee_array_c1 and P_ee_array_c2
if set_c == set_c1.union(set_c2):
    print("P_ee_array_c is a valid combination of P_ee_array_c1 and P_ee_array_c2, regardless of order.")
else:
    print("P_ee_array_c is NOT a valid combination of P_ee_array_c1 and P_ee_array_c2.")

with open('output.txt', 'w') as f:
    f.write("P_ee_array_c sample contents:\n")
    f.write(str(list(set_c)[:5]) + '\n')
    f.write("P_ee_array_c1 sample contents:\n")
    f.write(str(list(set_c1)[:5]) + '\n')
    f.write("P_ee_array_c2 sample contents:\n")
    f.write(str(list(set_c2)[:5]) + '\n')


# In[ ]:





# In[104]:


#make plots
print("making plots")
#Thorium
plt.plot(energy_array, Th_integral_values_c, marker='o', linestyle='-', color='b')
plt.plot(energy_array, Th_integral_values_c1 + Th_integral_values_c2, marker='o', linestyle='-', color='g')
plt.xlabel('Energy [MeV]')
plt.ylabel('Integral value')
plt.title('Integral over emission position, Thorium')
#plt.savefig('Integral over emission position, Thorium, standard oscillation parameters.pdf', format='pdf')
plt.show()

#Uranium
plt.plot(energy_array, U_integral_values_c, marker='o', linestyle='-', color='b')
plt.plot(energy_array, U_integral_values_c1 + U_integral_values_c2, marker='o', linestyle='-', color='g')
plt.xlabel('Energy [MeV]')
plt.ylabel('Integral value')
plt.title('Integral over emission position, Uranium')
#plt.savefig('Integral over emission position, Uranium, standard oscillation parameters.pdf', format='pdf')
plt.show()

#Both together

plt.plot(energy_array, Th_integral_values_c, marker='o', linestyle='-', color='b', label='Thorium')
plt.plot(energy_array, U_integral_values_c, marker='o', linestyle='-', color='b', label='Uranium')
plt.plot(energy_array, Th_integral_values_c1 + Th_integral_values_c2, marker='o', linestyle='-', color='g')
plt.plot(energy_array, U_integral_values_c1 + U_integral_values_c2, marker='o', linestyle='-', color='g')
plt.xlabel('Energy [MeV]')
plt.ylabel('Integral value')
plt.title('Integral over emission position')
plt.legend()
#plt.savefig('Integral over emission position, standard oscillation parameters.pdf', format='pdf')
plt.show()


# 

# In[105]:


# Ok so I think the issue isn't here either. Moving on:


# In[106]:


def integral_over_positions_Th_U_constant_P_ee(points_array, energy_array, grid_1d_size, A_Th, A_U, rho):
    dV = grid_1d_size**3
    
    relative_distance_array = calc_relative_dist(points_array)
    P_ee_mid = 0.521
    P_ee_stdev = 0.072
    P_ee_low = P_ee_mid - P_ee_stdev
    P_ee_high = P_ee_mid + P_ee_stdev
    P_ee_array_mid = np.full((len(energy_array), len(points_array)), P_ee_mid)
    P_ee_array_low = np.full((len(energy_array), len(points_array)), P_ee_low)
    P_ee_array_high = np.full((len(energy_array), len(points_array)), P_ee_high)
    
    
    # Compute sum_Th
    sum_mid_Th = np.sum(P_ee_array_mid * ((A_Th * rho) / (4 * np.pi * (relative_distance_array**2)))[np.newaxis, :] * dV, axis=1)
    sum_low_Th = np.sum(P_ee_array_low * ((A_Th * rho) / (4 * np.pi * (relative_distance_array**2)))[np.newaxis, :] * dV, axis=1)
    sum_high_Th = np.sum(P_ee_array_high * ((A_Th * rho) / (4 * np.pi * (relative_distance_array**2)))[np.newaxis, :] * dV, axis=1)
    
    # Compute sum_U
    sum_mid_U = np.sum(P_ee_array_mid * ((A_U * rho) / (4 * np.pi * (relative_distance_array**2)))[np.newaxis, :] * dV, axis=1)
    sum_low_U = np.sum(P_ee_array_low * ((A_U * rho) / (4 * np.pi * (relative_distance_array**2)))[np.newaxis, :] * dV, axis=1)
    sum_high_U = np.sum(P_ee_array_high * ((A_U * rho) / (4 * np.pi * (relative_distance_array**2)))[np.newaxis, :] * dV, axis=1)

    return sum_mid_Th, sum_mid_U, sum_low_Th, sum_low_U, sum_high_Th, sum_high_U


# In[107]:


print("computing integral values for standard P_ee (crust contribution only")
Th_integral_values_constant_P_ee_mid_c, U_integral_values_constant_P_ee_mid_c, Th_integral_values_constant_P_ee_low_c, U_integral_values_constant_P_ee_low_c, Th_integral_values_constant_P_ee_high_c, U_integral_values_constant_P_ee_high_c = integral_over_positions_Th_U_constant_P_ee (crust_grid, energy_array, grid_1d_size_crust, A_Th_c, A_U_c, rho_c)
print("done! moving on")


# In[108]:


print('   ')
print("computing integral values for standard P_ee (crust 1 contribution only")
Th_integral_values_constant_P_ee_mid_c1, U_integral_values_constant_P_ee_mid_c1, Th_integral_values_constant_P_ee_low_c1, U_integral_values_constant_P_ee_low_c1, Th_integral_values_constant_P_ee_high_c1, U_integral_values_constant_P_ee_high_c1 = integral_over_positions_Th_U_constant_P_ee (crust_grid_1, energy_array, grid_1d_size_crust, A_Th_c, A_U_c, rho_c)
print("done! moving on")


# In[109]:


print('   ')
print("computing integral values for standard P_ee (crust 1 contribution only")
Th_integral_values_constant_P_ee_mid_c2, U_integral_values_constant_P_ee_mid_c2, Th_integral_values_constant_P_ee_low_c2, U_integral_values_constant_P_ee_low_c2, Th_integral_values_constant_P_ee_high_c2, U_integral_values_constant_P_ee_high_c2 = integral_over_positions_Th_U_constant_P_ee (crust_grid_2, energy_array, grid_1d_size_crust, A_Th_c, A_U_c, rho_c)
print("done! moving on")


# In[110]:


#make plots, just for integrals with constant P_ee
print("making plots")
#Thorium
plt.plot(energy_array, Th_integral_values_constant_P_ee_mid_c, marker='o', linestyle='-', color='b')
plt.plot(energy_array, Th_integral_values_constant_P_ee_mid_c1 + Th_integral_values_constant_P_ee_mid_c2, marker='o', linestyle='-', color='g')
plt.xlabel('Energy [MeV]')
plt.ylabel('Integral value')
plt.title('Integral over emission position, Thorium')
plt.show()

#Uranium
plt.plot(energy_array, U_integral_values_constant_P_ee_mid_c, marker='o', linestyle='-', color='b')
plt.plot(energy_array, U_integral_values_constant_P_ee_mid_c1 + U_integral_values_constant_P_ee_mid_c2, marker='o', linestyle='-', color='g')
plt.xlabel('Energy [MeV]')
plt.ylabel('Integral value')
plt.title('Integral over emission position, Uranium')
plt.show()

#Both together

plt.plot(energy_array, Th_integral_values_constant_P_ee_mid_c, marker='o', linestyle='-', color='b', label='Thorium')
plt.plot(energy_array, U_integral_values_constant_P_ee_mid_c, marker='o', linestyle='-', color='b', label='Uranium')

plt.plot(energy_array, Th_integral_values_constant_P_ee_mid_c1 + Th_integral_values_constant_P_ee_mid_c2, marker='o', linestyle='-', color='g', label='Thorium')
plt.plot(energy_array, U_integral_values_constant_P_ee_mid_c1 + U_integral_values_constant_P_ee_mid_c2, marker='o', linestyle='-', color='g', label='Uranium')
plt.xlabel('Energy [MeV]')
plt.ylabel('Integral value')
plt.title('Integral over emission position')
plt.legend()
plt.show()
#make plots to compare integral for constant P_ee and non-constant(dependent on energy and oscillation parameters)

#Thorium
plt.plot(energy_array, Th_integral_values_constant_P_ee_mid_c, marker='o', linestyle='-', color='b', label='constant P_ee')
plt.plot(energy_array, Th_integral_values_c, marker='*', linestyle='-', color='b', label='non-constant P_ee, standard osc params')

plt.plot(energy_array, Th_integral_values_constant_P_ee_mid_c1 + Th_integral_values_constant_P_ee_mid_c2, marker='o', linestyle='-', color='g', label='constant P_ee')
plt.plot(energy_array, Th_integral_values_c1 + Th_integral_values_c2, marker='*', linestyle='-', color='g', label='non-constant P_ee, standard osc params')
plt.xlabel('Energy [MeV]')
plt.ylabel('Integral value')
plt.title('Integral over emission position, Thorium')
plt.legend()
#plt.savefig('Thorium Integral, standard osc params vs constant P_ee.pdf', format='pdf')
plt.show()

#Uranium
plt.plot(energy_array, U_integral_values_constant_P_ee_mid_c, marker='o', linestyle='-', color='b')
plt.plot(energy_array, U_integral_values_c, marker='*', linestyle='-', color='b', label = 'non-constant P_ee, standard osc params')

plt.plot(energy_array, U_integral_values_constant_P_ee_mid_c1 + U_integral_values_constant_P_ee_mid_c2, marker='o', linestyle='-', color='g')
plt.plot(energy_array, U_integral_values_c1 + U_integral_values_c2, marker='*', linestyle='-', color='g', label = 'non-constant P_ee, standard osc params')
plt.xlabel('Energy [MeV]')
plt.ylabel('Integral value')
plt.title('Integral over emission position, Uranium')
plt.legend()
#plt.savefig('Uranium Integral, standard osc params vs constant P_ee.pdf', format='pdf')
plt.show()


# In[111]:


#Still good ... 


# In[112]:


def sigma_IBD(energy_array):
    # optimized with numpy
    
    m_e = 0.511 #MeV
    m_p = 938 #MeV
    m_n = 941 #MeV
    E_th = 1.8 #MeV
    
    sigma = ((energy_array - E_th + m_e)**2) * ((1 - (m_e)**2/((energy_array - E_th + m_e)**2))**(1/2))
    
    print('   ')
    print('sigma computed for this energy array! It doesn\'t depend on position of courde :)' )
    print(f'size : {np.size(sigma)}')
    
    return sigma
#all these + sources in Earth model more bins -3 (for example)


# In[113]:


print("setting lambda and mu values")

lambda_U = 4.916
lambda_Th = 1.563
mu_U = 235
mu_Th = 232


# In[114]:


#get fluxes
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

plt.plot(energy_array_U, dn_dE_U)
plt.xlabel('E_nu [MeV]')
plt.yscale('log')
plt.ylabel('Intensity (arbitrary units)')
plt.title('Energy spectrum of geonus from U238 decay')
# 
# Add shaded region between 1.8 MeV and 3.2 MeV
plt.axvspan(1.8, 3.2, alpha=0.3, color='gray')
# Enable minor ticks on x-axis
plt.minorticks_on()
plt.savefig("Emission.pdf", format='pdf')
plt.show()

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

plt.plot(energy_array_Th, dn_dE_Th)
plt.xlabel('$E_{\nu}$ [MeV]')
plt.yscale('log')
plt.ylabel('Intensity (some units, doesn\'t matter)')
plt.title('Energy spectrum of geonus from Th232 decay')
# # Add shaded region between 1.8 MeV and 3.2 MeV
plt.axvspan(1.8, 3.2, alpha=0.3, color='gray')
# # Enable minor ticks on x-axis
plt.minorticks_on()
plt.show()
# Plot U238 decay data (blue line)
plt.plot(energy_array_U, dn_dE_U, label='U238 decays', color='blue')
# 
# # Plot Th232 decay data (red line)
plt.plot(energy_array_Th, dn_dE_Th, label='Th232 decays', color='red')
# 
plt.xlabel('E_nu [MeV]')
plt.yscale('log')
plt.ylabel('Intensity (arbitrary units)')
plt.title('Energy spectrum of geonus')
# 
# # Add shaded region between 1.8 MeV and 3.3 MeV
plt.axvspan(1.8, 3.3, alpha=0.3, color='gray')
# 
# # Enable minor ticks on x-axis
plt.minorticks_on()
# 
plt.legend(loc='upper right')
# 
#plt.show()

plt.savefig("Emission.pdf", format='pdf')
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
    counts_in_final_bins, _ = np.histogram(initial_bins, bins=np.concatenate([initial_bins, [2 * initial_bins[-1] - initial_bins[-2]]]), weights=counts_in_initial_bins)

    # Interpolate the counts to the final bin midpoints
    counts_in_final_bins = np.interp(final_bin_midpoints, bin_midpoints, counts_in_final_bins[:-1])

    return counts_in_final_bins

print("rebin to match energy array")

dn_dE_rebinned_U = rebin_counts(energy_array_U, dn_dE_U, energy_array)
dn_dE_rebinned_Th = rebin_counts(energy_array_Th, dn_dE_Th, energy_array)


# In[115]:


print("computing IBD cross section")
sigma = sigma_IBD(energy_array)


# In[116]:


print(sigma)


# In[117]:


print('    ')
print("computing fluxes at detector (standard oscillation parameters, crust contribution only")
N_Th_c = ((lambda_Th)/(mu_Th)) * sigma * dn_dE_rebinned_Th * Th_integral_values_c
N_U_c = ((lambda_U)/(mu_U)) * sigma * dn_dE_rebinned_U * U_integral_values_c
print(f'size N_Th_c : {np.size(N_Th_c)}')
print(f'size N_U_c : {np.size(N_U_c)}')


# In[118]:


print('    ')
print("computing fluxes at detector (standard oscillation parameters, crust 1 contribution only")
N_Th_c1 = ((lambda_Th)/(mu_Th)) * sigma * dn_dE_rebinned_Th * Th_integral_values_c1
N_U_c1 = ((lambda_U)/(mu_U)) * sigma * dn_dE_rebinned_U * U_integral_values_c1
print(f'size N_Th_c1 : {np.size(N_Th_c1)}')
print(f'size N_U_c1 : {np.size(N_U_c1)}')


# In[119]:


print('    ')
print("computing fluxes at detector (standard oscillation parameters, crust 2 contribution only")
N_Th_c2 = ((lambda_Th)/(mu_Th)) * sigma * dn_dE_rebinned_Th * Th_integral_values_c2
N_U_c2 = ((lambda_U)/(mu_U)) * sigma * dn_dE_rebinned_U * U_integral_values_c2
print(f'size N_Th_c2 : {np.size(N_Th_c2)}')
print(f'size N_U_c2 : {np.size(N_U_c2)}')


# In[120]:


print(N_U_c - N_U_c1 - N_U_c2)
print('    ')
print(N_Th_c - N_Th_c1 - N_U_c2)


# In[121]:


# Plot histogram for both
plt.step(energy_array, N_U_c + N_Th_c, where='mid', label='Total', color='green')
plt.step(energy_array, N_U_c, where='mid', label='U238 decays', color='blue')
plt.step(energy_array, N_Th_c, where='mid', label='Th232 decays', color='red')

plt.step(energy_array, N_U_c1 + N_U_c2 + N_Th_c1 + N_Th_c2, where='mid', label='Total', color='purple')
plt.step(energy_array, N_U_c1 + N_U_c2, where='mid', label='U238 decays', color='purple')
plt.step(energy_array, N_Th_c1 + N_Th_c2, where='mid', label='Th232 decays', color='purple')

plt.xlabel('E_nu [MeV]')
plt.yscale('log')
plt.ylabel('Expected counts (some units, doesn\'t matter)')
plt.title('Expected geonus')

# Set lower limit for y-axis to 10^-7
plt.ylim(bottom=6e-11)

# Enable minor ticks on x-axis
plt.minorticks_on()

# Add legend
plt.legend(loc='upper right')

plt.show()


# In[122]:


# Oh my god, it is himself!!! The issue!!!
# Try the steps above again and print everything

# Nevermind, I had a typo .... 


# In[ ]:





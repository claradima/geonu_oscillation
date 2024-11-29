# imports

import numpy as np
import matplotlib.pyplot as plt
import sys
import gc
import psutil
import os
import argparse

from functions import *
from test_functions import *

# Now test with actual crust

gridcount_c = 796

outer_rad_real = 6371

# Create crust points; one layer by default
print('creating 3d grid for crust')
crust_grid_specs, crust_grid_1d_size = create_3d_grid(grid_counts=gridcount_c)
print('cutting crust shell from 3d grid')

crust_points = cut_shell(inner_rad=6350 - np.sqrt(3)*crust_grid_1d_size/2, outer_rad=6371 + np.sqrt(3)*crust_grid_1d_size/2, sublayers=1, grid_specs=crust_grid_specs)
print(crust_points)
get_memory_usage()
print('deleting full crust grid')
del crust_grid_specs
print('deleted full crust grid')
get_memory_usage()
print(' ')
print(' ')
print(' ')

print(crust_points)
print(' ')
print(' ')
print(' ')

enlarged_points_c = crust_points[0]

print(f'len of enlarged crust point set is : {len(enlarged_points_c)}')
get_memory_usage()
print(' ')
print(' ')
print('computing raw set of midpoints')
unique_added_points_c_1 = create_midpoints(points_array = enlarged_points_c, spacing_1d = crust_grid_1d_size, acc = None )
print('computed')
print(f'len of added_points_raw_c_1 is : {len(unique_added_points_c_1)}')
get_memory_usage()
print(' ')
print(' ')


enlarged_points_c = np.concatenate((enlarged_points_c, unique_added_points_c_1), axis=0)
print(f'shape of enlarged_points_1 : {np.shape(enlarged_points_c)}')
print(f'shape of unique_added_points_1 : {np.shape(unique_added_points_c_1)}')
get_memory_usage()
print('removing unique_added_points_c_1')
del unique_added_points_c_1
get_memory_usage()

print(' ')
print('decrease enlarged shell')
crust_grid_1d_size_2 = crust_grid_1d_size/2

enlarged_points_c = trim_shell(inner_rad=6350, outer_rad=6371,  points_array=enlarged_points_c)
get_memory_usage()
gc.collect()
get_memory_usage()

csv_filename = f"crust_points_dist_{crust_grid_1d_size_2}.csv"
np.savetxt(csv_filename, enlarged_points_c, delimiter=",", header="x,y,z", comments="")
print(f"Saved trimmed points to {csv_filename}")
del csv_filename
gc.collect()


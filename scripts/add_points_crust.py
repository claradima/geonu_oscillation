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

gridcount_c = 396

outer_rad_real = 6371

# Create crust points; one layer by default
print('creating 3d grid for crust')
crust_grid_specs, crust_grid_1d_size = create_3d_grid(grid_counts=gridcount_c)
print('cutting crust shell from 3d grid')


crust_points = cut_shell(inner_rad=6350 - (np.sqrt(3)/2 + 1/10)*crust_grid_1d_size, outer_rad=6371 + (np.sqrt(3)/2 + 1/10)*crust_grid_1d_size, sublayers=1, grid_specs=crust_grid_specs)
print(crust_points)
get_memory_usage()
print('deleting full crust grid')
del crust_grid_specs
print('deleted full crust grid')
get_memory_usage()
print('gc collectiong')
gc.collect()
print('gc collected')
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
enlarged_points_new = remove_duplicates(enlarged_points_c, 0.1)
print(f'len of new enlarged crust point set is : {len(enlarged_points_new)}')
get_memory_usage()
print(' ')
print(' ')
print('computing raw set of midpoints')
unique_added_points_c_1 = create_midpoints(points_array = enlarged_points_c, spacing_1d = crust_grid_1d_size, acc = None )
print('computed')
print(f'len of added_points_raw_c_1 is : {len(unique_added_points_c_1)}')
added_points_new = remove_duplicates(unique_added_points_c_1, 0.1)
print(f'len(unique_added_points_c_1) : {len(unique_added_points_c_1)}')
print(f'len(added_points_new) : {len(added_points_new)}')
get_memory_usage()
print(' ')
print(' ')

enlarged_points_new = remove_duplicates(enlarged_points_c, 0.1)
print(f'len(enlarged_points_c) : {len(enlarged_points_c)}')
print(f'len(enlarged_points_new) before concatenation : {len(enlarged_points_new)}')
enlarged_points_c = np.concatenate((enlarged_points_c, unique_added_points_c_1), axis=0)
enlarged_points_new = remove_duplicates(enlarged_points_c, 0.1)
print(f'len(enlarged_points_c) : {len(enlarged_points_c)}')
print(f'len(enlarged_points_new) after concatenation : {len(enlarged_points_new)}')

print(f'shape of enlarged_points_1 : {np.shape(enlarged_points_c)}')
print(f'shape of unique_added_points_1 : {np.shape(unique_added_points_c_1)}')
get_memory_usage()
print('removing unique_added_points_c_1')
del unique_added_points_c_1
print('removed unique_added_points_c_1')
get_memory_usage()
print('gc collectiong')
gc.collect()
print('gc collected')
get_memory_usage()

print(' ')
print('decrease enlarged shell')
crust_grid_1d_size_2 = crust_grid_1d_size/2
print(f'crust_grid_1d_size_2 : {crust_grid_1d_size_2}')

enlarged_points_c = trim_shell(inner_rad=6350, outer_rad=6371,  points_array=enlarged_points_c)
get_memory_usage()
gc.collect()
get_memory_usage()

enlarged_points_new = remove_duplicates(enlarged_points_c, 0.1)
print(f'len(enlarged_points_c) : {len(enlarged_points_c)}')
print(f'len(enlarged_points_new) : {len(enlarged_points_new)}')

csv_filename = f"crust_points_dist_{crust_grid_1d_size_2}.csv"
np.savetxt(csv_filename, enlarged_points_c, delimiter=",", header="x,y,z", comments="")
print(f"Saved trimmed points to {csv_filename}")

# Add the file name to .gitignore (one folder back)
gitignore_path = os.path.join("..", ".gitignore")

# Ensure .gitignore exists
if not os.path.exists(gitignore_path):
    print(f".gitignore file not found at {gitignore_path}, creating a new one.")
    with open(gitignore_path, "w") as gitignore_file:
        pass  # Create an empty .gitignore file

# Append the file name to .gitignore if not already present
with open(gitignore_path, "r+") as gitignore_file:
    gitignore_content = gitignore_file.readlines()
    if csv_filename + "\n" not in gitignore_content:
        gitignore_file.write(csv_filename + "\n")
        print(f"Added {csv_filename} to {gitignore_path}")
    else:
        print(f"{csv_filename} is already in {gitignore_path}")

del csv_filename
gc.collect()


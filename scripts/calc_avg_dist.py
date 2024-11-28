import numpy as np
import matplotlib.pyplot as plt
import sys
import gc
import psutil
import os
import argparse

from functions import *

### NOTES
###
### You can also plot the unscaled spectra by uncommenting some lines

### TO DO :  add description of script

def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description = "Set global parameters.")

    # Add the -snor argument
    parser.add_argument('-snor', type = float, nargs = 3,
                        default=[0, 0, 6369],
                        help="Specify position of SNO detector SNO_r as three float values (default: [0, 0, 6369])")

    parser.add_argument('-cgridcount', type = int, nargs = 1,
                        default = 635,
                        help = 'Specify crust grid count (1d)') #690 for 18km
                                                                #820 for 15.5 km
                                                                #635 for 20km
    parser.add_argument('-mgridcount', type = int, nargs = 1,
                        default = 212, # 30 km for 420
                                       # 28 km for 450
                                       # 60 km for 212
                        help = 'Specify crust grid count (1d)')
    parser.add_argument('-Cshells', type = int, nargs = 1,
                        default = 1, help='Specify number of crust sublayers')
    parser.add_argument('-CLMshells', type = int, nargs = 1,
                        default = 1, help = 'Specify number of crust sublayers')
    parser.add_argument('-DMshells', type = int, nargs = 1,
                        default = 20, help = 'Specify number of DM sublayers')
    parser.add_argument('-EMshells', type = int, nargs = 1,
                        default = 10, help = 'Specify number of EM sublayers')

    # Parse the arguments
    args = parser.parse_args()

    # Create the SNO_r parameter as a NumPy array
    SNO_r = np.array(args.snor)

    # Print selected values for arguments
    print(f"SNO_r: {SNO_r}")


    # Generate Earth
    # Layered Earth Model : crust 6350 - 6371 km
    #                       CLM 6196 - 6350 km
    #                       DM 4202 - 6196 km
    #                       EM 3480 - 4202 km
    # Earth model stays the same in this script, except grid spacing

    # Create crust points; one layer by default
    print('creating 3d grid for crust')
    crust_grid_specs, crust_grid_1d_size = create_3d_grid(grid_counts = args.cgridcount)
    print('cutting crust shell from 3d grid')

    crust_points = cut_shell(inner_rad = 6350, outer_rad = 6371, sublayers = 1, grid_specs = crust_grid_specs)
    print(crust_points)
    get_memory_usage()
    print('deleting full crust grid')
    del crust_grid_specs
    print('deleted full crust grid')
    get_memory_usage()
    print(' ')

    print(crust_points)

    # Create mantle points
    print('creating 3d grid for mantle')
    mantle_grid_specs, mantle_grid_1d_size = create_3d_grid(grid_counts = args.mgridcount)

    # Create CLM shell; one layer by default
    print('cutting mantle CLM from 3d grid')
    CLM_points = cut_shell(inner_rad = 6196, outer_rad = 6350, sublayers = args.CLMshells, grid_specs = mantle_grid_specs)
    # Create DM shell;  25 layers by default
    print('cutting DM shell form 3d grid')
    DM_points = cut_shell(inner_rad = 4202, outer_rad = 6196, sublayers = args.DMshells, grid_specs = mantle_grid_specs)
    # Create EM shell;  10 layers by default
    print('cutting EM shell form 3d grid')
    EM_points = cut_shell(inner_rad = 3480, outer_rad = 4202, sublayers = args.EMshells, grid_specs=mantle_grid_specs)
    get_memory_usage()
    print('deleting full mantle grid')
    del mantle_grid_specs
    print('deleted full mantle grid')
    get_memory_usage()
    print(' ')

    # We have the Earth points now!
    # Can calculate average distance

    average_dist = avg_dist(crust_grid = crust_points, CLM_grid = CLM_points, DM_grids = DM_points, EM_grids = EM_points, grid_1d_size_crust = crust_grid_1d_size, grid_1d_size_mantle = mantle_grid_1d_size)

    print(f'average distance to emitting point for Layered Mantle model \n'
          f'crust spacing {crust_grid_1d_size}, mantle spacing {mantle_grid_1d_size}\n'
          f'DM shell number {args.DMshells}, EM shell number {args.EMshells}\n'
          f'is  : {average_dist}')

if __name__ == "__main__":
    main()
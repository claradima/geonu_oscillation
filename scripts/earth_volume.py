import numpy as np
import matplotlib.pyplot as plt
import sys
import gc
import psutil
import os
import argparse

from functions import *

def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description = "Set global parameters.")

    # Add the -snor argument
    parser.add_argument('-snor', type = float, nargs = 3,
                        default=[0, 0, 6369],
                        help="Specify position of SNO detector SNO_r as three float values (default: [0, 0, 6369])")

    parser.add_argument('-Ebins', type = int, nargs = 1,
                        default = 100, help = 'Specify number of energy bins')

    parser.add_argument('-cgridcount', type = int, nargs = 1,
                        default = 785,
                        help = 'Specify crust grid count (1d)') #690 for 18km
                                                                #820 for 15.5 km
                                                                #635 for 20km

    parser.add_argument('-mgridcount', type = int, nargs = 1,
                        default = 159, # 30 km for 420
                                       # 28 km for 450
                                       # 60 km for 212
                        help = 'Specify crust grid count (1d)')

    parser.add_argument('-Cshells', type = int, nargs = 1,
                        default = 1, help='Specify number of crust sublayers')

    parser.add_argument('-CLMshells', type = int, nargs = 1,
                        default = 1, help = 'Specify number of crust sublayers')

    parser.add_argument('-DMshells', type = int, nargs = 1,
                        default = 30, help = 'Specify number of DM sublayers')

    parser.add_argument('-EMshells', type = int, nargs = 1,
                        default = 10, help = 'Specify number of EM sublayers')

    # Parse the arguments
    args = parser.parse_args()

    # Create crust points; one layer by default

    '''
    print('creating 3d grid for crust')
    crust_grid_specs, crust_grid_1d_size = create_3d_grid(grid_counts=args.cgridcount)
    print('cutting crust shell from 3d grid')

    crust_points = cut_shell(inner_rad=6350, outer_rad=6371, sublayers=1, grid_specs=crust_grid_specs)
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
    mantle_grid_specs, mantle_grid_1d_size = create_3d_grid(grid_counts=args.mgridcount)

    # Create CLM shell; one layer by default
    print('cutting mantle CLM from 3d grid')
    CLM_points = cut_shell(inner_rad=6196, outer_rad=6350, sublayers=args.CLMshells, grid_specs=mantle_grid_specs)
    # Create DM shell;  25 layers by default
    print('cutting DM shell form 3d grid')
    DM_points = cut_shell(inner_rad=4202, outer_rad=6196, sublayers=args.DMshells, grid_specs=mantle_grid_specs)
    # Create EM shell;  10 layers by default
    print('cutting EM shell form 3d grid')
    EM_points = cut_shell(inner_rad=3480, outer_rad=4202, sublayers=args.EMshells, grid_specs=mantle_grid_specs)
    get_memory_usage()
    print('deleting full mantle grid')
    del mantle_grid_specs
    print('deleted full mantle grid')
    get_memory_usage()
    print(' ')
    '''

    # I want to calculate the volume if we're using the crust array with spacing 8

    # first import the content of the file
    # Specify the file name (replace with your actual value of crust_grid_1d_size_2)
    #csv_filename = f"crust_points_dist_8.csv"
    csv_filename = f"crust_points_dist_16.04785894206543.csv"
    # Read the CSV file into the crust_points variable
    crust_points_0 = np.loadtxt(csv_filename, delimiter=",", skiprows=1)  # Skip the header
    print(f'crust_points_0 : {crust_points_0}')

    print(f"Loaded points from {csv_filename}")
    get_memory_usage()
    #crust_grid_1d_size = 8.013836477987297
    crust_grid_1d_size = 16.04785894206543

    get_memory_usage()
    print('ha!')
    crust_points = np.array([crust_points_0])
    get_memory_usage()
    print('haha!')
    del crust_points_0
    get_memory_usage()
    print('hahaha!')
    gc.collect()
    get_memory_usage()


    # We have the Earth points now!


    # added function in functions.py to compute volume of a shell

    print('computing crust volume')
    no_points_c, volume_c = get_shell_volume(crust_points, crust_grid_1d_size)
    print(f'number of points in crust is : {no_points_c}')
    print(f'total volume of crust is : {volume_c}')
    print(' ')
    '''
    print('computing CLM volume')
    no_points_CLM, volume_CLM = get_shell_volume(CLM_points, mantle_grid_1d_size)
    print(f'number of points in CLM is : {no_points_CLM}')
    print(f'total volume of crust is : {volume_CLM}')
    print(' ')

    print('computing DM volume')
    no_points_DM, volume_DM = get_shell_volume(DM_points, mantle_grid_1d_size)
    print(f'number of points in DM is : {no_points_DM}')
    print(f'total volume of DM is : {volume_DM}')
    print(' ')

    print('computing EM volume')
    no_points_EM, volume_EM = get_shell_volume(EM_points, mantle_grid_1d_size)
    print(f'number of points in EM is : {no_points_EM}')
    print(f'total volume of EM is : {volume_EM}')
    print(' ')

    no_points_m = no_points_CLM + no_points_DM + no_points_EM
    print(f'number of points in mantle is : {no_points_m}')
    volume_m = volume_CLM + volume_DM + volume_EM
    print(f'total volume of mantle is : {volume_m}')

    no_points_tot = no_points_c + no_points_m
    print(f'number of points in total is : {no_points_tot}')
    volume_tot = volume_c + volume_m
    print(f'total volume of total is : {volume_tot}')
    '''

    print(' ')
    print(' ')
    print(' ')
    print(' ')

    print(f'SO!!! total crust points : {no_points_c}')
    #print(f'total mantle points : {no_points_m}')
    #print(f'total total points : {no_points_tot}')
    #print(' ')
    print(f'total crust volumes : {volume_c}')
    #print(f'total mantle volumes : {volume_m}')
    #print(f'total total volumes : {volume_tot}')


if __name__ == "__main__":
    main()

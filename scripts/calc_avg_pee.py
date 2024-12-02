import numpy as np
import matplotlib.pyplot as plt
import sys
import gc
import psutil
import os
import argparse

from functions import *
#from scripts.test_add_points_crust import csv_filename


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
    parser.add_argument('-abd', type = str, nargs = 1,
                        default = 'mid',
                        help = 'Specify abundances. Pick from \'mid\',\'low\', \'high\'')
    parser.add_argument('-Ebins', type = int, nargs = 1,
                        default = 100, help = 'Specify number of energy bins')
    parser.add_argument('-cgridcount', type = int, nargs = 1,
                        default = 799,
                        help = 'Specify crust grid count (1d)') #690 for 18km
                                                                #820 for 15.5 km
                                                                #635 for 20km
    parser.add_argument('-mgridcount', type = int, nargs = 1,
                        default = 231, # 30 km for 420
                                       # 28 km for 450
                                       # 60 km for 212
                        help = 'Specify crust grid count (1d)')
    parser.add_argument('-Cshells', type = int, nargs = 1,
                        default = 10, help='Specify number of crust sublayers')
    parser.add_argument('-CLMshells', type = int, nargs = 1,
                        default = 1, help = 'Specify number of crust sublayers')
    parser.add_argument('-DMshells', type = int, nargs = 1,
                        default = 30, help = 'Specify number of DM sublayers')
    parser.add_argument('-EMshells', type = int, nargs = 1,
                        default = 15, help = 'Specify number of EM sublayers')
    parser.add_argument('-theta12mid', type = float, nargs = 1,
                        default = 0.5882,
                        help = 'Specify standard value of theta_12 in rad; default = 0.5882 from James/Tony constrained fit')
    parser.add_argument('-dm21mid', type=float, nargs=1,
                        default = 7.58 * 10**(-5),
                        help = 'Specify standard value of delta_m_21^2 in eV^2; default = 7.58 * 10**(-5) from James/Tony constrained fit')
    parser.add_argument('-theta12err', type=float, nargs=1,
                        default = 0.014,
                        help = 'Specify error of theta_12 in rad; default = 0.0.014 from James/Tony constrained fit')
    parser.add_argument('-dm21err', type=float, nargs=1,
                        default = 0.18,
                        help = 'Specify error of delta_m_21^2 in eV^2; default = 0.18 from James/Tony constrained fit')
    parser.add_argument('-specsave', type = bool, nargs = 1,
                        default = False,
                        help = 'Specify whether to save the spectra (data in csv files and plots')
    parser.add_argument('-livetime', type = float, nargs = 1,
                        default = 134,
                        help = 'Specify detector livetime in days; default 100')
    parser.add_argument('-plotshow', type = bool, nargs = 1,
                        default = True,
                        help = 'Specify whether you want the plots shown while running the script; default no')
    # Parse the arguments
    args = parser.parse_args()

    # Create the SNO_r parameter as a NumPy array
    SNO_r = np.array(args.snor)

    # Print selected values for arguments
    print(f"SNO_r: {SNO_r}")
    print(f'Abundance set : {args.abd}')



    # Set abundances
    # Note : don't want these defined globally in case we want to compare
    #        spectra obtained from different abundances
    A_Th_c, A_U_c, A_Th_CLM, A_U_CLM, A_Th_DM, A_U_DM, A_Th_EM, A_U_EM, rho_c, rho_CLM, rho_DM, rho_EM = set_abund(args.abd)
    print(f'rho_c : {rho_c}')
    print(A_Th_CLM)
    # Setting fixed oscillation parameters
    set_fixed_params()
    # Set values for other 2 parameters

    theta_12_mid = args.theta12mid
    dm_21_sq_mid = args.dm21mid
    theta_12err = args.theta12err
    dm_21_sq_err = args.dm21err

    theta_12_low = theta_12_mid - theta_12err
    dm_21_sq_low = dm_21_sq_mid - dm_21_sq_err
    theta_12_high = theta_12_mid + theta_12err
    dm_21_sq_high = dm_21_sq_mid + dm_21_sq_err

    # Set e, n, p masses and E_th_IBD
    set_masses()
    # Set lambda and mu values for U and Th
    set_lambda_mu()


    # Set energy array
    # energy_array set globally!
    set_energy_array(args.Ebins)
    # Compute sigma_IBD for these E bins
    # sigma set globally!!!
    calc_sigma_IBD()

    # Get emission fluxes
    #  dn_dE_rebinned_U and dn_dE_rebinned_Th defined globally
    get_emission_fluxes(plot_spectrum=False)

    # Generate Earth
    # Layered Earth Model : crust 6350 - 6371 km
    #                       CLM 6196 - 6350 km
    #                       DM 4202 - 6196 km
    #                       EM 3480 - 4202 km
    # Earth model stays the same in this script, except grid spacing

    # Create crust points; one layer by default
    print('creating 3d grid for crust, for some reason')


    crust_grid_specs, crust_grid_1d_size = create_3d_grid(grid_counts = args.cgridcount)
    print('cutting crust shell from 3d grid')

    crust_points = cut_shell(inner_rad = 6350, outer_rad = 6371, sublayers = args.Cshells, grid_specs = crust_grid_specs)
    print(crust_points)
    get_memory_usage()
    print('deleting full crust grid')
    del crust_grid_specs
    print('deleted full crust grid')
    get_memory_usage()
    print(' ')
    '''
    print('hello')
    #csv_filename = f"crust_points_dist_8.csv"
    #csv_filename = 'crust_points_dist_10.064770932069678.csv'
    csv_filename = f'crust_points_dist_16.04785894206543.csv'

    # Read the CSV file into the crust_points variable
    crust_points_0 = np.loadtxt(csv_filename, delimiter=",", skiprows=1)  # Skip the header
    print(f'crust_points_0 : {crust_points_0}')

    print(f"Loaded points from {csv_filename}")
    get_memory_usage()
    crust_grid_1d_size = 16.04785894206543

    get_memory_usage()
    print('ha!')
    print(f'there are {len(crust_points_0)} crust points')
    print(f'splitting crust points shell into {args.Cshells} shells')
    crust_points = split_shell(inner_rad=6350, outer_rad=6371, sublayers=args.Cshells,
                               points_array_unwrapped=crust_points_0)
    # crust_points = np.array([crust_points_0])
    print('deleting initial crust array (unsplit)')
    get_memory_usage()
    del crust_points_0
    print('deleted; gc collecting')
    gc.collect()
    print('gc collected')

'''
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
    
    f_E_Th, f_E_U = set_energy_weights()

    print('COMPUTING NUMS AND DENOMS FOR C')

    numerator_Th_c, numerator_U_c, denominator_Th_c, denominator_U_c = calc_avg_P_ee_num_den(grids=crust_points,
                                                                                             theta_12=theta_12_mid,
                                                                                             delta_m_21_squared=dm_21_sq_mid,
                                                                                             grid_1d_size=crust_grid_1d_size,
                                                                                             A_Th=A_Th_c,
                                                                                             A_U=A_U_c,
                                                                                             rho=rho_c,
                                                                                             f_E_Th=f_E_Th,
                                                                                             f_E_U=f_E_U)
    print('COMPUTING NUMS AND DENOMS FOR CLM')
    numerator_Th_CLM, numerator_U_CLM, denominator_Th_CLM, denominator_U_CLM = calc_avg_P_ee_num_den(
        grids=CLM_points,
        theta_12=theta_12_mid,
        delta_m_21_squared=dm_21_sq_mid,
        grid_1d_size=mantle_grid_1d_size,
        A_Th=A_Th_CLM,
        A_U=A_U_CLM,
        rho=rho_CLM,
        f_E_Th=f_E_Th,
        f_E_U=f_E_U
    )
    print('COMPUTING NUMS AND DENOMS FOR DM')
    numerator_Th_DM, numerator_U_DM, denominator_Th_DM, denominator_U_DM = calc_avg_P_ee_num_den(
        grids=DM_points,
        theta_12=theta_12_mid,
        delta_m_21_squared=dm_21_sq_mid,
        grid_1d_size=mantle_grid_1d_size,
        A_Th=A_Th_DM,
        A_U=A_U_DM,
        rho=rho_DM,
        f_E_Th=f_E_Th,
        f_E_U=f_E_U
    )
    print('COMPUTING NUMS AND DENOMS FOR EM')
    numerator_Th_EM, numerator_U_EM, denominator_Th_EM, denominator_U_EM = calc_avg_P_ee_num_den(
        grids=EM_points,
        theta_12=theta_12_mid,
        delta_m_21_squared=dm_21_sq_mid,
        grid_1d_size=mantle_grid_1d_size,
        A_Th=A_Th_EM,
        A_U=A_U_EM,
        rho=rho_EM,
        f_E_Th=f_E_Th,
        f_E_U=f_E_U
    )
    get_memory_usage()
    average_P_ee = calc_P_ee_avg(numerator_Th_c, numerator_U_c, numerator_Th_CLM, numerator_U_CLM, numerator_Th_DM, numerator_U_DM, numerator_Th_EM, numerator_U_EM, denominator_Th_c, denominator_U_c, denominator_Th_CLM, denominator_U_CLM, denominator_Th_DM, denominator_U_DM, denominator_Th_EM, denominator_U_EM)
    get_memory_usage()
    print('deleting numerators')

    del numerator_Th_c
    del numerator_U_c
    del numerator_Th_CLM
    del numerator_U_CLM
    del numerator_Th_DM
    del numerator_U_DM
    del numerator_Th_EM
    del numerator_U_EM

    print('numerators deleted')
    get_memory_usage()

    print('COMPUTING VAR NUMS FOR C')
    numerator_var_Th_c, numerator_var_U_c = calc_var_P_ee_num(grids=crust_points,
                                                              theta_12=theta_12_mid,
                                                              delta_m_21_squared=dm_21_sq_mid,
                                                              grid_1d_size=crust_grid_1d_size,
                                                              A_Th=A_Th_c,
                                                              A_U=A_U_c,
                                                              rho=rho_c,
                                                              f_E_Th=f_E_Th,
                                                              f_E_U=f_E_U,
                                                              P_ee_average=average_P_ee)

    get_memory_usage()
    print(f'numerator_var_Th_c : {numerator_var_Th_c}')
    print(f'numerator_var_Th_c : {numerator_var_Th_c}')
    print(f'numerator_var_U_c : {numerator_var_U_c}')
    print('don\'t need crust shell anymore; delete')
    del crust_points
    print('deleted')
    get_memory_usage()
    print('COMPUTING VAR NUMS FOR CLM')
    numerator_var_Th_CLM, numerator_var_U_CLM = calc_var_P_ee_num(grids=CLM_points,
                                                                  theta_12=theta_12_mid,
                                                                  delta_m_21_squared=dm_21_sq_mid,
                                                                  grid_1d_size=mantle_grid_1d_size,
                                                                  A_Th=A_Th_CLM,
                                                                  A_U=A_U_CLM,
                                                                  rho=rho_CLM,
                                                                  f_E_Th=f_E_Th,
                                                                  f_E_U=f_E_U,
                                                                  P_ee_average=average_P_ee)

    get_memory_usage()
    print(f'numerator_var_Th_CLM : {numerator_var_Th_CLM}')
    print(f'numerator_var_U_CLM : {numerator_var_U_CLM}')
    print('don\'t need CLM shell anymore; delete')
    del CLM_points
    print('deleted')
    get_memory_usage()

    print('COMPUTING VAR NUMS FOR DM')
    numerator_var_Th_DM, numerator_var_U_DM = calc_var_P_ee_num(grids=DM_points,
                                                                theta_12=theta_12_mid,
                                                                delta_m_21_squared=dm_21_sq_mid,
                                                                grid_1d_size=mantle_grid_1d_size,
                                                                A_Th=A_Th_DM,
                                                                A_U=A_U_DM,
                                                                rho=rho_DM,
                                                                f_E_Th=f_E_Th,
                                                                f_E_U=f_E_U,
                                                                P_ee_average=average_P_ee)

    get_memory_usage()
    print(f'numerator_var_Th_DM : {numerator_var_Th_DM}')
    print(f'numerator_var_U_DM : {numerator_var_U_DM}')
    print('don\'t need DM shell anymore; delete')
    del DM_points
    print('deleted')
    get_memory_usage()

    print('COMPUTING VAR NUMS FOR EM')
    numerator_var_Th_EM, numerator_var_U_EM = calc_var_P_ee_num(grids=EM_points,
                                                                theta_12=theta_12_mid,
                                                                delta_m_21_squared=dm_21_sq_mid,
                                                                grid_1d_size=mantle_grid_1d_size,
                                                                A_Th=A_Th_EM,
                                                                A_U=A_U_EM,
                                                                rho=rho_EM,
                                                                f_E_Th=f_E_Th,
                                                                f_E_U=f_E_U,
                                                                P_ee_average=average_P_ee)

    get_memory_usage()
    print(f'numerator_var_Th_EM : {numerator_var_Th_EM}')
    print(f'numerator_var_U_EM : {numerator_var_U_EM}')
    print('don\'t need EM shell anymore; delete')
    del EM_points
    print('deleted')
    get_memory_usage()

    print('calculating stdev')

    stdev_P_ee = calc_P_ee_stdev(numerator_var_Th_c, numerator_var_U_c,
                                 numerator_var_Th_CLM, numerator_var_U_CLM,
                                 numerator_var_Th_DM, numerator_var_U_DM,
                                 numerator_var_Th_EM, numerator_var_U_EM,
                                 denominator_Th_c, denominator_U_c,
                                 denominator_Th_CLM, denominator_U_CLM,
                                 denominator_Th_DM, denominator_U_DM,
                                 denominator_Th_EM, denominator_U_EM)

    print(' ')
    print(' ')
    print(' ')
    print(f'average P_ee is : {average_P_ee}')
    print(f'P_ee stdev is : {stdev_P_ee}')


if __name__ == "__main__":
    main()

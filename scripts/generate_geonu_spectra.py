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
                        default = 820,
                        help = 'Specify crust grid count (1d)') #690 for 18km
                                                                #820 for 15.5 km
    parser.add_argument('-mgridcount', type = int, nargs = 1,
                        default = 420, # 30 km for 420
                                       # 28 km for 450
                        help = 'Specify crust grid count (1d)')
    parser.add_argument('-Cshells', type = int, nargs = 1,
                        default = 1, help='Specify number of crust sublayers')
    parser.add_argument('-CLMshells', type = int, nargs = 1,
                        default = 1, help = 'Specify number of crust sublayers')
    parser.add_argument('-DMshells', type = int, nargs = 1,
                        default = 25, help = 'Specify number of DM sublayers')
    parser.add_argument('-EMshells', type = int, nargs = 1,
                        default = 10, help = 'Specify number of EM sublayers')
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
                        default = True,
                        help = 'Specify whether to save the spectra (data in csv files and plots')
    parser.add_argument('-livetime', type = float, nargs = 1,
                        default = 100,
                        help = 'Specify detector livetime in days; default 100')
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
    # Can start calculating for each shell - standard osc params

    print('computing volume integrals for each layer')
    print(f'rho_c : {rho_c}')
    int_Th_c, int_U_c = add_vol_integrals(grids_array = crust_points,
                                          grid_1d_size = crust_grid_1d_size,
                                          theta_12 = theta_12_mid,
                                          delta_m_21_squared = dm_21_sq_mid,
                                          A_Th = A_Th_c, A_U = A_U_c, rho = rho_c)
    print(f'int_Th_c : {int_Th_c}')
    print(f'int_U_c : {int_U_c}')
    print(' ')
    print(' ')
    print('CRUST DONE')
    print(' ')
    print(' ')
    int_Th_CLM, int_U_CLM = add_vol_integrals(grids_array = CLM_points,
                                          grid_1d_size = mantle_grid_1d_size,
                                          theta_12 = theta_12_mid,
                                          delta_m_21_squared = dm_21_sq_mid,
                                          A_Th = A_Th_CLM, A_U = A_U_CLM, rho = rho_CLM)
    print(f'int_Th_CLM : {int_Th_CLM}')
    print(f'int_U_CLM : {int_U_CLM}')
    print(' ')
    print(' ')
    print('CLM DONE')
    print(' ')
    print(' ')
    int_Th_DM, int_U_DM = add_vol_integrals(grids_array=DM_points,
                                              grid_1d_size=mantle_grid_1d_size,
                                              theta_12=theta_12_mid,
                                              delta_m_21_squared=dm_21_sq_mid,
                                              A_Th=A_Th_DM, A_U=A_U_DM, rho=rho_DM)
    print(f'int_Th_DM : {int_Th_DM}')
    print(f'int_U_DM : {int_U_DM}')
    print(' ')
    print(' ')
    print('DM DONE')
    print(' ')
    print(' ')
    int_Th_EM, int_U_EM = add_vol_integrals(grids_array=EM_points,
                                              grid_1d_size=mantle_grid_1d_size,
                                              theta_12=theta_12_mid,
                                              delta_m_21_squared=dm_21_sq_mid,
                                              A_Th=A_Th_EM, A_U=A_U_EM, rho=rho_EM)
    print(f'int_Th_EM : {int_Th_EM}')
    print(f'int_U_EM : {int_U_EM}')
    print(' ')
    print(' ')
    print('EM DONE')
    print(' ')
    print(' ')
    print('ALL DONE')
    print('computing total volume integrals')
    # Computed all volume integrals with standard osc params
    # Now compute total volume integral over all layers
    int_Th = int_Th_c + int_Th_CLM + int_Th_DM + int_Th_EM
    int_U = int_U_c + int_U_CLM + int_U_DM + int_U_EM
    print(' ')
    print('computing total flux')

    N_Th_unsc, N_U_unsc = calc_exp_spec(U_vol_int=int_U, Th_vol_int=int_Th)
    print('computed unscaled spectrum')
    #plot_spec(N_Th = N_Th_unsc, N_U = N_U_unsc, spec_save = args.specsave,
    #          grid_1d_size_crust=crust_grid_1d_size,
    #          grid_1d_size_mantle=mantle_grid_1d_size,
    #          abd_set=args.abd,
    #          title_prefix = 'spec_unscaled_standard_osc_params' )

    print(' ')
    print('computing scaled spectrum')
    N_Th_scaled, N_U_scaled, geonus_tot = calc_exp_spec_scale(N_Th = N_Th_unsc, N_U = N_U_unsc, livetime = days_to_seconds(args.livetime) )
    print('computed scaled spectrum. now plotting ...')
    plot_spec(N_Th=N_Th_scaled, N_U=N_U_scaled, spec_save=args.specsave,
              grid_1d_size_crust=crust_grid_1d_size,
              grid_1d_size_mantle=mantle_grid_1d_size,
              abd_set=args.abd,
              title_prefix='spec_scaled_standard_osc_params')
    print('NICE!! We computed spectrum for standard osc params')
    print(' ')
    print(' ')

    # TO DO (VERY IMPORTANT) : scaling is messed up!!!!
    # Compute spectrum with constant P_ee
    #  Using P_ee 0.5278 ± 0.0721

    ### TO DO : Recalculate P_ee for various abundances and grids
    ###
    print('computing spectrum and for constant P_ee and +- err limits')
    print('using P_ee 0.5278 ± 0.0721')

    int_ctP_mid_Th_c, int_ctP_mid_U_c, int_ctP_low_Th_c, int_ctP_low_U_c, int_ctP_high_Th_c, int_ctP_high_U_c \
        = add_vol_integrals_const_P_ee(grids_array = crust_points,
                                       grid_1d_size=crust_grid_1d_size,
                                       A_Th = A_Th_c, A_U = A_U_c, rho = rho_c,
                                       P_ee_mid = 0.5278, P_ee_stdev = 0.0721)
    print(' ')
    print(' ')
    print('CRUST DONE')
    print(' ')
    print(' ')
    int_ctP_mid_Th_CLM, int_ctP_mid_U_CLM, int_ctP_low_Th_CLM, int_ctP_low_U_CLM, int_ctP_high_Th_CLM, int_ctP_high_U_CLM \
        = add_vol_integrals_const_P_ee(grids_array=CLM_points,
                                       grid_1d_size=mantle_grid_1d_size,
                                       A_Th=A_Th_CLM, A_U=A_U_CLM, rho=rho_CLM,
                                       P_ee_mid = 0.5278, P_ee_stdev = 0.0721)
    print(' ')
    print(' ')
    print('CLM DONE')
    print(' ')
    print(' ')
    int_ctP_mid_Th_DM, int_ctP_mid_U_DM, int_ctP_low_Th_DM, int_ctP_low_U_DM, int_ctP_high_Th_DM, int_ctP_high_U_DM \
        = add_vol_integrals_const_P_ee(grids_array=DM_points,
                                       grid_1d_size=mantle_grid_1d_size,
                                       A_Th=A_Th_DM, A_U=A_U_DM, rho=rho_DM,
                                       P_ee_mid = 0.5278, P_ee_stdev = 0.0721)
    print(' ')
    print(' ')
    print('DM DONE')
    print(' ')
    print(' ')
    int_ctP_mid_Th_EM, int_ctP_mid_U_EM, int_ctP_low_Th_EM, int_ctP_low_U_EM, int_ctP_high_Th_EM, int_ctP_high_U_EM \
        = add_vol_integrals_const_P_ee(grids_array=EM_points,
                                       grid_1d_size=mantle_grid_1d_size,
                                       A_Th=A_Th_EM, A_U=A_U_EM, rho=rho_EM,
                                       P_ee_mid = 0.5278, P_ee_stdev = 0.0721)
    print(' ')
    print(' ')
    print('EM DONE')
    print(' ')
    print(' ')
    print('all layers done, now adding up contributions')
    int_ctP_mid_U = int_ctP_mid_U_c + int_ctP_mid_U_CLM + int_ctP_mid_U_DM + int_ctP_mid_U_EM
    int_ctP_mid_Th = int_ctP_mid_Th_c + int_ctP_mid_Th_CLM + int_ctP_mid_Th_DM + int_ctP_mid_Th_EM
    print('integrals over whole Earth with constant P_ee computed')

    ### TO DO : only looking at plots with constant P_ee for now with no 'error bars'
    ###         so should add error bars later!

    print('computing unscaled fluxes, for midpoint values only')
    N_Th_ctP_mid_unsc, N_U_ctP_mid_unsc = calc_exp_spec(U_vol_int=int_ctP_mid_U, Th_vol_int=int_ctP_mid_Th)
    print('computing scaled fluxes, for midpoint values only')
    N_Th_ctP_mid_scaled, N_U_ctP_mid_scaled, geonus_ctP_mid_tot = calc_exp_spec_scale(N_Th=N_Th_ctP_mid_unsc, N_U=N_U_ctP_mid_unsc,
                                                              livetime=days_to_seconds(args.livetime))
    print('computed scaled spectrum. now plotting ...')

    print('all done! now plotting comparison between spectrum with constant P_ee and full oscillation formula')
    print('used standard osc params')

    plot_rat( N_Th_2 = N_Th_scaled, N_U_2 = N_U_scaled,
              N_Th_1 = N_Th_ctP_mid_scaled , N_U_1 = N_U_ctP_mid_scaled,
              spec_save = args.specsave,
              grid_1d_size_crust = crust_grid_1d_size,
              grid_1d_size_mantle = mantle_grid_1d_size,
              abd_set_1 = args.abd,
              abd_set_2 = args.abd,
              title_prefix_2="Scaled_standard_osc_params",
              title_prefix_1="Scaled_const_Pee")

    ##################################################################################################

    # Computing spectra for alternate oscillation parameters
    # Alternate parameters already calculated

    # Low theta_12

    print('low theta_12 : computing volume integrals for each layer, low theta_12')
    int_Th_low_theta_c, int_U_low_theta_c = add_vol_integrals(grids_array=crust_points,
                                          grid_1d_size=crust_grid_1d_size,
                                          theta_12=theta_12_low,
                                          delta_m_21_squared=dm_21_sq_mid,
                                          A_Th=A_Th_c, A_U=A_U_c, rho=rho_c)
    print(' ')
    print(' ')
    print('low theta_12 : CRUST DONE')
    print(' ')
    print(' ')
    int_Th_low_theta_CLM, int_U_low_theta_CLM = add_vol_integrals(grids_array=CLM_points,
                                              grid_1d_size=mantle_grid_1d_size,
                                              theta_12=theta_12_low,
                                              delta_m_21_squared=dm_21_sq_mid,
                                              A_Th=A_Th_CLM, A_U=A_U_CLM, rho=rho_CLM)
    print(' ')
    print(' ')
    print('low theta_12 : CLM DONE')
    print(' ')
    print(' ')
    int_Th_low_theta_DM, int_U_low_theta_DM = add_vol_integrals(grids_array=DM_points,
                                            grid_1d_size=mantle_grid_1d_size,
                                            theta_12=theta_12_low,
                                            delta_m_21_squared=dm_21_sq_mid,
                                            A_Th=A_Th_DM, A_U=A_U_DM, rho=rho_DM)
    print(' ')
    print(' ')
    print('low theta_12 : DM DONE')
    print(' ')
    print(' ')
    int_Th_low_theta_EM, int_U_low_theta_EM = add_vol_integrals(grids_array=EM_points,
                                            grid_1d_size=mantle_grid_1d_size,
                                            theta_12=theta_12_low,
                                            delta_m_21_squared=dm_21_sq_mid,
                                            A_Th=A_Th_EM, A_U=A_U_EM, rho=rho_EM)
    print(' ')
    print(' ')
    print('low theta_12 : EM DONE')
    print(' ')
    print(' ')
    print('low theta_12 : ALL DONE')
    print('low theta_12 : computing total volume integrals')
    # Computed all volume integrals with standard osc params
    # Now compute total volume integral over all layers
    int_Th_low_theta = int_Th_low_theta_c + int_Th_low_theta_CLM + int_Th_low_theta_DM + int_Th_low_theta_EM
    int_U_low_theta = int_U_low_theta_c + int_U_low_theta_CLM + int_U_low_theta_DM + int_U_low_theta_EM
    print(' ')
    print('low theta_12 : computing total flux')

    N_Th_low_theta_unsc, N_U_low_theta_unsc = calc_exp_spec(U_vol_int=int_U_low_theta, Th_vol_int=int_Th_low_theta)
    print('low theta_12 : computed unscaled spectrum.')
    #plot_spec(N_Th=N_Th_low_theta_unsc, N_U=N_U_low_theta_unsc, spec_save=args.specsave,
    #         grid_1d_size_crust=crust_grid_1d_size,
    #          grid_1d_size_mantle=mantle_grid_1d_size,
    #          abd_set=args.abd,
    #          title_prefix='spec_unscaled_low_theta_12')

    print(' ')
    print('low theta_12 : computing scaled spectrum')
    N_Th_low_theta_scaled, N_U_low_theta_scaled, geonus_low_theta_tot = calc_exp_spec_scale(N_Th=N_Th_low_theta_unsc,
                                                                                            N_U=N_U_low_theta_unsc,
                                                                                            livetime=days_to_seconds(args.livetime))
    print('low theta_12 : computed scaled spectrum. now plotting ...')
    plot_spec(N_Th=N_Th_low_theta_scaled, N_U=N_U_low_theta_scaled, spec_save=args.specsave,
              grid_1d_size_crust=crust_grid_1d_size,
              grid_1d_size_mantle=mantle_grid_1d_size,
              abd_set=args.abd,
              title_prefix='spec_scaled_low_theta_12')
    print('low theta_12 : NICE!! We computed spectrum for low theta')
    print(' ')
    print(' ')

    print('low theta_12 : now plotting comparison between spectrum with constant low and standard theta_12')

    plot_rat(N_Th_2=N_Th_scaled, N_U_2=N_U_scaled,
             N_Th_1=N_Th_low_theta_scaled, N_U_1=N_U_low_theta_scaled,
             spec_save=args.specsave,
             grid_1d_size_crust=crust_grid_1d_size,
             grid_1d_size_mantle=mantle_grid_1d_size,
             abd_set_1=args.abd,
             abd_set_2=args.abd,
             title_prefix_1="Scaled_low_theta_12",
             title_prefix_2="Scaled_standard_osc_params")

    # High theta_12

    print('high theta_12 : computing volume integrals for each layer, high theta_12')
    int_Th_high_theta_c, int_U_high_theta_c = add_vol_integrals(grids_array=crust_points,
                                          grid_1d_size=crust_grid_1d_size,
                                          theta_12=theta_12_high,
                                          delta_m_21_squared=dm_21_sq_mid,
                                          A_Th=A_Th_c, A_U=A_U_c, rho=rho_c)
    print(' ')
    print(' ')
    print('high theta_12 : CRUST DONE')
    print(' ')
    print(' ')
    int_Th_high_theta_CLM, int_U_high_theta_CLM = add_vol_integrals(grids_array=CLM_points,
                                              grid_1d_size=mantle_grid_1d_size,
                                              theta_12=theta_12_high,
                                              delta_m_21_squared=dm_21_sq_mid,
                                              A_Th=A_Th_CLM, A_U=A_U_CLM, rho=rho_CLM)
    print(' ')
    print(' ')
    print('high theta_12 : CLM DONE')
    print(' ')
    print(' ')
    int_Th_high_theta_DM, int_U_high_theta_DM = add_vol_integrals(grids_array=DM_points,
                                            grid_1d_size=mantle_grid_1d_size,
                                            theta_12=theta_12_high,
                                            delta_m_21_squared=dm_21_sq_mid,
                                            A_Th=A_Th_DM, A_U=A_U_DM, rho=rho_DM)
    print(' ')
    print(' ')
    print('high theta_12 : DM DONE')
    print(' ')
    print(' ')
    int_Th_high_theta_EM, int_U_high_theta_EM = add_vol_integrals(grids_array=EM_points,
                                            grid_1d_size=mantle_grid_1d_size,
                                            theta_12=theta_12_high,
                                            delta_m_21_squared=dm_21_sq_mid,
                                            A_Th=A_Th_EM, A_U=A_U_EM, rho=rho_EM)
    print(' ')
    print(' ')
    print('high theta_12 : EM DONE')
    print(' ')
    print(' ')
    print('high theta_12 : ALL DONE')
    print('high theta_12 : computing total volume integrals')
    # Computed all volume integrals with standard osc params
    # Now compute total volume integral over all layers
    int_Th_high_theta = int_Th_high_theta_c + int_Th_high_theta_CLM + int_Th_high_theta_DM + int_Th_high_theta_EM
    int_U_high_theta = int_U_high_theta_c + int_U_high_theta_CLM + int_U_high_theta_DM + int_U_high_theta_EM
    print(' ')
    print('high theta_12 : computing total flux')

    N_Th_high_theta_unsc, N_U_high_theta_unsc = calc_exp_spec(U_vol_int=int_U_high_theta, Th_vol_int=int_Th_high_theta)
    print('high theta_12 : computed unscaled spectrum')
    #plot_spec(N_Th=N_Th_high_theta_unsc, N_U=N_U_high_theta_unsc, spec_save=args.specsave,
    #          grid_1d_size_crust=crust_grid_1d_size,
    #          grid_1d_size_mantle=mantle_grid_1d_size,
    #          abd_set=args.abd,
    #          title_prefix='spec_unscaled_high_theta_12')

    print(' ')
    print('high theta_12 : computing scaled spectrum')
    N_Th_high_theta_scaled, N_U_high_theta_scaled, geonus_high_theta_tot = calc_exp_spec_scale(N_Th=N_Th_high_theta_unsc,
                                                                                            N_U=N_U_high_theta_unsc,
                                                                                            livetime=days_to_seconds(args.livetime))
    print('high theta_12 : computed scaled spectrum. now plotting ...')
    plot_spec(N_Th=N_Th_high_theta_scaled, N_U=N_U_high_theta_scaled, spec_save=args.specsave,
              grid_1d_size_crust=crust_grid_1d_size,
              grid_1d_size_mantle=mantle_grid_1d_size,
              abd_set=args.abd,
              title_prefix='spec_scaled_high_theta_12')
    print('high theta_12 : NICE!! We computed spectrum for high theta')
    print(' ')
    print(' ')

    print('high theta_12 : now plotting comparison between spectrum with constant high and standard theta_12')

    plot_rat(N_Th_2=N_Th_scaled, N_U_2=N_U_scaled,
             N_Th_1=N_Th_high_theta_scaled, N_U_1=N_U_high_theta_scaled,
             spec_save=args.specsave,
             grid_1d_size_crust=crust_grid_1d_size,
             grid_1d_size_mantle=mantle_grid_1d_size,
             abd_set_1=args.abd,
             abd_set_2=args.abd,
             title_prefix_1="Scaled_high_theta_12",
             title_prefix_2="Scaled_standard_osc_params")

    # Low delta_m_21^2

    print('low delta_m_21^2 : computing volume integrals for each layer, low delta_m_21^2')
    int_Th_low_dm_c, int_U_low_dm_c = add_vol_integrals(grids_array=crust_points,
                                                              grid_1d_size=crust_grid_1d_size,
                                                              theta_12=theta_12_mid,
                                                              delta_m_21_squared=dm_21_sq_low,
                                                              A_Th=A_Th_c, A_U=A_U_c, rho=rho_c)
    print(' ')
    print(' ')
    print('low delta_m_21^2 : CRUST DONE')
    print(' ')
    print(' ')
    int_Th_low_dm_CLM, int_U_low_dm_CLM = add_vol_integrals(grids_array=CLM_points,
                                                                  grid_1d_size=mantle_grid_1d_size,
                                                                  theta_12=theta_12_mid,
                                                                  delta_m_21_squared=dm_21_sq_low,
                                                                  A_Th=A_Th_CLM, A_U=A_U_CLM, rho=rho_CLM)
    print(' ')
    print(' ')
    print('low delta_m_21^2 : CLM DONE')
    print(' ')
    print(' ')
    int_Th_low_dm_DM, int_U_low_dm_DM = add_vol_integrals(grids_array=DM_points,
                                                                grid_1d_size=mantle_grid_1d_size,
                                                                theta_12=theta_12_mid,
                                                                delta_m_21_squared=dm_21_sq_low,
                                                                A_Th=A_Th_DM, A_U=A_U_DM, rho=rho_DM)
    print(' ')
    print(' ')
    print('low delta_m_21^2 : DM DONE')
    print(' ')
    print(' ')
    int_Th_low_dm_EM, int_U_low_dm_EM = add_vol_integrals(grids_array=EM_points,
                                                                grid_1d_size=mantle_grid_1d_size,
                                                                theta_12=theta_12_mid,
                                                                delta_m_21_squared=dm_21_sq_low,
                                                                A_Th=A_Th_EM, A_U=A_U_EM, rho=rho_EM)
    print(' ')
    print(' ')
    print('low delta_m_21^2 : EM DONE')
    print(' ')
    print(' ')
    print('low delta_m_21^2 : ALL DONE')
    print('low delta_m_21^2: computing total volume integrals')
    # Computed all volume integrals with standard osc params
    # Now compute total volume integral over all layers
    int_Th_low_dm = int_Th_low_dm_c + int_Th_low_dm_CLM + int_Th_low_dm_DM + int_Th_low_dm_EM
    int_U_low_dm = int_U_low_dm_c + int_U_low_dm_CLM + int_U_low_dm_DM + int_U_low_dm_EM
    print(' ')
    print('low delta_m_21^2 : computing total flux')

    N_Th_low_dm_unsc, N_U_low_dm_unsc = calc_exp_spec(U_vol_int=int_U_low_dm, Th_vol_int=int_Th_low_dm)
    print('low delta_m_21^2 : computed unscaled spectrum')
    #plot_spec(N_Th=N_Th_low_dm_unsc, N_U=N_U_low_dm_unsc, spec_save=args.specsave,
    #          grid_1d_size_crust=crust_grid_1d_size,
    #          grid_1d_size_mantle=mantle_grid_1d_size,
    #          abd_set=args.abd,
    #          title_prefix='spec_unscaled_low_dm_12')

    print(' ')
    print('low delta_m_21^2 : computing scaled spectrum')
    N_Th_low_dm_scaled, N_U_low_dm_scaled, geonus_low_dm_tot = calc_exp_spec_scale(N_Th=N_Th_low_dm_unsc,
                                                                                            N_U=N_U_low_dm_unsc,
                                                                                            livetime=days_to_seconds(
                                                                                                args.livetime))
    print('low delta_m_21^2 : computed scaled spectrum. now plotting ...')
    plot_spec(N_Th=N_Th_low_dm_scaled, N_U=N_U_low_dm_scaled, spec_save=args.specsave,
              grid_1d_size_crust=crust_grid_1d_size,
              grid_1d_size_mantle=mantle_grid_1d_size,
              abd_set=args.abd,
              title_prefix='spec_scaled_low_dm_12')
    print('low delta_m_21^2 : NICE!! We computed spectrum for low delta_m_21^2')
    print(' ')
    print(' ')

    print('low delta_m_21^2 : now plotting comparison between spectrum with constant low and standard delta_m_21^2')

    plot_rat(N_Th_2=N_Th_scaled, N_U_2=N_U_scaled,
             N_Th_1=N_Th_low_dm_scaled, N_U_1=N_U_low_dm_scaled,
             spec_save=args.specsave,
             grid_1d_size_crust=crust_grid_1d_size,
             grid_1d_size_mantle=mantle_grid_1d_size,
             abd_set_1=args.abd,
             abd_set_2=args.abd,
             title_prefix_1="Scaled_low_dm",
             title_prefix_2="Scaled_standard_osc_params")

    # High delta_m_21^2

    print('high delta_m_21^2 : computing volume integrals for each layer, high delta_m_21^2')
    int_Th_high_dm_c, int_U_high_dm_c = add_vol_integrals(grids_array=crust_points,
                                                          grid_1d_size=crust_grid_1d_size,
                                                          theta_12=theta_12_mid,
                                                          delta_m_21_squared=dm_21_sq_high,
                                                          A_Th=A_Th_c, A_U=A_U_c, rho=rho_c)
    print(' ')
    print(' ')
    print('high delta_m_21^2 : CRUST DONE')
    print(' ')
    print(' ')
    int_Th_high_dm_CLM, int_U_high_dm_CLM = add_vol_integrals(grids_array=CLM_points,
                                                              grid_1d_size=mantle_grid_1d_size,
                                                              theta_12=theta_12_mid,
                                                              delta_m_21_squared=dm_21_sq_high,
                                                              A_Th=A_Th_CLM, A_U=A_U_CLM, rho=rho_CLM)
    print(' ')
    print(' ')
    print('high delta_m_21^2 : CLM DONE')
    print(' ')
    print(' ')
    int_Th_high_dm_DM, int_U_high_dm_DM = add_vol_integrals(grids_array=DM_points,
                                                            grid_1d_size=mantle_grid_1d_size,
                                                            theta_12=theta_12_mid,
                                                            delta_m_21_squared=dm_21_sq_high,
                                                            A_Th=A_Th_DM, A_U=A_U_DM, rho=rho_DM)
    print(' ')
    print(' ')
    print('high delta_m_21^2 : DM DONE')
    print(' ')
    print(' ')
    int_Th_high_dm_EM, int_U_high_dm_EM = add_vol_integrals(grids_array=EM_points,
                                                            grid_1d_size=mantle_grid_1d_size,
                                                            theta_12=theta_12_mid,
                                                            delta_m_21_squared=dm_21_sq_high,
                                                            A_Th=A_Th_EM, A_U=A_U_EM, rho=rho_EM)
    print(' ')
    print(' ')
    print('high delta_m_21^2 : EM DONE')
    print(' ')
    print(' ')
    print('high delta_m_21^2 : ALL DONE')
    print('high delta_m_21^2: computing total volume integrals')
    # Computed all volume integrals with standard osc params
    # Now compute total volume integral over all layers
    int_Th_high_dm = int_Th_high_dm_c + int_Th_high_dm_CLM + int_Th_high_dm_DM + int_Th_high_dm_EM
    int_U_high_dm = int_U_high_dm_c + int_U_high_dm_CLM + int_U_high_dm_DM + int_U_high_dm_EM
    print(' ')
    print('high delta_m_21^2 : computing total flux')

    N_Th_high_dm_unsc, N_U_high_dm_unsc = calc_exp_spec(U_vol_int=int_U_high_dm, Th_vol_int=int_Th_high_dm)
    print('high delta_m_21^2 : computed unscaled spectrum')
    #plot_spec(N_Th=N_Th_high_dm_unsc, N_U=N_U_high_dm_unsc, spec_save=args.specsave,
    #          grid_1d_size_crust=crust_grid_1d_size,
    #          grid_1d_size_mantle=mantle_grid_1d_size,
    #          abd_set=args.abd,
    #          title_prefix='spec_unscaled_high_dm_12')

    print(' ')
    print('high delta_m_21^2 : computing scaled spectrum')
    N_Th_high_dm_scaled, N_U_high_dm_scaled, geonus_high_dm_tot = calc_exp_spec_scale(N_Th=N_Th_high_dm_unsc,
                                                                                      N_U=N_U_high_dm_unsc,
                                                                                      livetime=days_to_seconds(
                                                                                          args.livetime))
    print('high delta_m_21^2 : computed scaled spectrum. now plotting ...')
    plot_spec(N_Th=N_Th_high_dm_scaled, N_U=N_U_high_dm_scaled, spec_save=args.specsave,
              grid_1d_size_crust=crust_grid_1d_size,
              grid_1d_size_mantle=mantle_grid_1d_size,
              abd_set=args.abd,
              title_prefix='spec_scaled_high_dm_12')
    print('high delta_m_21^2 : NICE!! We computed spectrum for high delta_m_21^2')
    print(' ')
    print(' ')

    print('high delta_m_21^2 : now plotting comparison between spectrum with constant high and standard delta_m_21^2')

    plot_rat(N_Th_2=N_Th_scaled, N_U_2=N_U_scaled,
             N_Th_1=N_Th_high_dm_scaled, N_U_1=N_U_high_dm_scaled,
             spec_save=args.specsave,
             grid_1d_size_crust=crust_grid_1d_size,
             grid_1d_size_mantle=mantle_grid_1d_size,
             abd_set_1=args.abd,
             abd_set_2=args.abd,
             title_prefix_1="Scaled_high_dm",
             title_prefix_2="Scaled_standard_osc_params")


if __name__ == "__main__":
    main()
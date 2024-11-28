import numpy as np
import matplotlib.pyplot as plt
import sys
import gc
import psutil
import os
import argparse
import pandas as pd

from functions import *

### Compute and plot the ratio between two spectra provided .csv files
### Provided spectra must have been generated with the same energy bins !

def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description = "Set global parameters.")

    # Add the -snor argument

    parser.add_argument('-autopath1', type = bool, nargs = 1,
                        default = True,
                        help = 'Specify whether to use autopath for first spectrum')
    parser.add_argument('-autopath2', type=bool, nargs=1,
                        default=True,
                        help='Specify whether to use autopath for second spectrum')
    parser.add_argument('-csvtitle1', type = str, nargs = 1,
                        default = 'spec_scaled_standard_osc_params_spec_data_mid_100E20C60M.csv',
                        help = 'Specify the title of the first csv file')
                        # default is dumb here; change to '' TO DO
    parser.add_argument('-csvtitle2', type=str, nargs=1,
                        default='spec_scaled_standard_osc_params_spec_data_mid_100E20C60M.csv',
                        help='Specify the title of the second csv file')
                        # default is dumb here; change to '' TO DO
    parser.add_argument('-spec_save', type=bool, nargs=1,
                        default = True,
                        help = 'Specify whether to save the spectra')
    parser.add_argument('-plot_show', type=str, nargs=1,
                        default = False,)

    # Parse the arguments
    args = parser.parse_args()




    ### TO DO : take input from arg parser for path; hardcode what I need right now
    energy_array1, N_tot1 = extract_columns(args.csvtitle1, path='/home/claramariadima/SNO/geonu_oscillation/plots/sublayer_comparison/20C60M_10DM5EM/mid_100E20C60M')
    energy_array2, N_tot2 = extract_columns(args.csvtitle2, path='/home/claramariadima/SNO/geonu_oscillation/plots/sublayer_comparison/20C60M_2DM1EM/mid_100E20C60M')

    if np.array_equal(energy_array2, energy_array1):
        print('energy arrays are the same; all good!')
    else :
        print('energy arrays are not the same. boo! goodbye')
        sys.exit()

    ### NOTE : I tried using plot_rat but it's taking too many arguments to make
    ###        the whole title; I have to think about how to generalize that,
    ###        because I might want to plot various different things

    ### TO DO : make plot into simplified function, naming scheme

    plt.step(energy_array1, N_tot1 / N_tot2, where='mid', label='total', color='green')

    # Add grid
    plt.grid(True)  # Show grid
    plt.grid(which='both', color='lightgray', linestyle='--', linewidth=0.5)  # Customize grid appearance

    plt.xlabel('E_nu [MeV]')
    plt.ylabel('Expected geonu count ratio')
    plt.title(f'Ratio of expected geonus')

    plt.ylim(bottom=0.92, top=1.15)
    plt.minorticks_on()
    plt.legend(loc='lower left')

    if args.plot_show:
        print('showing plot')
        plt.show()

    if args.spec_save:
        plt.savefig('Ratio plot.pdf', format='pdf')

if __name__ == "__main__":
    main()
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
    parser = argparse.ArgumentParser(description="Set global parameters.")

    # Add the -snor argument
    parser.add_argument('-snor', type=float, nargs=3,
                        default=[0, 0, 6369],
                        help="Specify position of SNO detector SNO_r as three float values (default: [0, 0, 6369])")

    # Parse the arguments
    args = parser.parse_args()

    # Create the SNO_r parameter as a NumPy array
    SNO_r = np.array(args.snor)

    # Print the result
    print(f"SNO_r: {SNO_r}")


if __name__ == "__main__":
    main()
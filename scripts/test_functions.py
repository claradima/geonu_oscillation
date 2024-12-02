import numpy as np
import matplotlib.pyplot as plt
import sys
import gc
import psutil
import os
import argparse
from scipy.spatial import cKDTree
from functions import *

### Calculate distance between two points
###
### inputs : point1 : numpy array of len 3 with coords x1, y1, z1
###          point2 : ... x2, y2, z2

def get_dist(p1, p2):
    if len(p1) != 3 or len(p2) != 3:
        print('invalid input; points must be arrays of len 3 (coords x, y, z); exiting ...')
        sys.exit()
    # else:
    # print('points are valid')

    # print('computing distance between points')
    return np.sum((p1 - p2) ** 2)


### Check if a numpy array that contains cartesian coordinates
### of points contains duplicates (specified accuracy)
###
### inputs : points_array : a 2d numpy array, each element is a point of
###                         the form np.array([x, y, z])
###          acc : accuracy of check
###
### output : cleansed_points : a 2d numpy array with remaining points from
###                            points_array that are not duplicates

def remove_duplicates(points_array, acc):
    # Ensure points_array is a 2D numpy array
    if points_array.ndim != 2 or points_array.shape[1] != 3:
        print('invalid input; points_array must be a 2D array of shape (n, 3); exiting ...')
        sys.exit()

    # Build a KDTree for efficient nearest-neighbor search
    tree = cKDTree(points_array)

    # Find points within the specified accuracy range
    duplicates = tree.query_ball_tree(tree, acc)

    unique_points = []
    seen = set()

    print(f"Starting duplicate removal for {len(points_array)} points.")

    # Iterate over the points and keep track of the ones that are unique
    for i, neighbors in enumerate(duplicates):
        #print(f"Processing point {i+1}/{len(points_array)}...")

        # Check if any of the neighbors are already seen
        if not any(j in seen for j in neighbors):  # Ensure no duplicates are added
            unique_points.append(points_array[i])
            seen.update(neighbors)  # Mark these points as seen
            #print(f"Point {i+1}/{len(points_array)} is unique; added to the list.")
        #else:
            #print(f"Point {i+1}/{len(points_array)} is a duplicate; skipped.")

    print(f"Finished processing points. Found {len(points_array) - len(unique_points)} duplicates.")

    # Return the filtered points as a numpy array
    return np.array(unique_points)

def remove_duplicates_old(points_array, acc):
    get_memory_usage()
    print('copying array')
    cleansed_points = points_array.copy()
    print('copied')
    get_memory_usage()

    print('checking array for duplicates')
    duplicate_count = 0

    i = 0
    while i < len(cleansed_points):
        removed = False
        for j in range(i + 1, len(cleansed_points)):
            point1 = cleansed_points[i]
            point2 = cleansed_points[j]
            xdist = np.abs(point1[0] - point2[0])
            ydist = np.abs(point1[1] - point2[1])
            zdist = np.abs(point1[2] - point2[2])

            if (xdist < acc) and (ydist < acc) and (zdist < acc):
                print('found a duplicate!')
                print(f'looking at points of indices {i} and {j}')  # add
                print(f'xdistance between them is : {xdist}')  # add
                print(f'ydistance between them is : {ydist}')
                print(f'zdistance between them is : {zdist}')
                print(f'reminder: accuracy is {acc}')  # add
                cleansed_points = np.delete(cleansed_points, j, axis=0)
                duplicate_count += 1
                removed = True
                break  # Break from the inner loop since we already removed a point

        if not removed:
            i += 1  # Only move to the next element if no element was removed

    #print(' ')
    print(f'original length of the array : {len(points_array)}')
    print(f'length of new array : {len(cleansed_points)}')
    print(f'we found {duplicate_count} duplicate points')
    print(
        f'check that this is correct : len(points_array) - len(cleansed_points) = {len(points_array) - len(cleansed_points)}')

    return cleansed_points


### Given an array of points, create array of all the midpoints
### between 2 adjacent points;
### This includes points that share a unit side, 2d or 3d unit diagonal
###
### inputs : points_array : a 2d numpy array, each element is a point of
###                         the form np.array([x, y, z])
###          spacing_1d : dist between 2 adjacent points that share a unit
###                       edge
###          acc : accuracy of check : can check custom distances if this
###                value is not None
###
### outputs : added_points_raw : a 2d numpy array, each element is a point
###                              of the form np.array([x, y, z])

def create_midpoints(points_array, spacing_1d=None, acc=None):
    # Validate inputs
    if spacing_1d is None and acc is None:
        raise ValueError("Both spacing_1d and acc cannot be None")

    if acc is None:
        acc = spacing_1d * (np.sqrt(3) + (1 / 10))
        #acc = spacing_1d * (np.sqrt(3)) * 4

    # Build a KDTree for efficient nearest-neighbor search
    tree = cKDTree(points_array)

    # List to store unique midpoints
    added_points_raw_list = []
    midpoints_set = set()  # To track unique midpoints

    # Iterate over all points and find neighbors within the accuracy range
    for i, point1 in enumerate(points_array):

        # Find all neighbors within the accuracy range
        neighbors = tree.query_ball_point(point1, acc)

        for j in neighbors:
            if i < j:  # Ensure unique pairs (i, j)
                point2 = points_array[j]
                # Calculate midpoint
                midpoint = (point1 + point2) / 2

                # Check if the midpoint is already in the set of added points (using tuple for immutability)
                midpoint_tuple = tuple(midpoint)
                if midpoint_tuple not in midpoints_set:
                    added_points_raw_list.append(midpoint)
                    midpoints_set.add(midpoint_tuple)  # Mark the midpoint as added

        print(f"Processed point {i + 1}/{len(points_array)}...")

    print('outside of loop! converting to array')
    # Convert list of midpoints to a numpy array and return
    added_points_raw_list = np.array(added_points_raw_list)
    return added_points_raw_list

def create_midpoints_with_duplicates(points_array, spacing_1d=None, acc=None):
    # Validate inputs
    if spacing_1d is None and acc is None:
        raise ValueError("Both spacing_1d and acc cannot be None")

    if acc is None:
        acc = spacing_1d * (np.sqrt(3) + (1 / 5))

    # Build a KDTree for efficient nearest-neighbor search
    tree = cKDTree(points_array)

    added_points_raw_list = []

    # Iterate over all points and find neighbors within the accuracy range
    for i, point1 in enumerate(points_array):
        # Find all neighbors within the accuracy range
        neighbors = tree.query_ball_point(point1, acc)

        for j in neighbors:
            if i < j:  # Ensure unique pairs (i, j)
                point2 = points_array[j]
                # Calculate midpoint
                midpoint = (point1 + point2) / 2
                added_points_raw_list.append(midpoint)

    # Convert list of midpoints to a numpy array and return
    added_points_raw_array = np.array(added_points_raw_list)
    return added_points_raw_array

def create_midpoints_old(points_array, spacing_1d=None, acc=None):
    if spacing_1d is None and acc is None:
        print('invalid inputs; spacing_1d and acc can\'t both be None')
        sys.exit()

    if acc is None:
        acc = spacing_1d * (np.sqrt(3) + (1 / 5))

    print(f'spacing is : {spacing_1d}')  # add
    print(f'accuracy is : {acc}')  # add

    print(f'len(points_array) : {len(points_array)}')
    print('creating list of midpoints (will convert to array at the end')
    get_memory_usage()

    added_points_raw_list = []

    for i in range(len(points_array)):
        for j in range(i + 1, len(points_array)):
            #print(' ')  # add

            point1 = points_array[i]
            point2 = points_array[j]

            dist = get_dist(point1, point2)
            print(f'looking at points of indices {i} and {j} out of {len(points_array)}')#add
            # print(f'distance between them is : {dist}')#add
            # print(f'reminder: accuracy is {acc}')#add

            if dist < acc:
                print('adjacent points found')
                get_memory_usage()
                print('creating midpoint')

                x_mid = (point1[0] + point2[0]) / 2
                y_mid = (point1[1] + point2[1]) / 2
                z_mid = (point1[2] + point2[2]) / 2

                new_point = np.array([x_mid, y_mid, z_mid])
                print(f'len(added_points): {len(added_points_raw_list)}; appending midpoint')
                added_points_raw_list.append(new_point)
                print(f'appended; len(added_points): {len(added_points_raw_list)}')
                get_memory_usage()
                print('deleting intermediate bits')
                del x_mid
                del y_mid
                del z_mid
                del new_point
                print('deleted')
                get_memory_usage()
                print('gc collecting')
                gc.collect()
                print('gc collected')
                get_memory_usage()
                print('')

    print('done appending points')
    print(f'len(added_points): {len(added_points_raw_list)}')
    print(f'len(points_array) = {len(points_array)}')

    print('converting list to array')
    added_points_raw_array = np.array(added_points_raw_list)
    print('done')
    get_memory_usage()
    print('deleting list')
    del added_points_raw_list
    print('gc collecting')
    gc.collect()
    print('gc collected')
    get_memory_usage()
    print('done; nice!')

    return added_points_raw_array

### TO DO : these work for one layer, with no sublayers; need to adapt for actual data structure, but test first



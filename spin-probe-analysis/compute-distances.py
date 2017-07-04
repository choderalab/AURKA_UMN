#!/bin/env python

"""
Compute spin probe distances.

John D. Chodera
3 July 2017
"""

import numpy as np
import mdtraj as md
import os, os.path
import glob

fahdata_path = '/cbio/jclab/projects/fah/fah-data/munged4/11431/'
output_path = 'data' # location for output data
nthreads = 16

if not os.path.exists(output_path):
    os.makedirs(output_path)

def process_clone(clone_path):
    """
    Process a clone HDF5 trajectory file, writing a numpy file containing inter-spin-probe distances. 

    Parameters
    ----------
    clone_path : str
        File path to clone HDF5 file

    """
    print('Processing %s...' % clone_path)

    # Read trajectory
    traj = md.load(clone_path)

    # Determine spin probe oxygens
    oxygens = traj.top.select('resn CYR and name ON')

    # Compute distances
    distances = md.compute_distances(traj, [oxygens])
    distances = distances[:,0] # squeeze out unnecessary coordinates

    # Clean up
    del traj
    
    return distances

if __name__ == '__main__':
    from multiprocessing import Pool
    pool = Pool(nthreads)

    nruns = 12 # TODO: Automatically determine
    for run in range(nruns):
        # Find all CLONEs for this RUN
        h5_filenames = glob.glob(os.path.join(fahdata_path, 'run%d-*.h5'% run))
        print('RUN %5d : There are %d clones to process.' % (run, len(h5_filenames)))

        # Process these clones
        print('Processing on pool of %d threads...' % nthreads)
        distances = pool.map(process_clone, h5_filenames)

        # Save
        output_filename = os.path.join(output_path, 'run%d.npy' % run)
        np.save(output_filename, distances)



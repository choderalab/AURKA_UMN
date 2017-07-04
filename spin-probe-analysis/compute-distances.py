#!/bin/env python

"""
Compute spin probe distances.

John D. Chodera
3 July 2017
"""

import numpy as np
import mdtraj as md
import os, os.path
import re
import glob

fahdata_path = '/cbio/jclab/projects/fah/fah-data/munged4/11431/'
output_path = 'data' # location for output data
nthreads = 16

# offset_index[run] is the ACTUAL first residue number in the AurA protein sequence of the first residue of the topology PDB
# THIS IS A MASSIVE HACK DUE TO A FAILURE TO PRESERVE PROVENANCE INFORMATION IN PROTEIN SEQID
offset_index = {
    0 : 123, # 1OL5
    1 : 123,
    2 : 123,
    3 : 123,
    4 : 127, # 1OL7
    5 : 127,
    6 : 127,
    7 : 127,
    8 : 126, # 5L8K
    9 : 126,
    10 : 126,
    11 : 126,
}

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

    # Determine offset for first residue
    filename = os.path.basename(clone_path)
    [prefix, extension] = os.path.splitext(filename)
    match = re.match('run(\d+)-clone(\d+)', prefix)
    run = int(match[1])
    offset = offset_index[run]

    # Read trajectory
    traj = md.load(clone_path)

    # Determine spin probe oxygen atom indices
    oxygens = traj.top.select('resn CYR and name ON')
    if len(oxygens) != 2:
        raise Exception('Selection for O-O distance did not return exactly two atoms')

    # Determine R255 CZ - T288 CA distance
    RT = traj.top.select('(resSeq %d and resname ARG and name CZ) or (resSeq %d and resname THR and name CA)' % (255 - offset + 1, 288 - offset + 1))
    if len(RT) != 2:
        raise Exception('Selection for R255 CZ - T288 CA distance did not return exactly two atoms')

    # Compute distances
    distances = md.compute_distances(traj, [oxygens, RT])

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



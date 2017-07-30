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

    # Determine spin probe NO oxygen atom indices
    oxygens = traj.top.select('resn CYR and name NN')
    if len(oxygens) != 2:
        raise Exception('%s: Selection for spin probe NN-NN distance did not return exactly two atoms' % clone_path)

    # Determine spin probe CA distances
    alpha_carbons = traj.top.select('resn CYR and name CA')
    if len(alpha_carbons) != 2:
        raise Exception('%s: Selection for spin probe CA-CA distance did not return exactly two atoms' % clone_path)

    # Determine R255 CZ - T288 CA distance
    RT = traj.top.select('(resSeq %d and resname ARG and name CZ) or (resSeq %d and resname THR and name CA)' % (255 - offset + 1, 288 - offset + 1))
    if len(RT) != 2:
        raise Exception('%s: Selection for R255 CZ - T288 CA distance did not return exactly two atoms' % clone_path)

    # Determine F275 CZ - I193 CG2 distance
    DFG_PHE_ILE = traj.top.select('(resSeq %d and resname PHE and name CZ) or (resSeq %d and resname ILE and name CG2)' % (275 - offset + 1, 193 - offset + 1))
    if len(DFG_PHE_ILE) != 2:
        raise Exception('%s: Selection for F275 CZ - I193 CG2 distance did not return exactly two atoms' % clone_path)

    # Determine W277 CE2 - I193 CG2 distance
    DFG_TRP_ILE = traj.top.select('(resSeq %d and resname TRP and name CE2) or (resSeq %d and resname ILE and name CG2)' % (277 - offset + 1, 193 - offset + 1))
    if len(DFG_TRP_ILE) != 2:
        raise Exception('%s: Selection for W277 CE2 - I193 CG2 distance did not return exactly two atoms' % clone_path)

    # Determine P282 O - R285 H distance (activation loop helical contact)
    helix_hbonds_1 = traj.top.select('(resSeq %d and resname PRO and name O) or (resSeq %d and resname ARG and name H)' % (282 - offset + 1, 285 - offset + 1))
    if len(helix_hbonds_1) != 2:
        raise Exception('%s: Selection for P282 O - R285 HN distance did not return exactly two atoms (%d instead)' % (clone_path, len(helix_hbonds_1)))

    # Determine S283 O - R286 H distance (activation loop helical contact)
    helix_hbonds_2 = traj.top.select('(resSeq %d and resname SER and name O) or (resSeq %d and resname ARG and name H)' % (283 - offset + 1, 286 - offset + 1))
    if len(helix_hbonds_2) != 2:
        raise Exception('%s: Selection for S283 O - R286 HN distance did not return exactly two atoms (%d instead)' % (clone_path, len(helix_hbonds_2)))

    # Determine L225 CA - S284 CA distance (CYR-CYR here)
    #L225_S284 = traj.top.select('(resSeq %d and resname LEU and name CA) or (resSeq %d and resname SER and name CA)' % (225 - offset + 1, 284 - offset + 1))
    L225_S284 = traj.top.select('(resSeq %d and resname CYR and name CA) or (resSeq %d and resname CYR and name CA)' % (225 - offset + 1, 284 - offset + 1))
    if len(L225_S284) != 2:
        raise Exception('%s: Selection for L225 CA -  S284 CA distance did not return exactly two atoms (%d instead)' % (clone_path, len(L225_S284)))

    # Compute distances
    distances = md.compute_distances(traj, [oxygens, alpha_carbons, RT, DFG_PHE_ILE, DFG_TRP_ILE, helix_hbonds_1, helix_hbonds_2, L225_S284])

    # Determine pseudotorsion for residues 282-285
    r282ca = traj.top.select('(resSeq %d and name CA)' % (282 - offset + 1,))[0]
    r283ca = traj.top.select('(resSeq %d and name CA)' % (283 - offset + 1,))[0]
    r284ca = traj.top.select('(resSeq %d and name CA)' % (284 - offset + 1,))[0]
    r285ca = traj.top.select('(resSeq %d and name CA)' % (285 - offset + 1,))[0]
    r286ca = traj.top.select('(resSeq %d and name CA)' % (286 - offset + 1,))[0]
    torsions = md.compute_dihedrals(traj, np.array([
                [r282ca, r283ca, r284ca, r285ca],
                [r283ca, r284ca, r285ca, r286ca]]))                                                                     

    data = np.concatenate((distances, torsions), axis=1)

    # Clean up
    del traj
                                    
    return data

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



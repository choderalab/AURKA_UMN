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

fahdata_path = '/cbio/jclab/projects/fah/fah-data/munged3/no-solvent'
output_path = 'data' # location for output data
nthreads = 12

runs = {
    11414 : [0], # AurA +TPX2, pdb 1OL5
    11419 : [0, 1, 2, 3], # AurA +TPX2, pdb 1OL5
    11418 : [0, 1, 2, 3, 4], # AurA -TPX2, pdb 1OL5
    11423 : [], # AurA -TPX2, pdb 1OL5
    11428 : [0], # AurA +TPX2, pdb 1OL5 - pT288
    11429 : [0], # AurA -TPX2, pdb 1OL5 - pT288
}
projects = runs.keys()

if not os.path.exists(output_path):
    os.makedirs(output_path)

"""
  * index 0 : R255 CZ - T288 CA distance
  * index 1 : F275 CZ - I193 CG2 distance (low in DFG-in)
  * index 2 : W277 CE2 - I193 CG2 distance (low in DFG-out)
  * index 3 : P282 O - R285 H distance (activation loop helical contact)
  * index 4 : S283 O - R286 H distance (activation loop helical contact)
  * index 5 : L225 CA - S284 CA distance
  * index 6 : Pseudotorsion for residues 282-285
  * index 7 : Pseudotorsion for residues 283-286
"""

def get_atom_index(clone_path, traj, selection):
    indices = traj.top.select(selection)
    if len(indices) != 1:
        msg = '%s : Selection "%s" did not match a unique atom.' % (clone_path, selection)
        raise Exception(msg)

    return indices[0]

def compute_distance(clone_path, traj, *args):
    """
    Compute the specified distance.
    """
    indices = [ get_atom_index(clone_path, traj, selection) for selection in args ]

    # Compute distances in Angstroms
    distances = md.compute_distances(traj, [indices]).squeeze() * 10.0
    
    return distances

def compute_torsion(clone_path, traj, *args):
    """
    Compute the specified torsion.
    """
    indices = [ get_atom_index(clone_path, traj, selection) for selection in args ]

    # Compute torsion in degrees
    torsions = md.compute_dihedrals(traj, [indices]).squeeze() * (180.0 / np.pi)
    
    return torsions

def get_residue_index(clone_path, traj, selection):
    atom_indices = traj.top.select(selection)
    residue_indices = set([ atom.residue.index for atom in traj.top.atoms if atom.index in atom_indices ])
    if len(residue_indices) != 1:
        msg = '%s : Selection "%s" did not match a unique residue. ' % (clone_path, selection)
        msg += 'Matched atom indices: %s ' % (atom_indices)
        msg += 'Matched residue_indices: %s ' % (residue_indices)
        #traj[0].save('error.pdb')
        raise Exception(msg)
    residue_indices = list(residue_indices)
    
    return residue_indices[0]

def compute_contact_distance(clone_path, traj, *args):
    """
    Compute the specified sidechain closest heavy atom distance.
    """    
    residue_indices = [ get_residue_index(clone_path, traj, selection) for selection in args ]

    # Compute distances in Angstroms
    [distances, residue_pairs] = md.compute_contacts(traj, contacts=[residue_indices], scheme='sidechain-heavy')

    return distances.squeeze()

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

    # Compute observables
    data = dict()

    data['T288 - R180 distance (A)'] = compute_contact_distance(clone_path, traj, '(resSeq 288)', '(resSeq 180 and resname ARG)')
    data['T288 - R255 distance (A)'] = compute_contact_distance(clone_path, traj, '(resSeq 288)', '(resSeq 255 and resname ARG)')
    data['T288 - R286 distance (A)'] = compute_contact_distance(clone_path, traj, '(resSeq 288)', '(resSeq 286 and resname ARG)')

    data['R255 CZ - T288 CA distance (A)'] = compute_distance(clone_path, traj, '(resSeq 255 and resname ARG and name CZ)', '(resSeq 288 and name CA)')
    data['F275 CZ - I193 CG2 distance (A)'] = compute_distance(clone_path, traj, '(resSeq 275 and resname PHE and name CZ)', '(resSeq 193 and resname ILE and name CG2)')
    data['W277 CE2 - I193 CG2 distance (A)'] = compute_distance(clone_path, traj, '(resSeq 277 and resname TRP and name CE2)', '(resSeq 193 and resname ILE and name CG2)')
    data['P282 O - R285 H distance (A)'] = compute_distance(clone_path, traj, '(resSeq 282 and resname PRO and name O)', '(resSeq 285 and resname ARG and name H)')
    data['S283 O - R286 H distance (A)'] = compute_distance(clone_path, traj, '(resSeq 283 and resname SER and name O)', '(resSeq 286 and resname ARG and name H)')
    data['L225 CA - S284 CA distance (A)'] = compute_distance(clone_path, traj, '(resSeq 225 and resname LEU and name CA)', '(resSeq 284 and resname SER and name CA)')

    data['282-285 pseudodihedral (degrees)'] = compute_torsion(clone_path, traj, *['(resSeq %d and name CA)' % resSeq for resSeq in (282, 283, 284, 285)])
    data['283-286 pseudodihedral (degrees)'] = compute_torsion(clone_path, traj, *['(resSeq %d and name CA)' % resSeq for resSeq in (283, 284, 285, 286)])

    # Clean up
    del traj
                                    
    return data

if __name__ == '__main__':
    from multiprocessing import Pool
    pool = Pool(nthreads)

    nruns = 12 # TODO: Automatically determine
    for project in projects:
        # Find all WT RUNs for this project
        for run in runs[project]:
            # Find all CLONEs for this RUN
            h5_filenames = glob.glob(os.path.join(fahdata_path, '%s/run%d-*.h5'% (project, run)))
            print('RUN %5d : There are %d clones to process.' % (run, len(h5_filenames)))

            # Process these clones
            print('Processing on pool of %d threads...' % nthreads)
            data = pool.map(process_clone, h5_filenames)

            # Save
            output_filename = os.path.join(output_path, 'proj%d-run%d.npy' % (project, run))
            np.save(output_filename, data)



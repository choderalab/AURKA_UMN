#!/usr/bin/python

import numpy as np
import sys
import math
import os
import time
import mdtraj as md

"""
Extract all trajectories and associated waters meeting W1/W2 definition for that trajectory.
Requires /cbio/jclab/conditions/behrj/AURKA_UMN/output-{condition#}/data/{condition#}_{run#}_{residue#}_distHBonds.npy files

"""

from mpi4py import MPI
rank = MPI.COMM_WORLD.rank
size = MPI.COMM_WORLD.size

if rank==0: print('rank = %d, size = %d' % (rank, size))

project_basepath = '/cbio/jclab/conditions/fah/fah-data/munged3/all-atoms' # location of FAH all-atom trajectories
hbond_basepath = '/cbio/jclab/conditions/behrj/AURKA_UMN/'
output_basepath = 'trajectories-with-key-waters'

if not os.path.exists(output_basepath):
    os.makedirs(output_basepath)

def get_molecule(molecules, index):
    """Return the molecule containing the specified atom index"""
    for m in molecules:
        atom_indices = set([ a.index for a in m ])
        if len(atom_indices.intersection(W_set)) > 0:
            return m
    return None

def accumulate_frame_count(molecules, molecule_frames, frame_set):
    """Accumulate count of all molecules that atoms in specified frame_set belong to.
    
    Parameters
    ----------
    molecules : set of molecules
       The set of molecules in the system
    molecule_frames : dict of molecule : int
       molecule_frames[m] is the number of frames molecule m appears in
    frame_set : set of int
       Set of atom indices corresponding to water oxygens
    """

    if W1_frame_set is not None:
        for oxygen_index in frame_set:
            m = get_molecule(molecules, oxygen_index)
            molecule_frames[m] += 1

def split_major_minor(molecule_frames, MIN_FRAME_COUNT=20):
    """Extract atom indices of major and minor molecule frame counts.
    
    Parameters
    ----------
    molecule_frames : dict of molecule : count
       molecule_frames[m] is the number of frames molecule m appears in

    Returns
    -------
    major_atom_indices : list of int
       Atom indices of atoms that appear >= MIN_FRAME_COUNT times
    minor_atom_indices : list of int
       Atom indices of atoms that appear < MIN_FRAME_COUNT times

    """
    major_atom_indices = list()
    minor_atom_indices = list()
    for m in molecule_frames:
        if molecule_frames[m] >= MIN_FRAME_COUNT:
            major_atom_indices += [ a.index for a in m ]
        else:
            minor_atom_indices += [ a.index for a in m ]
    return (major_atom_indices, minor_atom_indices)

def write_trajectory(t, prefix='trajectory'):
    """Write PDB and DCD files"""
    t[0].save(prefix + '.pdb')
    t.save(prefix + '.dcd')

def molecule_contains_atom_in_set(m, atom_indices):
    """Return True if molecule m contains an atom index in set 'atom_indices'."""
    return len(set(atom_indices).intersection([a.index for a in m])) > 0
            
def atom_mapping(atomslice):
    return { old : new for (new, old) in enumerate(atomslice) }

nclones = 50 # number of CLONEs per RUN
projects = ['11414', '11419', '11418', '11423']
nruns = 7 # number of runs per condition
resid = 185 # residue number 
for project in projects:
    for run in range(nruns):
        if rank==0: print('PROJECT %s RUN %d' % (project, run))

        # Load hbond info
        try:
            W1_filename = os.path.join(hbond_basepath, 'output-%(project)s/data/RUN%(run)d-274N-oxygen-indices.npy' % vars())
            W2_filename = os.path.join(hbond_basepath, 'output-%(project)s/data/RUN%(run)d-W2-181185-or275-oxygen-indices.npy' % vars())
            W1 = np.load(W1_filename)
            W2 = np.load(W2_filename)
        except Exception as e:
            print(e)
            continue

        # Process trajectories
        for clone in range(rank,nclones,size):
            time.sleep(rank*0.2)

            # WX_trajectory[frame] is the set of water oxygen indices satisfying criteria for water X at frame 'frame'
            W1_trajectory = W1[clone,:]
            W2_trajectory = W2[clone,:]

            # Get all waters that satisfy W1 or W2 in this trajectory
            W1_set = set.union(*[s for s in W1_trajectory if s is not None])
            W2_set = set.union(*[s for s in W2_trajectory if s is not None])
            W_set = set.union(W1_set, W2_set)
            print('  CLONE %5d : W1 %8d W2 %8d : total %8d' % (clone, len(W1_set), len(W2_set), len(W_set)))

            # Read trajectory
            if rank==0: print('      Reading trajectory...')
            initial_time = time.time()
            h5_filename = os.path.join(project_basepath, '%s/run%d-clone%d.h5' % (project, run, clone))
            t = md.load(h5_filename)
            if rank==0: print('          completed in %8.3f s' % (time.time() - initial_time))
            
            # Report statistics
            components = {
                'AurA' : "(resSeq >= 123) and (resSeq <= 388)",
                'Tpx2' : "(resSeq < 123)",
                'ADP' : "resname == MOL",
                'ions' : "resname == MG",
                'water' : "water"
                }
            components = { name : t.topology.select(selection) for (name, selection) in components.iteritems() }
            for (name, indices) in components.iteritems():
                print('%12s : %12d' % (name, len(indices)))

            # Purge waters that aren't involved in W1 or W2
            if rank==0: print('      Stripping inessential waters...')
            solute_atom_indices = list(t.topology.select("not water"))
            # Filter waters
            molecules = t.topology.find_molecules()
            water_molecules = [ m for m in molecules if molecule_contains_atom_in_set(m, W_set) ] # waters involved in W1 or W2 at some point in trajectory
            water_atom_indices = [ a.index for m in water_molecules for a in m] # atom indices for waters
            # Slice out inessential waters
            atomslice = solute_atom_indices + water_atom_indices
            #t.atom_slice(atomslice, inplace=True)
            # Update water index trajectories
            atom_mapping = atom_mapping(atomslice)
            W1_trajectory = [ set([atom_mapping[old_index] for old_index in old_set]) for old_set in W1_trajectory ]
            W2_trajectory = [ set([atom_mapping[old_index] for old_index in old_set]) for old_set in W2_trajectory ]
            if rank==0: print('          completed in %8.3f s' % (time.time() - initial_time))
                        
            # Align the trajectory to AurA
            if rank==0: print('      Aligning trajectory...')
            initial_time = time.time()
            AurA = t.topology.select("(resSeq >= 123) and (resSeq <= 388)")
            #t.superpose(t, frame=0, atom_indices=AurA, parallel=True)
            if rank==0: print('          completed in %8.3f s' % (time.time() - initial_time))

            # Image the trajectory around the protein and peptide.
            if rank==0: print('      Imaging trajectory...')
            initial_time = time.time()
            molecules = t.topology.find_molecules()
            anchor_molecules = molecules[0:2]
            #t.image_molecules(inplace=True, anchor_molecules=anchor_molecules)
            if rank==0: print('          completed in %8.3f s' % (time.time() - initial_time))

            # Extract atom subsets of interest
            if rank==0: print('      Slicing trajectory for protein...')
            initial_time = time.time()
            solute = t.topology.select('not water')            
            write_trajectory(t.atom_slice(solute, inplace=False), prefix=os.path.join(output_basepath, 'condition%s-run%d-clone%d-protein' % (project, run, clone)))
            if rank==0: print('          completed in %8.3f s' % (time.time() - initial_time))

            # Determine major and minor waters for W1 and W2
            if rank==0: print('      Slicing trajectory for waters...')
            initial_time = time.time()
            molecules = t.topology.find_molecules()
            W1_molecule_frames = { m : 0 for m in molecules } # molecule_frames[molecule] is the number of frames 'molecule' appears in
            W2_molecule_frames = { m : 0 for m in molecules } # molecule_frames[molecule] is the number of frames 'molecule' appears in
            for (W1_frame_set, W2_frame_set) in zip(W1[clone,:], W2[clone,:]):
                accumulate_frame_count(molecules, W1_molecule_frames, W1_frame_set)
                accumulate_frame_count(molecules, W2_molecule_frames, W2_frame_set)
            (W1_major_waters, W1_minor_waters) = split_major_minor(W1_molecule_frames)
            (W2_major_waters, W2_minor_waters) = split_major_minor(W2_molecule_frames)            
            # Slice out major and minor water populations
            write_trajectory(t.atom_slice(W1_major_waters, inplace=False), prefix=os.path.join(output_basepath, 'condition%s-run%d-clone%d-W1-major' % (project, run, clone)))
            write_trajectory(t.atom_slice(W1_minor_waters, inplace=False), prefix=os.path.join(output_basepath, 'condition%s-run%d-clone%d-W1-minor' % (project, run, clone)))
            write_trajectory(t.atom_slice(W2_major_waters, inplace=False), prefix=os.path.join(output_basepath, 'condition%s-run%d-clone%d-W2-major' % (project, run, clone)))
            write_trajectory(t.atom_slice(W2_minor_waters, inplace=False), prefix=os.path.join(output_basepath, 'condition%s-run%d-clone%d-W2-minor' % (project, run, clone)))
            if rank==0: print('          completed in %8.3f s' % (time.time() - initial_time))
            
            MPI.COMM_WORLD.Barrier()

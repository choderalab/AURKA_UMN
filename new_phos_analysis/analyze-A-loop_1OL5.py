#!/usr/bin/python

import numpy as np
import sys
import math
import os
import time
import mdtraj as md

"""
Analyze C-helix

"""

from mpi4py import MPI
rank = MPI.COMM_WORLD.rank
size = MPI.COMM_WORLD.size

if rank==0: print('rank = %d, size = %d' % (rank, size))

project_basepath = '/cbio/jclab/projects/fah/fah-data/munged4' # location of FAH trajectories
output_basepath = './data/aloop'

if not os.path.exists(output_basepath):
    os.makedirs(output_basepath)

# Load reference structure for comparison of alphaC RMSD
reference_pdbfile = '../fah-setup/1OL5-WT-pdbfixer.pdb'
reference = md.load(reference_pdbfile)
alignment_selection_dsl = '(resSeq >= 123) and (resSeq <= 387) and (name CA)'
alignment_reference_indices = reference.topology.select(alignment_selection_dsl)

rmsd_selection_dsl = '(resSeq >= 286) and (resSeq <= 293) and (name CA)'
rmsd_reference_indices = reference.topology.select(rmsd_selection_dsl)

nclones = 250 # number of CLONEs per RUN
nframes = 4000 # max frames / trajectory
projects = ['11428', '11429']
#conditions = ['11428', '11429']
nruns = 7 # number of runs per condition
for project in projects:
    for run in range(nruns):
        h5_filename = os.path.join(project_basepath, '%s/run%d-clone%d.h5' % (project, run, 0))
        if not os.path.exists(h5_filename):
            continue
        
        if rank==0: print('PROJECT %s RUN %d' % (project, run))

        # Process trajectories
        rmsds = list()
        for clone in range(rank,nclones,size):
            # Read trajectory
            initial_time = time.time()
            h5_filename = os.path.join(project_basepath, '%s/run%d-clone%d.h5' % (project, run, clone))
            t = md.load(h5_filename)
            print('  CLONE %5d : %5d frames' % (clone, t.n_frames))

            # Align to reference using AurA CA atoms only
            alignment_trajectory_indices = t.topology.select(alignment_selection_dsl)
            t.superpose(reference, atom_indices=alignment_trajectory_indices, ref_atom_indices=alignment_reference_indices, parallel=False)
            # Compute the RMSD manually without additional alignment
            rmsd_trajectory_indices = t.topology.select(rmsd_selection_dsl)
            rmsd = np.sqrt(3*np.mean((t.xyz[:, rmsd_trajectory_indices, :] - reference.xyz[:, rmsd_reference_indices, :])**2, axis=(1,2)))
            rmsds.append(rmsd)

        # Gather data to root and write it
        gathered_rmsds = MPI.COMM_WORLD.gather(rmsds, root=0)
        if rank==0:
            rmsds = np.zeros([nclones, nframes], np.float32) - 1
            for clone in range(nclones):
                rmsd = gathered_rmsds[clone % size][clone / size]
                N = len(rmsd)
                rmsds[clone,0:N] = rmsd[0:N]
            output_filename = os.path.join(output_basepath, '%s-run%d-aloop-rmsd_1OL5.npy' % (project, run))
            np.save(output_filename, rmsds)

        


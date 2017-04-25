
"""
Analyze C-helix

"""

from mpi4py import MPI
import mdtraj as md
import numpy as np
import sys
import math
import os
import time
import mdtraj as md
rank = MPI.COMM_WORLD.rank
size = MPI.COMM_WORLD.size

if rank == 0: print('rank = %d, size = %d' % (rank, size))

project_basepath = '/cbio/jclab/home/albaness/trajectories/AURKA'  # location of FAH trajectories
output_basepath = '../data/e-fret'

nclones = 50  # number of CLONEs per RUN
nframes = 2040  # max frames / trajectory
conditions = ['AURKA_nophos_notpx2/11418', 'AURKA_phos_notpx2/11429']
save_file_label = ['NoPhos', 'Phos']
nruns = 1
for i, condition in enumerate(conditions):
    for run in range(nruns):
        if rank == 0: print('PROJECT %s RUN %d' % (condition, run))
        # Process trajectories
        distances = list()
        for clone in range(rank, nclones, size):
            # Read trajectory

            h5_filename = os.path.join(project_basepath, '%s/run%d-clone%d.h5' % (condition, run, clone))
            t = md.load(h5_filename)
            min_frame = 400
            end_frame = t.n_frames
            print('  CLONE %5d : %5d frames' % (clone, t.n_frames))

            short_traj = t.slice(range(min_frame, end_frame), copy=False)
            [distance, res_list] = md.compute_contacts(short_traj, scheme='ca')
            distances.append(distance)

        # Gather data to root and write it
        gathered_dist = MPI.COMM_WORLD.gather(distances, root=0)
        if rank == 0:
            output_filename = os.path.join(output_basepath, '%s-run%d-contacts.npy' % (save_file_label[i], run))
            np.save(output_filename, gathered_dist)



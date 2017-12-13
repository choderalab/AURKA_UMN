#!/usr/bin/python

import numpy as np
import sys
import math
import os
import time
import mdtraj as md

"""
Analyze A Loop

"""

project_basepath = '/cbio/jclab/projects/fah/fah-data/munged4' # location of FAH trajectories
output_basepath = '../data/rmsd'

#if not os.path.exists(output_basepath):
#    os.makedirs(output_basepath)

# Residue numbering corrected by offset for the reference structure
offset = 0
start_alignment = 123 - offset
end_alignment = 387 - offset
start_rmsd = 286 - offset
end_rmsd = 293 - offset

# Load reference structure for comparison of alphaC RMSD
reference_pdbfile = '../../wt_charmm_setup/input/1ol5-prepped.pdb'
reference = md.load(reference_pdbfile)
alignment_selection_dsl = '(resSeq >= %s) and (resSeq <= %s) and (name CA)' % (start_alignment, end_alignment)
alignment_reference_indices = reference.topology.select(alignment_selection_dsl)
rmsd_selection_dsl = '(resSeq >= %s) and (resSeq <= %s) and (name CA)' % (start_rmsd, end_rmsd)
rmsd_reference_indices = reference.topology.select(rmsd_selection_dsl)

nclones = 100 # number of CLONEs per RUN
nframes = 4000 # max frames / trajectory
projects = ['11432']

nruns = 4 # number of runs per condition
for project in projects:
    for run in range(nruns):
        h5_filename = os.path.join(project_basepath, '%s/run%d-clone%d.h5' % (project, run, 0))
        if not os.path.exists(h5_filename):
            continue

        # Process trajectories
        rmsds = list()
        for clone in range(nclones):
            rmsds = list()
            # Read trajectory
            initial_time = time.time()
            h5_filename = os.path.join(project_basepath, '%s/run%d-clone%d.h5' % (project, run, clone))
            t = md.load(h5_filename)
            print('  CLONE %5d : %5d frames' % (clone, t.n_frames))

            # Align to reference using AurA CA atoms only
            # Residue numbering corrected by the offset for the trajectories
            offset = 122
            start_alignment = 123 - offset
            end_alignment = 387 - offset
            start_rmsd = 286 - offset
            end_rmsd = 293 - offset
            alignment_selection_dsl = '(resSeq >= %s) and (resSeq <= %s) and (name CA)' % (start_alignment,
                                                                                           end_alignment)
            rmsd_selection_dsl = '(resSeq >= %s) and (resSeq <= %s) and (name CA)' % (start_rmsd, end_rmsd)
            alignment_trajectory_indices = t.topology.select(alignment_selection_dsl)
            t.superpose(reference, atom_indices=alignment_trajectory_indices,
                        ref_atom_indices=alignment_reference_indices, parallel=False)
            # Compute the RMSD manually without additional alignment
            rmsd_trajectory_indices = t.topology.select(rmsd_selection_dsl)
            rmsd = np.sqrt(
                3 * np.mean((t.xyz[:, rmsd_trajectory_indices, :] - reference.xyz[:, rmsd_reference_indices, :]) ** 2,
                            axis=(1, 2)))
            rmsds.append(rmsd)

            # Gather data to root and write it
            output_filename = os.path.join(output_basepath, '%s-run%d-fullaloop-rmsd_1OL5.npy' % (project, run))
            np.save(output_filename, rmsds)


        


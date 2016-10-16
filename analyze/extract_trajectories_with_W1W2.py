#!/usr/bin/python

import numpy as np
import sys
import math
import os
import time
import mdtraj as md

"""
Extract all trajectories and associated waters meeting W1/W2 definition for that trajectory.
Requires /cbio/jclab/projects/behrj/AURKA_UMN/output-{project#}/data/{project#}_{run#}_{residue#}_distHBonds.npy files

"""

from mpi4py import MPI
rank = MPI.COMM_WORLD.rank
size = MPI.COMM_WORLD.size

if rank==0: print('rank = %d, size = %d' % (rank, size))

project_basepath = '/cbio/jclab/projects/fah/fah-data/munged3/all-atoms' # location of FAH all-atom trajectories
hbond_basepath = '/cbio/jclab/projects/behrj/AURKA_UMN/'
output_basepath = 'trajectories-with-key-waters'

if not os.path.exists(output_basepath):
    os.makedirs(output_basepath)

nclones = 50 # number of CLONEs per RUN
projects = ['11414', '11419', '11418', '11423']
nruns = 7 # number of runs per project
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
            # Get all waters that satisfy W1 or W2 in this trajectory
            W1_set = set.union(*[s for s in W1[clone,:] if s is not None])
            W2_set = set.union(*[s for s in W2[clone,:] if s is not None])
            W_set = set.union(W1_set, W2_set)
            print('  CLONE %5d : W1 %8d W2 %8d : total %8d' % (clone, len(W1_set), len(W2_set), len(W_set)))

            # Read trajectory
            h5_filename = os.path.join(project_basepath, '%s/run%d-clone%d.h5' % (project, run, clone))
            t = md.load(h5_filename)

            # Slice the trajectory to retain only key waters
            molecules = t.topology.find_molecules()
            atoms_to_keep = [ a.index for a in t.topology.atoms if a.residue.is_protein ]
            for m in molecules:
                atom_indices = set([ a.index for a in m ])
                if len(atom_indices.intersection(W_set)) > 0:
                    atoms_to_keep += atom_indices
                    break
            t.restrict_atoms(atoms_to_keep)

            # Align the trajectory to AurA
            AurA_atoms = [ a.index for a in t.topology.atoms if ((a.residue.resSeq >= 123) and (a.residue.resSeq <= 388)) ]
            t.superpose(t, frame=0, atom_indices=AurA_atoms, parallel=True)

            # Image the trajectory around the protein.
            molecules = t.topology.find_molecules()
            t.image_molecules(inplace=True, anchor_molecules=molecules[0])

            # Write PDB and DCD
            pdb_filename = os.path.join(output_basepath, 'project%s-run%d-clone%d.pdb' % (project, run, clone))
            t[0].save(pdb_filename)
            dcd_filename = os.path.join(output_basepath, 'project%s-run%d-clone%d.dcd' % (project, run, clone))
            t.save(dcd_filename)
            
            MPI.COMM_WORLD.Barrier()

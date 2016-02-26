from msmbuilder import dataset
import mdtraj as md
import numpy as np
import math
import os
import copy

projects = ['11410','11411']
project_dirs = {'11410':'/cbio/jclab/projects/behrj/AURKA_UMN/output-1OL5','11411':'/cbio/jclab/projects/behrj/AURKA_UMN/output-1OL7'}
system = {'11410':'with TPX2','11411':'without TPX2'}

with open('/cbio/jclab/projects/behrj/AURKA_UMN/output-1OL7/run-index.txt','r') as fi:
    run_index = fi.read()

mutant = dict()
for entry in run_index.split('\n'):
    try:
        mutant[entry.split(' ')[0]] = entry.split(' ')[1]
    except:
        pass

for project in projects:
    project_dir = project_dirs[project]
    for run in range(5):
        HB_total = np.load('%s/data/%s_%s_HBonds.npy' % (project_dir, project, run))

        if os.path.exists('%s/data/%s_RUN%s_water_indices.npy' % (project_dir, project, run)):
            identified_hydrogens = np.load('%s/data/%s_RUN%s_water_indices.npy' % (project_dir, project, run))
            identified_hydrogens = identified_hydrogens.tolist()
            identified_hydrogens = set(identified_hydrogens)
        else:
            identified_hydrogens = set()
            for clone_id, hbonds_traj in enumerate(HB_total):
                for frame in hbonds_traj:
                    num_hbonds_in_frame_in_traj = frame.shape[0]
                    if num_hbonds_in_frame_in_traj > 0:
                        for bond_in_frame in frame:
                            identified_hydrogens.add(bond_in_frame[1])
            np.save('%s/data/%s_RUN%s_water_indices.npy' % (project_dir, project, run), identified_hydrogens)

        trajectories = dataset.MDTrajDataset("/cbio/jclab/projects/fah/fah-data/munged2/all-atoms/%s/run%d-clone*.h5" % (project, run))
        for clone, traj in enumerate(trajectories):
            if clone != 0 and os.path.exists('/cbio/jclab/projects/AURKA_UMN/%s_RUN%s_clone%s.dcd' % (project, run, clone)):
                continue
            if clone==0:
                important_atoms = copy.deepcopy(identified_hydrogens)
                [important_atoms.add(i) for i in traj.top.select("protein or resname MOL").tolist()]
                important_atoms = list(important_atoms)
            traj.atom_slice(important_atoms, inplace=True)
            if clone==0:
                if os.path.exists('/cbio/jclab/projects/AURKA_UMN/%s_RUN%s_clone%s.dcd' % (project, run, clone)):
                    continue
                first_frame = traj[0]
                try:
                    first_frame.save_pdb('/cbio/jclab/projects/AURKA_UMN/%s_RUN%s.pdb' % (project, run))
                except:
                    pass
            traj.save_dcd('/cbio/jclab/projects/AURKA_UMN/%s_RUN%s_clone%s.dcd' % (project, run, clone))

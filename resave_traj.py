from msmbuilder import dataset
import mdtraj as md
import numpy as np
import math
import os
import copy
import sys

projects = ['11410','11411']
project_dirs = {'11410':'/cbio/jclab/projects/behrj/AURKA_UMN/output-1OL5','11411':'/cbio/jclab/projects/behrj/AURKA_UMN/output-1OL7'}
system = {'11410':'with TPX2','11411':'without TPX2'}
runs = range(5)

with open('/cbio/jclab/projects/behrj/AURKA_UMN/output-1OL7/run-index.txt','r') as fi:
    run_index = fi.read()

mutant = dict()
for entry in run_index.split('\n'):
    try:
        mutant[entry.split(' ')[0]] = entry.split(' ')[1]
    except:
        pass

try:
    this_run = int(sys.argv[1]) # 1 - 10
except:
    this_run = None

if this_run is not None:
    projects = [projects[this_run%2]]
    runs = [runs[this_run%5]]

for project in projects:
    project_dir = project_dirs[project]
    for run in runs:
        HB_total = np.load('%s/data/%s_%s_HBonds.npy' % (project_dir, project, run))

        if os.path.exists('%s/data/%s_RUN%s_water_indices.npy' % (project_dir, project, run)):
            identified_waters = np.load('%s/data/%s_RUN%s_water_indices.npy' % (project_dir, project, run))
            identified_waters = identified_waters.tolist()
            identified_waters = set(identified_waters)
        else:
            identified_waters = set()
            for clone_id, hbonds_traj in enumerate(HB_total):
                for frame in hbonds_traj:
                    num_hbonds_in_frame_in_traj = frame.shape[0]
                    if num_hbonds_in_frame_in_traj > 0:
                        for bond_in_frame in frame:
                            identified_waters.add(bond_in_frame[0])
                            identified_waters.add(bond_in_frame[1])
                            identified_waters.add(bond_in_frame[2])
            np.save('%s/data/%s_RUN%s_water_indices.npy' % (project_dir, project, run), identified_waters)

        trajectories = dataset.MDTrajDataset("/cbio/jclab/projects/fah/fah-data/munged2/all-atoms/%s/run%d-clone*.h5" % (project, run))
        for clone, traj in enumerate(trajectories):
            if clone != 0 and os.path.exists('/cbio/jclab/projects/AURKA_UMN/%s_RUN%s_clone%s.dcd' % (project, run, clone)):
                continue
            if clone==0:
                important_atoms = copy.deepcopy(identified_waters)
                [important_atoms.add(i) for i in traj.top.select("protein or resname MOL").tolist()]
                important_atoms = list(important_atoms)
            traj.atom_slice(important_atoms, inplace=True)
            if clone==0:
                first_frame = traj[0]
                try:
                    first_frame.save_pdb('/cbio/jclab/projects/AURKA_UMN/%s_RUN%s.pdb' % (project, run))
                except:
                    pass
            traj.save_dcd('/cbio/jclab/projects/AURKA_UMN/%s_RUN%s_clone%s.dcd' % (project, run, clone))

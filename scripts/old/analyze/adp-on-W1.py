"""
Create true/false filter for whether ADP is bound in a given frame
"""
import numpy as np
import sys
import math
import os
import mdtraj as md
from msmbuilder import dataset
from itertools import chain

local_path = os.path.dirname(os.path.realpath(__file__))

projects = ['11414','11418','11419','11423','11424','11425']
project_dirs = {
    '11410':'%s/../output-1OL5' % local_path,
    '11411':'%s/../output-1OL7' % local_path,
    '11414':'%s/../output-11414' % local_path,
    '11418':'%s/../output-11418' % local_path,
    '11419':'%s/../output-11419' % local_path,
    '11423':'%s/../output-11423' % local_path,
    '11424':'%s/../output-11424' % local_path,
    '11425':'%s/../output-11425' % local_path,
}

system = {
    '11410':'with TPX2',
    '11411':'without TPX2',
    '11414':'with TPX2',
    '11418':'with TPX2 removed',
    '11419':'with TPX2',
    '11423':'with TPX2 removed',
    '11424':'with TPX2; charmm',
    '11425':'with TPX2 removed; charmm'
}

run_guide = list()
num_runs = dict()
for project in projects:
    num_runs[project] = 0
    filename = project_dirs[project]+'/run-index.txt'
    with open(filename, 'r') as fi:
        project_run_index = fi.read()
    for entry in project_run_index.split('\n'):
        try:
            run = entry.split(' ')[0]
            mutant = entry.split(' ')[1]
            run = int(run[3:])
            run_guide.append((project, run))
            num_runs[project] += 1
        except:
            pass
print(run_guide)
print(len(run_guide))


try:
    this_run = int(sys.argv[1]) # 1-48
except:
    this_run = None
if this_run is not None:
    run_guide = run_guide[this_run%48]
    projects = [run_guide[0]]
    runs = [run_guide[1]]
else:
    runs = None

def find_hbonds_for_this_traj(traj, residue, haystack, sidechain=False, backbone=False):
    if sidechain:
        residue_atoms = [atom.index for atom in residue.atoms if atom.is_sidechain]
    elif backbone:
        residue_atoms = [atom.index for atom in residue.atoms if atom.is_backbone]
    else:
        residue_atoms = [atom.index for atom in residue.atoms]
    neighbor_set = find_neighbor_set(traj, residue_atoms, haystack)
    if verbose:
        print('Finding hbonds for run%s clone%s' % (run, i))
    hbonds0 = md.wernet_nilsson(traj, exclude_water=False, proposed_donor_indices=residue_atoms, proposed_acceptor_indices=neighbor_set)
    hbonds1 = md.wernet_nilsson(traj, exclude_water=False, proposed_donor_indices=neighbor_set, proposed_acceptor_indices=residue_atoms)
    hbonds = list()
    for frame, bondlist in enumerate(hbonds0):
        try:
            hbonds.append(np.concatenate((bondlist,hbonds1[frame])))
        except Exception as e:
            print('hbonds0')
            print(bondlist)
            print(bondlist.shape)
            print('hbonds1')
            print(hbonds1[frame])
            print(hbonds1[frame].shape)
            raise(e)
    return hbonds

def find_neighbor_set(traj, residue, haystack):
    neighbors = md.compute_neighbors(traj, 0.4, residue, haystack_indices=haystack)
    neighbor_set = set(chain.from_iterable(neighbors))
    return list(neighbor_set)

verbose = True
for project in projects:
    project_dir = project_dirs[project]
    if this_run is None:
        runs = range(num_runs[project])
    for run in runs:
        HB_adp_water = list()
        if verbose:
            print("Loading Project %s RUN%s..." % (project, run))
        trajectories = dataset.MDTrajDataset("/cbio/jclab/projects/fah/fah-data/munged3/all-atoms/%s/run%s-clone*.h5" % (project, run))
        for i,traj in enumerate(trajectories):
            if i == 0:
                for residue in traj.topology.residues:
                    if str(residue).startswith('MOL') or str(residue).startswith('ADP'):
                        adp = residue
                        break
                haystack = traj.top.select("water")
            HB_adp_water.append(find_hbonds_for_this_traj(traj, adp, haystack))
        np.save('%s/data/%s_%s_adp-hbond-water.npy' % (project_dir, project, run), HB_adp_water)


if verbose:
    print('Complete!')

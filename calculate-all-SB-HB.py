import numpy as np
import sys
import math
import os
import mdtraj as md
from msmbuilder import dataset
from itertools import chain

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


verbose = False
for project in projects:
    project_dir = project_dirs[project]
    for run in runs:
        if verbose:
            print("Loading RUN%s..." % run)
        SB_eq_total = []
        SB_ek_total = []
        HB_185_total = []
        HB_181_total = []
        HB_275_total = []
        trajectories = dataset.MDTrajDataset("/cbio/jclab/projects/fah/fah-data/munged2/all-atoms/%s/run%d-clone*.h5" % (project, run))
        for i,traj in enumerate(trajectories):
            if i == 0:
                for residue in traj.topology.residues:
                    if str(residue) == 'GLU181':
                        e181 = residue
                    if str(residue) == 'GLN185':
                        q185 = residue
                    if str(residue) == 'LYS162':
                        k162 = residue
                    if str(residue) == 'PHE275':
                        f275 = residue
            res185 = traj.topology.residue(q185.index)

            distances_181_185, residue_pairs = md.compute_contacts(traj, contacts=[[e181.index,q185.index]])
            distances_181_162, residue_pairs = md.compute_contacts(traj, contacts=[[e181.index,k162.index]])
            SB_eq_total.append(distances_181_185[:,0])
            SB_ek_total.append(distances_181_162[:,0])

            res185atoms = [atom.index for atom in res185.atoms]
            haystack = traj.top.select("water")
            neighbors = md.compute_neighbors(traj, 0.3, res185atoms, haystack_indices=haystack)

            neighbor_set = set(chain.from_iterable(neighbors))
            neighbor_set = list(neighbor_set)
            np.save('%s/data/%s_RUN%s_clone%s_water_indices.npy' % (project_dir, project, run, i), neighbor_set)

            # using wernet_nilsson because the output is hbonds per frame
            hbonds0 = md.wernet_nilsson(traj, exclude_water=False, proposed_donor_indices=res185atoms, proposed_acceptor_indices=neighbor_set)
            hbonds1 = md.wernet_nilsson(traj, exclude_water=False, proposed_donor_indices=neighbor_set, proposed_acceptor_indices=res185atoms)
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
            HB_185_total.append(hbonds)

            res181atoms = [atom.index for atom in e181.atoms]
            hbonds0 = md.wernet_nilsson(traj, exclude_water=False, proposed_donor_indices=res181atoms, proposed_acceptor_indices=neighbor_set)
            hbonds1 = md.wernet_nilsson(traj, exclude_water=False, proposed_donor_indices=neighbor_set, proposed_acceptor_indices=res181atoms)
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
            HB_181_total.append(hbonds)

            res275atoms = [atom.index for atom in f275.atoms]
            hbonds0 = md.wernet_nilsson(traj, exclude_water=False, proposed_donor_indices=res275atoms, proposed_acceptor_indices=neighbor_set)
            hbonds1 = md.wernet_nilsson(traj, exclude_water=False, proposed_donor_indices=neighbor_set, proposed_acceptor_indices=res275atoms)
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
            HB_275_total.append(hbonds)


        np.save('%s/data/%s_%s_181_HBonds.npy' % (project_dir, project, run), HB_181_total)
        np.save('%s/data/%s_%s_185_HBonds.npy' % (project_dir, project, run), HB_185_total)
        np.save('%s/data/%s_%s_275_HBonds.npy' % (project_dir, project, run), HB_275_total)
        np.save('%s/data/%s_%s_181-185_SB_total.npy' % (project_dir, project, run), SB_eq_total)
        np.save('%s/data/%s_%s_181-162_SB_total.npy' % (project_dir, project, run), SB_ek_total)



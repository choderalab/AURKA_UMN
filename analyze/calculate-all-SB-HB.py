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

DIST_FOR_HBOND = True

def find_hbondsdists_for_this_traj(traj, residue, haystack, sidechain=False, backbone=False):
    if sidechain:
        residue_atoms = [atom.index for atom in residue.atoms if (atom.is_sidechain and (str(atom.element)=='oxygen' or str(atom.element)=='nitrogen'))]
    elif backbone:
        residue_atoms = [atom.index for atom in residue.atoms if (atom.is_backbone and (str(atom.element)=='oxygen' or str(atom.element)=='nitrogen'))]
    else:
        residue_atoms = [atom.index for atom in residue.atoms if (str(atom.element)=='oxygen' or str(atom.element)=='nitrogen')]

    contact_size = (len(haystack)*len(residue_atoms),2)
    contacts = np.zeros(contact_size)
    index = 0
    for resatom in residue_atoms:
        for water in haystack:
            contacts[index][0] = resatom
            contacts[index][1] = water
            index+=1
    distances = md.compute_distances(traj, contacts)
    hbonds = list()
    for frame in distances:
        num_bonds = 0
        for index, little_dist in enumerate(frame):
            if little_dist < 0.35:
                num_bonds+=1
        local_hbonds = np.zeros((num_bonds,3))
        bond_index = 0
        for index, little_dist in enumerate(frame):
            if little_dist < 0.35:
                local_hbonds[bond_index][0] = contacts[index][0]
                local_hbonds[bond_index][1] = little_dist
                local_hbonds[bond_index][2] = contacts[index][1]
                bond_index+=1
        hbonds.append(local_hbonds)
    return hbonds

def find_hbonds_for_this_traj(traj, residue, haystack, sidechain=False, backbone=False):
    if sidechain:
        residue_atoms = [atom.index for atom in residue.atoms if atom.is_sidechain]
    elif backbone:
        residue_atoms = [atom.index for atom in residue.atoms if atom.is_backbone]
    else:
        residue_atoms = [atom.index for atom in residue.atoms]
    neighbor_set = find_neighbor_set(traj, residue_atoms, haystack)
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

try:
    this_run = int(sys.argv[1]) # 1 - 10
except:
    this_run = None

if this_run is not None:
    projects = [projects[this_run%2]]
    runs = [runs[this_run%5]]


verbose = True
for project in projects:
    project_dir = project_dirs[project]
    for run in runs:
        if verbose:
            print("Loading Project %s RUN%s..." % (project, run))
        SB_eq_total = []
        SB_ek_total = []
        HB_185_total = []
        HB_181_total = []
        HB_274_total = []
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
                    if str(residue) == 'ASP274':
                        d274 = residue
                    if str(residue) == 'PHE275':
                        f275 = residue
            res185 = traj.topology.residue(q185.index)

            distances_181_185, residue_pairs = md.compute_contacts(traj, contacts=[[e181.index,q185.index]],scheme='sidechain-heavy')
            distances_181_162, residue_pairs = md.compute_contacts(traj, contacts=[[e181.index,k162.index]],scheme='sidechain-heavy')
            SB_eq_total.append(distances_181_185[:,0])
            SB_ek_total.append(distances_181_162[:,0])

            if DIST_FOR_HBOND:
                haystack = traj.top.select("water and name O")
            else:
                haystack = traj.top.select("water")

            if DIST_FOR_HBOND:
                print('185')
                HB_185_total.append(find_hbondsdists_for_this_traj(traj, res185, haystack, sidechain=True))
                print('181')
                HB_181_total.append(find_hbondsdists_for_this_traj(traj, e181, haystack, sidechain=True))
                print('274')
                HB_274_total.append(find_hbondsdists_for_this_traj(traj, d274, haystack, backbone=True))
                print('275')
                HB_275_total.append(find_hbondsdists_for_this_traj(traj, f275, haystack, backbone=True))
            else:
                HB_185_total.append(find_hbonds_for_this_traj(traj, res185, haystack, sidechain=True))
                HB_181_total.append(find_hbonds_for_this_traj(traj, e181, haystack, sidechain=True))
                HB_274_total.append(find_hbonds_for_this_traj(traj, d274, haystack, backbone=True))
                HB_275_total.append(find_hbonds_for_this_traj(traj, f275, haystack, backbone=True))
        if DIST_FOR_HBOND:
            np.save('%s/data/%s_%s_181_distHBonds.npy' % (project_dir, project, run), HB_181_total)
            np.save('%s/data/%s_%s_185_distHBonds.npy' % (project_dir, project, run), HB_185_total)
            np.save('%s/data/%s_%s_274_distHBonds.npy' % (project_dir, project, run), HB_274_total)
            np.save('%s/data/%s_%s_275_distHBonds.npy' % (project_dir, project, run), HB_275_total)
        else:
            np.save('%s/data/%s_%s_181_HBonds.npy' % (project_dir, project, run), HB_181_total)
            np.save('%s/data/%s_%s_185_HBonds.npy' % (project_dir, project, run), HB_185_total)
            np.save('%s/data/%s_%s_274_HBonds.npy' % (project_dir, project, run), HB_274_total)
            np.save('%s/data/%s_%s_275_HBonds.npy' % (project_dir, project, run), HB_275_total)
        np.save('%s/data/%s_%s_181-185_SB_total.npy' % (project_dir, project, run), SB_eq_total)
        np.save('%s/data/%s_%s_181-162_SB_total.npy' % (project_dir, project, run), SB_ek_total)

if verbose:
    print('Complete!')

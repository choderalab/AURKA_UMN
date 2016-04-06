"""
Create true/false filter for whether ADP is bound in a given frame
"""
import numpy as np
import sys
import math
import os
import mdtraj as md
from itertools import chain
import plot_function

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

overwrite = False

which_pdb = 'system'

OFFSET = plot_function.OFFSET
OFFSET = 0

def find_hbonds_for_this_traj(traj, residue, haystack, sidechain=False, backbone=False):
    if sidechain:
        residue_atoms = [atom.index for atom in residue.atoms if atom.is_sidechain]
    elif backbone:
        residue_atoms = [atom.index for atom in residue.atoms if atom.is_backbone]
    else:
        residue_atoms = [atom.index for atom in residue.atoms]
    neighbor_set = set()
    for atom in haystack.atoms:
        if str(atom.element) == 'nitrogen' or str(atom.element) == 'hydrogen':
            neighbor_set.add(atom.index)   
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

def save_adp_status(little_distances, big_distances, hbonds, project_dir):
    adp_active = np.empty(len(little_distances),dtype=bool)
    for run, little_distance in enumerate(little_distances):
        print('Finding if ADP is bound in initial RUN%s' % run)
        this_lil_dist = little_distance
        this_big_dist = big_distances[run]
                #hbond_count = hbonds[clone][index].shape[0]
        hbond_dist = hbonds[run]
        if this_lil_dist < this_big_dist and hbond_dist < .40:
            adp_active = True
        else:
            adp_active = False
        print(adp_active)

try:
    this_project = int(sys.argv[1]) # 1 - 10
except:
    this_project = None

if this_project is not None:
    projects = [projects[this_project%2]]

project_adp_active = dict()

verbose = True

for project in projects:
    project_dir = project_dirs[project]
    a213n_total = list()
    a213c_total = list()
    e211nh2_total = list()
    for run in runs:
        if verbose:
            print("Loading Project %s RUN%s..." % (project, run))
        traj = md.load("%s/RUN%s/%s.pdb" % (project_dir, run, which_pdb))
        for residue in traj.topology.residues:
            if str(residue) == 'GLU211':
                e211 = residue
            if str(residue) == 'ALA213':
                a213 = residue
            if str(residue).startswith('MOL') or str(residue).startswith('ADP'):
                adp = residue
                break
        found = False
        for atom in e211.atoms:
            if atom.is_backbone and str(atom.element) == 'oxygen':
                e211_backbone_carbonyl = atom
                found = True
        assert found
        found = False
        for atom in a213.atoms:
            if atom.is_backbone and str(atom.element) == 'nitrogen':
                a213_backbone_amide = atom
                found = True
        assert found
        adp_nitrogens = dict()
        for atom in adp.atoms:
            if str(atom.element) == 'nitrogen':
                adp_nitrogens[atom] = list()
        for bond in traj.topology.bonds:
            if bond[0] in adp_nitrogens.keys():
                nitrogen = bond[0]
                atom = bond[1]
            elif bond[1] in adp_nitrogens.keys():
                nitrogen = bond[1]
                atom = bond[0]
            else:
                continue
            adp_nitrogens[nitrogen].append(atom)
        foundn4 = False
        foundn3 = False
        foundc7 = False
        for key, value in adp_nitrogens.items():
            if len(value) == 3 and 'hydrogen' in [str(atom.element) for atom in value]:
                nh2_group = key
                foundn3 = True
        for key, value in adp_nitrogens.items():
            if key == nh2_group:
                continue
            for atom in value:
                if atom in adp_nitrogens[nh2_group]:
                    n4_adp = key
                    foundn4 = True
                    for bond in traj.topology.bonds:
                        if atom == bond[0] and n4_adp not in bond and nh2_group not in bond:
                            c7_adp = bond[1]
                            foundc7 = True
                        elif atom == bond[1] and n4_adp not in bond and nh2_group not in bond:
                            c7_adp = bond[0]
                            foundc7 = True
                        else:
                            continue
        assert foundc7
        assert foundn4
        assert foundn3
#            distances, residue_pairs = md.compute_contacts(traj, contacts=[[a213.index,adp.index]],scheme='active-adp')
#            SB_a213_total.append(distances[:,0])
#            distances, residue_pairs = md.compute_contacts(traj, contacts=[[e211.index,adp.index]],scheme='active-adp')
#            HB_e211_total.append(distances[:,0])
        atom_pair_AN = np.zeros((1,2))
        atom_pair_AN[0,0] = a213_backbone_amide.index
        atom_pair_AN[0,1] = n4_adp.index
        atom_pair_AC = np.zeros((1,2))
        atom_pair_AC[0,0] = a213_backbone_amide.index
        atom_pair_AC[0,1] = c7_adp.index
        atom_pair_E = np.zeros((1,2))
        atom_pair_E[0,0] = e211_backbone_carbonyl.index
        atom_pair_E[0,1] = nh2_group.index
        a213n_total.append(md.compute_distances(traj, atom_pair_AN))
        a213c_total.append(md.compute_distances(traj, atom_pair_AC))
        e211nh2_total.append(md.compute_distances(traj, atom_pair_E))
    save_adp_status(a213n_total, a213c_total, e211nh2_total, project_dir)

if verbose:
    print('Complete!')

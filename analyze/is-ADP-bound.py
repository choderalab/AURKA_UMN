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
    neighbors = md.compute_neighbors(traj, 0.3, residue, haystack_indices=haystack)
    neighbor_set = set(chain.from_iterable(neighbors))
    return list(neighbor_set)

def save_adp_status(distances, hbonds, project_dir):
    adp_active = np.empty((5*50,2000),dtype=bool)
    for clone, distance in enumerate(distances):
        count = 0.0
        length = 0
        print('Finding how many ADPs are bound in RUN%s clone%s' % (clone/50, clone%50))
        for index in range(2000):
            try:
                this_dist = distance[index]
                #hbond_count = hbonds[clone][index].shape[0]
                hbond_dist = hbonds[clone][index]
                length = index + 1.0
            except:
                adp_active[clone][index] = False
                continue
            #if this_dist < 4.5 and hbond_count > 0:
            if this_dist < .45 and hbond_dist < .35:
                adp_active[clone][index] = True
                count += 1.0
            else:
                adp_active[clone][index] = False
        print(count/length)
    np.save('%s/is-ADP-bound.npy' % project_dir, adp_active)

try:
    this_project = int(sys.argv[1]) # 1 - 10
except:
    this_project = None

if this_project is not None:
    projects = [projects[this_project%2]]

verbose = True
for project in projects:
    project_dir = project_dirs[project]
    SB_a213_total = []
    HB_e211_total = []
    for run in runs:
        if verbose:
            print("Loading Project %s RUN%s..." % (project, run))
        trajectories = dataset.MDTrajDataset("/cbio/jclab/projects/fah/fah-data/munged2/all-atoms/%s/run%d-clone*.h5" % (project, run))
        for i,traj in enumerate(trajectories):
            if i == 0:
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
                for key, value in adp_nitrogens.items():
                    if len(value) == 3 and 'hydrogen' in [str(atom.element) for atom in value]:
                        nh2_group = key
                        foundn3 = True
                for key, value in adp_nitrogens.items():
                    if key == nh2_group:
                        continue
                    if any([atom in adp_nitrogens[nh2_group] for atom in value]):
                        n4_adp = key
                        foundn4 = True
                assert foundn4
                assert foundn3
#            distances, residue_pairs = md.compute_contacts(traj, contacts=[[a213.index,adp.index]],scheme='active-adp')
#            SB_a213_total.append(distances[:,0])
#            distances, residue_pairs = md.compute_contacts(traj, contacts=[[e211.index,adp.index]],scheme='active-adp')
#            HB_e211_total.append(distances[:,0])
            atom_pair_A = np.zeros((1,2))
            atom_pair_A[0,0] = a213_backbone_amide.index
            atom_pair_A[0,1] = n4_adp.index
            atom_pair_E = np.zeros((1,2))
            atom_pair_E[0,0] = e211_backbone_carbonyl.index
            atom_pair_E[0,1] = nh2_group.index
            SB_a213_total.append(md.compute_distances(traj, atom_pair_A))
            HB_e211_total.append(md.compute_distances(traj, atom_pair_E))
    save_adp_status(SB_a213_total, HB_e211_total, project_dir)

if verbose:
    print('Complete!')

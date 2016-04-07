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
import plot_function

local_path = os.path.dirname(os.path.realpath(__file__))

projects = ['11410','11411','11418']
projects = ['11410','11411']
project_dirs = {'11410':'%s/../output-1OL5' % local_path,'11411':'%s/../output-1OL7' % local_path,'11418':'%s/../output1OL5-TPX2' % local_path}
system = {'11410':'with TPX2','11411':'without TPX2','11418': 'with TPX2 removed'}
runs = range(5)

with open('/cbio/jclab/projects/behrj/AURKA_UMN/output-1OL7/run-index.txt','r') as fi:
    run_index = fi.read()

mutant = dict()
for entry in run_index.split('\n'):
    try:
        mutant[entry.split(' ')[0]] = entry.split(' ')[1]
    except:
        pass

overwrite = True

OFFSET = 0

class PO(object):
    def __init__(self, dist, index):
        self.distance = dist
        self.index = index

def plot_2dhist(key, x_axis, adp_count, weights, title, filename):
    ylabel = 'ADP mlc bound'
    plot_function.plot_2dhist(key, x_axis, adp_count, weights, title, ylabel, filename)

def plot_initial_adp(project_adp_active):
    bin_x = np.arange(6) - 0.25
    for project in project_adp_active.keys():
        ADPs = np.zeros(5*50) - 1
        x_axis = np.zeros(5*50)
        weights = np.zeros(5*50)
        project_dir = project_dirs[project]
        adp_active = project_adp_active[project]
        for clone, traj in enumerate(adp_active):
            x_axis[clone] = clone/50
            if traj[0]:
                ADPs[clone] = 1
            else:
                ADPs[clone] = 0
            weights[clone] = 1.00
        title = 'Is ADP bound in initial frame in WT AURKA %s' % system[project]
        filename = "/cbio/jclab/projects/behrj/AURKA_UMN/plots/ADP-frame0-binding-%s.png" % project
        key = 'ADP0'
        plot_2dhist(key, x_axis, ADPs, weights, title, filename)

def find_po_dist(traj, k162_sidechain_amineN, value):
    k162_index = k162_sidechain_amineN.index
#    hbonds = md.wernet_nilsson(traj, proposed_donor_indices=k162_indices, proposed_acceptor_indices=oxygen_indices)
    # construct a thing that still looks like hbonds? or that has 2 entries at least ok fine like a tuple (?)
    # dist, atom.index for chosen O
    atom_pairs_PO = np.zeros((len(value),2))
    for i, oxygen in enumerate(value):
        atom_pairs_PO[i,0] = k162_index
        atom_pairs_PO[i,1] = oxygen.index
    find_O = md.compute_distances(traj[0], atom_pairs_PO)[0]
    oxygen = value[find_O.argmin()]
    atom_pair_PO = np.zeros((1,2))
    atom_pair_PO[0,0] = k162_index
    atom_pair_PO[0,1] = oxygen.index
    distances = md.compute_distances(traj, atom_pair_PO) # ASSUME Os DON'T FLIP
    return PO(distances, oxygen.index)

def find_o_to_o(traj, k162op_alpha, k162op_beta):
    atom_pair_OO = np.zeros((1,2))
    atom_pair_OO[0,0] = k162op_alpha.index
    atom_pair_OO[0,1] = k162op_beta.index
    return md.compute_distances(traj, atom_pair_OO)

def plot_adp_active(project_adp_active):
    bin_x = np.arange(OFFSET/4,510,10) - 0.25
    for project in project_adp_active.keys():
        ADPs = np.zeros((5*50,2000-OFFSET)) - 1
        x_axis = np.zeros((5*50,2000-OFFSET))
        weights = np.zeros((5*50,2000-OFFSET))
        column_count = np.zeros(bin_x.shape)
        project_dir = project_dirs[project]
        adp_active = project_adp_active[project]
        residue = 185
        for run in range(5):
            if run == 0:
                frame_reference = np.load('%s/data/%s_%s_%s_HBonds.npy' % (project_dir, project, run, residue))
            else:
                new_run = np.load('%s/data/%s_%s_%s_HBonds.npy' % (project_dir, project, run, residue))
                frame_reference = np.concatenate((frame_reference, new_run))
        for clone, traj in enumerate(adp_active):
            for index in range(OFFSET,2000):
                x_axis[clone][index-OFFSET] = index*0.25
                try:
                    frame_exists = frame_reference[clone][index]
                    column_count[(index-OFFSET-0.25)/40] += 1
                except:
                    continue
                if traj[index]:
                    ADPs[clone][index-OFFSET] = 1
                else:
                    ADPs[clone][index-OFFSET] = 0
        print(column_count)
        for clone, traj in enumerate(adp_active):
            for index in range(OFFSET,2000):
                weights[clone][index-OFFSET] = 1.00 / column_count[(index-OFFSET-0.25)/40]
        x_axis = x_axis.flatten()
        ADPs = ADPs.flatten()
        weights = weights.flatten()
        title = 'Frames with ADP bound in WT AURKA %s' % system[project]
        filename = "/cbio/jclab/projects/behrj/AURKA_UMN/plots/ADP-binding-%s.png" % project
        key = 'ADP'
        plot_2dhist(key, x_axis, ADPs, weights, title, filename)

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

def save_adp_status(little_distances, big_distances, hbonds, k162op, project_dir):
    adp_active = np.empty((len(little_distances),2000),dtype=bool)
    for clone, little_distance in enumerate(little_distances):
        count = 0.0
        length = 0
        print('Finding how many ADPs are bound in RUN%s clone%s' % (clone/50, clone%50))
        for index in range(2000):
            try:
                this_lil_dist = little_distance[index]
                this_big_dist = big_distances[clone][index]
                #hbond_count = hbonds[clone][index].shape[0]
                hbond_dist = hbonds[clone][index]
                length = index + 1.0
            except:
                adp_active[clone][index] = False
                continue
            #if this_dist < 4.5 and hbond_count > 0:
            k162p1 = min([k162op[i][clone].distance[index] for i in range(2)])
            k162p2 = max([k162op[i][clone].distance[index] for i in range(2)])
            o_to_o = k162op[2][clone][index]
            allowed_p2 = (o_to_o**2 + k162p1**2)**(1.0/2.0)
            if (this_lil_dist < this_big_dist and
                hbond_dist < .40 and
                k162p1 < .33 and
                k162p2 < allowed_p2):
                adp_active[clone][index] = True
                count += 1.0
            else:
                adp_active[clone][index] = False
        print(count/length)
    np.save('%s/is-ADP-bound.npy' % project_dir, adp_active)
    project_adp_active[project] = adp_active

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
    if not overwrite and os.path.exists('%s/is-ADP-bound.npy' % project_dir):
        project_adp_active[project] = np.load('%s/is-ADP-bound.npy' % project_dir)
        continue
    a213n_total = []
    a213c_total = []
    e211nh2_total = []
    k162op_total = [list(),list(), list()] # N-Oa, N-Ob, Oa-Ob
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
                    if str(residue) == 'LYS162':
                        k162 = residue
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
                found = False
                for atom in k162.atoms:
                    if atom.is_sidechain and str(atom.element) == 'nitrogen':
                        k162_sidechain_amineN = atom
                        k162_sidechain_aminegroup = [atom] # in case we switch back to hbond from dist
                        found = True
                assert found
                adp_oxygens = dict()
                adp_nitrogens = dict()
                oxygens = 0
                for atom in adp.atoms:
                    if str(atom.element) == 'nitrogen':
                        adp_nitrogens[atom] = list()
                    if str(atom.element) == 'phosphorus':
                        adp_oxygens[atom] = list()
                for bond in traj.topology.bonds:
                    if bond[0] in adp_nitrogens.keys():
                        adp_nitrogens[bond[0]].append(bond[1])
                    elif bond[1] in adp_nitrogens.keys():
                        adp_nitrogens[bond[1]].append(bond[0])
                    elif bond[0] in adp_oxygens.keys():
                        adp_oxygens[bond[0]].append(bond[1])
                        oxygens+=1
                    elif bond[1] in adp_oxygens.keys():
                        adp_oxygens[bond[1]].append(bond[0])
                        oxygens+=1
                    elif k162_sidechain_amineN in bond:
                        [k162_sidechain_aminegroup.append(atom) for atom in bond if str(atom.element) == 'hydrogen']
                assert oxygens == 8
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
            for p, value in enumerate(adp_oxygens.values()):
                k162op_total[p].append(find_po_dist(traj, k162_sidechain_amineN, value))
            k162op_total[2].append(find_o_to_o(traj, k162op_total[0][-1], k162op_total[1][-1]))
    save_adp_status(a213n_total, a213c_total, e211nh2_total, k162op_total, project_dir)
plot_initial_adp(project_adp_active)
plot_adp_active(project_adp_active)

if verbose:
    print('Complete!')

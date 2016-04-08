"""
Create true/false filter for whether ADP is bound in a given frame
"""
import numpy as np
import sys
import math
import os
import mdtraj as md
from itertools import chain

local_path = os.path.dirname(os.path.realpath(__file__))

projects = ['11410','11411','11418']
project_dirs = {'11410':'%s/../output-1OL5' % local_path,'11411':'%s/../output-1OL7' % local_path,'11418':'%s/../output1OL5-TPX2' % local_path}
system = {'11410':'with TPX2','11411':'without TPX2','11418': 'with TPX2 removed'}
runs = range(5)

overwrite = False

which_pdb = 'system'

class PO(object):
    def __init__(self, dist, index):
        self.distance = dist
        self.index = index

def find_po_dist(traj, k162_sidechain_amineN, value):
    k162_index = k162_sidechain_amineN.index
#    hbonds = md.wernet_nilsson(traj, proposed_donor_indices=k162_indices, proposed_acceptor_indices=oxygen_indices)
    # construct a thing that still looks like hbonds? or that has 2 entries at least ok fine like a tuple (?)
    # dist, atom.index for chosen O
    atom_pairs_PO = np.zeros((len(value),2))
    for i, oxygen in enumerate(value):
        atom_pairs_PO[i,0] = k162_index
        atom_pairs_PO[i,1] = oxygen.index
    distances = md.compute_distances(traj, atom_pairs_PO)[0]
    dist = min(distances) # only 1 frame
    oxygen = value[distances.argmin()]
    return PO(dist, oxygen.index)

def find_o_to_o(traj, k162op_alpha, k162op_beta):
    atom_pair_OO = np.zeros((1,2))
    atom_pair_OO[0,0] = k162op_alpha.index
    atom_pair_OO[0,1] = k162op_beta.index
    return md.compute_distances(traj, atom_pair_OO)

def save_adp_status(little_distances, big_distances, hbonds, k162op, project_dir):
    adp_active = np.empty(len(little_distances),dtype=bool)
    for run, little_distance in enumerate(little_distances):
        print('Finding if ADP is bound in initial RUN%s' % run)
        this_lil_dist = little_distance
        this_big_dist = big_distances[run]
        hbond_dist = hbonds[run]
        k162p1 = min([k162op[i][run].distance for i in range(2)])
        k162p2 = max([k162op[i][run].distance for i in range(2)])
        o_to_o = k162op[2][run]
        allowed_p2 = (o_to_o**2 + k162p1**2)**(1.0/2.0)
        if (this_lil_dist < this_big_dist and
            hbond_dist < .40 and
            k162p1 < .33 and
            k162p2 < allowed_p2):
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
    k162op_total = [list(),list(), list()] # N-Oa, N-Ob, Oa-Ob
    for run in runs:
        if verbose:
            print("Loading Project %s RUN%s..." % (project, run))
        traj = md.load("%s/RUN%s/%s.pdb" % (project_dir, run, which_pdb))
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
                k162_sidechain_aminegroup = [atom]
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

if verbose:
    print('Complete!')

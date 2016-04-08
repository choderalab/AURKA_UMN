import numpy as np
import sys
import math
import os
import mdtraj as md
import plot_function
from itertools import chain

        # make plots of all data past t = 250ns
            # quantify how much P(salt bridge) and (1-P)
                # P = Ntraj_frames_bound

        # look at trajectories : how long do waters stay in place
            # may need to account for the waters exchanging

#OFFSET = plot_function.OFFSET
OFFSET = 0

residues_with_H = [185,181,274,275]
reference = 185

local_path = os.path.dirname(os.path.realpath(__file__))

projects = ['11410','11411','11418']
projects = ['11410','11411']
project_dirs = {'11410':'%s/../output-1OL5' % local_path,'11411':'%s/../output-1OL7' % local_path,'11418':'%s/../output1OL5-TPX2' % local_path}
system = {'11410':'with TPX2','11411':'without TPX2','11418': 'with TPX2 removed'}

try:
    this_project = int(sys.argv[1]) # 1 - 10
except:
    this_project = None

if this_project is not None:
    projects = [projects[this_project%2]]

class bondlife(object):
    def __init__(self, oxygen1_index, oxygen2_index, birth_frame, clone):
        self.oxygens = [oxygen1_index, oxygen2_index]
        o1 = min([oxygen1_index, oxygen2_index])
        o2 = max([oxygen1_index, oxygen2_index])
        self.name = 'CLONE'+str(clone)+'-'+str(o1)+'-'+str(o2)
        self.frames = [birth_frame]
        self.clone = clone
        self._longest_life = len(self.frames)
        self._current_life = 1

    def add_frame(self, frame_id):
        if frame_id in self.frames:
            return
        if self.frames[-1] == frame_id - 1:
            self._current_life += 1
        else:
            self._current_life = 1
        if self._current_life > self._longest_life:
            self._longest_life = self._current_life
        self.frames.append(frame_id)

    @property
    def age(self):
        return self._longest_life

    @property
    def lifetime(self):
        return len(self.frames)

def water_set(topology, frame, hydrogens=False):
    waters = set()

    [[waters.add(int(atom)) for atom in [bond[0],bond[2]]
      if topology.top.atom(int(atom)).residue.is_water]
     for bond in frame]
    if not hydrogens:
        return waters
    oxygens = list(waters)
    [[waters.add(atom.index) for atom in topology.top.atom(oatom).residue.atoms]
     for oatom in oxygens]
    return waters

def find_other_waters(frame, first_waters):
    waters = set()

    for bond in frame:
        if int(bond[0]) in first_waters:
            waters.add(int(bond[2]))
        if int(bond[2]) in first_waters:
            waters.add(int(bond[0]))
    return waters

def build_water_cloud(waters, traj_slice):
    neighbors = md.compute_neighbors(traj_slice, 0.4, waters)
    water_cloud = set(chain.from_iterable(neighbors))
    return list(water_cloud)

def find_hbonds_between_waters(HB_total):
    if os.path.exists('%s/data/%s_%s_tryagain-IntraWater.npy' % (project_dir, project, 0)):
        WB_total = np.load('%s/data/%s_%s_tryagain-IntraWater.npy' % (project_dir, project, 0))
        return WB_total
    HB_res_total = HB_total[274]
    WB_total = list()

    for clone, traj in enumerate(HB_res_total):
        trajectory = md.load("/cbio/jclab/projects/fah/fah-data/munged2/all-atoms/%s/run%s-clone%s.h5" % (project, clone/50, clone%50))
        if clone%50 == 0:
            print('Now loading trajectories for RUN%s' % str(clone/50))
            topology = md.load('/cbio/jclab/projects/AURKA_UMN/trajectories/%s_RUN%s.pdb' % (project, 0))
        old_chunk = 0
        hbonds = None
        print('Finding RUN%s clone%s hbonds' % (str(clone/50),str(clone%50)))
        for chunk in range(100,2100,100):
            if chunk > 2000:
                break
            waters = set()
            chunk_frames = range(old_chunk,chunk)
            for index in chunk_frames:
                try:
                    frame_274 = traj[index]
                    frame_185 = HB_total[185][clone][index]
                    frame_181 = HB_total[181][clone][index]
                except:
                    chunk_frames = range(old_chunk,index)
                    break
                waters = waters.union(water_set(topology, frame_274, hydrogens=True))
                waters = waters.union(water_set(topology, frame_185, hydrogens=True))
                waters = waters.union(water_set(topology, frame_181, hydrogens=True))
            traj_slice = trajectory.slice(chunk_frames, copy=False)
            water_cloud = build_water_cloud(waters, traj_slice)
            hbond_chunk = md.wernet_nilsson(traj_slice, exclude_water=False, proposed_donor_indices=water_cloud, proposed_acceptor_indices=water_cloud)
            if hbonds is None:
                hbonds = hbond_chunk
            else:
                try:
                    hbonds.extend(hbond_chunk)
                except Exception as e:
                    print(len(hbonds))
                    print(len(hbond_chunk))
                    raise(e)
            old_chunk = chunk
        WB_total.append(hbonds)
        del(trajectory)
    np.save('%s/data/%s_%s_tryagain-IntraWater.npy' % (project_dir, project, 0),WB_total)
    return WB_total

def find_anchors(topology):
    anchors = [atom.index 
               for atom in topology.atoms 
               if (atom.residue.is_protein or 
                   atom.residue.name in ['MOL','ADP'])]
    return anchors

WB_dict = dict()
for i, project in enumerate(projects):
    project_dir = project_dirs[project]
    ADP_bound = np.load('%s/is-ADP-bound.npy' % project_dir)
    HB_total = dict()
    if os.path.exists('%s/data/%s_%s_tryagain-IntraWater.npy' % (project_dir, project, 0)):
        WB_total = find_hbonds_between_waters(None)
        WB_dict[project] = WB_total
        continue
    for residue in residues_with_H:
        for run in range(5):
            if not os.path.exists('%s/data/%s_%s_%s_distHBonds.npy' % (project_dir, project, run, residue)):
                continue
            if not HB_total.has_key(residue):
                HB_total[residue] = np.load('%s/data/%s_%s_%s_distHBonds.npy' % (project_dir, project, run, residue))
            else:
                new_run = np.load('%s/data/%s_%s_%s_distHBonds.npy' % (project_dir, project, run, residue))
                HB_total[residue] = np.concatenate((HB_total[residue], new_run))
    WB_total = find_hbonds_between_waters(HB_total)
    WB_dict[project] = WB_total

anchors = dict()
for project, WB_total in WB_dict.values():
    static_positions = dict()
    position_definitions = list()
    if os.path.exists('%s/data/%s_%s_1234-IntraWater-sig.npy' % (project_dir, project, 0)):
        long_bonds = np.load('%s/data/%s_%s_1234-IntraWater-sig.npy' % (project_dir, project, 0))
    else:
        print('Tracking bonds in project %s' % project)
        bond_lifes = dict()
        for clone, traj in enumerate(WB_total):
            for frame_id, frame in enumerate(traj):
                for bond in frame:
                    min_O = min([bond[0],bond[2]])
                    max_O = max([bond[0],bond[2]])
                    name = 'CLONE'+str(clone)+'-'+str(min_O)+'-'+str(max_O)
                    if name in bond_lifes.keys():
                        this_bond = bond_lifes[name]
                        this_bond.add_frame(frame_id)
                    else:
                        this_bond = bondlife(min_O, max_O, frame_id, clone)
                        bond_lifes[this_bond.name] = this_bond
        long_bonds = [bond for bond in bond_lifes.values() if bond.age >= 100]
        np.save('%s/data/%s_%s_1234-IntraWater-sig.npy' % (project_dir, project, 0), long_bonds)
    print('Bonds of considerable length:')
    print([(bond.name, bond.age) for bond in long_bonds])

    for bond in long_bonds:
        if not static_positions.has_key(bond.clone):
            static_positions[bond.clone] = dict()
        if not anchors.has_key(project):
            protein_and_adp = find_anchors(traj.top)
            anchors[project] = protein_and_adp
        traj = md.load("/cbio/jclab/projects/fah/fah-data/munged2/all-atoms/%s/run%s-clone%s.h5" % (project, bond.clone/50, bond.clone%50))
        for water in bond.oxygens:
            if not static_positions[bond.clone].has_key(water):
                static_positions[bond.clone][water] = dict()
            neighbors = list()
            for frame_id in bond.frames:
                if not static_positions[bond.clone][water].has_key(frame_id):
                    instant_neighbors = md.compute_neighbors(traj.slice(frame_id, copy=False), 0.45, water, haystack_indices = protein_and_adp)
                    static_positions[bond.clone][water][frame_id] = instant_neighbors
                    neighbors.append(instant_neighbors)
                else:
                    neighbors.append(static_positions[bond.clone][water][frame_id])
            neighbors = chain.from_iterable(neighbors)
            unique_neighbors = set(neighbors)
            significant_neighbors = [neighbor for neighbor in unique_neighbors if neighbors.count(neighbor) > bond.life/2.0]
            if len(significant_neighbors) > 0:
                static_atoms = [str(traj.top.atom(index)) for index in significant_neighbors if index in anchors[project]]
                static_atoms.sort()
                [static_atoms.append('water') for index in significant_neighbors if traj.top.atom(index).residue.is_water]
                position_definitions.append(static_atoms)
        del(traj)
    unique_positions = set(position_definitions)
    to_save_position_info = list()
    print('All identified significantly maintained positions for project %s' % project)
    print('Atoms involved, how many times this position was seen')
    print('(may occur multiple times within one clone)')
    for position in unique_positions:
        local = (position, position_definitions.count(position))
        print(local)
        to_save_position_info.append(local)
    np.save('%s/data/%s_%s_1234-sig-positions-defined.npy' % (project_dir, project, 0),to_save_position_info)

print('Complete!')


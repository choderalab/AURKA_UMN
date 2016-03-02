import numpy as np
import sys
import math
import os
import mdtraj as md
from msmbuilder import dataset
from itertools import chain
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import seaborn as sns
sns.set_style("white")
sns.set_context("poster")


class AURKAnalysis(object):
    """
    calculate-all-SB-HB.py >> calculate_all_SB_HB

    RUN0_hbonds_2dhists.py >> plot_hbonds_on_all_residues
    RUN0_search_for_W1_W2.py >> plot_hbonds_on_all_residues

    RUN0_salt_bridge_2dhists.py >> plot_salt_bridges

    RUN0_275hbond_EQsaltbridge.py >> correlate_hbonds_bridges
    """
    def __init__(self, combine=True, store_arrays=False):
        self.projects = ['11410','11411']
        self.project_dirs = {'11410':'/cbio/jclab/projects/behrj/AURKA_UMN/output-1OL5','11411':'/cbio/jclab/projects/behrj/AURKA_UMN/output-1OL7'}
        self.system = {'11410':'with TPX2','11411':'without TPX2'}
        self.runs = range(5)

        self.residues_with_H = [185,181,274,275]
        self.reference = 185
        self.bridges = ['181-185','181-162']

        with open('/cbio/jclab/projects/behrj/AURKA_UMN/output-1OL7/run-index.txt','r') as fi:
            run_index = fi.read()
        self.mutant = dict()
        for entry in run_index.split('\n'):
            try:
                self.mutant[entry.split(' ')[0]] = entry.split(' ')[1]
            except:
                pass

        self.combine = combine

        self.store_arrays = store_arrays
        if self.store_arrays:
            self.HBonds = dict()
            self.SaltBridges = dict()

        self.offset = 400

        self.bin_x = np.arange(self.offset/4,500,10) - 0.25
        self.bond_bin_y = np.arange(7) - 0.5
        self.bridge_bin_y = np.arange(20) * 0.02 + 0.25

        self.bond_axis = [self.offset/4,500,-0.5,6.5]
        self.bridge_axis = [self.offset/4,500,0.25,0.65]

        self.xlabel = 't (nanoseconds)'


    def _find_hbonds_for_this_traj(self, traj, residue, haystack, sidechain=False, backbone=False):
        if sidechain:
            residue_atoms = [atom.index for atom in residue.atoms if atom.is_sidechain]
        elif backbone:
            residue_atoms = [atom.index for atom in residue.atoms if atom.is_backbone]
        else:
            residue_atoms = [atom.index for atom in residue.atoms]
        neighbor_set = self._find_neighbor_set(traj, residue_atoms, haystack)
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


    def _find_neighbor_set(self, traj, residue, haystack):
        neighbors = md.compute_neighbors(traj, 0.3, residue, haystack_indices=haystack)
        neighbor_set = set(chain.from_iterable(neighbors))
        return list(neighbor_set)

    def calculate_all_SB_HB(self, this_project=None, this_run=None, verbose=False, overwrite=False, save_files=False):
        if this_project is not None:
            projects = [this_project]
        else:
            projects = self.projects
        if this_run is not None:
            runs = [this_run]
        else:
            runs = self.runs
        for project in projects:
            project_dir = self.project_dirs[project]
            for run in runs:
                if verbose:
                    print("Loading Project %s RUN%s..." % (project, run))
                if not overwrite:
                    if (os.path.exists('%s/data/%s_%s_181_HBonds.npy' % (project_dir, project, run)) and
                        os.path.exists('%s/data/%s_%s_185_HBonds.npy' % (project_dir, project, run)) and
                        os.path.exists('%s/data/%s_%s_274_HBonds.npy' % (project_dir, project, run)) and
                        os.path.exists('%s/data/%s_%s_275_HBonds.npy' % (project_dir, project, run)) and
                        os.path.exists('%s/data/%s_%s_181-185_SB_total.npy' % (project_dir, project, run)) and
                        os.path.exists('%s/data/%s_%s_181-162_SB_total.npy' % (project_dir, project, run))):
                        continue
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

                    haystack = traj.top.select("water")

                    HB_185_total.append(self._find_hbonds_for_this_traj(traj, res185, haystack, sidechain=True))
                    HB_181_total.append(self._find_hbonds_for_this_traj(traj, e181, haystack, sidechain=True))
                    HB_274_total.append(self._find_hbonds_for_this_traj(traj, d274, haystack, backbone=True))
                    HB_275_total.append(self._find_hbonds_for_this_traj(traj, f275, haystack, backbone=True))
                if save_files:
                    np.save('%s/data/%s_%s_181_HBonds.npy' % (project_dir, project, run), HB_181_total)
                    np.save('%s/data/%s_%s_185_HBonds.npy' % (project_dir, project, run), HB_185_total)
                    np.save('%s/data/%s_%s_274_HBonds.npy' % (project_dir, project, run), HB_274_total)
                    np.save('%s/data/%s_%s_275_HBonds.npy' % (project_dir, project, run), HB_275_total)
                    np.save('%s/data/%s_%s_181-185_SB_total.npy' % (project_dir, project, run), SB_eq_total)
                    np.save('%s/data/%s_%s_181-162_SB_total.npy' % (project_dir, project, run), SB_ek_total)

        if verbose:
            print('Complete!')
        if not self.store_arrays:
            return
        self.HBonds['%s-RUN%s-181' % (project, run)] = HB_181_total
        self.HBonds['%s-RUN%s-185' % (project, run)] = HB_185_total
        self.HBonds['%s-RUN%s-274' % (project, run)] = HB_274_total
        self.HBonds['%s-RUN%s-275' % (project, run)] = HB_275_total
        self.SaltBridges['%s-RUN%s-181-185' % (project, run)] = SB_eq_total
        self.SaltBridges['%s-RUN%s-181-162' % (project, run)] = SB_ek_total

    def _plot_2dhist(self, x_axis, y_axis, bins, run, project, title, xlabel, ylabel, axis, filename, weights=None):
        fig1 = plt.figure()
        plt.hist2d(x_axis,y_axis,bins=bins,weights=weights,cmap=plt.get_cmap('jet'))
        plt.title(title)
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        plt.colorbar()
        plt.axis(axis)
        plt.savefig(filename,dpi=300)
        plt.close(fig1)
        print('Saved %s' % filename)

    def _plot_hbond_2dhist(self, residue, x_axis, hbond_count, weights, run, project, compare_to, count_waters):
        x_axis = x_axis[hbond_count > -1]
        hbond_count = hbond_count[hbond_count > -1]
        bins = [self.bin_x, self.bond_bin_y]
        weights = weights[hbond_count > -1]
        xlabel = 't (nanoseconds)'
        if count_waters:
            title = 'Waters interacting with AURKA %s residue %s over time %s' % (self.mutant['RUN%s' % run], residue, self.system[project])
            ylabel = 'number of water molecules'
            if compare_to is not None:
                residue = str(residue)+'-possible-W1-W2'
            filename = '"./plots/AURKA-%s-waters-hist2d-entire-traj-%s-combined-RUN%s" % (residue, project, run)'
        else:
            title = 'AURKA %s number of hydrogen bonds on residue %s over time %s' % (self.mutant['RUN%s' % run], residue, self.system[project])
            ylabel = 'number of hydrogen bonds'
            if compare_to is not None:
                residue = str(residue)+'-possible-W1-W2'
        axis = self.bond_axis
        filename = "./plots/AURKA-%s-hbonds-hist2d-entire-traj-%s-combined-RUN%s.png" % (residue, project, run)
        self._plot_2dhist(x_axis, hbond_count, bins, run, project, title, xlabel, ylabel, axis, filename, weights=weights)

    def _plot_bridge_2dhist(self, bridge, x_axis, minimum_distance, weights, run, project):
        x_axis = x_axis[minimum_distance > 0]
        minimum_distance = minimum_distance[minimum_distance > 0]
        bins = [self.bin_x, self.bridge_bin_y]
        weights = weights[minimum_distance > 0]
        title = 'AURKA %s minimum %s salt bridge distance over time %s' % (self.mutant['RUN%s' % run], bridge, self.system[project])
        xlabel = 't (nanoseconds)'
        ylabel = 'distance r (nanometers) between residues %s and %s' % (bridge.split('-')[0], bridge.split('-')[1])
        axis = self.bridge_axis
        filename = "./plots/AURKA-salt-bridge-%s-hist2d-entire-traj-%s-combined-RUN%s.png" % (bridge, project, run)
        self._plot_2dhist(x_axis, minimum_distance, bins, run, project, title, xlabel, ylabel, axis, filename, weights=weights)

    def _plot_correlated_bridge_2dhist(self, bridge, x_axis, minimum_distance, weights, run, project, has_bonds):
        x_axis = x_axis[minimum_distance > 0]
        minimum_distance = minimum_distance[minimum_distance > 0]
        bins = [self.bin_x, self.bridge_bin_y]
        weights = weights[minimum_distance > 0]
        if has_bonds:
            title = 'AURKA %s, with 275 Hbond,  minimum %s salt bridge distance over time %s' % (self.mutant['RUN%s' % run], bridge, self.system[project])
            filename = "./plots/AURKA-275bonds-salt-bridge-%s-hist2d-entire-traj-%s-combined-RUN%s" % (bridge, project, run)
        else:
            title = 'AURKA %s, no 275 Hbonds,  minimum %s salt bridge distance over time %s' % (self.mutant['RUN%s' % run], bridge, self.system[project])
            filename = "./plots/AURKA-no-bonds-salt-bridge-%s-hist2d-entire-traj-%s-combined-RUN%s" % (bridge, project, run)
        xlabel = 't (nanoseconds)'
        ylabel = 'distance r (nanometers) between residues %s and %s' % (bridge.split('-')[0], bridge.split('-')[1])
        axis = self.bridge_axis
        self._plot_2dhist(x_axis, minimum_distance, bins, run, project, title, xlabel, ylabel, axis, filename, weights=weights)

    def _count_and_plot_res_bonds(self, residue, HB_res_total, run, project, compare_to=None, correlate=False, count_waters=False):
        total_clones = len(HB_res_total)
        offset = self.offset
        bin_x = self.bin_x

        if compare_to is not None and np.array_equal(HB_res_total, compare_to):
            compare_to = None

        hbond_count = np.zeros((total_clones,2000-offset)) - 1
        x_axis = np.zeros((total_clones,2000-offset))
        weights = np.zeros((total_clones,2000-offset))
        column_count = np.zeros(bin_x.shape)
        for clone, traj in enumerate(HB_res_total):
            for index in range(offset,2000):
                x_axis[clone][index-offset] = index*0.25
                try:
                    hbonds_frame = traj[index]
                except:
                    continue
                if compare_to is None:
                    if count_waters:
                        hbond_count[clone][index-offset] = len(self._find_waters(hbonds_frame))
                    else:
                        hbond_count[clone][index-offset] = hbonds_frame.shape[0]
                    column_count[(index-offset-0.25)/40] += 1
                else: # you can only count waters if not comparing to something else
                    reference_donors = [bond[0] for bond in compare_to[clone][index]]
                    reference_acceptors = [bond[2] for bond in compare_to[clone][index]]
                    count = 0
                    if count_waters:
                        for water in self._find_waters(hbonds_frame):
                            if water in reference_donors or water in reference_acceptors:
                                count += 1
                    else:
                        for bond in hbonds_frame:
                            if (bond[2] in reference_donors or bond[0] in reference_donors or
                                bond[2] in reference_acceptors or bond[0] in reference_acceptors):
                                count += 1
                    hbond_count[clone][index-offset] = count
                    column_count[(index-offset-0.25)/40] += 1
        if correlate:
            self._build_bridges(hbond_count, x_axis, run, project)
        else:
            for clone, traj in enumerate(HB_res_total):
                for index in range(offset,2000):
                    weights[clone][index-offset] = 1.00 / column_count[(index-offset-0.25)/40]
            x_axis = x_axis.flatten()
            hbond_count = hbond_count.flatten()
            weights = weights.flatten()
            self._plot_hbond_2dhist(residue, x_axis, hbond_count, weights, run, project, compare_to, count_waters)

    def _find_waters(self, hbonds_frame):
        waters = set()
        for bond in hbonds_frame:
            if abs(bond[0] - bond[1]) <= 2:
                waters.add(bond[0])
            else:
                waters.add(bond[2])
        return waters

    def plot_hbonds_on_all_residues(self, this_project=None, this_run=None, use_reference=False, count_waters=False):
        # combine = True
        run = 0
        if this_project is not None:
            projects = [this_project]
        else:
            projects = self.projects
        if this_run is not None:
            runs = [this_run]
        else:
            runs = self.runs
        for i, project in enumerate(projects):
            HB_total = self._load_hbonds(project, this_run=this_run)
            for key in HB_total.keys():
                if use_reference:
                    self._count_and_plot_res_bonds(key, HB_total[key], run, project, compare_to=HB_total[self.reference],count_waters=count_waters)
                else:
                    self._count_and_plot_res_bonds(key, HB_total[key], run, project, count_waters=count_waters)

    def correlate_hbonds_bridges(self, this_project=None, this_run=None):
        run = 0
        correlate_to = [275]
        if this_project is not None:
            projects = [this_project]
        else:
            projects = self.projects
        if this_run is not None:
            runs = [this_run]
        else:
            runs = self.runs
        for i, project in enumerate(projects):
            HB_total = self._load_hbonds(project, this_run=this_run)
            for key in correlate_to:
                self._count_and_plot_res_bonds(key, HB_total[key], run, project, correlate=True)

    def _load_hbonds(self, project, this_run=None):
        if this_run is not None:
            runs = [this_run]
        else:
            runs = self.runs
        project_dir = self.project_dirs[project]
        HB_total = dict()
        for residue in self.residues_with_H:
            for run in runs:
                if self.store_arrays and self.HBonds.has_key('%s-RUN%s-%s' % (project, run, residue)):
                    new_run = self.HBonds['%s-RUN%s-%s' % (project, run, residue)]
                else:
                    try:
                        new_run = np.load('%s/data/%s_%s_%s_HBonds.npy' % (project_dir, project, run, residue))
                        if self.store_arrays:
                            self.HBonds['%s-RUN%s-%s' % (project, run, residue)] = new_run
                    except:
                        continue
                if not HB_total.has_key(residue):
                    HB_total[residue] = new_run
                else:
                    HB_total[residue] = np.concatenate((HB_total[residue], new_run))
        return HB_total

    def plot_salt_bridges(self, this_project=None, this_run=None):
        # only going to work with combine right now
        if this_project is not None:
            projects = [this_project]
        else:
            projects = self.projects
        if this_run is not None:
            runs = [this_run]
        else:
            runs = self.runs
        if self.combine:
            total_clones = 250
        else:
            total_clones = 50
        offset = self.offset
        for i, project in enumerate(projects):
            project_dir = self.project_dirs[project]
            for bridge in self.bridges:
                minimum_distance = np.zeros((total_clones, 2000-self.offset))
                x_axis = np.zeros((total_clones, 2000-self.offset))
                weights = np.zeros((total_clones, 2000-self.offset))
                column_count = np.zeros(self.bin_x.shape)
                for run in runs:
                    if self.store_arrays and self.SaltBridges.has_key('%s-RUN%s-%s' % (project, run, bridge)):
                        SB_total = self.SaltBridges['%s-RUN%s-%s' % (project, run, bridge)]
                    else:
                        try:
                            SB_total = np.load('%s/data/%s_%s_%s_SB_total.npy' % (project_dir, project, run, bridge))
                            if self.store_arrays:
                                self.SaltBridges['%s-RUN%s-%s' % (project, run, bridge)] = SB_total
                        except:
                            continue
                    for clone, traj in enumerate(SB_total):
                        for index in range(offset,2000):
                            try:
                                minimum_distance[run*50+clone][index-offset] = traj[index]
                                column_count[(index-offset-0.25)/40] += 1
                            except:
                                pass
                            x_axis[run*50+clone][index-offset] = index*0.25
                for clone, traj in enumerate(SB_total):
                    for index in range(offset,2000):
                        weights[clone][index-offset] = 1.00 / column_count[(index-offset-0.25)/40]
                x_axis = x_axis.flatten()
                minimum_distance = minimum_distance.flatten()
                weights = weights.flatten()
                if self.combine:
                    run = 0
                self._plot_bridge_2dhist(bridge, x_axis, minimum_distance, weights, run, project)

    def _build_bridges(self, hbond_count, x_axis, run, project):
        correlated_bridges = ['181-185']
        project_dir = self.project_dirs[project]
        if self.combine:
            total_clones = 250
        else:
            total_clones = 50
        offset = self.offset
        for bridge in correlated_bridges:
            minimum_distance_0bond = np.zeros((total_clones, 2000-offset))
            minimum_distance_1bond = np.zeros((total_clones, 2000-offset))
            weights_0bond = np.zeros((total_clones, 2000-offset))
            weights_1bond = np.zeros((total_clones, 2000-offset))
            column_count_0bond = np.zeros(self.bin_x.shape)
            column_count_1bond = np.zeros(self.bin_x.shape)
            for run in self.runs:
                if self.store_arrays and self.SaltBridges.has_key('%s-RUN%s-%s' % (project, run, bridge)):
                    SB_total = self.SaltBridges['%s-RUN%s-%s' % (project, run, bridge)]
                else:
                    try:
                        SB_total = np.load('%s/data/%s_%s_%s_SB_total.npy' % (project_dir, project, run, bridge))
                        if self.store_arrays:
                            self.SaltBridges['%s-RUN%s-%s' % (project, run, bridge)] = SB_total
                    except:
                        continue
                for clone, traj in enumerate(SB_total):
                    for index in range(offset, 2000):
                        try:
                            frame_distance = traj[index]
                        except:
                            continue
                        if hbond_count[run*50+clone][index-offset] > 0: # COMBINE ONLY
                            minimum_distance_1bond[run*50+clone][index-offset] = frame_distance
                            column_count_1bond[(index-offset-0.25)/40] += 1
                        elif hbond_count[run*50+clone][index-offset] == 0:
                            minimum_distance_0bond[run*50+clone][index-offset] = frame_distance
                            column_count_0bond[(index-offset-0.25)/40] += 1
            for clone, traj in enumerate(SB_total):
                for index in range(offset,2000):
                    weights_0bond[clone][index-offset] = 1.00 / column_count_0bond[(index-offset-0.25)/40]
                    weights_1bond[clone][index-offset] = 1.00 / column_count_1bond[(index-offset-0.25)/40]
            x_axis = x_axis.flatten()
            minimum_distance_0bond = minimum_distance_0bond.flatten()
            minimum_distance_1bond = minimum_distance_1bond.flatten()
            weights_0bond = weights_0bond.flatten()
            weights_1bond = weights_1bond.flatten()
            self._plot_correlated_bridge_2dhist(bridge, x_axis, minimum_distance_0bond, weights_0bond, 0, project, False)
            self._plot_correlated_bridge_2dhist(bridge, x_axis, minimum_distance_1bond, weights_1bond, 0, project, True)
















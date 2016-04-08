"""
Creates two separate plots of salt bridge distances between residues 181 and 185
based on whether or not a water is bound to F275 in that frame
(salt bridge is more likely to form when F275 is bound elsewhere)
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
import math
from matplotlib.pyplot import cm
import seaborn as sns
import os

sns.set_style("white")
sns.set_context("poster")

        # make plots of all data past t = 250ns
            # quantify how much P(salt bridge) and (1-P)
                # P = Ntraj_frames_bound

        # look at trajectories : how long do waters stay in place
            # may need to account for the waters exchanging

offset = 400

residues_with_H = [275]
bridges = ['181-185']

local_path = os.path.dirname(os.path.realpath(__file__))

projects = ['11410','11411','11418']
project_dirs = {'11410':'%s/../output-1OL5' % local_path,'11411':'%s/../output-1OL7' % local_path,'11418':'%s/../output1OL5-TPX2' % local_path}
system = {'11410':'with TPX2','11411':'without TPX2','11418': 'with TPX2 removed'}

with open('../output-1OL7/run-index.txt','r') as fi:
    run_index = fi.read()
mutant = dict()
for entry in run_index.split('\n'):
    try:
        mutant[entry.split(' ')[0]] = entry.split(' ')[1]
    except:
        pass

offset = 400
bin_x = np.arange(offset/4,510,10) - 0.25
bin_y = {'181-185': np.arange(21) * 0.04 + 0.25, '181-162': np.arange(20) * 0.02 + 0.25}
axis = {'181-185': [offset/4,500,0.25,1.05], '181-162':[offset/4,500,0.25,0.63]}

def plot_bridge_2dhist(bridge, x_axis, minimum_distance, weights, run, project, has_bonds):
    fig1 = plt.figure()
    plt.hist2d(x_axis[minimum_distance > 0],minimum_distance[minimum_distance > 0],bins=[bin_x,bin_y[bridge]],weights=weights[minimum_distance > 0],cmap=plt.get_cmap('jet'))
    if has_bonds:
        plt.title('AURKA %s, with 275 Hbond,  minimum %s salt bridge distance over time %s' % (mutant['RUN%s' % run], bridge, system[project]))
    else:
        plt.title('AURKA %s, no 275 Hbonds,  minimum %s salt bridge distance over time %s' % (mutant['RUN%s' % run], bridge, system[project]))
    plt.ylabel('distance r (nanometers) between residues %s and %s' % (bridge.split('-')[0], bridge.split('-')[1]))
    plt.xlabel('t (nanoseconds)')
    plt.colorbar()
    plt.axis(axis[bridge])
    if has_bonds:
        filename = "../plots/AURKA-275bonds-salt-bridge-%s-hist2d-entire-traj-%s-combined-RUN%s" % (bridge, project, run)
    else:
        filename = "../plots/AURKA-no-bonds-salt-bridge-%s-hist2d-entire-traj-%s-combined-RUN%s" % (bridge, project, run)
    plt.savefig(filename,dpi=300)
    plt.close(fig1)
    print('Saved %s' % filename)


def count_res_bonds_plot_bridges(residue, HB_res_total, project, ADP_bound, compare_to=None):
    hbond_count = np.zeros((5*50,2000-offset)) - 1
    x_axis = np.zeros((5*50,2000-offset))
    weights = np.zeros((5*50,2000-offset))
    column_count = np.zeros(bin_x.shape)
    for clone, traj in enumerate(HB_res_total):
        for index in range(offset,2000):
            x_axis[clone][index-offset] = index*0.25
            if not ADP_bound[clone][index]:
                continue
            if compare_to is None:
                try:
                    hbond_count[clone][index-offset] = traj[index].shape[0]
                    column_count[(index-offset-0.25)/40] += 1
                except:
                    pass
            else:
                try:
                    reference_donors = [bond[0] for bond in compare_to[clone][index]]
                    reference_acceptors = [bond[2] for bond in compare_to[clone][index]]
                    count = 0
                    for bond in traj[index]:
                        if (bond[2] in reference_donors or bond[0] in reference_donors or
                            bond[2] in reference_acceptors or bond[0] in reference_acceptors):
                            count += 1
                    hbond_count[clone][index-offset] = count
                    column_count[(index-offset-0.25)/40] += 1
                except:
                    pass
    for clone, traj in enumerate(HB_res_total):
        for index in range(offset,2000):
            weights[clone][index-offset] = 1.00 / column_count[(index-offset-0.25)/40]
    build_bridges(hbond_count, project, ADP_bound)

def build_bridges(hbond_count, project, ADP_bound):
    project_dir = project_dirs[project]
    for bridge in bridges:
        x_axis = np.zeros((5*50,2000-offset))
        minimum_distance_0bond = np.zeros((5*50,2000-offset))
        minimum_distance_1bond = np.zeros((5*50,2000-offset))
        weights_0bond = np.zeros((5*50,2000-offset))
        weights_1bond = np.zeros((5*50,2000-offset))
        column_count_0bond = np.zeros(bin_x.shape)
        column_count_1bond = np.zeros(bin_x.shape)
        for run in range(5):
            if not os.path.exists('%s/data/%s_%s_%s_SB_total.npy' % (project_dir, project, run, bridge)):
                continue
            SB_total = np.load('%s/data/%s_%s_%s_SB_total.npy' % (project_dir, project, run, bridge))

            for clone, traj in enumerate(SB_total):
                for index in range(offset,2000):
                    x_axis[run*50+clone][index-offset] = index*0.25
                    if not ADP_bound[clone][index]:
                        continue
                    try:
                        frame_distance = traj[index]
                    except:
                        continue
                    if hbond_count[run*50+clone][index-offset] > 0:
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
        weights_0bond = weights_0bond.flatten()
        minimum_distance_1bond = minimum_distance_1bond.flatten()
        weights_1bond = weights_1bond.flatten()
        plot_bridge_2dhist(bridge, x_axis, minimum_distance_0bond, weights_0bond, 0, project, False)
        plot_bridge_2dhist(bridge, x_axis, minimum_distance_1bond, weights_1bond, 0, project, True)


for i, project in enumerate(projects):
    project_dir = project_dirs[project]
    HB_total = dict()
    ADP_bound = np.load('%s/is-ADP-bound.npy' % project_dir)
    for residue in residues_with_H:
        for run in range(5):
            if not os.path.exists('%s/data/%s_%s_%s_HBonds.npy' % (project_dir, project, run, residue)):
                continue
            if not HB_total.has_key(residue):
                HB_total[residue] = np.load('%s/data/%s_%s_%s_HBonds.npy' % (project_dir, project, run, residue))
            else:
                new_run = np.load('%s/data/%s_%s_%s_HBonds.npy' % (project_dir, project, run, residue))
                HB_total[residue] = np.concatenate((HB_total[residue], new_run))

    for key in HB_total.keys():
        count_res_bonds_plot_bridges(key, HB_total[key], project, ADP_bound)

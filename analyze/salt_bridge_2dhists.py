import numpy as np
import sys
import math
import os
import plot_function

        # make plots of all data past t = 250ns
            # quantify how much P(salt bridge) and (1-P)
                # P = Ntraj_frames_bound

        # look at trajectories : how long do waters stay in place
            # may need to account for the waters exchanging

local_path = os.path.dirname(os.path.realpath(__file__))

projects = ['11410','11411','11418']
project_dirs = {'11410':'%s/../output-1OL5' % local_path,'11411':'%s/../output-1OL7' % local_path,'11418':'%s/../output1OL5-TPX2' % local_path}
system = {'11410':'with TPX2','11411':'without TPX2','11418': 'with TPX2 removed'}

with open('./output-1OL7/run-index.txt','r') as fi:
    run_index = fi.read()
mutant = dict()
for entry in run_index.split('\n'):
    try:
        mutant[entry.split(' ')[0]] = entry.split(' ')[1]
    except:
        pass

offset = plot_function.OFFSET
bin_x = plot_function.BIN_X

def plot_2dhist(bridge, x_axis, minimum_distance, weights, run, project):
    title = 'AURKA %s minimum %s salt bridge distance over time %s' % (mutant['RUN%s' % run], bridge, system[project])
    filename = "./plots/AURKA-salt-bridge-%s-hist2d-entire-traj-%s-RUN%s" % (bridge, project, run)
    ylabel = 'distance r (nanometers) between residues %s and %s' % (bridge.split('-')[0], bridge.split('-')[1])
    plot_function.plot_2dhist(bridge, x_axis, minimum_distance, weights, title, ylabel, filename)

for i, project in enumerate(projects):
    project_dir = project_dirs[project]
    for run in range(5):
        for bridge in ['181-185','181-162']:
            if not os.path.exists('%s/data/%s_%s_%s_SB_total.npy' % (project_dir, project, run, bridge)):
                continue
            SB_total = np.load('%s/data/%s_%s_%s_SB_total.npy' % (project_dir, project, run, bridge))

            minimum_distance = np.zeros((50,2000))
            x_axis = np.zeros((50,2000))
            weights = np.zeros((50,2000-offset))
            column_count = np.zeros(bin_x.shape)
            for clone, traj in enumerate(SB_total):
                for index in range(offset, 2000):
                    try:
                        minimum_distance[clone][index-offset] = traj[index]
                        column_count[(index-offset-0.25)/40] += 1
                    except:
                        pass
                    x_axis[clone][index] = index*0.25
            for clone, traj in enumerate(SB_total):
                for index in range(offset,2000):
                    weights[clone][index-offset] = 1.00 / column_count[(index-offset-0.25)/40]
            x_axis = x_axis.flatten()
            minimum_distance = minimum_distance.flatten()
            plot_2dhist(bridge, x_axis, minimum_distance, run, project)

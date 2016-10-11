"""
Plot minimum distances between 181-185 and 181-162
"""
import os
import numpy as np
import sys
import math
import plot_function

        # make plots of all data past t = 250ns
            # quantify how much P(salt bridge) and (1-P)
                # P = Ntraj_frames_bound

        # look at trajectories : how long do waters stay in place
            # may need to account for the waters exchanging

local_path = os.path.dirname(os.path.realpath(__file__))

projects = ['11414','11419','11418','11423']
project_dirs = {
    '11410':'%s/../output-1OL5' % local_path,
    '11411':'%s/../output-1OL7' % local_path,
    '11414':'%s/../output-11414' % local_path,
    '11418':'%s/../output-11418' % local_path,
    '11419':'%s/../output-11419' % local_path,
    '11423':'%s/../output-11423' % local_path,
    '11424':'%s/../output-11424' % local_path,
    '11425':'%s/../output-11425' % local_path,
}

system = {
    '11410':'with TPX2',
    '11411':'without TPX2',
    '11414':'with TPX2',
    '11418':'with TPX2 removed',
    '11419':'with TPX2',
    '11423':'with TPX2 removed',
    '11424':'with TPX2; charmm',
    '11425':'with TPX2 removed; charmm'
}
run_guide = dict()
mutants = dict()
for project in projects:
    run_guide[project] = 0
    filename = project_dirs[project]+'/run-index.txt'
    with open(filename, 'r') as fi:
        project_run_index = fi.read()
    for entry in project_run_index.split('\n'):
        try:
            run = entry.split(' ')[0]
            mutant = entry.split(' ')[1]
            run = int(run[3:])
            run_guide[project] += 1
            mutants[(project, run)] = mutant
        except:
            pass
print(run_guide)

offset = plot_function.OFFSET
bin_x = plot_function.BIN_X

def plot_2dhist(bridge, x_axis, minimum_distance, weights, run, project):
    title = 'AURKA %s minimum %s salt bridge distance over time %s' % (mutants[(project, run)], bridge, system[project])
    filename = "../plots/AURKA-salt-bridge-%s-hist2d-entire-traj-%s-combined-RUN%s" % (bridge, project, run)
    ylabel = 'distance r (nanometers) between residues %s and %s' % (bridge.split('-')[0], bridge.split('-')[1])
    plot_function.plot_2dhist(bridge, x_axis, minimum_distance, weights, title, ylabel, filename)

for i, project in enumerate(projects):
    project_dir = project_dirs[project]
    for bridge in ['181-162']:
        runs = run_guide[project]
        for run in range(runs):
            if project == '11419' and run == 4:
                continue
            minimum_distance = np.zeros((50,2000-offset)) - 1
            x_axis = np.zeros((50,2000-offset))
            weights = np.zeros((50,2000-offset))
            column_count = np.zeros(bin_x.shape)
            if not os.path.exists('%s/data/%s_%s_%s_SB_total.npy' % (project_dir, project, run, bridge)):
                print('BRIDGE DATA MISSING FOR %s RUN %s' % (project, run))
                continue
            SB_total = np.load('%s/data/%s_%s_%s_SB_total.npy' % (project_dir, project, run, bridge))

            for clone, traj in enumerate(SB_total):
                for index in range(offset,2000):
                    x_axis[clone][index-offset] = index*0.25
                    try:
                        minimum_distance[clone][index-offset] = traj[index]
                        column_count[(index-offset-0.25)/40] += 1
                    except:
                        pass
            for clone, traj in enumerate(SB_total):
                for index in range(offset,2000):
                    weights[clone][index-offset] = 1.00 / column_count[(index-offset-0.25)/40]
            x_axis = x_axis.flatten()
            minimum_distance = minimum_distance.flatten()
            weights = weights.flatten()
            plot_2dhist(bridge, x_axis, minimum_distance, weights, run, project)

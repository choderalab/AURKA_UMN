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


projects = ['11414','11419']
#conditions = ['11418','11423']

offset = plot_function.OFFSET
bin_x = plot_function.BIN_X

def plot_2dhist(bridge, x_axis, minimum_distance, weights, mutant, project):
    title = 'AURKA %s minimum %s CB-CB distance over time %s' % (mutant, bridge, system[project])
    filename = "../plots/AURKA-CB-%s-hist2d-entire-traj-%s-combined-%s.pdf" % (bridge, project, mutant)
    ylabel = 'distance r (nanometers) between residues %s and %s' % (bridge.split('-')[0], bridge.split('-')[1])
    plot_function.plot_2dhist(bridge, x_axis, minimum_distance, weights, title, ylabel, filename)

for bridge in ['284-225','287-225']:
    minimum_distance_WT = np.zeros((50*5,2000-offset)) - 1
    x_axis_WT = np.zeros((50*5,2000-offset))
    weights_WT = np.zeros((50*5,2000-offset))
    minimum_distance_MUT = np.zeros((50*4,2000-offset)) - 1
    x_axis_MUT = np.zeros((50*4,2000-offset))
    weights_MUT = np.zeros((50*4,2000-offset))
    wt_index = 0
    mut_index = 0

    for i, project in enumerate(projects):
        project_dir = project_dirs[project]
        runs = run_guide[project]
        for run in range(runs):
            if project == '11419' and run == 4:
                continue
            if mutants[(project, run)] == 'WT':
                minimum_distance = minimum_distance_WT
                x_axis = x_axis_WT
                weights = weights_WT
                run_count = wt_index
                wt_index += 1
            elif mutants[(project, run)] in ['Q185N','Q185M','Q185C','Q185L']:
                minimum_distance = minimum_distance_MUT
                x_axis = x_axis_MUT
                weights = weights_MUT
                run_count = mut_index
                mut_index += 1
            else:
                continue
            if not os.path.exists('%s/data/%s_%s_%s_SB_total.npy' % (project_dir, project, run, bridge)):
                print('BRIDGE DATA MISSING FOR %s RUN %s' % (project, run))
                continue
            SB_total = np.load('%s/data/%s_%s_%s_SB_total.npy' % (project_dir, project, run, bridge))

            column_count = np.zeros(bin_x.shape)
            for clone, traj in enumerate(SB_total):
                for index in range(offset,2000):
                    x_axis[50*run_count+clone][index-offset] = index*0.25
                    try:
                        minimum_distance[50*run_count+clone][index-offset] = traj[index]
                        column_count[(index-offset-0.25)/40] += 1
                    except:
                        pass
            for clone, traj in enumerate(SB_total):
                for index in range(offset,2000):
                    weights[50*run_count+clone][index-offset] = 1.00 / column_count[(index-offset-0.25)/40]

    x_axis_WT = x_axis_WT.flatten()
    minimum_distance_WT = minimum_distance_WT.flatten()
    weights_WT = weights_WT.flatten()
    plot_2dhist(bridge, x_axis_WT, minimum_distance_WT, weights_WT, 'WT', project)

    x_axis_MUT = x_axis_MUT.flatten()
    minimum_distance_MUT = minimum_distance_MUT.flatten()
    weights_MUT = weights_MUT.flatten()
    plot_2dhist(bridge, x_axis_MUT, minimum_distance_MUT, weights_MUT, 'Q185mutants', project)


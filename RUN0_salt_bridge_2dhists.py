import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
import math
from matplotlib.pyplot import cm
import seaborn as sns
import os

sns.set_style("whitegrid")
sns.set_context("poster")

        # make plots of all data past t = 250ns
            # quantify how much P(salt bridge) and (1-P)
                # P = Ntraj_frames_bound

        # look at trajectories : how long do waters stay in place
            # may need to account for the waters exchanging

projects = ['11410','11411']
project_dirs = {'11410':'./output-1OL5','11411':'./output-1OL7'}
system = {'11410':'with TPX2','11411':'without TPX2'}

with open('./output-1OL7/run-index.txt','r') as fi:
    run_index = fi.read()
mutant = dict()
for entry in run_index.split('\n'):
    try:
        mutant[entry.split(' ')[0]] = entry.split(' ')[1]
    except:
        pass

offset = 400
bin_x = np.arange(offset/4,500,10) - 0.25
bin_y = np.arange(20) * 0.02 + 0.25

def plot_2dhist(bridge, x_axis, minimum_distance, weights, run, project):
    fig1 = plt.figure()
    plt.hist2d(x_axis[minimum_distance > 0],minimum_distance[minimum_distance > 0],bins=[bin_x,bin_y],weights=weights[minimum_distance > 0],cmap=plt.get_cmap('jet'))
    plt.title('AURKA %s minimum %s salt bridge distance over time %s' % (mutant['RUN%s' % run], bridge, system[project]))
    plt.ylabel('distance r (nanometers) between residues %s and %s' % (bridge.split('-')[0], bridge.split('-')[1]))
    plt.xlabel('t (nanoseconds)')
    plt.colorbar()
    plt.axis([offset/4,500,0.25,0.65])
    plt.savefig("./plots/AURKA-salt-bridge-%s-hist2d-entire-traj-%s-combined-RUN%s" % (bridge, project, run),dpi=300)
    plt.close(fig1)
    print('Saved ./plots/AURKA-salt-bridge-%s-hist2d-entire-traj-%s-combined-RUN%s.png' % (bridge, project, run))

for i, project in enumerate(projects):
    project_dir = project_dirs[project]
    for bridge in ['181-185','181-162']:
        minimum_distance = np.zeros((5*50,2000-offset))
        x_axis = np.zeros((5*50,2000-offset))
        weights = np.zeros((5*50,2000-offset))
        column_count = np.zeros(bin_x.shape)
        for run in range(5):
            if not os.path.exists('%s/data/%s_%s_%s_SB_total.npy' % (project_dir, project, run, bridge)):
                continue
            SB_total = np.load('%s/data/%s_%s_%s_SB_total.npy' % (project_dir, project, run, bridge))

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
        plot_2dhist(bridge, x_axis, minimum_distance, weights, 0, project)

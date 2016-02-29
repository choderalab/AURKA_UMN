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

def plot_2dhist(bridge, x_axis, minimum_distance, run, project):
    fig1 = plt.figure()
    plt.hist2d(x_axis[minimum_distance > 0],minimum_distance[minimum_distance > 0],bins=[100,20],cmap=plt.get_cmap('jet'))
    plt.title('AURKA %s minimum %s salt bridge distance over time %s' % (mutant['RUN%s' % run], bridge, system[project]))
    plt.ylabel('distance r (nanometers) between residues %s and %s' % (bridge.split('-')[0], bridge.split('-')[1]))
    plt.xlabel('t (nanoseconds)')
    plt.colorbar()
    plt.axis([0,500,0.25,0.65])
    plt.savefig("./plots/AURKA-salt-bridge-%s-hist2d-entire-traj-%s-RUN%s" % (bridge, project, run),dpi=300)
    plt.close(fig1)

for i, project in enumerate(projects):
    project_dir = project_dirs[project]
    for bridge in ['181-185','181-162']:
        minimum_distance = np.zeros((5*50,2000))
        x_axis = np.zeros((5*50,2000))
        for run in range(5):
            if not os.path.exists('%s/data/%s_%s_%s_SB_total.npy' % (project_dir, project, run, bridge)):
                continue
            SB_total = np.load('%s/data/%s_%s_%s_SB_total.npy' % (project_dir, project, run, bridge))

            for clone, traj in enumerate(SB_total):
                for index in range(2000):
                    try:
                        minimum_distance[run*50+clone][index] = traj[index]
                    except:
                        pass
                    x_axis[run*50+clone][index] = index*0.25
        x_axis = x_axis.flatten()
        minimum_distance = minimum_distance.flatten()
        plot_2dhist(bridge, x_axis, minimum_distance, 0, project)

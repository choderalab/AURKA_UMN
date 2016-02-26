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

offset = 400

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

bin_x = np.arange(offset/4,500,10) - 0.25
bin_y = np.arange(7) - 0.5

def plot_2dhist(residue, x_axis, hbond_count, run, project):
    fig1 = plt.figure()
    plt.hist2d(x_axis[hbond_count > -1],hbond_count[hbond_count > -1],bins=[bin_x,bin_y],cmap=plt.get_cmap('jet'))
    plt.title('AURKA %s number of hydrogen bonds on residue %s over time %s' % (mutant['RUN%s' % run], residue, system[project]))
    plt.ylabel('number of hydrogen bonds')
    plt.xlabel('t (nanoseconds)')
    plt.colorbar()
    plt.axis([offset/4,500,-0.5,6.5])
    plt.savefig("./plots/AURKA-%s-hbonds-hist2d-entire-traj-%s-RUN%s" % (residue, project, run),dpi=300)
    plt.close(fig1)

for i, project in enumerate(projects):
    project_dir = project_dirs[project]
    for run in range(5):
        for residue in [181,185,274,275]:
            if not os.path.exists('%s/data/%s_%s_%s_HBonds.npy' % (project_dir, project, run, residue)):
                continue
            HB_total = np.load('%s/data/%s_%s_%s_HBonds.npy' % (project_dir, project, run, residue))

            hbond_count = np.zeros((50,2000-offset)) - 1
            x_axis = np.zeros((50,2000-offset))
            for clone, traj in enumerate(HB_total):
                for index in range(offset,2000):
                    try:
                        hbond_count[clone][index-offset] = traj[index].shape[0]
                    except:
                        pass
                    x_axis[clone][index-offset] = index*0.25
            x_axis = x_axis.flatten()
            hbond_count = hbond_count.flatten()
            plot_2dhist(residue, x_axis, hbond_count, run, project)

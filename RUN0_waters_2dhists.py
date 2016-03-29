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

#residues_with_H = [185,181,274,275]
residues_with_H = [181]
reference = 185

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

bin_x = np.arange(offset/4,510,10) - 0.25
bin_y = np.arange(11) - 0.5

def plot_2dhist(residue, x_axis, hbond_count, weights, run, project):
    fig1 = plt.figure()
    plt.hist2d(x_axis[hbond_count > -1],hbond_count[hbond_count > -1],bins=[bin_x,bin_y],weights=weights[hbond_count > -1],cmap=plt.get_cmap('jet'))
    plt.title('Waters interacting with AURKA %s residue %s over time %s' % (mutant['RUN%s' % run], residue, system[project]))
    plt.ylabel('number of water molecules')
    plt.xlabel('t (nanoseconds)')
    plt.colorbar()
    plt.axis([offset/4,500,-0.5,9.5])
    plt.savefig("./plots/AURKA-%s-waters-hist2d-entire-traj-%s-combined-RUN%s" % (residue, project, run),dpi=300)
    plt.close(fig1)
    print('Saved ./plots/AURKA-%s-waters-hist2d-entire-traj-%s-combined-RUN%s.png' % (residue, project, run))

def count_and_plot_res_waters(residue, HB_res_total,compare_to=None):
    hbond_count = np.zeros((5*50,2000-offset)) - 1
    x_axis = np.zeros((5*50,2000-offset))
    weights = np.zeros((5*50,2000-offset))
    column_count = np.zeros(bin_x.shape)
    for clone, traj in enumerate(HB_res_total):
        for index in range(offset,2000):
            x_axis[clone][index-offset] = index*0.25
            try:
                hbonds_frame = traj[index]
            except:
                continue
            if compare_to is None:
                count = set()
                for bond in hbonds_frame:
                    if abs(bond[0] - bond[1]) <= 2:
                        count.add(bond[0])
                    else:
                        count.add(bond[2])
                hbond_count[clone][index-offset] = len(count)
                column_count[(index-offset-0.25)/40] += 1
            else:
                reference_donors = [bond[0] for bond in compare_to[clone][index]]
                reference_acceptors = [bond[2] for bond in compare_to[clone][index]]
                count = 0
                for bond in hbonds_frame:
                    if (bond[2] in reference_donors or bond[0] in reference_donors or
                        bond[2] in reference_acceptors or bond[0] in reference_acceptors):
                        count += 1
                hbond_count[clone][index-offset] = count
                column_count[(index-offset-0.25)/40] += 1
    for clone, traj in enumerate(HB_res_total):
        for index in range(offset,2000):
            weights[clone][index-offset] = 1.00 / column_count[(index-offset-0.25)/40]
    x_axis = x_axis.flatten()
    hbond_count = hbond_count.flatten()
    weights = weights.flatten()
    plot_2dhist(residue, x_axis, hbond_count, weights, 0, project)

for i, project in enumerate(projects):
    project_dir = project_dirs[project]
    HB_total = dict()
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
        count_and_plot_res_waters(key, HB_total[key])

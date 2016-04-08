"""
PREFER TO USE RUN0_waters_2dhists.py 

Identify how many hydrogen bonds are forming on residues 181, 185, 274, and 275
in each frame.
Creates one 2D hist per residue per project

Can be modified to consider only hydrogen bonds with molecules also bound to 185
"""
import numpy as np
import sys
import math
import plot_function
        # make plots of all data past t = 250ns
            # quantify how much P(salt bridge) and (1-P)
                # P = Ntraj_frames_bound

        # look at trajectories : how long do waters stay in place
            # may need to account for the waters exchanging

offset = plot_function.OFFSET

residues_with_H = [185,181,274,275]
#residues_with_H = [181]
reference = 185

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

def plot_2dhist(residue, x_axis, hbond_count, weights, run, project):
    title = 'AURKA %s number of hydrogen bonds on residue %s over time %s' % (mutant['RUN%s' % run], residue, system[project])
    filename = "../plots/AURKA-%s-hbonds-hist2d-entire-traj-%s-combined-RUN%s" % (residue, project, run)
    ylabel = 'number of hydrogen bonds'
    plot_function.plot_2dhist(residue, x_axis, hbond_count, weights, title, ylabel, filename)

def count_and_plot_res_bonds(residue, HB_res_total, ADP_bound, compare_to=None):
    bin_x = plot_function.BIN_X
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
    x_axis = x_axis.flatten()
    hbond_count = hbond_count.flatten()
    weights = weights.flatten()
    plot_2dhist(residue, x_axis, hbond_count, weights, 0, project)

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
        count_and_plot_res_bonds(key, HB_total[key], ADP_bound)

"""
Plot minimum distances between 181-185 and 181-162
"""
import os
import numpy as np
import sys
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import seaborn as sns

sns.set_style("white")
sns.set_context("poster")

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

clumps = ['+','-','+muts','-muts']
runs_in_clumps = {
    '+'    : [(11419, 0),(11419, 1),(11419, 2),(11419, 3),(11414, 0)],
    '-'    : [(11418, 0),(11418, 1),(11418, 2),(11418, 3),(11418, 4)],
    '+muts' : [(11414, 1),(11414, 2),(11414, 3),(11414, 4),(11419, 5),(11419, 6)], 
    '-muts' : [(11423, 0),(11423, 1),(11423, 2),(11423, 3),(11423, 4),(11423, 5),(11423, 6)],
}

def plot_scatter(bridge, group, bridge_total, dihedral, mutant):
    colors = matplotlib.cm.hsv(np.linspace(0, 1, len(bridge_total)))
    title = 'AURKA %s dihedral vs %s CB-CB distance %sTPX2' % (mutant, bridge, group[0])
    filename = "../plots/AURKA-dihedral-CB-%s-combined-%s%sTPX2" % (bridge, mutant, group[0])
    xlabel = 'distance r (nanometers) between residues %s and %s' % (bridge.split('-')[0], bridge.split('-')[1])
    ylabel = 'Dihedral (radians)'
    fig1 = plt.figure()
    plt.title(title)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.axis([2.3,5.0,-5,3])
    for (color, clone_dihedral,clone_bridge) in zip(colors, dihedral, bridge_total):
#    for clone_dihedral, clone_bridge in zip(dihedral, bridge_total):
        clone_dihedral = np.ravel(clone_dihedral)
        clone_bridge = np.ravel(clone_bridge)
        minlen = min(len(clone_dihedral), len(clone_bridge))
        try:
            plt.plot(clone_bridge[:minlen], clone_dihedral[:minlen],'.',color=color, alpha=0.15,ms=4)
        except:
            try:
                clone_bridge = clone_bridge[0]
                plt.plot(clone_bridge[:minlen], clone_dihedral[:minlen],'.',color=color, alpha=0.15,ms=4)
            except Exception as e:
                print(clone_dihedral)
                raise(e)

    plt.savefig(filename,dpi=300)
    plt.close(fig1)
    print('Saved %s' % filename)

for bridge in ['284-225','287-225']:
    for group in clumps:
        di_filename = "../kinalysis/results/AURKA"+group+"/dihedral_by_something.npy"
        dihedrals = np.load(di_filename)
        for i, pro_run in enumerate(runs_in_clumps[group]):
            project = str(pro_run[0])
            run = pro_run[1]
            project_dir = project_dirs[project]
            bridge_run = np.load('%s/data/%s_%s_%s_SB_total.npy' % (project_dir, project, run, bridge)).tolist()
            if i == 0:
                bridge_total = bridge_run
            else:
                try:
                    bridge_total.append(bridge_run)
                except Exception as e:
                    print(group, project, run)
                    print(type(bridge_total))
                    print(type(bridge_run))
                    raise(e)
        if len(group) == 1:
            mutant = 'WT'
        else:
            mutant = 'Q185mutants'
        plot_scatter(bridge, group, bridge_total, dihedrals, mutant)



"""
Identifies how many unique water molecules are interacting with E181 in every frame
saves 2D histogram (1 per residue per project)

Can be modified to only consider certain trajectory of waters (for instance, 
identify how many unique waters that interact with Q185 also interact with E181)

Can be modified to filter out ADP-unbound frames or not
"""
import numpy as np
import sys
import math
import os
import mdtraj as md
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

USE_ADP = True

def water_set(topology, frame, hydrogens=False):
    waters = set()

    [[waters.add(atom) for atom in [bond[0],bond[2]]
              if topology.top.atom(atom).residue.is_water]
             for bond in frame]
    if not hydrogens:
        return waters
    oxygens = list(waters)
    [[waters.add(atom.index) for atom in topology.top.atom(oatom).residue.atoms]
     for oatom in oxygens]
    return waters

def find_other_waters(frame, first_waters):
    waters = set()

    for bond in frame:
        if bond[0] in first_waters:
            waters.add(bond[2])
        if bond[2] in first_waters:
            waters.add(bond[0])
    return waters

def plot_2dhist(residue, x_axis, hbond_count, weights, run, project):
    title = 'Waters interacting with AURKA %s residue %s over time %s' % (mutant['RUN%s' % run], residue, system[project])
    if USE_ADP:
        filename = "../plots/AURKA-%s-waters-hist2d-entire-traj-%s-combined-RUN%s_ADPfiltered" % (residue, project, run)
    else:
        filename = "../plots/AURKA-%s-waters-hist2d-entire-traj-%s-combined-RUN%s_nofilter" % (residue, project, run)
    ylabel = 'number of water molecules'
    plot_function.plot_2dhist(residue, x_axis, hbond_count, weights, title, ylabel, filename)

def count_and_plot_res_waters(residue, HB_res_total, ADP_bound = None, compare_to=None):
    bin_x = plot_function.BIN_X
    hbond_count = np.zeros((5*50,2000-offset)) - 1
    x_axis = np.zeros((5*50,2000-offset))
    weights = np.zeros((5*50,2000-offset))
    column_count = np.zeros(bin_x.shape)
    unbound = 0
    bound = 0
    for clone, traj in enumerate(HB_res_total):
        if clone == 0:
            topology = md.load('/cbio/jclab/projects/AURKA_UMN/trajectories/%s_RUN%s.pdb' % (project, run))
        for index in range(offset,2000):
            x_axis[clone][index-offset] = index*0.25
            if USE_ADP and not ADP_bound[clone][index]:
                unbound += 1
                continue
            try:
                hbonds_frame = traj[index]
                bound += 1
            except:
                continue
            if compare_to is None:
                waters = water_set(topology, hbonds_frame)
                hbond_count[clone][index-offset] = len(waters)
                column_count[(index-offset-0.25)/40] += 1
            else:
                reference = water_set(topology, compare_to[clone][index])
                waters = find_other_waters(hbonds_frame, reference)
                hbond_count[clone][index-offset] = len(waters)
                column_count[(index-offset-0.25)/40] += 1
            #if residue == 181 and index%10 == 0 and index-offset < 200 and hbond_count[clone][index-offset] > 4:
            #    print('PROJECT %s: FOUND A CLONE WITH TOO MANY WATERS: CLONE %s AT FRAME %s' % (project, clone, index))
            #    return
    print('%s frames found unbound' % unbound)
    print('%s frames found bound' % bound)
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
    if USE_ADP:
        ADP_bound = np.load('%s/is-ADP-bound.npy' % project_dir)
    else:
        ADP_bound = None
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
        count_and_plot_res_waters(key, HB_total[key], ADP_bound)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
import math
from matplotlib.pyplot import cm
import seaborn as sns
import os
import mdtraj as md

sns.set_style("white")
sns.set_context("poster")

        # make plots of all data past t = 250ns
            # quantify how much P(salt bridge) and (1-P)
                # P = Ntraj_frames_bound

        # look at trajectories : how long do waters stay in place
            # may need to account for the waters exchanging

OFFSET = 400

residues_with_H = [185,181,274,275]
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

bin_x = np.arange(OFFSET/4,510,10) - 0.25
bin_y = np.arange(8) - 0.5

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

def plot_2dhist(x_axis, hbond_count, weights, title, filename):
    fig1 = plt.figure()
    plt.hist2d(x_axis[hbond_count > -1],hbond_count[hbond_count > -1],bins=[bin_x,bin_y],weights=weights[hbond_count > -1],cmap=plt.get_cmap('jet'))
    plt.title(title)
    plt.ylabel('number of waters found')
    plt.xlabel('t (nanoseconds)')
    plt.colorbar()
    plt.axis([OFFSET/4,500,-0.5,6.5])
    plt.savefig(filename,dpi=300)
    plt.close(fig1)
    print('Saved %s' % filename)

def find_W1(HB_total, WB_total):
    HB_res_total = HB_total[274]
    W1s = np.empty((5*50,2000-OFFSET),dtype=set)
    hbond_count = np.zeros((5*50,2000-OFFSET)) - 1
    x_axis = np.zeros((5*50,2000-OFFSET))
    weights = np.zeros((5*50,2000-OFFSET))
    column_count = np.zeros(bin_x.shape)
    for clone, traj in enumerate(HB_res_total):
        if clone == 0:
            topology = md.load('/cbio/jclab/projects/AURKA_UMN/trajectories/%s_RUN%s.pdb' % (project, 0))
        for index in range(OFFSET,2000):
            x_axis[clone][index-OFFSET] = index*0.25
            try:
                this_frame = traj[index]
                column_count[(index-OFFSET-0.25)/40] += 1
            except:
                continue
            waters = water_set(topology, this_frame)
            W1s[clone][index-OFFSET] = waters
            hbond_count[clone][index-OFFSET] = len(waters)
    for clone, traj in enumerate(HB_res_total):
        for index in range(OFFSET,2000):
            weights[clone][index-OFFSET] = 1.00 / column_count[(index-OFFSET-0.25)/40]
    hbond_count = find_W2(HB_total, WB_total, hbond_count, W1s)
    x_axis = x_axis.flatten()
    hbond_count = hbond_count.flatten()
    weights = weights.flatten()
    title = 'Possible W1 identified on AURKA %s over time %s' % (mutant['RUN%s' % 0], system[project])
    filename = "./plots/W1-AURKA-hist2d-entire-traj-%s-combined-RUN%s.png" % (project, 0)
    plot_2dhist(x_axis, hbond_count, weights, title, filename)

def find_W2(HB_total, WB_total, W1_hbond_count, W1s):
    W2s = np.empty((5*50,2000-OFFSET),dtype=set)
    hbond_count = np.zeros((5*50,2000-OFFSET)) - 1
    x_axis = np.zeros((5*50,2000-OFFSET))
    weights = np.zeros((5*50,2000-OFFSET))
    column_count = np.zeros(bin_x.shape)

    # Make W2 anything that is bound to a W1
    for clone, traj in enumerate(WB_total):
        if clone == 0:
            topology = md.load('/cbio/jclab/projects/AURKA_UMN/trajectories/%s_RUN%s.pdb' % (project, 0))
        for index in range(OFFSET,2000):
            x_axis[clone][index-OFFSET] = index*0.25
            try:
                this_frame = traj[index]
                column_count[(index-OFFSET-0.25)/40] += 1
            except:
                continue
            if this_frame.shape[0] == 0:
                hbond_count[clone][index-OFFSET] = 0
                W1_hbond_count[clone][index-OFFSET] = 0
                continue
            waters = water_set(topology, this_frame)
            reference = W1s[clone][index-OFFSET]
            waters = reference.intersection(waters)
            W1_hbond_count[clone][index-OFFSET] = len(waters)
#            frame_185 = HB_total[185][clone][index]
#            waters_185 = water_set(topology, frame_185)
            frame_181 = HB_total[181][clone][index]
            waters_181 = water_set(topology, frame_181)
#            waters = waters.intersection(waters_185)
            waters = waters.intersection(waters_181)
            W2s[clone][index-OFFSET] = waters            
            hbond_count[clone][index-OFFSET] = len(W2s[clone][index-OFFSET])
    for clone, traj in enumerate(WB_total):
        for index in range(OFFSET,2000):
            weights[clone][index-OFFSET] = 1.00 / column_count[(index-OFFSET-0.25)/40]
    x_axis = x_axis.flatten()
    hbond_count = hbond_count.flatten()
    weights = weights.flatten()
    title = 'Possible W2 identified on AURKA %s over time %s' % (mutant['RUN%s' % 0], system[project])
    filename = "./plots/W2-AURKA-hist2d-entire-traj-%s-combined-RUN%s.png" % (project, 0)
    plot_2dhist(x_axis, hbond_count, weights, title, filename)
    return W1_hbond_count

def find_hbonds_between_waters(HB_total):
    from msmbuilder import dataset
    HB_res_total = HB_total[274]
    WB_total = list()

    for clone, traj in enumerate(HB_res_total):
        if clone%50 == 0:
            print('Now loading trajectories for RUN%s' % str(clone/50))
            trajectories = dataset.MDTrajDataset("/cbio/jclab/projects/fah/fah-data/munged2/all-atoms/%s/run%d-clone*.h5" % (project, clone/50))
            topology = md.load('/cbio/jclab/projects/AURKA_UMN/trajectories/%s_RUN%s.pdb' % (project, 0))
        trajectory = trajectories[clone%50]
        old_chunk = 0
        hbonds = None
        print('Finding RUN%s clone%s hbonds' % (str(clone/50),str(clone%50)))
        for chunk in range(100,2100,100):
            if chunk > 2000:
                break
            waters = set()
            chunk_frames = range(old_chunk,chunk)
            for index in chunk_frames:
                try:
                    frame_274 = traj[index]
                    frame_185 = HB_total[185][clone][index]
                    frame_181 = HB_total[181][clone][index]
                except:
                    chunk_frames = range(old_chunk,index)
                    break
                waters = waters.union(water_set(topology, frame_274, hydrogens=True))
                waters = waters.union(water_set(topology, frame_185, hydrogens=True))
                waters = waters.union(water_set(topology, frame_181, hydrogens=True))
            hbond_chunk = md.wernet_nilsson(trajectory[chunk_frames], exclude_water=False, proposed_donor_indices=waters, proposed_acceptor_indices=waters)
            if hbonds is None:
                hbonds = hbond_chunk
            else:
                try:
                    hbonds.extend(hbond_chunk)
                except Exception as e:
                    print(len(hbonds))
                    print(len(hbond_chunk))
                    raise(e)
            old_chunk = chunk
        WB_total.append(hbonds)
    np.save('%s/data/%s_%s_IntraWater.npy' % (project_dir, project, 0),WB_total)
    return WB_total

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
    if not os.path.exists('%s/data/%s_%s_IntraWater.npy' % (project_dir, project, 0)):
        WB_total = find_hbonds_between_waters(HB_total)
    else:
        WB_total = np.load('%s/data/%s_%s_IntraWater.npy' % (project_dir, project, 0))
    find_W1(HB_total, WB_total)

"""
DIFFERENT FROM SEARCH: ALL WATERS, CHANGED DEFS

Attempt to identify W1 and W2 via definitions
W1 is bound to 274 and ADP
W2 is bound to 185 and 181

DESIRED NEW DEFINITIONS:
W1 is bound to 274 backbone N and W2
W2 is bound to 275 backbone N and W1

Plot 2D hists: W1 per frame, W2 per frame, per project
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

#OFFSET = plot_function.OFFSET
OFFSET = 0

residues_with_H = [185,181,274,275]
reference = 185

local_path = os.path.dirname(os.path.realpath(__file__))

projects = ['11414','11418','11419','11423','11424','11425']
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
            run = run[3:]
            run_guide[project] += 1
            mutants[(project, run)] = mutant
        except:
            pass
print(run_guide)
print(len(run_guide))
# remove condition 11414 because it will be analyzed with condition 11419
# remove condition 11423 because there are no WT runs in it
projects = ['11418','11419','11424','11425']

try:
    this_run = int(sys.argv[1]) # 1-4
except:
    this_run = None
if this_run is not None:
    projects = [projects[this_run%4]]
if '11419' in projects:
    projects.append('11414')

def water_set(topology, frame, hydrogens=False):
    waters = set()

    [[waters.add(int(atom)) for atom in [bond[0],bond[2]]
      if topology.top.atom(int(atom)).residue.is_water]
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
        if int(bond[0]) in first_waters:
            waters.add(int(bond[2]))
        if int(bond[2]) in first_waters:
            waters.add(int(bond[0]))
    return waters

def plot_2dhist(x_axis, hbond_count, weights, title, filename):
    key = 'W1W2'
    ylabel = 'number of waters found'
    plot_function.plot_2dhist(key, x_axis, hbond_count, weights, title, ylabel, filename)

def find_W1(HB_total, ADP_bound):
    bin_x = np.arange(OFFSET/4,510,10) - 0.25
    HB_res_total = HB_total[274]
    W1s = np.empty((5*50,2000-OFFSET),dtype=set)
    hbond_count = np.zeros((5*50,2000-OFFSET)) - 1
    x_axis = np.zeros((5*50,2000-OFFSET))
    weights = np.zeros((5*50,2000-OFFSET))
    column_count = np.zeros(bin_x.shape)
    for clone, traj in enumerate(HB_res_total):
        if clone%50 == 0:
            if project == '11419' and clone/50 == 4:
                project = '11414'
                project_dir = project_dirs[project]
            water_to_adp = np.load('%s/data/%s_%s_adp-hbond-water.npy' % (project_dir, project, clone/50))
        if clone == 0:
            topology = md.load('/cbio/jclab/conditions/AURKA_UMN/trajectories/%s_RUN%s.pdb' % (project, 0))
        for index in range(OFFSET,2000):
            x_axis[clone][index-OFFSET] = index*0.25
#            if not ADP_bound[clone][index]:
            if not ADP_bound:
                continue
            try:
                this_frame = traj[index]
                column_count[(index-OFFSET-0.25)/40] += 1
            except:
                continue
            waters = water_set(topology, this_frame)
            adp_waters = water_set(topology, water_to_adp[clone%50][index])
            waters = waters.intersection(adp_waters)
            W1s[clone][index-OFFSET] = waters
            hbond_count[clone][index-OFFSET] = len(waters)
    for clone, traj in enumerate(HB_res_total):
        for index in range(OFFSET,2000):
            weights[clone][index-OFFSET] = 1.00 / column_count[(index-OFFSET-0.25)/40]
    x_axis = x_axis.flatten()
    hbond_count = hbond_count.flatten()
    weights = weights.flatten()
    title = 'Possible W1 identified on AURKA %s over time %s' % (mutants[(project, 0)], system[project])
    filename = "/cbio/jclab/conditions/behrj/AURKA_UMN/plots/W1-AURKA-hist2d-entire-traj-%s-combined-RUN%s.png" % (project, 0)
    plot_2dhist(x_axis, hbond_count, weights, title, filename)
    np.save("%s/data/W1-oxygen-indices.npy" % project_dirs[project], W1s)

def find_W2(HB_total, ADP_bound):
    bin_x = np.arange(OFFSET/4,510,10) - 0.25
    #bin_x = plot_function.BIN_X
    W2s = np.empty((5*50,2000-OFFSET),dtype=set)
    hbond_count = np.zeros((5*50,2000-OFFSET)) - 1
    x_axis = np.zeros((5*50,2000-OFFSET))
    weights = np.zeros((5*50,2000-OFFSET))
    column_count = np.zeros(bin_x.shape)

    for clone, traj in enumerate(HB_total[185]):
        if clone == 0:
            topology = md.load('/cbio/jclab/conditions/AURKA_UMN/trajectories/%s_RUN%s.pdb' % (project, 0))
        for index in range(OFFSET,2000):
            x_axis[clone][index-OFFSET] = index*0.25
            #if not ADP_bound[clone][index]:
            if not ADP_bound:
                continue
            try:
                this_frame = traj[index]
                column_count[(index-OFFSET-0.25)/40] += 1
            except:
                continue
            if this_frame.shape[0] == 0:
                hbond_count[clone][index-OFFSET] = 0
                continue
            waters = water_set(topology, this_frame)
            frame_181 = HB_total[181][clone][index]
            waters_181 = water_set(topology, frame_181)
            waters = waters.intersection(waters_181)
            W2s[clone][index-OFFSET] = waters            
            hbond_count[clone][index-OFFSET] = len(W2s[clone][index-OFFSET])
    for clone, traj in enumerate(HB_total[185]):
        for index in range(OFFSET,2000):
            weights[clone][index-OFFSET] = 1.00 / column_count[(index-OFFSET-0.25)/40]
    x_axis = x_axis.flatten()
    hbond_count = hbond_count.flatten()
    weights = weights.flatten()
    title = 'Possible W2 identified on AURKA %s over time %s' % (mutants[(project, 0)], system[project])
    filename = "/cbio/jclab/conditions/behrj/AURKA_UMN/plots/W2-AURKA-hist2d-entire-traj-%s-combined-RUN%s.png" % (project, 0)
    plot_2dhist(x_axis, hbond_count, weights, title, filename)
    np.save("%s/data/W2-oxygen-indices.npy" % project_dirs[project], W2s)

def find_hbonds_between_waters(HB_total):
    from msmbuilder import dataset
    HB_res_total = HB_total[274]
    WB_total = list()

    for clone, traj in enumerate(HB_res_total):
        if clone%50 == 0:
            print('Now loading trajectories for RUN%s' % str(clone/50))
            trajectories = dataset.MDTrajDataset("/cbio/jclab/conditions/fah/fah-data/munged3/all-atoms/%s/run%d-clone*.h5" % (project, clone/50))
            topology = md.load('/cbio/jclab/conditions/AURKA_UMN/trajectories/%s_RUN%s.pdb' % (project, 0))
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

HB_total = dict()
for i, project in enumerate(projects):
    project_dir = project_dirs[project]
    #ADP_bound = np.load('%s/is-ADP-bound.npy' % project_dir) #--> what to do about this
    ADP_bound = True
    for run in range(run_guide[project]):
        if not mutants[(project, run)] == 'WT': continue
        for residue in residues_with_H:
            if not os.path.exists('%s/data/%s_%s_%s_distHBonds.npy' % (project_dir, project, run, residue)):
                continue
            if not HB_total.has_key(residue):
                HB_total[residue] = np.load('%s/data/%s_%s_%s_distHBonds.npy' % (project_dir, project, run, residue))
            else:
                new_run = np.load('%s/data/%s_%s_%s_distHBonds.npy' % (project_dir, project, run, residue))
                HB_total[residue] = np.concatenate((HB_total[residue], new_run))
project = projects[0]
find_W1(HB_total, ADP_bound)
find_W2(HB_total, ADP_bound)

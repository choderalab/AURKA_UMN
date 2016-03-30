import matplotlib
matplotlib.use('Agg')
import numpy as np
import sys
import math
import os
import mdtraj as md

        # make plots of all data past t = 250ns
            # quantify how much P(salt bridge) and (1-P)
                # P = Ntraj_frames_bound

        # look at trajectories : how long do waters stay in place
            # may need to account for the waters exchanging

OFFSET = 400

residues_with_H = [185,181,274,275]
reference = 185

projects = ['11410','11411']
project_dirs = {'11410':'/cbio/jclab/projects/behrj/AURKA_UMN/output-1OL5','11411':'/cbio/jclab/projects/behrj/AURKA_UMN/output-1OL7'}
system = {'11410':'with TPX2','11411':'without TPX2'}

with open('/cbio/jclab/projects/behrj/AURKA_UMN/output-1OL7/run-index.txt','r') as fi:
    run_index = fi.read()
mutant = dict()
for entry in run_index.split('\n'):
    try:
        mutant[entry.split(' ')[0]] = entry.split(' ')[1]
    except:
        pass

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
    print(project)
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
    WB_total = find_hbonds_between_waters(HB_total)



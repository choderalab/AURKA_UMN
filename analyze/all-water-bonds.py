import numpy as np
import sys
import math
import os
import mdtraj as md
from msmbuilder import dataset
from itertools import chain

projects = ['11410','11411']
project_dirs = {'11410':'/cbio/jclab/projects/behrj/AURKA_UMN/output-1OL5','11411':'/cbio/jclab/projects/behrj/AURKA_UMN/output-1OL7'}
system = {'11410':'with TPX2','11411':'without TPX2'}
runs = range(5)

with open('/cbio/jclab/projects/behrj/AURKA_UMN/output-1OL7/run-index.txt','r') as fi:
    run_index = fi.read()

mutant = dict()
for entry in run_index.split('\n'):
    try:
        mutant[entry.split(' ')[0]] = entry.split(' ')[1]
    except:
        pass

DIST_FOR_HBOND = False

verbose = True
overwrite = False

def find_hbonds_for_this_traj(traj, haystack, traj_has_frames, sidechain=False, backbone=False):
    protein_atoms = find_protein_atoms(traj, sidechain=sidechain, backbone=backbone)
    if verbose:
        print('Starting search for RUN%s clone%s hbonds...' % (run, i))
    hbonds = list()
    old_chunk = 0
    keep_going = True
    for chunk in range(100,2100,100):
        if chunk > 2000:
            break
        chunk_frames = range(old_chunk,chunk)
        for index in chunk_frames:
            try:
                is_frame = traj_has_frames[index]
            except:
                chunk_frames = range(old_chunk,index)
                keep_going = False
                break
        if len(chunk_frames) == 0:
            break
        if verbose:
            print('Calculating bonds in frames %s-%s' % (old_chunk, chunk_frames[-1]))
        hbonds0 = md.wernet_nilsson(traj[chunk_frames], exclude_water=False, proposed_donor_indices=protein_atoms, proposed_acceptor_indices=haystack)
        hbonds1 = md.wernet_nilsson(traj[chunk_frames], exclude_water=False, proposed_donor_indices=haystack, proposed_acceptor_indices=protein_atoms)
        old_chunk = chunk
        for frame, bondlist in enumerate(hbonds0):
            try:
                hbonds.append(np.concatenate((bondlist,hbonds1[frame])))
            except Exception as e:
                print('hbonds0')
                print(bondlist)
                print(bondlist.shape)
                print('hbonds1')
                print(hbonds1[frame])
                print(hbonds1[frame].shape)
                raise(e)
        if not keep_going:
            break
    if verbose:
        print('Found all hbonds!')
    return hbonds, protein_atoms

def find_protein_atoms(traj, sidechain=False, backbone=False):
    top = traj.topology
    if sidechain:
        protein_atoms = [atom.index for atom in top.atoms if atom.is_sidechain]
    elif backbone:
        protein_atoms = [atom.index for atom in top.atoms if atom.is_backbone]
    else:
        protein_atoms = [atom.index for atom in top.atoms if atom.residue.is_protein]
    return protein_atoms

try:
    this_run = int(sys.argv[1]) # 1 - 10
except:
    this_run = None

if this_run is not None:
    projects = [projects[this_run%2]]
    runs = [runs[this_run%5]]

run_protein_atoms = dict()
run_water_atoms = dict()
run_all_water_hbonds = dict()

for project in projects:
    project_dir = project_dirs[project]
    for run in runs:
        traj_has_frame = np.load('%s/data/%s_%s_185_HBonds.npy' % (project_dir, project, run))
        if not overwrite and os.path.exists('%s/data/%s_%s_all-water-bonds_total.npy' % (project_dir, project, run)):
            all_water_hbonds = np.load('%s/data/%s_%s_all-water-bonds_total.npy' % (project_dir, project, run))
            trajectories = dataset.MDTrajDataset("/cbio/jclab/projects/fah/fah-data/munged2/all-atoms/%s/run%d-clone0.h5" % (project, run))
            run_protein_atoms[str(project)+str(run)] = find_protein_atoms(trajectories[0])
            run_water_atoms[str(project)+str(run)] = trajectories[0].top.select("water")
            run_all_water_hbonds[str(project)+str(run)] = all_water_hbonds
            continue
        if verbose:
            print("Loading Project %s RUN%s..." % (project, run))
        all_water_hbonds = []
        trajectories = dataset.MDTrajDataset("/cbio/jclab/projects/fah/fah-data/munged2/all-atoms/%s/run%d-clone*.h5" % (project, run))
        for i,traj in enumerate(trajectories):
            if i == 0:
                haystack = traj.top.select("water")
                run_water_atoms[str(project)+str(run)] = haystack
            hbonds, run_protein_atoms[str(project)+str(run)] = find_hbonds_for_this_traj(traj, haystack, traj_has_frame[i])
            all_water_hbonds.append(hbonds)
        np.save('%s/data/%s_%s_all-water-bonds_total.npy' % (project_dir, project, run), all_water_hbonds)
        run_all_water_hbonds[str(project)+str(run)] = all_water_hbonds

if verbose:
    print('Found and saved bonds!')

for project in projects:
    project_dir = project_dirs[project]
    for run in runs:
        if not overwrite and os.path.exists('%s/data/%s_%s_all-protein-significant-water-bonds.npy' % (project_dir, project, run)):
            significant_bonds = np.load('%s/data/%s_%s_all-protein-significant-water-bonds.npy' % (project_dir, project, run))
            continue
        hbonds = np.empty(len(run_protein_atoms[run]), dtype=dict)
        significant_bonds = dict()
        for frame_id, frame in enumerate(run_all_water_hbonds[str(project)+str(run)]):
            for bond in frame:
                index = None
                key = None
                index = [atom for atom in [bond[0],bond[2]] if atom in run_protein_atoms[str(project)+str(run)]][0]
                key = [atom for atom in [bond[0],bond[2]] if atom in run_water_atoms[str(project)+str(run)]][0]
                if index is None or key is None:
                    continue
                if hbonds[index].has_key(key):
                    if len(hbonds[index][key]) == 100:
                        if significant_bonds.has_key(key):
                            significant_bonds[key].append(index)
                        else:
                            significant_bonds[key] = [index]
                    hbonds[index][key].append(frame_id)
                else:
                    hbonds[index][key] = [frame_id]
        for key in significant_bonds.keys():
            print('Significant bonds found at water %s:' % key)
            print(significant_bonds[key])
        np.save('%s/data/%s_%s_all-protein-significant-water-bonds.npy' % (project_dir, project, run), significant_bonds)

if verbose:
    print('Complete!')

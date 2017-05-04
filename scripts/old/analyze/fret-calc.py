import numpy as np
import sys
import math
import os
import mdtraj as md
from msmbuilder import dataset
from itertools import chain

local_path = os.path.dirname(os.path.realpath(__file__))

projects = ['11410','11411','11418']
project_dirs = {'11410':'%s/../output-1OL5' % local_path,'11411':'%s/../output-1OL7' % local_path,'11418':'%s/../output1OL5-TPX2' % local_path}
system = {'11410':'with TPX2','11411':'without TPX2','11418': 'with TPX2 removed'}
projects = ['11410','11411','11418']
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
    '11424':'with TPX2',
    '11425':'with TPX2 removed'
}

projects = ['11414','11418','11419','11423']

run_guide = list()
#mutants = dict()
for project in projects:
    filename = project_dirs[project]+'/run-index.txt'
    with open(filename, 'r') as fi:
        project_run_index = fi.read()
    for entry in project_run_index.split('\n'):
        try:
            run = entry.split(' ')[0]
            mutant = entry.split(' ')[1]
            run = run[3:]
            run_guide.append([project, run])
            #mutants[(condition, run)] = mutant
        except:
            pass
print(run_guide)
print(len(run_guide))


try:
    this_run = int(sys.argv[1]) # 1 - 48
except:
    this_run = None

if this_run is not None:
    projects = [run_guide[this_run][0]]
    runs = [run_guide[this_run][1]]
else: runs = [] # so basically, never do this

verbose = True
for project in projects:
    project_dir = project_dirs[project]
    if project in ['11424','11425']:
        charmm = True
    else:
        charmm = False
    for run in runs:
        if verbose:
            print("Loading Project %s RUN%s..." % (project, run))
        SB_sl_total = []
        SB_tl_total = []
        trajectories = dataset.MDTrajDataset("/cbio/jclab/conditions/fah/fah-data/munged3/no-solvent/%s/run%s-clone*.h5" % (project, run))
        for i,traj in enumerate(trajectories):
            print(i)
            if charmm and i == 0:
                for residue in traj.topology.residues:
                    if str(residue) == 'SER162':
                        s284 = residue
                    if str(residue) == 'LEU103':
                        l225 = residue
                    if str(residue) == 'THR165':
                        t287 = residue

            elif i == 0:
                for residue in traj.topology.residues:
                    if str(residue) == 'SER284':
                        s284 = residue
                    if str(residue) == 'LEU225':
                        l225 = residue
                    if str(residue) == 'THR287':
                        t287 = residue

            if i == 0:
                for atom in s284.atoms:
                    if atom.name == 'CB':
                        s284CB = atom
                for atom in l225.atoms:
                    if atom.name == 'CB':
                        l225CB = atom
                for atom in t287.atoms:
                    if atom.name == 'CB':
                        t287CB = atom
                pairssl = np.zeros((1,2),dtype=int)
                pairstl = np.zeros((1,2),dtype=int)
                pairssl[0][0] = s284CB.index
                pairssl[0][1] = l225CB.index
                pairstl[0][0] = t287CB.index
                pairstl[0][1] = l225CB.index
                print(pairstl)
                print(pairssl)

            distances_284_225 = md.compute_distances(traj, pairssl)
            distances_287_225 = md.compute_distances(traj, pairstl)
            print(distances_284_225[36])
            print(distances_287_225[51])

            SB_sl_total.append(distances_284_225)
            SB_tl_total.append(distances_287_225)
        np.save('%s/data/%s_%s_284-225_SB_total.npy' % (project_dir, project, run), SB_sl_total)
        np.save('%s/data/%s_%s_287-225_SB_total.npy' % (project_dir, project, run), SB_tl_total)

if verbose:
    print('Complete!')

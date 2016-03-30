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

sns.set_style("whitegrid")
sns.set_context("poster")

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

residues = [181,185,274,275]
res_names = {181:'GLU',185:'GLN',274:'ASP',275:'PHE'}
residue_atoms = dict()

for project in projects:
    project_dir = project_dirs[project]
    run = 0
    for residue in residues:
        Save_waters = list()
        HB_total = np.load('%s/data/%s_%s_%s_HBonds.npy' % (project_dir, project, run, residue))

        for clone_id, hbonds_traj in enumerate(HB_total):
            num_frames = len(hbonds_traj)
            if clone_id == 0:
                topology = md.load('/cbio/jclab/projects/AURKA_UMN/trajectories/%s_RUN%s.pdb' % (project, run))
                for res in topology.top.residues:
                    if str(res) == res_names[residue]+str(residue):
                        print(str(res))
                        break
                residue_atoms[residue] = [atom.index for atom in res.atoms]
                print('Indices of atoms in this residue:')
                print(residue_atoms[residue])
            #if os.path.exists("%s/hbonds/AURKA-%s-RUN%s-clone%s-%s-significant-waters.png" % (project_dir, project, run, clone_id, residue)):
            #    continue
            tracking_waters = dict()
            for frame_id, frame in enumerate(hbonds_traj):
                if frame.shape[0] > 0:
                    for bond_num, bond_in_frame in enumerate(frame):
                        this_water = [atom for atom in [bond_in_frame[0],bond_in_frame[2]]
                                           if topology.top.atom(atom).residue.is_water]
                        assert len(this_water) == 1
                        this_water = this_water[0]
                        if not tracking_waters.has_key(this_water):
                            tracking_waters[this_water] = [frame_id]
                        else:
                            tracking_waters[this_water].append(frame_id)

            print('%s waters are tracked' % len(tracking_waters))

            for key in tracking_waters.keys():
                tracking_waters[key] = np.asarray(tracking_waters[key])

            plt.figure()
            significant_waters = dict()
            for water_id, (this_water, water_bonds) in enumerate(tracking_waters.items()):
                try:
                    y_axis = np.ones(len(water_bonds)) * (water_id + 1)
                    plt.plot(water_bonds*0.25, y_axis, '.')
                except Exception as e:
                    print(clone_id)
                    print(this_water)
                    print('Failed to plot %s RUN%s clone%s water %s' % (project, run, clone_id, water_id))
                    #raise(e)
                    print(str(e))
                    continue
                if len(water_bonds) >= num_frames * 0.01:
                    significant_waters[this_water] = water_bonds
            plt.xlabel('t (nanoseconds)')
            plt.ylabel('water molecule')
            plt.title('Hydrogen bonds in AURKA %s %s clone%s' % (mutant['RUN%s' % run], system[project], clone_id))
            plt.savefig("%s/hbonds/AURKA-%s-RUN%s-clone%s-%s-ALL-waters.png" % (project_dir, project, run, clone_id, residue),dpi=300)
            plt.close()

            print('Number of significant waters: %s' % len(significant_waters))
            plt.figure()
            for water_id, (this_water, water_bonds) in enumerate(significant_waters.items()):
                y_axis = np.ones(len(water_bonds)) * (water_id + 1)
                plt.plot(water_bonds*0.25, y_axis, '.')
            plt.xlabel('t (nanoseconds)')
            plt.ylabel('water molecule')
            plt.title('Significant (more than 1percent) hydrogen bonds in %s AURKA %s %s clone%s' % (mutant['RUN%s' % run], system[project], clone_id, residue))
            plt.savefig("%s/hbonds/AURKA-%s-RUN%s-clone%s-%s-significant-waters.png" % (project_dir, project, run, clone_id, residue),dpi=300)
            plt.close()

            Save_waters.append([significant_waters.keys()])
        np.save("%s/hbonds/AURKA-%s-RUN%s-%s-significant-waters.npy" % (project_dir, project, run, residue), Save_waters)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
import math
from matplotlib.pyplot import cm
import seaborn as sns
import os

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

for project in projects:
    project_dir = project_dirs[project]
    for run in range(5):
        HB_total = np.load('%s/data/%s_%s_HBonds.npy' % (project_dir, project, run))

        for clone_id, hbonds_traj in enumerate(HB_total):
            identified_hydrogens = set()
            for frame in hbonds_traj:
                num_hbonds_in_frame_in_traj = frame.shape[0]
                if num_hbonds_in_frame_in_traj > 0:
                    for bond_in_frame in frame:
                        identified_hydrogens.add(bond_in_frame[1])

            identified_hydrogens = list(identified_hydrogens)
            print('%s total hydrogens are found' % len(identified_hydrogens))

            tracking_hydrogens = np.zeros((len(identified_hydrogens),len(hbonds_traj)))
            for frame_id, frame in enumerate(hbonds_traj):
                num_hbonds_in_frame_in_traj = frame.shape[0]
                if num_hbonds_in_frame_in_traj > 0:
                    for bond_in_frame in frame:
                        try:
                            assert bond_in_frame[1] in identified_hydrogens
                            assert identified_hydrogens.index(bond_in_frame[1]) in range(tracking_hydrogens.shape[0])
                            tracking_hydrogens[identified_hydrogens.index(bond_in_frame[1])][frame_id] = identified_hydrogens.index(bond_in_frame[1]) + 1
                        except Exception as e:
                            print('Dimensions of tracking_hydrogens: (%s,%s)' % (tracking_hydrogens.shape[0],tracking_hydrogens.shape[1]))
                            print('Index of this Hydrogen: %s' % identified_hydrogens.index(bond_in_frame[1]))
                            print('Index of this frame: %s' % frame_id)
                            raise(e)

            plt.figure()
            for h_id, hydrogen in enumerate(tracking_hydrogens):
                try:
                    plt.plot(hydrogen, '.')
                except Exception as e:
                    print(hydrogen)
                    print('Failed on %s RUN%s clone%s hydrogen %s' % (project, run, clone_id, h_id))
                    #raise(e)
                    print(str(e))
                    continue
            plt.xlabel('frame')
            plt.ylabel('hydrogen atom')
            plt.title('Hydrogen bonds in AURKA %s %s clone%s' % (mutant['RUN%s' % run], system[project], clone_id))
            plt.savefig("%s/hbonds/AURKA-%s-RUN%s-clone%s-ALL-HBonds.png" % (project_dir, project, run, clone_id),dpi=300)
            plt.close()

            significant_count = 0
            for index, hydrogen in enumerate(tracking_hydrogens):
                marker_value = index+1
                if hydrogen.mean() > marker_value * 0.01:
                    significant_count += 1
            significant_hydrogens = np.zeros((significant_count,len(hbonds_traj)))
            print('Number of significant hydrogens: %s' % significant_count)

            significant_found = 0
            for index, hydrogen in enumerate(tracking_hydrogens):
                marker_value = index+1
                if hydrogen.mean() > marker_value * 0.01:
                    significant_hydrogens[significant_found] = hydrogen
                    significant_found+=1
            for index, sig_hydrogen in enumerate(significant_hydrogens):
                for frame_id, frame in enumerate(sig_hydrogen):
                    if frame > 0:
                        significant_hydrogens[index][frame_id] = index + 1

            plt.figure()
            for hydrogen in  significant_hydrogens:
                plt.plot(hydrogen, '.')
            plt.xlabel('frame')
            plt.ylabel('hydrogen atom')
            plt.title('Significant (more than 1percent) hydrogen bonds in AURKA %s %s clone%s' % (mutant['RUN%s' % run], system[project], clone_id))
            plt.savefig("%s/hbonds/AURKA-%s-RUN%s-clone%s-significant-HBonds.png" % (project_dir, project, run, clone_id),dpi=300)
            plt.close()

            np.save("%s/hbonds/AURKA-%s-RUN%s-clone%s-significant-hydrogens" % (project_dir, project, run, clone_id), significant_hydrogens)

            del(hydrogen)
            del(significant_hydrogens)
            del(index)
            del(frame_id)
            del(sig_hydrogen)
            del(significant_found)
            del(marker_value)
            del(significant_count)
            del(tracking_hydrogens)
            del(identified_hydrogens)
            try:
                del(identified_hydrogens)
            except:
                pass
            del(bond_in_frame)
            del(frame)
            del(num_hbonds_in_frame_in_traj)
            del(hbonds_traj)

            # only do clone 0 until everthing is right
            #break

        del(HB_total)

    # only do 11410 until 11411 is ready
    #break

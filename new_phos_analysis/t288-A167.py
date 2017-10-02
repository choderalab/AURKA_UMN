# A script to compute the fraction of a given contact in all of your frames + the standard error

import matplotlib
import sys
import math

matplotlib.use('Agg')
from matplotlib.pyplot import cm
import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np
from msmbuilder import dataset
import seaborn as sns

sns.set_style("whitegrid")
sns.set_context("poster")

# Define condition.
condition = sys.argv[2]

# Define kinase.
kinase_definition = sys.argv[1]



# Define Alpha carbon coordinates (1-indexed), give me residue numbers here
res_pairs = {'AURKA': [[288, 167]]}




def alpha_distances(traj, residue_pair):
    """
    :param traj: the trajectory to be analyzed
    :param spine: the residues involved
    :return: two flattened numpy arrays
    """
    min_frame = 400
    end_frame = len(traj)
    short_traj = traj.slice(range(min_frame, end_frame), copy=False)
    atom1 = short_traj.topology.select("residue %s and name == 'CA'" % residue_pair[0])
    atom2 = short_traj.topology.select("residue %s and name == 'CA'" % residue_pair[1])
    list_of_atoms = [atom1, atom2]
    atom_array = np.asanyarray(list_of_atoms)
    atom_array = atom_array.reshape(1,2)
    dists = md.compute_distances(short_traj, atom_array)



    # Append difference and individual distances
    dists = np.multiply(dists, 10)
    dists = dists.flatten()

    # flatten list of arrays
    return dists


if __name__ == "__main__":
    trajectories = dataset.MDTrajDataset(
        '/cbio/jclab/home/albaness/trajectories2/AURKA/%s/*/*.h5' % condition)
    master_dist_list = []
    for pair in range(len(res_pairs[kinase_definition])):
        dist_list = []
        pair_list = res_pairs[kinase_definition][pair]
        res1 = pair_list[0]
        res2 = pair_list[1]
        print(pair_list)
        for traj_in in trajectories:
            distance1 = alpha_distances(traj_in, pair_list)
            dist_list.extend(distance1)
        np.save('./data/distances/distances_%s_%s-pair%s-%s.npy' % (kinase_definition, condition, res1, res2), dist_list)
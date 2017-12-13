# A script to compute the fraction of a given contact in all of your frames + the standard error

import sys
import mdtraj as md
import numpy as np
from msmbuilder import dataset


# Define condition.
condition = sys.argv[1]

# Define kinase offset.
offset = int(sys.argv[2])

# Define Alpha carbon coordinates (1-indexed), give me residue numbers here
res_pairs = {'AURKA': [[288, 180]]}




def alpha_distances(traj, residue_pair):
    """
    :param traj: the trajectory to be analyzed
    :param spine: the residues involved
    :return: two flattened numpy arrays
    """
    min_frame = 400
    end_frame = len(traj)
    frames_to_slice = list(range(min_frame, end_frame))
    short_traj = traj.slice(frames_to_slice, copy=False)
    atom1 = short_traj.topology.select("residue %s and name == 'CB'" % residue_pair[0])
    atom2 = short_traj.topology.select("residue %s and name == 'CZ'" % residue_pair[1])
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
        '/cbio/jclab/home/albaness/trajectories2/AURKA/CHARMM/%s/*.h5' % condition)

    master_dist_list = []
    for pair in range(len(res_pairs['AURKA'])):
        dist_list = []
        pair_list = res_pairs['AURKA'][pair]
        res1 = pair_list[0] - offset
        res2 = pair_list[1] - offset
        pair_list = [res1, res2]
        pair_list_offset = [res1, res2]
        for traj_in in trajectories:
            distance1 = alpha_distances(traj_in, pair_list_offset)
            dist_list.extend(distance1)
            np.save(
                '../data/distances/distances_CHARMM_AURKA_%s-pair%s-%s.npy' % (condition, pair_list[0], pair_list[1]),
                dist_list)
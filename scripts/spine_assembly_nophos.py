# A script to compute the fraction of a given contact in all of your frames + the standard error

# import libraries
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

# Define project.
project_num = sys.argv[1]

# Define kinase.
kinase_definition = sys.argv[2]



# Define hydrogen bond coordinates (0-indexed)
spine = {'AURKA': [[73, 62], [62, 152], [152, 131]]}




def spine_distances(traj, spine):
    """
    :param traj: the trajectory to be analyzed
    :param spine: the residues involved
    :return: two flattened numpy arrays
    """
    min_frame = 0
    end_frame = len(traj)
    flat_spine = [item for sublist in spine for item in sublist]
    sidechains = []
    short_traj = traj.slice(range(min_frame, end_frame), copy=False)

    for resid in flat_spine:
        atoms = short_traj.topology.select("resid %s and sidechain" % resid)
        sidechains.extend(atoms.tolist())

    short_traj.atom_slice(sidechains, inplace=True)

    [d2358_2202, res_list_one] = md.compute_contacts(short_traj, [[1,0]])
    [d2202_2222, res_list_two] = md.compute_contacts(short_traj, [[0,3]])
    [d2222_2326, res_list_two] = md.compute_contacts(short_traj, [[3,2]])
    #[dist4_nm, res_list_three] = md.compute_contacts(short_traj, [spine[3]])

    # Append difference and individual distances
    dist1= np.multiply(d2358_2202, 10)
    dist2 = np.multiply(d2202_2222, 10)
    dist3 = np.multiply(d2222_2326, 10)
    #dist4 = np.multiply(dist4_nm, 10)

    # flatten list of arrays
    return [dist1, dist2, dist3]


def stat_analyze(distances, window, cutoff):
    tmax = 0
    ntraj = len(distances)  # number of trajectories

    # Compute the maximum time
    for n in range(ntraj):
        if len(distances[n]) > tmax:
            tmax = len(distances[n])
            print('T max for run %s : %d' % (n, len(distances[n])))
    # Compute the contact fraction
    contact_fraction = np.zeros([tmax - window], np.float64)
    contact_fraction_stderr = np.zeros([tmax - window], np.float64)
    for t in range(0, tmax - window):
        # Count the number of trajectoties that are 't+sliding window' in length
        ntraj_t = 0
        for n in range(ntraj):
            if len(distances[n]) >= t + window:
                ntraj_t += 1
        contact_fraction_n = np.zeros(ntraj_t)
        index = 0
        for n in range(ntraj):
            if len(distances[n]) >= t + window:
                contact_fraction_n[index] = (distances[n][t:(t + window)] < cutoff).mean()
                index += 1
        contact_fraction[t] = contact_fraction_n.mean()
        contact_fraction_stderr[t] = contact_fraction_n.std() / np.sqrt(ntraj_t)

    return contact_fraction, contact_fraction_stderr


if __name__ == "__main__":

    sliding_window = 40
    cutoff_dist = 4
    list_of_runs = [0]

    dist_list1 = []
    dist_list2 = []
    dist_list3 = []
    #dist_list4 = []
    trajectories = dataset.MDTrajDataset(
        "/cbio/jclab/home/albaness/trajectories/%s/*/*.h5" % (project_num, run))
    for traj_in in trajectories:
        [distance1, distance2, distance3] = spine_distances(traj_in, spine[kinase_definition])
        dist_list1.append(distance1[:, 0])
        dist_list2.append(distance2[:, 0])
        dist_list3.append(distance3[:, 0])
        #dist_list4.append(distance4[:, 0])
    [dist1_fraction, dist1_stderr] = stat_analyze(dist_list1, sliding_window, cutoff_dist)
    [dist2_fraction, dist2_stderr] = stat_analyze(dist_list2, sliding_window, cutoff_dist)
    [dist3_fraction, dist3_stderr] = stat_analyze(dist_list3, sliding_window, cutoff_dist)
    #[dist4_fraction, dist4_stderr] = stat_analyze(dist_list4, sliding_window, cutoff_dist)
    np.save('../data/%s_%s_dist1_fraction_%s.npy' % (project_num, mutations[run], cutoff_dist), dist1_fraction)
    np.save('../data/%s_%s_dist1_stderr_%s.npy' % (project_num, mutations[run], cutoff_dist), dist1_stderr)
    np.save('../data/%s_%s_dist2_fraction_%s.npy' % (project_num, mutations[run], cutoff_dist), dist2_fraction)
    np.save('../data/%s_%s_dist2_stderr_%s.npy' % (project_num, mutations[run], cutoff_dist), dist2_stderr)
    np.save('../data/%s_%s_dist3_fraction_%s.npy' % (project_num, mutations[run], cutoff_dist), dist3_fraction)
    np.save('../data/%s_%s_dist3_stderr_%s.npy' % (project_num, mutations[run], cutoff_dist), dist3_stderr)
    #np.save('../data/%s_%s_dist4_fraction_%s.npy' % (project_num, mutations[run], cutoff_dist), dist4_fraction)
    #np.save('../data/%s_%s_dist4_stderr_%s.npy' % (project_num, mutations[run], cutoff_dist), dist4_stderr)

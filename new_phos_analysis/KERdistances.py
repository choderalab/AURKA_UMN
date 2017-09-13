# This script allows you to make plots according to the coordinates in Fig 2 of Shukla et al Nature Communications 2014.
#
# Right now the syntax is pretty silly but for now this is how it works:
#
#     python plotting_Shukla.py <condition> <kinase> <ref_kinase>
#  Example: python plotting_Shukla.py 11401 'SRC'
#
# The only kinases currently availabe are 'SRC', 'ABL', and 'DDR1', but it should be simple to add your own.
#
# Made by Sonya Hanson Jan 2016, for better or worse.

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

# Define condition.
project_num = sys.argv[1]

# Define kinase.
kinase_definition = sys.argv[2]


# Define hydrogen bond coordinates (0-indexed)
KER_hbond = {'SRC': [[28, 43], [43, 142]],
             'ABL': [[29, 44], [44, 144]],
             'DDR1': [[51, 68], [68, 185]],
             'EGFR': [[49, 66], [66, 164]],
             'MTOR': [[5, 8], [8, 248]],
             'MTOR_long': [[811, 814], [814, 1054]],
             'AURKA' : [[39, 58], [58, 157]],
             'AURKA_phos': [[57, 165], [57, 166]]
             }


# Define mutational groups
Mut_group = {'Single_nonactivating': [0, 4, 7],
             'Single_activating': [1, 2, 3, 5, 6, 8, 9],
             'Double_nonactivating': [12, 15, 19, 22, 25, 28, 31, 32, 33, 34, 35, 37, 40, 43, 44],
             'Double_activating': [10, 11, 13, 14, 16, 17, 18, 20, 21, 23, 24, 26, 27, 29, 30, 36, 38, 39, 41, 42, 45]}


def shukla_coords(traj, KER):
    """
    :param traj: the trajectory to be analyzed
    :param KER: the residues involved
    :return: two flattened numpy arrays
    """

    min_frame = 0
    end_frame = len(traj)
    flat_KER  = [item for sublist in KER for item in sublist]
    sidechains = []
    short_traj = traj.slice(range(min_frame, end_frame), copy=False)

    for resid in flat_KER:
        atoms = short_traj.topology.select("resid %s and sidechain" % resid)
        sidechains.extend(atoms.tolist())


    short_traj.atom_slice(sidechains, inplace=True)

    [k2187e2195, res_list_one] = md.compute_contacts(short_traj, [[0,1]])
    [e2195r2430, res_list_two] = md.compute_contacts(short_traj, [[1,2]])

    # Append difference and individual distances
    K_E = np.multiply(k2187e2195, 10)
    R_E = np.multiply(e2195r2430, 10)

    # flatten list of arrays
    return [R_E, K_E]


def plot(kinase, mutgroup, project):
    RE_graph = []
    KE_graph = []
    for index in Mut_group[mutgroup]:
        trajectories = dataset.MDTrajDataset(
            "/cbio/jclab/conditions/fah/fah-data/munged2/no-solvent/%s/run%d-clone*.h5" % (project, index))
        for i, traj in enumerate(trajectories):
            if i == 0:
                [RE, KE] = shukla_coords(traj, KER_hbond[kinase])
                RE_graph = list(RE[:, 0])
                KE_graph = list(KE[:, 0])
            else:
                [RE, KE] = shukla_coords(traj, KER_hbond[kinase])
                np.hstack((RE_graph, RE[:, 0]))
                np.hstack((KE_graph, KE[:, 0]))

    # y_max = max(KE)
    # x_max = max(RE)
    ax = plt.gca()
    ax.set_xlim(0, 20)
    ax.set_ylim(0, 20)
    plt.hexbin(RE_graph, KE_graph, gridsize=36, cmap='jet')
    plt.xlabel('d(E2195-R2430) ($\AA$)')
    plt.ylabel('d(K2187-E2195) ($\AA$)')
    plt.colorbar()
    plt.title('%s sims - %s' % (mutgroup, project))

    plt.savefig('Hexbin_%s_%s_%s.png' % (kinase, project, mutgroup), dpi=500)
    plt.close()


def stat_analyze(distances, window, cutoff):
    tmax = 0
    ntraj = len(distances)  # number of trajectories

    # Compute the maximum time
    for n in range(ntraj):
        if len(distances[n]) > tmax:
            tmax = len(distances[n])
            print('The maximum time is: %d' % len(distances[n]))
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
    cutoff_dist = 5
    list_of_runs = [0]
    for run in list_of_runs:
        RE_total = []
        KE_total = []
        trajectories = dataset.MDTrajDataset(
            "/cbio/jclab/conditions/fah/fah-data/munged4/%s/run%s-clone*.h5" % (project_num, run))
        for traj_in in trajectories:
            [RE_dist, KE_dist] = shukla_coords(traj_in, KER_hbond[kinase_definition])
            RE_total.append(RE_dist[:, 0])
            KE_total.append(KE_dist[:, 0])
        [RE_fraction, RE_stderr] = stat_analyze(RE_total, sliding_window, cutoff_dist)
        [KE_fraction, KE_stderr] = stat_analyze(KE_total, sliding_window, cutoff_dist)
        np.save('./data/%s_%s_RE_fraction_%s.npy' % (project_num, kinase_definition, cutoff_dist), RE_fraction)
        np.save('./data/%s_%s_RE_stderr_%s.npy' % (project_num, kinase_definition, cutoff_dist), RE_stderr)
        np.save('./data/%s_%s_KE_fraction_%s.npy' % (project_num, kinase_definition, cutoff_dist), KE_fraction)
        np.save('./data/%s_%s_KE_stderr_%s.npy' % (project_num, kinase_definition, cutoff_dist), KE_stderr)

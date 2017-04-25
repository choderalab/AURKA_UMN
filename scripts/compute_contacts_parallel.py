# Compute Contacts Script
# Steven Albanese, Chodera Lab

# Global Imports

import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
from mpi4py import MPI
from msmbuilder import dataset
import argparse
import os


def analyze(trajectories):
    """ Main Analysis function

    :param trajectory: trajectory contact map will be calculated from
    :return: a symmetrical square numpy.array of distances and residue pairs

    """

    min_frame = 400
    distances = list()
    for i,traj in enumerate(trajectories):
        # Various variables to be used later on in script
        end_frame = len(traj)

        short_traj = traj.slice(range(min_frame, end_frame), copy=False)

        # Calculate the contact map
        [distance_trajectory, contact_list] = md.compute_contacts(short_traj, contacts='all', scheme='ca')
        distances.append(distance_trajectory)
        distances = np.asarray(distances)

    return distances


################################
#        Data Analysis         #
################################

if __name__ == "__main__":

    n_runs = 1
    project_basepath = '/cbio/jclab/home/albaness/trajectories/AURKA'
    conditions = ['AURKA_nophos_notpx2/11418', 'AURKA_phos_notpx2/11429']
    save_file_label = ['NoPhos', 'Phos']
    for i,condition in enumerate(conditions):
        trajectories = dataset.MDTrajDataset(os.path.join(project_basepath, '%s/run0-clone*.h5' % condition))
        data = analyze(trajectories)
        np.save("../data/e-fret/%s.npy" % save_file_label[i], data)


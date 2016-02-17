import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
from msmbuilder import dataset
from itertools import chain
import sys
import math
from matplotlib.pyplot import cm
import seaborn as sns

sns.set_style("whitegrid")
sns.set_context("poster")



def stat_analyze(distances, window, cutoff):
    tmax = 0
    ntraj = len(distances)  # number of trajectories

    # Compute the maximum time
    for n in range(ntraj):
        if len(distances[n]) > tmax:
            tmax = len(distances[n])
            print 'The value you are looking for is: %d' % len(distances[n])
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


project = '11411'

sliding_window = 40
cutoff_dist = 3

for index in range(5):
    SB_total = []
    trajectories = dataset.MDTrajDataset("/cbio/jclab/projects/fah/fah-data/munged2/all-atoms/%s/run%d-clone*.h5" % (project, index))
    for i,traj in enumerate(trajectories):
        if i == 0:
            for residue in traj.topology.residues:
                if str(residue) == 'GLU181':
                    e181 = residue
                if str(residue) == 'GLN185':
                    q185 = residue
        res185 = traj.topology.residue(q185.index)

        distances, residue_pairs = md.compute_contacts(traj, contacts=[[e181.index,q185.index]])
        SB_total.append(distances[:,0])

        res185atoms = [atom.index for atom in res185.atoms]
        haystack = traj.top.select("water")
        neighbors = md.compute_neighbors(traj, 0.6, res185atoms, haystack_indices=haystack)

        neighbor_set = set(chain.from_iterable(neighbors))
        neighbor_set = list(neighbor_set)
        # using wernet_nilsson because the output is hbonds per frame
        hbonds0 = md.wernet_nilsson(traj, exclude_water=False, proposed_donor_indices=res185atoms, proposed_acceptor_indices=neighbor_set)
        hbonds1 = md.wernet_nilsson(traj, exclude_water=False, proposed_donor_indices=neighbor_set, proposed_acceptor_indices=res185atoms)
        hbonds = list()
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

        # this is only for funs
        label = lambda hbond : '%s -- %s' % (traj.topology.atom(hbond[0]), traj.topology.atom(hbond[2]))
        for index, frame in enumerate(hbonds):
            if len(frame) > 0:
                break
        for hbond in hbonds[index]:
            print label(hbond)

        # make plots of all data past t = 250ns
            # quantify how much P(salt bridge) and (1-P)
                # P = Ntraj_frames_bound

        # look at trajectories : how long do waters stay in place
            # may need to account for the waters exchanging
        
    SB_fraction, SB_stderr = stat_analyze(SB_total, sliding_window, cutoff_dist)
    np.save('./data/%s_%s_SB_fraction_3.npy' % (project, index), SB_fraction)
    np.save('./data/%s_%s_SB_stderr_3.npy' % (project, index), SB_stderr)
    break

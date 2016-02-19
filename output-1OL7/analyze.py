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

        # make plots of all data past t = 250ns
            # quantify how much P(salt bridge) and (1-P)
                # P = Ntraj_frames_bound

        # look at trajectories : how long do waters stay in place
            # may need to account for the waters exchanging


def bond_analyze(HBonds, window, verbose=False):

    tmax = 0
    ntraj = len(HBonds)  # number of trajectories
    if verbose:
        print("\nNumber of trajectories in this array of bonds: %s " % ntraj)
        print("Type of array of bonds: %s" % type(HBonds))

    # Compute the maximum time
    for n in range(ntraj):
        if len(HBonds[n]) > tmax:
            tmax = len(HBonds[n])
    if verbose:
        print("Maximum trajectory length: %s frames" % tmax)
    # Compute the contact fraction
    contact_fraction = np.zeros([tmax - window], np.float64)
    contact_fraction_stderr = np.zeros([tmax - window], np.float64)
    for frame in range(0, tmax - window): # frame is an index, not an object
        # Count the number of trajectoties that are 't+sliding window' in length
        ntraj_t = 0
        for n in range(ntraj):
            if len(HBonds[n]) >= frame + window:
                ntraj_t += 1
        contact_fraction_n = np.zeros(ntraj_t)
        frame_id = 0
        printed_type = False
        for traj in range(ntraj): # traj is an index, not trajectory
            hbonds_traj = HBonds[traj]
            if len(hbonds_traj) >= frame + window:
                num_hbonds_in_frame_in_traj = hbonds_traj[frame].shape[0]
                if verbose:
                    print("Number of hbonds in this frame of this traj: %s" % num_hbonds_in_frame_in_traj)
                if verbose and num_hbonds_in_frame_in_traj>0 and not printed_type:
                    print("Type of entry representing an hbond in traj of bonds: %s" % type(hbonds_traj[frame][0]))
                    printed_type = True
                local_frame_ids = range(frame,frame+window)
                shapes = list()
                for local_id in local_frame_ids:
                    shapes.append(hbonds_traj[local_id].shape[0])
                if verbose and frame < 5 and traj==0:
                    print("These are the number of rows that the thing is coming up with:")
                    print(shapes)
                    print("These are the frames it is checking:")
                    print(local_frame_ids)
                np_shapes = np.asarray(shapes, dtype=int)
                if verbose:
                    print("The ndarray version of shapes:")
                    print(np_shapes)
                contact_fraction_n[frame_id] = np_shapes.mean()
                if verbose:
                    print("Contact fraction frame %s for trajectory %s" % (frame_id,traj)) # why tf does index exist instead of using n
                    print(contact_fraction_n[frame_id])
                frame_id += 1
        contact_fraction[frame] = contact_fraction_n.mean()
        contact_fraction_stderr[frame] = contact_fraction_n.std() / np.sqrt(ntraj_t)

    return contact_fraction, contact_fraction_stderr


def distance_analyze(distances, window, verbose=False):
    tmax = 0
    ntraj = len(distances)  # number of trajectories

    # Compute the maximum time
    for n in range(ntraj):
        if len(distances[n]) > tmax:
            tmax = len(distances[n])
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
        frame_id = 0
        for n in range(ntraj):
            if verbose and t==0 and n==0:
                print("Entry 0 in trajectory %s of distances:" % n)
                print(distances[n][0])
                print("Type of trajectory of distances:")
                print(type(distances[n]))
                print("Type of an entry within trajectory of distances:")
                print(type(distances[n][0]))
            if len(distances[n]) >= t + window:
                contact_fraction_n[frame_id] = (distances[n][t:(t + window)]).mean()
                frame_id += 1
        contact_fraction[t] = contact_fraction_n.mean()
        contact_fraction_stderr[t] = contact_fraction_n.std() / np.sqrt(ntraj_t)

    return contact_fraction, contact_fraction_stderr

project = '11411'

verbose = True
sliding_window = 10
# or like this isn't quite what i want, i need a mean distance not a mean number of distances (?)

fig = plt.figure()

for index in range(5):
    if verbose:
        print("Loading RUN%s..." % index)
    SB_total = []
    HB_total = []
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
                num_hbonds0 = bondlist.shape[0]
                num_hbonds1 = hbonds1[frame].shape[0]
                num_hbonds = hbonds[frame].shape[0]
                assert num_hbonds == num_hbonds0 + num_hbonds1
            except Exception as e:
                print('hbonds0')
                print(bondlist)
                print(bondlist.shape)
                print('hbonds1')
                print(hbonds1[frame])
                print(hbonds1[frame].shape)
                raise(e)

        HB_total.append(hbonds)

        # this is for funs
        label = lambda hbond : '%s -- %s' % (traj.topology.atom(hbond[0]), traj.topology.atom(hbond[2]))
        for frame_id, frame in enumerate(hbonds):
            if len(frame) > 0:
                break
        for hbond in hbonds[frame_id]:
            print label(hbond)

        # end loop for each trajectory in a given run

    SB_fraction, SB_stderr = distance_analyze(SB_total, sliding_window)
    HB_fraction, HB_stderr = bond_analyze(HB_total, sliding_window)
    if verbose:
        print("Saving Salt Bridge fraction...")
    np.save('./data/%s_%s_SB_fraction_3.npy' % (project, index), SB_fraction)
    if verbose:
        print("Saving Salt Bridge standard error...")
    np.save('./data/%s_%s_SB_stderr_3.npy' % (project, index), SB_stderr)
    if verbose:
        print("Saving Hydrogen Bonding fraction...")
    np.save('./data/%s_%s_HB_fraction_3.npy' % (project, index), HB_fraction)
    if verbose:
        print("Saving Hydrogen Bonding standard error...")
    np.save('./data/%s_%s_HB_stderr_3.npy' % (project, index), HB_stderr)
    if verbose:
        print("Saving total Hydrogen Bonding data...")
    np.save('./data/%s_%s_HBonds.npy' % (project, index), HB_total)
    if verbose:
        print("Saving total Salt Bridge data...")
    np.save('./data/%s_%s_SB_total.npy' % (project, index), SB_total)
    if verbose:
        print("Adding RUN%s to Salt Bridge plot..." % index)
    plt.fill_between(range(len(SB_fraction)),SB_fraction-SB_stderr, SB_fraction+SB_stderr)
    plt.plot(SB_fraction)
    # end loop for each run


plt.legend(['RUN0','RUN1','RUN2','RUN3','RUN4'])
if verbose:
    print("Saving Plot...")
plt.savefig("salt-bridge-distances.png",dpi=300)
plt.close(fig)
if verbose:
    print("Complete!")

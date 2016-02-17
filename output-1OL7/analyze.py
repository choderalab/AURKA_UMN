import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
from msmbuilder import dataset
from itertools import chain

project = '11411'
fig = plt.figure()
for index in range(5):
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
        distances_to_plot = [d[0] for d in distances]
        plt.plot(distances)

        res185atoms = [atom.index for atom in res185.atoms]
        haystack = traj.top.select("water")
        neighbors = md.compute_neighbors(traj, 0.6, res185atoms, haystack_indices=haystack)

        neighbor_set = set(chain.from_iterable(neighbors))
        neighbor_set = list(neighbor_set)
        # using wernet_nilsson because the output is hbonds per frame
        hbonds0 = md.wernet_nilsson(traj, exclude_water=False, proposed_donor_indices=res185atoms, proposed_acceptor_indices=neighbor_set)
        hbonds1 = md.wernet_nilsson(traj, exclude_water=False, proposed_donor_indices=neighbor_set, proposed_acceptor_indices=res185atoms)

        # this is wrong now -- want to match up corresponding frames!
        hbonds = np.array(hbonds0 + hbonds1)

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

        break
    break
plt.legend(['RUN0','RUN1','RUN2','RUN3','RUN4'])
plt.savefig("salt-bridge-distances.png",dpi=300)
plt.close(fig)

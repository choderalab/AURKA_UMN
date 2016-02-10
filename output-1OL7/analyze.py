import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
from msmbuilder import dataset

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
        min_index = min([atom.index for atom in e181.atoms])
        max_index = min([atom.index for atom in traj.topology.residue(q185.index+1).atoms])
        atom_indices = [i for i in range(min_index,max_index)]
        # how to include relevant waters?
        salt_bridge = traj.atom_slice(atom_indices)

        distances, residue_pairs = md.compute_contacts(salt_bridge, contacts=[[0,4]])
#        distances, residue_pairs = md.compute_contacts(traj, contacts=[[e181.index,q185.index]])
        distances_to_plot = [d[0] for d in distances]
        print(len(distances))
        plt.plot(distances)


        possible_hbonds = traj.atom_slice(atom_indices)
        hbonds = md.baker_hubbard(possible_hbonds, exclude_water=False)
        label = lambda hbond : '%s -- %s' % (possible_hbonds.topology.atom(hbond[0]), possible_hbonds.topology.atom(hbond[2]))
        for hbond in hbonds:
            print label(hbond)

        # md.compute_contacts(traj, contacts=[[0-idx,0-idx]])
            # can plot with time, all on top of each other
            # plot distance over time
            # store in hdf5? make plot right away? each will be a different length (list of the numpy arrays or something and then store that or pickle it)(why)
        # hydrogen bonding: how to tell it one residue instead of looking at the whole thing
            # ask steve how to do the slice thing
        # how does hydrogen bonding change with and without TPX2
        # make plots of all data past t = 250ns
            # quantify how much P(salt bridge) and (1-P)
                # P = Ntraj_frames_bound

        # look at trajectories : how long do waters stay in place
            # may need to account for the waters exchanging
        # ask nick for which residues (why.)
        # make a cut-out traj of residues plus closest waters in each frame ** might need fanciness
        # john looks at code for specifying residues in hbonding 

#        hbonds = md.baker_hubbard(traj, exclude_water=False)
#        label = lambda hbond : '%s -- %s' % (traj.topology.atom(hbond[0]), traj.topology.atom(hbond[2]))
#        for hbond in hbonds:
#            print label(hbond)
        break
plt.legend(['RUN0','RUN1','RUN2','RUN3','RUN4'])
plt.savefig("salt-bridge-distances.png",dpi=300)
plt.close(fig)

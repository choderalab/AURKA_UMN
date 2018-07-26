
import sys
import mdtraj as md
import numpy as np
from msmbuilder import dataset




# Define kinase.
condition = sys.argv[1]


def get_atom_index(traj, selection):
    indices = traj.top.select(selection)
    if len(indices) != 1:
        msg = 'Selection "%s" did not match a unique atom.' % (selection)
        raise Exception(msg)

    return indices[0]


def compute_torsion(traj, *args):
    """
    Compute the specified torsion.
    """
    indices = [get_atom_index(traj, selection) for selection in args]
    min_frame = 400
    end_frame = len(traj)
    frames_to_slice = list(range(min_frame, end_frame))
    short_traj = traj.slice(frames_to_slice, copy=False)
    # Compute torsion in degrees
    torsions = md.compute_dihedrals(short_traj, [indices]).squeeze() * (180.0 / np.pi)

    return torsions


if __name__ == "__main__":
    trajectories = dataset.MDTrajDataset(
        '/cbio/jclab/home/albaness/trajectories2/AURKA/CHARMM/%s/*.h5' % condition)
    torsion1_list = []
    torsion2_list = []
    for traj_in in trajectories:
        if len(traj_in) > 400:
            torsion1 = compute_torsion(traj_in, *['(resSeq %d and name CA and chainid 0)' % resSeq for resSeq in (160, 161, 162, 163)])
            torsion1_list.extend(torsion1)
            torsion2 = compute_torsion(traj_in,
                                       *['(resSeq %d and name CA and chainid 0)' % resSeq for resSeq in (161, 162, 163, 164)])
            torsion2_list.extend(torsion2)
        else:
            pass
    np.save('../data/dihedral/dihedral_%s-%s-%s.npy' % (condition, 282, 285), torsion1_list)
    np.save('../data/dihedral/dihedral_%s-%s-%s.npy' % (condition, 283, 286),
            torsion2_list)
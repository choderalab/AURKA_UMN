
import mdtraj as md

base_path = '/cbio/jclab/projects/fah/fah-data/munged3/no-solvent/11418'

trajectories = ['run0-clone8.h5',
                'run0-clone27.h5',
                'run1-clone5.h5',
                'run1-clone25.h5',
                'run2-clone6.h5',
                'run3-clone11.h5',
                'run3-clone19.h5',
                'run4-clone18.h5',
                'run4-clone49.h5']

for traj in trajectories:
    name = traj.split('.')[0]
    t = md.load('%s/%s'%(base_path,traj), stride=10)
    t[0].save_pdb('%s.pdb'%name)
    t.save_xtc('%s.xtc'%name)

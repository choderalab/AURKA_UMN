import numpy as np

out_line = -0.5
for group in ['+', '-', '+muts', '-muts']:
    filename = "./results/AURKA"+group+"/dihedral_AURKA.npy"
    #filename = "./results/AURKA"+group+"/dihedral_by_something.npy"
    dihedrals = np.load(filename)
    print(group)
    print('Shape of dihedrals:')
    print(dihedrals.shape)
    dihedrals = dihedrals.flatten()
    outhedrals = dihedrals[dihedrals <= out_line]
    print('%s out of %s frames' % (len(outhedrals), len(dihedrals)))
    print(outhedrals)
    filename = "./results/AURKA"+group+"/dihedral_by_something.npy"
    dihedrals = np.load(filename)
    for clone, traj in enumerate(dihedrals):
        if any([angle <= out_line for angle in traj]):
            print(clone)
    del(outhedrals)
    del(dihedrals)



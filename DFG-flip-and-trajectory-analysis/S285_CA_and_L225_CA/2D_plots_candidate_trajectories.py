#This script plots our trajectory on the coordinates that allowed us to select in the first place (Steven's plot)

# Imports

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import seaborn
import mdtraj as md
import re
from glob import glob
import os

#Define our trajectories
base_path = '/cbio/jclab/projects/fah/fah-data/munged3/no-solvent/11418'

files = ['run0-clone8.h5',
         'run0-clone27.h5',
         'run1-clone5.h5',
         'run1-clone25.h5',
         'run2-clone6.h5',
         'run3-clone11.h5',
         'run3-clone19.h5',
         'run4-clone18.h5',
         'run4-clone49.h5']


#Define our coordinates
protein = 'AURA'

SL_distance = dict()

SL_distance['AURA'] = ['resid 161 and resname SER and name CA','resid 102 and resname LEU and name CA']

TR_distance = dict()

TR_distance['AURA'] = ['resid 164 and resname THR and name CA','resid 132 and resname ARG and name CZ']

for my_file in files:
     filename = '%s/%s'%(base_path, my_file)
     clone_name = my_file.split('.')[0]

     traj = md.load(filename)
     topology = traj.topology

     SL_atom_1 = topology.select(SL_distance['AURA'][0])
     SL_atom_2 = topology.select(SL_distance['AURA'][1])
     SL_atoms = [SL_atom_1[0], SL_atom_2[0]]

     SL = md.compute_distances(traj,[SL_atoms])

     TR_atom_1 = topology.select(TR_distance['AURA'][0])
     TR_atom_2 = topology.select(TR_distance['AURA'][1])
     TR_atoms = [TR_atom_1[0], TR_atom_2[0]]

     TR = md.compute_distances(traj,[TR_atoms])

     plt.figure(figsize=(5,4))

     plt.plot(TR*10,SL*10,'.', alpha=0.5, ms='3', label='all frames')
     plt.axhline(y=38,color='0.7',linestyle='--',label='cut-off')

     plt.plot(TR[0]*10,SL[0]*10,'.', alpha=0.5, ms='8', color = 'green', label='first frame')
     plt.plot(TR[-1]*10,SL[-1]*10,'.', alpha=0.5, ms='8', color = 'red', label='last frame')

     plt.ylabel('S285 CA - L225 CA distance ($\AA$)')
     plt.xlabel('T288 CA - R255 CZ distance ($\AA$)')

     #plt.ylim((28,46))
     #plt.xlim((7,20))

     plt.title('-Phos -Tpx2: Proj 11418:%s' %clone_name)
     plt.legend()

     plt.savefig('%s-2D.png'%clone_name,bbox_inches='tight',dpi=700)
     plt.close()

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import seaborn
import mdtraj as md
import re
from glob import glob
import os

# Pick which protein your analyzing
protein = 'AURA'

projects = dict()

projects['AURA'] = 11418

# These distances are meant to describe the DFG flip in AURA
# The first two were suggested by Nick Levinson
# The third had been found previously useful for Src and Abl structures
# And uses an alanine in the alphaE helix: ELAN_A_LSY in AURK_A

distances_dict = {}

distances_dict['F275 CZ - I193 CG2'] = ['resid 152 and resname PHE and name CZ','resid 70 and resname ILE and name CG2']
distances_dict['W277 CG2 - I193 CG2'] = ['resid 154 and resname TRP and name CE2','resid 70 and resname ILE and name CG2']
distances_dict['F275 CZ - A243 CA'] = ['resid 152 and resname PHE and name CZ','resid 120 and resname ALA and name CA']

# These are other distances of interest in AURA

distances_dict['S285 CA - L225 CA'] = ['resid 161 and resname SER and name CA','resid 102 and resname LEU and name CA']
distances_dict['T288 CA - R255 CZ'] = ['resid 164 and resname THR and name CA','resid 132 and resname ARG and name CZ']

# We want all the results files to go to different destinations depending on what distance they're for

destination_dict = {}

destination_dict['F275 CZ - I193 CG2'] = 'F275_CZ_and_I193_CG2'
destination_dict['W277 CG2 - I193 CG2'] = 'W277_CG2_and_I193_CG2'
destination_dict['F275 CZ - A243 CA'] = 'F275_CZ_and_A243_CA'
destination_dict['S285 CA - L225 CA'] = 'S285_CA_and_L225_CA'
destination_dict['T288 CA - R255 CZ'] = 'T288_CA_and_R255_CZ'

#Define plotting things
cutoff_dict = {}

cutoff_dict['F275 CZ - I193 CG2'] = 0.6
cutoff_dict['W277 CG2 - I193 CG2'] = 1.3
cutoff_dict['F275 CZ - A243 CA'] = 1.3
cutoff_dict['S285 CA - L225 CA'] = 3.8 

ylim_dict = {}

ylim_dict['F275 CZ - I193 CG2'] = (0.0,1.2)
ylim_dict['W277 CG2 - I193 CG2'] = (0.6,2.2)
ylim_dict['F275 CZ - A243 CA'] = (0.6,2.2)
ylim_dict['S285 CA - L225 CA'] = (3.2,4.6)

# Here we manually define the number of runs (same for all distances)
num_runs = 5

def DFG_distance_bytraj(files,def_DFG):

    distance_with_labels = []

    for filename in files:

        distance = []
        run = re.search('run([^-]+)',filename).group(1)
        clone = re.search('clone([^.]+)',filename).group(1)

        traj = md.load(filename)

        topology = traj.topology

        def_DFG_atom_1 = topology.select(def_DFG[0])
        def_DFG_atom_2 = topology.select(def_DFG[1])
   
        def_DFG_atoms = [def_DFG_atom_1[0], def_DFG_atom_2[0]]
        #print 'These correspond to atom numbers %s.' %def_DFG_atoms

        #print 'Atom distances computed between %s and %s' %(topology.atom(def_DFG_atom_1),topology.atom(def_DFG_atom_2))  

        distance.append(md.compute_distances(traj,[def_DFG_atoms]))

        short_filename = filename.split('/')[-1]

        distance_with_labels.append([[short_filename,run,clone],distance])

    return distance_with_labels

for my_distance in cutoff_dict:

    print 'Using %s' %my_distance

    destination = destination_dict[my_distance]
    cutoff = cutoff_dict[my_distance]
    ylim = ylim_dict[my_distance]

    try:
        os.stat(destination)
    except:
        os.mkdir(destination)

    for i in range(num_runs):

        print 'working on run%s ...'%i

        files = glob('/cbio/jclab/projects/fah/fah-data/munged3/no-solvent/%s/run%s-*.h5'%(projects[protein],i))
        distances = DFG_distance_bytraj(files,distances_dict[my_distance])
        np.save('%s/distance_%s_bytraj_run%s.npy'%(destination,protein,i),distances)

        colors = matplotlib.cm.summer(np.linspace(0, 1, len(distances)))

        plt.figure(figsize=(18,14))

        frame_add = [0]
        flipfinder = []
        for j in range(len(distances)):
            frame_add = [frame_add[-1] + a for a in range(len(distances[j][1][0]))]
            plt.plot(frame_add,distances[j][1][0],'o',color=colors[j],label='%s'%distances[j][0][0])
            if np.mean(distances[0][1][0:100]) >= cutoff:
                if any(distance<=(cutoff-0.2) for distance in distances[j][1][0]):
                    flipfinder.append(distances[j][0][0])
                    print 'FlipFinder says: you might have a flip in %s!' %distances[j][0][0]
                    plt.text(frame_add[0],min(distances[j][1][0]),'%s'%distances[j][0][0])
            else:
                if any((cutoff+0.1)<=distance for distance in distances[j][1][0]):
                    flipfinder.append(distances[j][0][0])
                    print 'FlipFinder says: you might have a flip in %s!' %distances[j][0][0]
                    plt.text(frame_add[0],max(distances[j][1][0]),'%s'%distances[j][0][0])

        plt.axhline(y=cutoff,color='red',label='cut-off')
        plt.title('Concatenated Trajectories: %s_%s_run%s'%(protein,projects[protein],i),fontsize=26)
        #if np.mean(distances[0][1][0:100]) >= cutoff:
        #    plt.title('Concatenated Trajectories %s_%s_run%s: DFG-out start'%(protein,projects[protein],i),fontsize=26)
        #else:
        #    plt.title('Concatenated Trajectories %s_%s_run%s: DFG-in start'%(protein,projects[protein],i),fontsize=26)
        #plt.text(1000,0.7,'DFG-in',fontsize=26,fontweight='bold')
        #plt.text(1000,2.1,'DFG-out',fontsize=26,fontweight='bold')

        plt.ylim(ylim)
        plt.ylabel('%s distance (nm)'%my_distance,fontsize=20)
        plt.xlabel('Frame',fontsize=20)
        plt.tick_params(axis='both',which='major',labelsize=20)
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.savefig('%s/%s_%s_DFG_by_traj_run%s.png'%(destination,protein,projects[protein],i),bbox_inches='tight')
        plt.close()

        for traj in flipfinder:
            base_name = os.path.splitext(traj)[0] 
            flip_file = ['/cbio/jclab/projects/fah/fah-data/munged3/no-solvent/%s/%s'%(projects[protein],traj)]
            distances = DFG_distance_bytraj(flip_file,distances_dict[my_distance])

            plt.clf()
            plt.figure(figsize=(8,6))
            plt.plot(distances[0][1][0],'o',label='%s'%distances[0][0][0])
            plt.axhline(y=cutoff,color='red',label='cut-off')

            plt.ylim(ylim)
            plt.title('%s_%s:%s'%(protein,projects[protein],base_name))
            plt.ylabel('%s distance (nm)'%my_distance,fontsize=20)
            plt.xlabel('Frame',fontsize=20)
            plt.legend(loc=0)
            plt.savefig('%s/%s_%s_DFG_by_clone-%s.png'%(destination,protein,projects[protein],base_name))
            plt.close()



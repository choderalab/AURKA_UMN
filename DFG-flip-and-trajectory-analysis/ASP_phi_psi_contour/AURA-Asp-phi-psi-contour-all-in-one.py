import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import mdtraj as md
import re
from glob import glob
import os

from msmbuilder import dataset
import re

plt.rc('font',family='serif')

sns.set_style("white")
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Pick which protein your analyzing

protein = 'AURA'

state_order = {}
state_order['AURKA_phos_tpx2'] = ['0','+Phos +TPX2']
state_order['AURKA_nophos_tpx2'] = ['1','-Phos +TPX2']
state_order['AURKA_phos_notpx2'] = ['2','+Phos -TPX2']
state_order['AURKA_nophos_notpx2'] = ['3','-Phos -TPX2']

ASP_phi = dict()
ASP_phi['AURA'] = ['resid 151 and resname ASP and name C','resid 151 and resname ASP and name CA','resid 151 and resname ASP and name N','resid 150 and resname ALA and name C']

ASP_psi = {}
ASP_psi['AURA'] =  ['resid 152 and resname PHE and name N','resid 151 and resname ASP and name C','resid 151 and resname ASP and name CA','resid 151 and resname ASP and name N']

def phi_psi_byrun(files,phi_def,psi_def):

    phi = []
    psi = []

    phi_combinetrajs = []
    psi_combinetrajs = []

    glob_files = glob(files)

    for file in glob_files:
    
        traj = md.load(file)

        topology = traj.topology

        phi_atom_1 = topology.select(phi_def[0])
        phi_atom_2 = topology.select(phi_def[1])
        phi_atom_3 = topology.select(phi_def[2])
        phi_atom_4 = topology.select(phi_def[3])

        #print 'The dihedral computed is between %s, %s, %s, and %s' %(topology.atom(phi_atom_1),topology.atom(phi_atom_2),topology.atom(phi_atom_3),topology.atom(phi_atom_4))       

        phi_atoms = [phi_atom_1[0], phi_atom_2[0],phi_atom_3[0], phi_atom_4[0]]
        #print 'These correspond to atom numbers %s.' %phi_atoms

        phi_combinetrajs.append(md.compute_dihedrals(traj,[phi_atoms]) * (180.0 / np.pi)) # this will be in degrees
            
        psi_atom_1 = topology.select(psi_def[0])
        psi_atom_2 = topology.select(psi_def[1])
        psi_atom_3 = topology.select(psi_def[2])
        psi_atom_4 = topology.select(psi_def[3])

        #print 'The dihedral computed is between %s, %s, %s, and %s' %(topology.atom(psi_atom_1),topology.atom(psi_atom_2),topology.atom(psi_atom_3),topology.atom(psi_atom_4))       

        psi_atoms = [psi_atom_1[0], psi_atom_2[0],psi_atom_3[0], psi_atom_4[0]]
        #print 'These correspond to atom numbers %s.' %psi_atoms

        psi_combinetrajs.append(md.compute_dihedrals(traj,[psi_atoms])* (180.0 / np.pi)) # this will be in degrees

    # flatten list of arrays
    phi_combinetrajs = np.asarray([val for sublist in phi_combinetrajs for val in sublist])
    psi_combinetrajs = np.asarray([val for sublist in psi_combinetrajs for val in sublist])

    phi.append(phi_combinetrajs)
    phi_combinetrajs = []

    psi.append(psi_combinetrajs)
    psi_combinetrajs = []

    return [phi, psi]

colors = ["amber","dark lavender","cerulean", "pale red"]
colorp = sns.xkcd_palette(colors)

plt.figure(figsize=(20,15), dpi=300)

for state in state_order:

    i = int(state_order[state][0])

    files = '/cbio/jclab/home/albaness/trajectories2/AURKA/%s/*/*.h5'%state

    [phi_separate,psi_separate] = phi_psi_byrun(files,ASP_phi[protein],ASP_psi[protein])

    #save rmsd and difference data, these are separated by run
    np.save('phi_separate_%s.npy'%state,phi_separate)
    np.save('psi_separate_%s.npy'%state,psi_separate)

    # here we combine all the runs
    phi = np.asarray([val for sublist in phi_separate for val in sublist])
    psi = np.asarray([val for sublist in psi_separate for val in sublist])

    plt.subplot(2,2,i+1)

    cmap = sns.light_palette(colorp[i], as_cmap=True)

    sns.kdeplot(phi.flatten(), psi.flatten(), shaded=True, cmap=cmap)

    plt.xlabel('ASP 274 $\phi$', fontsize = 18)
    plt.ylabel('ASP 274 $\psi$', fontsize = 18) 
    plt.yticks([-180,-90,0,90,180],fontsize = 16)
    plt.xticks([-180,-90,0,90,180],fontsize = 16)
    plt.ylim((-180,180))
    plt.xlim((-180,180))

    plt.title('%s' % (state_order[state][1]),fontsize = 20 )

plt.savefig('phi_psi_contour_all-in-one.pdf',bbox_inches='tight',dpi=300)


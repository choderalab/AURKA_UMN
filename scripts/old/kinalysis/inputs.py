
# This dictionary defines which condition you are running by your protein name.
# This will have to be hardcorded.

projects = dict()

projects['AURKA'] = 11414

# This dictionary defines which dihedrals define the DFG flip for each protein. Hopefully this will eventually
# be defined programmatically rather than by hand.

DFG = dict()

import mdtraj as md
aurora = md.load("../../../AURKA_UMN/trajectories/%s_RUN0.pdb" % projects['AURKA'])
for residue in aurora.topology.residues:
    if str(residue) == 'ALA273':
        a273 = residue
    if str(residue) == 'ASP274':
        d274 = residue
    if str(residue) == 'LYS162':
        k162 = residue
    if str(residue) == 'GLU181':
        e181 = residue
    if str(residue) == 'GLU299':
        e299 = residue
    if str(residue) == 'ARG285':
        r285 = residue
for atom in a273.atoms:
    if atom.name == 'CB':
        alacb = atom.index
    if atom.name == 'CA':
        alaca = atom.index
for atom in d274.atoms:
    if atom.name == 'CA':
        aspca = atom.index
    if atom.name == 'CG':
        aspcg = atom.index

DFG['AURKA'] = [alacb, alaca, aspca, aspcg]

print 'Projects:'
print projects
print 'DFG dihedral definitions:'
print DFG

# These two dictionary defines our Shukla coordinates

# Define hydrogen bond coordinates (0-indexed)
KER_hbond = { 'AURKA' : [[k162.index,e181.index],[e181.index,r285.index]]}

# Define Activation loop (resid)
Aloop_def = { 'AURKA': [d274.index,e299.index]}

print 'KER hydrogen bond definitions:'
print KER_hbond
print 'Activation loop definitions :'
print Aloop_def

# Here we store these dictionaries as JSON files.

import json
with open('conditions.json', 'w') as fp:
    json.dump(projects, fp)
with open('DFG.json', 'w') as fp:
    json.dump(DFG, fp)
with open('KER_hbond.json', 'w') as fp:
    json.dump(KER_hbond, fp)
with open('Aloop_def.json', 'w') as fp:
    json.dump(Aloop_def, fp)

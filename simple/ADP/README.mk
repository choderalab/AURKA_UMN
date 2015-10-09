Within this directory:

minimize.py: loads protein, adds ADP, attempts to minimize through 2 minimization steps 
(constraints=None and constraints=app.HBonds), runs MD steps, writes out to a pdb file
(ADP-and-protein-minimized.py) if successful

ADP-alone.py: loads ADP mol2, attempts to minimize and run MD steps, writes out to a pdb
file (ADP-minimizing.pdb) after each individual step, will print out forces upon failure

1OL5_pdb_setup.py: runs pdbfixer and generates 1OL5-pdbfixer.pdb, which is the input to minimize.py

ADP5.mol2: generated from openmoltools, extracted ADP molecule from 1OL5 pdb

ADP.xml: generated from processAmberForceField.py, manually corrected (atom types, loop closing bonds)

with-ADP.pdb: pdbfile written out during minimize.py after ADP is added via modeller, to 
make sure there are no steric clashes when the molecule is added



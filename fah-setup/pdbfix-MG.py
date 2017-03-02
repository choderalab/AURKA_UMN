import pdbfixer
from simtk.openmm.app import PDBFile

pdbfile = "1OL5-WT-pdbfixer.pdb"

fixer = pdbfixer.PDBFixer(pdbfile)
fixer.findNonstandardResidues()
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()

PDBFile.writeFile(fixer.topology, fixer.positions, open('1OL5-modifiedMG-WT.pdb', 'w'), keepIds=True)


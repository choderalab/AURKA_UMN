from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile

pdbid = '1OL5'

fixer = PDBFixer(pdbid=pdbid)

fixer.missingResidues = {}
fixer.findMissingAtoms()
fixer.addMissingAtoms()

fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()

#fixer.removeHeterogens(keepWater=False)

fixer.addMissingHydrogens(pH=7.4)

PDBFile.writeFile(fixer.topology, fixer.positions, open(pdbid+'-ions-pdbfixer.pdb', 'w'), keepIds=True)


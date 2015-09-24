from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile

pdbid = '1OL7'

fixer = PDBFixer(pdbid=pdbid)

fixer.missingResidues = {}
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.removeHeterogens()

fixer.removeChains(chainIndices=[1])

fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()
fixer.addMissingHydrogens(pH=7.4)

PDBFile.writeFile(fixer.topology, fixer.positions, open(pdbid+'-AUR-only.pdb', 'w'), keepIds=True)




from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile

RESIDUES_TO_REMOVE = ['SO4']

def delete_this_line(line):
    for resi in RESIDUES_TO_REMOVE:
        if line.find(resi)!=-1:
            return True
    return False

for pdbid in ['1OL5','1OL7']:

    fixer = PDBFixer(pdbid=pdbid)

    fixer.findMissingResidues()
    for key in fixer.missingResidues.keys():
        if key[0]==0:
            fixer.missingResidues.pop(key)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.addMissingHydrogens(pH=7.4)

    PDBFile.writeFile(fixer.topology, fixer.positions, open(pdbid+'-WT-int-pdbfixer.pdb', 'w'), keepIds=True)

    with open(pdbid+'-WT-int-pdbfixer.pdb', 'r') as fold:
        with open(pdbid+'-WT-pdbfixer.pdb', 'w') as fnew:
            for line in fold.readlines():
                if not delete_this_line(line):
                    fnew.write(line)






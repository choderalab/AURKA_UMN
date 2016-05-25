from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile, Modeller, ForceField

RESIDUES_TO_REMOVE = ['SO4']

def delete_this_line(line):
    for resi in RESIDUES_TO_REMOVE:
        if line.find(resi)!=-1:
            return True
    return False

for TPX2 in [True, False]:

    if TPX2:
        id = '+'
    else:
        id = '-'
    pdbid = '1OL5'
    fixer = PDBFixer(pdbid=pdbid)

    fixer.findMissingResidues()
    for key in fixer.missingResidues.keys():
        if key[0]==0:
            fixer.missingResidues.pop(key)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)

    if not TPX2:
        fixer.removeChains(chainIds=['B'])

    forcefield = ForceField('amber99sbildn.xml','TPO.xml','tip3p.xml','ions.xml','adp.xml')
    modeller = Modeller(fixer.topology, fixer.positions)
    modeller.addHydrogens(forcefield=forcefield,pH=7.4)

    PDBFile.writeFile(modeller.topology, modeller.positions, open(pdbid+'-WT-int-pdbfixer.pdb', 'w'), keepIds=True)
    with open(pdbid+'-WT-int-pdbfixer.pdb', 'r') as fold:
        with open(pdbid+'-WT'+id+'TPX2.pdb', 'w') as fnew:
            for line in fold.readlines():
                if not delete_this_line(line):
                    fnew.write(line)


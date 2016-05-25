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

    for residue in modeller.topology.residues():
        if residue.name == 'TPO':
            tpo1 = residue
            break
    for residue in modeller.topology.residues():
        if residue.index == tpo1.index-1:
            prev = residue
        elif residue.index == tpo1.index+1:
            tpo2 = residue
            break
    for atom in prev.atoms():
        if atom.name == 'C':
            prevC = atom
            break
    for atom in tpo1.atoms():
        if atom.name == 'N':
            tpo1N = atom
        elif atom.name == 'C':
            tpo1C = atom         
    for atom in tpo2.atoms():
        if atom.name == 'N':
            tpo2N = atom
            break
    modeller.topology.addBond(prevC, tpo1N)
    modeller.topology.addBond(tpo1C, tpo2N)

    modeller.addHydrogens(forcefield=forcefield,pH=7.4)

    PDBFile.writeFile(modeller.topology, modeller.positions, open(pdbid+'-WT-int-pdbfixer.pdb', 'w'), keepIds=True)

    with open(pdbid+'-WT-int-pdbfixer.pdb', 'r') as fold:
        with open(pdbid+'-WT'+id+'TPX2.pdb', 'w') as fnew:
            for line in fold.readlines():
                if not delete_this_line(line):
                    fnew.write(line)






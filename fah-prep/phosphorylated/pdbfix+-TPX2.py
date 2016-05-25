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

    tpo_ff = ForceField("TPO.xml")
    modeller = Modeller(fixer.topology, fixer.positions)

    new_atoms = list()
    for residue in modeller.topology.residues():
        if residue.name != 'TPO':
            continue
        template = tpo_ff._templates['T1P']
        standard_atoms = set(atom.name for atom in template.atoms)
        template_atoms = list(template.atoms)
        atom_names = set(atom.name for atom in residue.atoms())
        missing = list()
        for atom in template_atoms:
            if atom.name not in atom_names:
                missing.append(atom)
        print("This better be 7H...")
        print(missing)
        for atom in missing:
            new_atom = modeller.topology.addAtom(atom.name, atom.element, residue)
            new_atoms.append(new_atom)
        new_res_atoms = dict()
        for atom in residue.atoms():
            new_res_atoms[atom.name] = atom
        new_res_bonds = list()
        for bond in modeller.topology._bonds:
            if bond[0].residue == residue and bond[1].residue == residue:
                new_res_bonds.append((bond[0].name, bond[1].name))
        template_bonds = [(template.atoms[bond[0]].name, template.atoms[bond[1]].name) for bond in template.bonds]
        for bond in new_res_bonds:
            if bond not in template_bonds and (bond[1],bond[0]) not in template_bonds:
                bonded_0 = new_res_atoms[bond[0]]
                bonded_1 = new_res_atoms[bond[1]]
                modeller.topology._bonds.remove((bonded_0, bonded_1))
        for bond in template_bonds:
            if bond not in new_res_bonds and (bond[1],bond[0]) not in new_res_bonds:
                new_bonded_0 = new_res_atoms[bond[0]]
                new_bonded_1 = new_res_atoms[bond[1]]
                modeller.topology.addBond(new_bonded_0, new_bonded_1)
    modeller.topology._numAtoms = len(list(modeller.topology.atoms()))










    PDBFile.writeFile(modeller.topology, modeller.positions, open(pdbid+'-WT-int-pdbfixer.pdb', 'w'), keepIds=True)
    with open(pdbid+'-WT-int-pdbfixer.pdb', 'r') as fold:
        with open(pdbid+'-WT'+id+'TPX2.pdb', 'w') as fnew:
            for line in fold.readlines():
                if not delete_this_line(line):
                    fnew.write(line)


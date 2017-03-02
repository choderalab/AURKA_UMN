"""
Set up AurKA mutant simulations
@author John D. Chodera
@date 8 Aug 2014
"""

#
# IMPORTS
#

import os, os.path
import copy
import numpy
import shutil
import tempfile
import numpy as np

import pdbfixer
from simtk import openmm, unit
from simtk.openmm import app
import mdtraj as md

#

TPX2 = False

if TPX2:
    id = '+'
else:
    id = '-'
# Path to put all output in
output_path = "1OL5"+id+"TPX2-TPOs"

# Source PDB
pdbfilename = "1OL5-WT"+id+"TPX2.pdb"

adp_mol2 = "ADP5.mol2"

print "Source PDB filename: %s" % pdbfilename
print "Output directory for mutations: %s" % output_path

#
# PARAMETERS
#

chain_id_to_mutate = 'A' # chain to mutate
pH = 7.4 # pH to model
keep_crystallographic_water = False # keep crystallographic waters?

# Single point mutants
point_mutants = ['Q185C', 'Q185L','Q185M','Q185N','Q185H','C247A','C247L']

# Forcefield
ff_name = 'amber99sbildn'
water_name = 'tip3p'
ion_ff_name = 'ions'
ADP_ff_name = 'adp'
phos_res = 'TPO'

solvate = True # if True, will add water molecules using simtk.openm.app.modeller
padding = 11.0 * unit.angstroms
nonbonded_cutoff = 9.0 * unit.angstroms
nonbonded_method = app.PME
max_minimization_iterations = 5000
temperature = 300.0 * unit.kelvin
pressure = 1.0 * unit.atmospheres
collision_rate = 5.0 / unit.picoseconds
barostat_frequency = 50
timestep = 2.0 * unit.femtoseconds
nsteps = 5000 # number of steps to take for testing
ionicStrength = 300 * unit.millimolar

# Verbosity level
verbose = True

#===============================================================================
# DATA
#===============================================================================

three_letter_code = {
    'A' : 'ALA',
    'C' : 'CYS',
    'D' : 'ASP',
    'E' : 'GLU',
    'F' : 'PHE',
    'G' : 'GLY',
    'H' : 'HIS',
    'I' : 'ILE',
    'K' : 'LYS',
    'L' : 'LEU',
    'M' : 'MET',
    'N' : 'ASN',
    'P' : 'PRO',
    'Q' : 'GLN',
    'R' : 'ARG',
    'S' : 'SER',
    'T' : 'THR',
    'V' : 'VAL',
    'W' : 'TRP',
    'Y' : 'TYR'
}

one_letter_code = dict()
for one_letter in three_letter_code.keys():
    three_letter = three_letter_code[one_letter]
    one_letter_code[three_letter] = one_letter

def decompose_mutation(mutation):
    import re
    match = re.match('(\D)(\d+)(\D)', mutation)
    original_residue_name = three_letter_code[match.group(1)]
    residue_index = int(match.group(2))
    mutated_residue_name = three_letter_code[match.group(3)]
    return (original_residue_name, residue_index, mutated_residue_name)

def generate_pdbfixer_mutation_code(original_residue_name, residue_index, mutated_residue_name):
    return '%s-%d-%s' % (original_residue_name, residue_index, mutated_residue_name)

def write_file(filename, contents):
    with open(filename, 'w') as outfile:
        outfile.write(contents)

#
# Read reference PDB file to create a list of possible alterations.
#

print "Reading reference PDB file..."
fixer = pdbfixer.PDBFixer(filename=pdbfilename)
residues = list()
for chain in fixer.topology.chains():
    for residue in chain.residues():
        if residue.name != 'HOH':
            key = (chain.id, residue.id, residue.name)
            residues.append(key)
print residues
residues = set(residues)

#
# FILENAMES
#

exception_filename = os.path.join(output_path, 'exceptions.out') # to store exceptions
run_index_filename = os.path.join(output_path, 'run-index.txt') # to store index of which mutants are which

#
# Read the list of mutants already set up.
#

existing_mutants = list()
if os.path.exists(run_index_filename):
    infile = open(run_index_filename, 'r')
    lines = infile.readlines()
    for line in lines:
        [run_name, name] = line.strip().split()
        existing_mutants.append(name)
    infile.close()
print "Existing mutants:"
print existing_mutants

#
# TODO: Determine offset to apply to residue numbers
# TODO: We should look at first residue of source PDB file.
#

#residue_offset = 1 - first_residue # offset that must be added to desired mutation index to index into actual PDB file

#
# Generate list of mutants.
#

npoint_mutants = len(point_mutants)

mutant_names = list()
mutant_codes = list()

# Append wild type (no mutation).
mutant_names.append('WT')
mutant_codes.append([])

# Append point mutants.
for mutation in point_mutants:
    (original_residue_name, residue_index, mutated_residue_name) = decompose_mutation(mutation)
    #residue_index += residue_offset
    key = (chain_id_to_mutate, str(residue_index), original_residue_name)
    print key
    if key in residues:
        mutant_names.append(mutation)
        mutant_codes.append([generate_pdbfixer_mutation_code(original_residue_name, residue_index, mutated_residue_name)])

# Append all pairs of point mutants.
#for i in range(npoint_mutants):
#    for j in range(i+1, npoint_mutants):
#        mutation_i = point_mutants[i]
#        mutation_j = point_mutants[j]

#        (original_residue_name_i, residue_index_i, mutated_residue_name_i) = decompose_mutation(mutation_i)
#        (original_residue_name_j, residue_index_j, mutated_residue_name_j) = decompose_mutation(mutation_j)

        #residue_index_i += residue_offset
        #residue_index_j += residue_offset

#        key_i = (chain_id_to_mutate, str(residue_index_i), original_residue_name_i)
#        key_j = (chain_id_to_mutate, str(residue_index_j), original_residue_name_j)

#        if (key_i in residues) and (key_j in residues):
#            mutant_names.append(mutation_i + '+' + mutation_j)
#            mutant_codes.append([
#                    generate_pdbfixer_mutation_code(original_residue_name_i, residue_index_i, mutated_residue_name_i), 
#                    generate_pdbfixer_mutation_code(original_residue_name_j, residue_index_j, mutated_residue_name_j)
#                    ])

print ""
print "Feasible mutants:"
print mutant_names
print ""


#
# MAIN
#

# Create output directory.
if not os.path.exists(output_path):
    os.makedirs(output_path)

# Create temporary directory.
tmp_path = tempfile.mkdtemp()
print "Working in temporary directory: %s" % tmp_path
tmp_path = output_path

# Open file to write all exceptions that occur during execution.
exception_outfile = open(exception_filename, 'a')
run_index_outfile = open(run_index_filename, 'a')
runs = len(existing_mutants) # number for RUN to set up
for (name, mutant) in zip(mutant_names, mutant_codes):
    if name in existing_mutants:
        # Skip this.
        print "%s : %s already exists, skipping" % (name, str(mutant))        
        continue

    print "%s : %s" % (name, str(mutant))
    simulation = None
    if True:
#    try:

        # Create directory to store files in.
        workdir = os.path.join(tmp_path, name)
        if not os.path.exists(workdir):
            os.makedirs(workdir)
            print "Creating path %s" % workdir

        # Create PDBFixer, retrieving PDB template
        fixer = pdbfixer.PDBFixer(filename=pdbfilename)
#        pdb_filename = os.path.join(workdir, 'just_loaded.pdb')
#        outfile = open(pdb_filename, 'w')
#        app.PDBFile.writeFile(fixer.topology, fixer.positions, outfile, keepIds=True)
#        outfile.close()

        # Attempt to make mutations.
        if len(mutant) > 0:
            try:
                fixer.applyMutations(mutant, chain_id_to_mutate)
            except Exception as e:
                # Mutant could not be constructed.
                print e
                exception_outfile.write("%s : %s" % (name, str(mutant)) + '\n')
                exception_outfile.write(str(e) + '\n')
                continue

        modeller = app.Modeller(fixer.topology, fixer.positions)
        pdb_filename = os.path.join(workdir, 'complex.pdb')
        with open(pdb_filename, 'w') as fo:
            app.PDBFile.writeFile(fixer.topology, fixer.positions, fo, keepIds=True)

        # If everything worked, add this RUN.
        run_name = 'RUN%d' % runs
        run_dir = os.path.join(output_path, run_name)
        shutil.move(workdir, run_dir)
        run_index_outfile.write('%s %s\n' % (run_name, name))
        run_index_outfile.flush()
        runs += 1

    if False:
        e = 'No'
#    except Exception as e:
        print str(e)
        exception_outfile.write("%s : %s : %s\n" % (name, str(mutant), str(e)))
        exception_outfile.flush()

        # Clean up.
        if simulation:
            if simulation.context: del simulation.context        

exception_outfile.close()
run_index_outfile.close()








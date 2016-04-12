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

print("Input PDB structure: 1OL5")
pdbid = '1OL5'

if TPX2:
    output_identifier = '+'
else:
    output_identifier = '-'
# Path to put all output in
output_path = "output"+pdbid+output_identifier+"TPX2"

# Source PDB
pdbfilename = pdbid+"-WT-pdbfixer.pdb"
pdbfilename = pdbid+"-modifiedMG-WT.pdb"

adp_mol2 = "ADP"+pdbid[-1]+".mol2"

print "Source PDB filename: %s" % pdbfilename
print "Output directory for mutations: %s" % output_path

#
# PARAMETERS
#

chain_id_to_mutate = 'A' # chain to mutate
pH = 7.4 # pH to model
keep_crystallographic_water = True # not using

# Forcefield
ff_name = 'amber99sbildn'
water_name = 'tip3p'
ion_ff_name = 'mg'
ADP_ff_name = 'adp'

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

mutant_names = list()
mutant_codes = list()

# Append wild type (no mutation).
for i in range(5):
    mutant_names.append('WT')
    mutant_codes.append([])

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
for name, mutant in enumerate(mutant_codes):
    runs = name
    name = str(name)
    print "%s : %s" % (name, str(mutant))
    simulation = None
    try:
        # Create PDBFixer, retrieving PDB template
        print("creating Modeller...")
        pdb = app.PDBFile(pdbfilename, extraParticleIdentifier='')
        modeller = app.Modeller(pdb.topology, pdb.positions)
        print('removing ADP')
        for residue in modeller.topology.residues():
            if residue.name == 'ADP':
                modeller.delete([residue])
        print('fixing dummy particles')
        ion_forcefield = app.ForceField(ion_ff_name+'.xml', ff_name+'.xml', water_name+'.xml', ADP_ff_name+'.xml')
        modeller.addExtraParticles(ion_forcefield)
        if not TPX2:
            print('removing TPX2...')
            for chain in modeller.topology.chains():
                if chain.id == 'B':
                    to_delete = [chain]
            modeller.delete(to_delete)
        print("adding hydrogens...")
        modeller.addHydrogens(pH=pH)
        #fixer.addSolvent(fixer.topology.getUnitCellDimensions())

        # Create directory to store files in.
        workdir = os.path.join(tmp_path, 'RUN'+name)
        if not os.path.exists(workdir):
            os.makedirs(workdir)
            print "Creating path %s" % workdir

        # Solvate
        if verbose: print "Loading forcefield..."
        forcefield = app.ForceField(ff_name+'.xml',water_name+'.xml',ion_ff_name+'.xml',ADP_ff_name+'.xml')

        # Convert positions to numpy format and remove offset
        if verbose: print "Subtracting offset..."
        modeller.positions = unit.Quantity(np.array(modeller.positions / unit.nanometer), unit.nanometer)
        offset = modeller.positions[0,:]
        if verbose: print "Shifting model by %s" % str(offset)
        modeller.positions -= offset

        # Add correct ADP (with hydrogens)
        if verbose: print "Replacing ADP with bonded ADP..."
        adp = md.load_mol2(adp_mol2)
        adp.topology = adp.top.to_openmm()
        adp.positions = unit.Quantity(np.array([(x,y,z) for x,y,z in adp.xyz[0]]), unit.nanometer)
        if verbose: print "Shifting protonated ADP by %s" % str(offset)
        adp.positions -= offset
        modeller.add(adp.topology,adp.positions)

        # Write PDB file for solute only.
        if verbose: print "Writing pdbfixer output..."
        pdb_filename = os.path.join(workdir, 'pdbfixer.pdb')
        outfile = open(pdb_filename, 'w')
        app.PDBFile.writeFile(modeller.topology, modeller.positions, outfile, keepIds=True)
        outfile.close()

        # Create OpenMM system.
        if verbose: print "Creating OpenMM system..."
        #system = forcefield.createSystem(modeller.topology, nonbondedMethod=nonbonded_method, nonbondedCutoff=nonbonded_cutoff, constraints=app.HBonds)
        system = forcefield.createSystem(modeller.topology, nonbondedMethod=nonbonded_method, nonbondedCutoff=nonbonded_cutoff, constraints=None)
        if verbose: print "Adding barostat..."
        system.addForce(openmm.MonteCarloBarostat(pressure, temperature, barostat_frequency))

        # Create simulation.
        if verbose: print "Creating simulation..."
        integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
        #platform = openmm.Platform.getPlatformByName('CPU')
        platform = openmm.Platform.getPlatformByName('OpenCL')
        platform.setPropertyDefaultValue('OpenCLPrecision', 'double') # use double precision
        simulation = app.Simulation(modeller.topology, system, integrator, platform=platform)
        simulation.context.setPositions(modeller.positions)

        # Compute forces.
#        if verbose: print "Computing forces..."
#        forces = simulation.context.getState(getForces=True).getForces(asNumpy=True)
#        print forces
#        atoms = [ atom for atom in simulation.topology.atoms() ]
#        force_unit = unit.kilojoules_per_mole / unit.nanometers
#        for (index, atom) in enumerate(atoms):            
#            force_norm = np.sqrt(((forces[index,:] / force_unit)**2).sum())
#            if force_norm > 1000000.0:
#                print "%8d %8s %20s %8d %8s %5d : %24f kJ/nm" % (index, atom.name, str(atom.element), atom.index, atom.residue.name, atom.residue.index, force_norm)    
        
        # Write modeller positions.
#        if verbose: print "Writing modeller output..."
#        filename = os.path.join(workdir, 'modeller.pdb')
#        positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
#        print abs(positions / unit.nanometers).max()
#        app.PDBFile.writeFile(simulation.topology, positions, open(filename, 'w'), keepIds=True)

        # Minimize energy.
        if verbose: print "Minimizing energy..."
        potential_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
        if numpy.isnan(potential_energy / unit.kilocalories_per_mole):
            raise Exception("Potential energy is NaN before minimization.")
        if verbose: print "Initial potential energy : %10.3f kcal/mol" % (potential_energy / unit.kilocalories_per_mole)
        simulation.minimizeEnergy(maxIterations=max_minimization_iterations)
        potential_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
        if numpy.isnan(potential_energy / unit.kilocalories_per_mole):
            raise Exception("Potential energy is NaN after minimization.")
        if verbose: print "Final potential energy:  : %10.3f kcal/mol" % (potential_energy / unit.kilocalories_per_mole)

        del(modeller)
        positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
        modeller = app.Modeller(simulation.topology, positions)

        del(system)
        del(integrator)
        del(platform)
        del(simulation.context)
        del(simulation)
        simulation = None

        if verbose: print "Creating constrained OpenMM system..."
        system = forcefield.createSystem(modeller.topology, nonbondedMethod=nonbonded_method, nonbondedCutoff=nonbonded_cutoff, constraints=app.HBonds)
        if verbose: print "Adding barostat..."
        system.addForce(openmm.MonteCarloBarostat(pressure, temperature, barostat_frequency))

        # Create simulation.
        if verbose: print "Creating simulation..."
        integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
        #platform = openmm.Platform.getPlatformByName('CPU')
        platform = openmm.Platform.getPlatformByName('OpenCL')
        platform.setPropertyDefaultValue('OpenCLPrecision', 'double') # use double precision
        simulation = app.Simulation(modeller.topology, system, integrator, platform=platform)
        simulation.context.setPositions(modeller.positions)

        # Write modeller positions.
        if verbose: print "Writing modeller output..."
        filename = os.path.join(workdir, 'modeller.pdb')
        positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
        print abs(positions / unit.nanometers).max()
        app.PDBFile.writeFile(simulation.topology, positions, open(filename, 'w'), keepIds=True)

        # Minimize energy.
        if verbose: print "Minimizing energy..."
        potential_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
        if numpy.isnan(potential_energy / unit.kilocalories_per_mole):
            raise Exception("Potential energy is NaN before minimization.")
        if verbose: print "Initial potential energy : %10.3f kcal/mol" % (potential_energy / unit.kilocalories_per_mole)
        simulation.minimizeEnergy(maxIterations=max_minimization_iterations)
        potential_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
        if numpy.isnan(potential_energy / unit.kilocalories_per_mole):
            raise Exception("Potential energy is NaN after minimization.")
        if verbose: print "Final potential energy:  : %10.3f kcal/mol" % (potential_energy / unit.kilocalories_per_mole)


        if solvate:
            # Write initial positions.
            filename = os.path.join(workdir, 'minimized_vacuum.pdb')
            positions = simulation.context.getState(getPositions=True).getPositions()
            app.PDBFile.writeFile(simulation.topology, positions, open(filename, 'w'), keepIds=True)

            del(modeller)
            positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
            modeller = app.Modeller(simulation.topology, positions)

            del(system)
            del(integrator)
            del(platform)
            del(simulation.context)
            del(simulation)
            simulation = None

            # Solvate and create system.
            if verbose: print "Solvating with %s..." % water_name
            modeller.addSolvent(forcefield, padding=padding, model=water_name, ionicStrength=ionicStrength)

            # Separate waters into a separate 'W' chain with renumbered residues.
            # This is kind of a hack, as waters are numbered 1-9999 and then we repeat
            print "Renumbering waters..."
            for chain in modeller.topology.chains(): water_id = chain.id
            nwaters_and_ions = 0
            for residue in modeller.topology.residues():
                if residue.chain.id == water_id:
                    residue.id = (nwaters_and_ions % 9999) + 1
                    nwaters_and_ions += 1
            print "System contains %d waters and ions." % nwaters_and_ions

            # Create OpenMM system.
            if verbose: print "Creating solvated OpenMM system..."
            system = forcefield.createSystem(modeller.topology, nonbondedMethod=nonbonded_method, nonbondedCutoff=nonbonded_cutoff, constraints=app.HBonds)
            #system = forcefield.createSystem(modeller.topology, nonbondedMethod=nonbonded_method, nonbondedCutoff=nonbonded_cutoff, constraints=None)
            if verbose: print "Adding barostat..."
            system.addForce(openmm.MonteCarloBarostat(pressure, temperature, barostat_frequency))

            # Create simulation.
            if verbose: print "Creating solvated simulation..."
            integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
            #platform = openmm.Platform.getPlatformByName('CPU')
            platform = openmm.Platform.getPlatformByName('OpenCL')
            platform.setPropertyDefaultValue('OpenCLPrecision', 'double') # use double precision
            simulation = app.Simulation(modeller.topology, system, integrator, platform=platform)
            simulation.context.setPositions(modeller.positions)

            # Compute forces.
#            if verbose: print "Computing forces..."
#            forces = simulation.context.getState(getForces=True).getForces(asNumpy=True)
#            print forces
#            atoms = [ atom for atom in simulation.topology.atoms() ]
#            force_unit = unit.kilojoules_per_mole / unit.nanometers
#            for (index, atom) in enumerate(atoms):
#                force_norm = np.sqrt(((forces[index,:] / force_unit)**2).sum())
#                if force_norm > 10.0:
#                    print "%8d %8s %20s %8d %8s %5d : %24f kJ/nm" % (index, atom.name, str(atom.element), atom.index, atom.residue.name, atom.residue.index, force_norm)

            # Write modeller positions.
            if verbose: print "Writing modeller output..."
            filename = os.path.join(workdir, 'modeller_solvent.pdb')
            positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
            print abs(positions / unit.nanometers).max()
            app.PDBFile.writeFile(simulation.topology, positions, open(filename, 'w'), keepIds=True)

            # Minimize energy.
            if verbose: print "Minimizing energy..."
            potential_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
            if numpy.isnan(potential_energy / unit.kilocalories_per_mole):
                raise Exception("Potential energy is NaN before minimization.")
            if verbose: print "Initial solvated potential energy : %10.3f kcal/mol" % (potential_energy / unit.kilocalories_per_mole)
            simulation.minimizeEnergy(maxIterations=max_minimization_iterations)
            potential_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
            if numpy.isnan(potential_energy / unit.kilocalories_per_mole):
                raise Exception("Potential energy is NaN after minimization.")
            if verbose: print "Final solvated potential energy:  : %10.3f kcal/mol" % (potential_energy / unit.kilocalories_per_mole)


        # Write initial positions.
        filename = os.path.join(workdir, 'minimized.pdb')
        positions = simulation.context.getState(getPositions=True).getPositions()
        app.PDBFile.writeFile(simulation.topology, positions, open(filename, 'w'), keepIds=True)

        # Assign temperature
        simulation.context.setVelocitiesToTemperature(temperature)

        # Take a few steps to relax structure.
        if verbose: print "Taking a few steps..."
        simulation.step(nsteps)

        # Write initial positions.
        if verbose: print "Writing positions..."
        filename = os.path.join(workdir, 'system.pdb')
        positions = simulation.context.getState(getPositions=True).getPositions()
        app.PDBFile.writeFile(simulation.topology, positions, open(filename, 'w'), keepIds=True)

        # Write mutation.
        filename = os.path.join(workdir, 'mutation.txt')
        outfile = open(filename, 'w')
        outfile.write('%s\n' % name)
        outfile.flush()
        outfile.close()

        # Serialize to XML files.
        if verbose: print "Serializing to XML..."
        system_filename = os.path.join(workdir, 'system.xml')
        integrator_filename = os.path.join(workdir, 'integrator.xml')
        write_file(system_filename, openmm.XmlSerializer.serialize(system))
        write_file(integrator_filename, openmm.XmlSerializer.serialize(integrator))
        simulation.context.setVelocitiesToTemperature(temperature)
        state = simulation.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)
        state_filename = os.path.join(workdir, 'state.xml')
        serialized = openmm.XmlSerializer.serialize(state)
        write_file(state_filename, serialized)

        # If everything worked, add this RUN.
        run_name = 'RUN%d' % runs
        run_dir = os.path.join(output_path, run_name)
        shutil.move(workdir, run_dir)
        run_index_outfile.write('%s %s\n' % (run_name, name))
        run_index_outfile.flush()
        runs += 1

        # Clean up.
        del simulation.context
        del simulation
        del system
        del positions

    except Exception as e:
        print str(e)
        exception_outfile.write("%s : %s : %s\n" % (name, str(mutant), str(e)))
        exception_outfile.flush()

        # Clean up.
        if simulation:
            if simulation.context: del simulation.context        

exception_outfile.close()
run_index_outfile.close()








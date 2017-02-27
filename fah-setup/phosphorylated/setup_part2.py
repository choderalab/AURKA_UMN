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

TPX2 = True

if TPX2:
    id = '+'
else:
    id = '-'
# Path to put all output in
output_path = "1OL5"+id+"TPX2-TPOs"

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

# Verbosity level
verbose = True

def write_file(filename, contents):
    with open(filename, 'w') as outfile:
        outfile.write(contents)

exception_filename = os.path.join(output_path, 'exceptions.out') # to store exceptions
run_index_filename = os.path.join(output_path, 'run-index.txt') # to store index of which mutants are which

point_mutants = ['Q185C', 'Q185L','Q185M','Q185N','Q185H','C247A','C247L']

# Create temporary directory.
tmp_path = tempfile.mkdtemp()
print "Working in temporary directory: %s" % tmp_path
tmp_path = output_path

# Open file to write all exceptions that occur during execution.
exception_outfile = open(exception_filename, 'a')
run_index_outfile = open(run_index_filename, 'a')

for runs, name in enumerate(point_mutants):

    print('Building %s' % name)
    run_name = 'RUN%d' % runs
    run_dir = os.path.join(output_path, run_name)
    print(run_dir)
    inpcrdfilename = os.path.join(run_dir, 'complex.inpcrd')
    prmtopfilename = os.path.join(run_dir, 'complex.prmtop')

    simulation = None
    if True:

        if verbose: print("Loading inpcrd file")
        inpcrd = app.AmberInpcrdFile(inpcrdfilename)
        if verbose: print("Loading prmtop file")
        prmtop = app.AmberPrmtopFile(prmtopfilename)

        if verbose: print "Creating OpenMM system..."
        system = prmtop.createSystem(nonbondedMethod=nonbonded_method, nonbondedCutoff=nonbonded_cutoff, constraints=None, temperature=temperature)
        if verbose: print "Adding barostat..."
        system.addForce(openmm.MonteCarloBarostat(pressure, temperature, barostat_frequency))

        # Create simulation.
        if verbose: print "Creating simulation..."
        integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
        #platform = openmm.Platform.getPlatformByName('CPU')
        platform = openmm.Platform.getPlatformByName('OpenCL')
        platform.setPropertyDefaultValue('OpenCLPrecision', 'double') # use double precision
        simulation = app.Simulation(prmtop.topology, system, integrator, platform=platform)
        simulation.context.setPositions(inpcrd.positions)

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

        positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
        modeller = app.Modeller(simulation.topology, positions)

        del(system)
        del(integrator)
        del(platform)
        del(simulation.context)
        del(simulation)
        simulation = None

        if verbose: print "Creating constrained OpenMM system..."
        system = prmtop.createSystem(nonbondedMethod=nonbonded_method, nonbondedCutoff=nonbonded_cutoff, constraints=app.HBonds, temperature=temperature)
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
        filename = os.path.join(run_dir, 'modeller.pdb')
        positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
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

        # Write initial positions.
        filename = os.path.join(run_dir, 'minimized.pdb')
        positions = simulation.context.getState(getPositions=True).getPositions()
        app.PDBFile.writeFile(simulation.topology, positions, open(filename, 'w'), keepIds=True)

        # Assign temperature
        simulation.context.setVelocitiesToTemperature(temperature)

        # Take a few steps to relax structure.
        if verbose: print "Taking a few steps..."
        simulation.step(nsteps)


        # Write initial positions.
        if verbose: print "Writing positions..."
        filename = os.path.join(run_dir, 'system.pdb')
        positions = simulation.context.getState(getPositions=True).getPositions()
        app.PDBFile.writeFile(simulation.topology, positions, open(filename, 'w'), keepIds=True)

        # Write mutation.
        filename = os.path.join(run_dir, 'mutation.txt')
        outfile = open(filename, 'w')
        outfile.write('%s\n' % name)
        outfile.flush()
        outfile.close()

        # Serialize to XML files.
        if verbose: print "Serializing to XML..."
        system_filename = os.path.join(run_dir, 'system.xml')
        integrator_filename = os.path.join(run_dir, 'integrator.xml')
        write_file(system_filename, openmm.XmlSerializer.serialize(system))
        write_file(integrator_filename, openmm.XmlSerializer.serialize(integrator))
        simulation.context.setVelocitiesToTemperature(temperature)
        state = simulation.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)
        state_filename = os.path.join(run_dir, 'state.xml')
        serialized = openmm.XmlSerializer.serialize(state)
        write_file(state_filename, serialized)

        # If everything worked, add this RUN.
        run_index_outfile.write('%s %s\n' % (run_name, name))
        run_index_outfile.flush()

        # Clean up.
        del simulation.context
        del simulation
        del system
        del positions

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








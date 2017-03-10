"""
Set up AurKA mutant simulations
@author John D. Chodera
@date 8 Aug 2014
"""

#
# IMPORTS
#
from __future__ import division, print_function

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

import sys

# ParmEd Imports
#from parmed.charmm import CharmmPsfFile, CharmmCrdFile, CharmmParameterSet
from simtk.openmm.app import CharmmPsfFile, CharmmCrdFile, CharmmParameterSet

from parmed import unit as u ### ???

#


print("Input PDB structure: 1OL5")
pdbid = '1OL5'
runs = 0

output_path = '1OL5-phos-tpx2'

# Source PDB
pdbfilename = 'step2_solvator.pdb'

print("Source PDB filename: %s" % pdbfilename)
print("Output directory for mutations: %s" % output_path)

#
# PARAMETERS
#

universal_parameter_files = [
    'toppar/par_all36_prot.prm',
    'toppar/par_all36_na.prm',
    'toppar/top_all36_prot.rtf',
    'toppar/toppar_all36_na_nad_ppi.str',
    'toppar/toppar_water_ions.str',
    'step2_solvator.str',
    'step2.2_ions.prm',
    'step2.1_waterbox.prm',
    'step1_pdbreader.str',
    'step1_labelrot.str',
    'toppar/toppar_all36_label_spin.str',
    'step1_pdbreader.str',
    'toppar/top_all36_cgenff.rtf',
    'toppar/top_all36_na.rtf',
    'toppar/toppar_all36_na_rna_modified.str',
    'toppar/toppar_all36_prot_na_combined.str',
    'toppar/par_all36_cgenff.prm',
    'toppar/par_all36_na.prm',
]

psf_file = 'step2_solvator.psf'

padding = 11.0 * unit.angstroms
nonbonded_cutoff = 9.0 * unit.angstroms
nonbonded_method = app.PME
max_minimization_iterations = 10000
temperature = 300.0 * unit.kelvin
pressure = 1.0 * unit.atmospheres
collision_rate = 5.0 / unit.picoseconds
barostat_frequency = 50
timestep = 2.0 * unit.femtoseconds
nsteps = 50000 # number of steps to take for testing
ionicStrength = 20 * unit.millimolar

# Verbosity level
verbose = True

def write_file(filename, contents):
    with open(filename, 'w') as outfile:
        outfile.write(contents)

exception_filename = os.path.join(output_path, 'exceptions.out') # to store exceptions
run_index_filename = os.path.join(output_path, 'run-index.txt') # to store index of which mutants are which

# Create output directory.
if not os.path.exists(output_path):
    os.makedirs(output_path)

# Create temporary directory.
tmp_path = tempfile.mkdtemp()
print("Working in temporary directory: %s" % tmp_path)
tmp_path = output_path

# Open file to write all exceptions that occur during execution.
exception_outfile = open(exception_filename, 'a')
run_index_outfile = open(run_index_filename, 'a')


name = str(runs)
simulation = None


# Load the CHARMM files
print('Loading CHARMM files...')
param_files = universal_parameter_files #+ charmm_parameter_files
params = CharmmParameterSet(*param_files)
psf = CharmmPsfFile(psf_file)
pdb = app.PDBFile(pdbfilename)

coords = pdb.positions
min_crds = [coords[0][0], coords[0][1], coords[0][2]]
max_crds = [coords[0][0], coords[0][1], coords[0][2]]

for coord in coords:
    min_crds[0] = min(min_crds[0], coord[0])
    min_crds[1] = min(min_crds[1], coord[1])
    min_crds[2] = min(min_crds[2], coord[2])
    max_crds[0] = max(max_crds[0], coord[0])
    max_crds[1] = max(max_crds[1], coord[1])
    max_crds[2] = max(max_crds[2], coord[2])

psf.setBox(max_crds[0]-min_crds[0],
           max_crds[1]-min_crds[1],
           max_crds[2]-min_crds[2],
           )

# Create PDBFixer, retrieving PDB template
print("creating Modeller...")
modeller = app.Modeller(pdb.topology, pdb.positions)

# Create directory to store files in.
workdir = os.path.join(tmp_path, 'RUN'+name)
if not os.path.exists(workdir):
    os.makedirs(workdir)
    print("Creating path %s" % workdir)

# Write PDB file for solute only.
if verbose: print("Writing initial output...")
pdb_filename = os.path.join(workdir, 'initial.pdb')
outfile = open(pdb_filename, 'w')
app.PDBFile.writeFile(modeller.topology, modeller.positions, outfile)
outfile.close()

# Create OpenMM system.
if verbose: print("Creating OpenMM system...")
#system = psf.createSystem(params, nonbondedMethod=nonbonded_method, nonbondedCutoff=nonbonded_cutoff, constraints=app.HBonds)
system = psf.createSystem(params, nonbondedMethod=nonbonded_method, nonbondedCutoff=nonbonded_cutoff, constraints=None)
if verbose: print("Adding barostat...")
system.addForce(openmm.MonteCarloBarostat(pressure, temperature, barostat_frequency))

# Create simulation.
if verbose: print("Creating simulation...")
integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
#platform = openmm.Platform.getPlatformByName('CPU')
platform = openmm.Platform.getPlatformByName('OpenCL')
platform.setPropertyDefaultValue('OpenCLPrecision', 'double') # use double precision
simulation = app.Simulation(modeller.topology, system, integrator, platform=platform)
try:
    simulation.context.setPositions(modeller.positions)
except Exception as e:
    print(len(modeller.positions))
    print(simulation.context.getSystem().getNumParticles())
    raise(e)

# Minimize energy.
if verbose: print("Minimizing energy...")
potential_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
if numpy.isnan(potential_energy / unit.kilocalories_per_mole):
    raise Exception("Potential energy is NaN before minimization.")
if verbose: print("Initial potential energy : %10.3f kcal/mol" % (potential_energy / unit.kilocalories_per_mole))
simulation.minimizeEnergy(maxIterations=max_minimization_iterations)
potential_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
if numpy.isnan(potential_energy / unit.kilocalories_per_mole):
    raise Exception("Potential energy is NaN after minimization.")
if verbose: print("Final potential energy:  : %10.3f kcal/mol" % (potential_energy / unit.kilocalories_per_mole))

del(modeller)
positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
modeller = app.Modeller(simulation.topology, positions)

del(system)
del(integrator)
del(platform)
del(simulation.context)
del(simulation)
simulation = None

if verbose: print("Creating constrained OpenMM system...")
coords = modeller.positions
min_crds = [coords[0][0], coords[0][1], coords[0][2]]
max_crds = [coords[0][0], coords[0][1], coords[0][2]]

for coord in coords:
    min_crds[0] = min(min_crds[0], coord[0])
    min_crds[1] = min(min_crds[1], coord[1])
    min_crds[2] = min(min_crds[2], coord[2])
    max_crds[0] = max(max_crds[0], coord[0])
    max_crds[1] = max(max_crds[1], coord[1])
    max_crds[2] = max(max_crds[2], coord[2])

psf.setBox(max_crds[0]-min_crds[0],
           max_crds[1]-min_crds[1],
           max_crds[2]-min_crds[2],
)
system = psf.createSystem(params, nonbondedMethod=nonbonded_method, nonbondedCutoff=nonbonded_cutoff, constraints=app.HBonds)
if verbose: print("Adding barostat...")
system.addForce(openmm.MonteCarloBarostat(pressure, temperature, barostat_frequency))

# Create simulation.
if verbose: print("Creating simulation...")
integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
#platform = openmm.Platform.getPlatformByName('CPU')
platform = openmm.Platform.getPlatformByName('OpenCL')
platform.setPropertyDefaultValue('OpenCLPrecision', 'double') # use double precision
simulation = app.Simulation(modeller.topology, system, integrator, platform=platform)
simulation.context.setPositions(modeller.positions)

# Write modeller positions.
if verbose: print("Writing modeller output...")
filename = os.path.join(workdir, 'modeller.pdb')
positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
print(abs(positions / unit.nanometers).max())
app.PDBFile.writeFile(simulation.topology, positions, open(filename, 'w'))

# Minimize energy.
if verbose: print("Minimizing energy...")
potential_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
if numpy.isnan(potential_energy / unit.kilocalories_per_mole):
    raise Exception("Potential energy is NaN before minimization.")
if verbose: print("Initial potential energy : %10.3f kcal/mol" % (potential_energy / unit.kilocalories_per_mole))
simulation.minimizeEnergy(maxIterations=max_minimization_iterations)
potential_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
if numpy.isnan(potential_energy / unit.kilocalories_per_mole):
    raise Exception("Potential energy is NaN after minimization.")
if verbose: print("Final potential energy:  : %10.3f kcal/mol" % (potential_energy / unit.kilocalories_per_mole))


# Write initial positions.
filename = os.path.join(workdir, 'minimized.pdb')
positions = simulation.context.getState(getPositions=True).getPositions()
app.PDBFile.writeFile(simulation.topology, positions, open(filename, 'w'))

# Assign temperature
simulation.context.setVelocitiesToTemperature(temperature)

# Take a few steps to relax structure.
if verbose: print("Taking a few steps...")
simulation.step(nsteps)

# Write initial positions.
if verbose: print("Writing positions...")
filename = os.path.join(workdir, 'system.pdb')
positions = simulation.context.getState(getPositions=True).getPositions()
app.PDBFile.writeFile(simulation.topology, positions, open(filename, 'w'))

# Write mutation.
filename = os.path.join(workdir, 'mutation.txt')
outfile = open(filename, 'w')
outfile.write('%s\n' % name)
outfile.flush()
outfile.close()

# Serialize to XML files.
if verbose: print("Serializing to XML...")
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


exception_outfile.close()
run_index_outfile.close()








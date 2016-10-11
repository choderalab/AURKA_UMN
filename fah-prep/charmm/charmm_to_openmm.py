from __future__ import division, print_function

import sys

# OpenMM Imports
import simtk.openmm as mm
import simtk.openmm.app as app

# ParmEd Imports
#from parmed.charmm import CharmmPsfFile, CharmmCrdFile, CharmmParameterSet
from simtk.openmm.app import CharmmPsfFile, CharmmCrdFile, CharmmParameterSet
from parmed.openmm import StateDataReporter
from parmed import unit as u

param_files = ['par_all36_prot.prm','par_all36_na.prm','top_all36_prot.rtf','toppar_all36_na_nad_ppi.str','toppar_water_ions.str']

# Load the CHARMM files
print('Loading CHARMM files...')
params = CharmmParameterSet(*param_files)
ala5_gas = CharmmPsfFile('charmm_1OL5+TPX2.psf')
ala5_crds = app.PDBFile('charmm_1OL5+TPX2.pdb')

# NOTE NOTE
# The parameter set we used here is the CHARMM 36 force field, but this is
# strictly an example. It is important that you use the most accurate (typically
# most up-to-date) force fields for your own simulation. See the CHARMM
# parameter web page for updates:
# http://mackerell.umaryland.edu/CHARMM_ff_params.html
# END NOTE
coords = ala5_crds.positions
min_crds = [coords[0][0], coords[0][1], coords[0][2]]
max_crds = [coords[0][0], coords[0][1], coords[0][2]]

for coord in coords:
    min_crds[0] = min(min_crds[0], coord[0])
    min_crds[1] = min(min_crds[1], coord[1])
    min_crds[2] = min(min_crds[2], coord[2])
    max_crds[0] = max(max_crds[0], coord[0])
    max_crds[1] = max(max_crds[1], coord[1])
    max_crds[2] = max(max_crds[2], coord[2])

ala5_gas.setBox(max_crds[0]-min_crds[0],
                   max_crds[1]-min_crds[1],
                   max_crds[2]-min_crds[2],
)


# Create the OpenMM system
print('Creating OpenMM System')
system = ala5_gas.createSystem(params, nonbondedMethod=app.PME,
                               constraints=app.HBonds)

# Create the integrator to do Langevin dynamics
integrator = mm.LangevinIntegrator(
                        300*u.kelvin,       # Temperature of heat bath
                        1.0/u.picoseconds,  # Friction coefficient
                        2.0*u.femtoseconds, # Time step
)

# Define the platform to use; CUDA, OpenCL, CPU, or Reference. Or do not specify
# the platform to use the default (fastest) platform
platform = mm.Platform.getPlatformByName('CPU')

# Create the Simulation object
sim = app.Simulation(ala5_gas.topology, system, integrator, platform)

# Set the particle positions
sim.context.setPositions(ala5_crds.positions)

# Minimize the energy
print('Minimizing energy')
sim.minimizeEnergy(maxIterations=500)


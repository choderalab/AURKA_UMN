import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as u
import mdtraj as md
import numpy 
import pdbfixer

use_fixer = True
fixer_steps = True
solvate = True
include_ADP = False

if use_fixer:
    fixer = pdbfixer.PDBFixer("1OL5-pdbfixer.pdb")
else:
    pdb = app.PDBFile("1OL5-pdbfixer.pdb")

if solvate:
    ff = app.ForceField("amber99sbildn.xml", "tip3p.xml")
else:
    ff = app.ForceField("amber99sbildn.xml")
print("loaded topology and forcefield")

if fixer_steps:
    print "findMissingResidues..."
    fixer.missingResidues = {}
    print "findNonstandardResidues..."
    fixer.findNonstandardResidues()
    print "replaceNonstandardResidues..."
    fixer.replaceNonstandardResidues()
    print "findMissingAtoms..."
    fixer.findMissingAtoms()
    print "addMissingAtoms..."
    fixer.addMissingAtoms()
    print "addingmissinghydrogens..."
    fixer.addMissingHydrogens(7.4)

if use_fixer:
    modeller = app.Modeller(fixer.topology, fixer.positions)
else:
    modeller = app.Modeller(pdb.topology, pdb.positions)

if include_ADP:
    adp = md.load_mol2('ADP5.mol2')
    adp.topology = adp.top.to_openmm()
    adp.positions = unit.Quantity(np.array([(x,y,z) for x,y,z in adp.xyz[0]]), u.nanometer)
    modeller.add(adp.topology,adp.positions)
    

if solvate:
    print("solvating")
    padding = 11.0 * u.angstroms
    water_name = 'tip3p'
    ionicStrength = 300 * u.millimolar
    modeller.addSolvent(ff, padding=padding, model=water_name, ionicStrength=ionicStrength)

system = ff.createSystem(modeller.topology,nonbondedMethod=app.PME,constraints=app.HBonds)
integrator = mm.LangevinIntegrator(300.0*u.kelvin , 10.0 / u.picoseconds, 2.0*u.femtosecond)
system.addForce(mm.MonteCarloBarostat(1.0 * u.atmospheres, 300.0*u.kelvin, 25))
print("systemed")

simulation = app.Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)
print('minimizing')

potential_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
if numpy.isnan(potential_energy / u.kilocalories_per_mole):
    raise Exception("Potential energy is NaN before minimization.")
print "Initial potential energy : %10.3f kcal/mol" % (potential_energy / u.kilocalories_per_mole)
simulation.minimizeEnergy()
potential_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
if numpy.isnan(potential_energy / u.kilocalories_per_mole):
    raise Exception("Potential energy is NaN after minimization.")
print "Final potential energy:  : %10.3f kcal/mol" % (potential_energy / u.kilocalories_per_mole)

simulation.context.setVelocitiesToTemperature(300.0*u.kelvin)

print('production')
simulation.step(5000)
positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
if numpy.isnan(positions / u.nanometers).any() : print("Everything is broken and terrible")
print("Complete!")


import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as u
import mdtraj as md
import numpy 
import pdbfixer

solvate = False

ff = app.ForceField("amber99sbildn.xml")
ff.loadFile("ADP.xml")
if solvate:
    ff.loadFile("tip3p.xml")
print("loaded forcefield")


print("adding ADP")
adp = md.load_mol2('ADP5.mol2')
adp.topology = adp.top.to_openmm()
xcoord = [x for x,y,z in adp.xyz[0]]
ycoord = [y for x,y,z in adp.xyz[0]]
zcoord = [z for x,y,z in adp.xyz[0]]
xVec = mm.Vec3(max(xcoord)-min(xcoord)+2.00,0.00,0.00)
yVec = mm.Vec3(0.00,max(ycoord)-min(ycoord)+2.00,0.00)
zVec = mm.Vec3(0.00,0.00,max(zcoord)-min(zcoord)+2.00)
adp.positions = u.Quantity(numpy.array([(x,y,z) for x,y,z in adp.xyz[0]]), u.nanometer)
modeller = app.Modeller(adp.topology,adp.positions)
modeller.topology.setPeriodicBoxVectors((xVec, yVec, zVec))

    
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
positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
if numpy.isnan(potential_energy / u.kilocalories_per_mole):
    raise Exception("Potential energy is NaN before minimization.")
if numpy.isnan(positions / u.nanometers).any() : raise Exception("Positions are NaN before minimization.")
print "Initial potential energy : %10.3f kcal/mol" % (potential_energy / u.kilocalories_per_mole)
with open('minimizing.pdb','w') as filename:
    app.PDBFile.writeHeader(simulation.topology, filename)
    for i in range(200):
        simulation.minimizeEnergy(maxIterations=1)
        positions = simulation.context.getState(getPositions=True).getPositions()
        app.PDBFile.writeModel(simulation.topology, positions, file=filename, modelIndex=i, keepIds=True)
#    app.PDBFile.writeFooter(simulation.topology, file=filename)
    potential_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    if numpy.isnan(potential_energy / u.kilocalories_per_mole):
        raise Exception("Potential energy is NaN after minimization.")
    print "Final potential energy:  : %10.3f kcal/mol" % (potential_energy / u.kilocalories_per_mole)
    positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
    if numpy.isnan(positions / u.nanometers).any() : print("Everything is broken and terrible")

    print('Setting velocities to temp')
    simulation.context.setVelocitiesToTemperature(300.0*u.kelvin)

    print('production')
    for i in range(100):
        simulation.step(1)
        positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
        if numpy.isnan(positions / u.nanometers).any(): 
            print("Everything is broken and terrible")
            break
        try:
            app.PDBFile.writeModel(simulation.topology, positions, file=filename, modelIndex=i+200, keepIds=True)
        except:
            print("The thing broke on step "+str(i))
            break
    app.PDBFile.writeFooter(simulation.topology, file=filename)

print("Complete!")


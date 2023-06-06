from openmm import *
from openmm.app import *
from openmm.unit import *

# Monkey-patch the DCD.writeModel method to avoid writing dummy-atom positions.
writeModel = DCDFile.writeModel
def writeModelPatched(self, positions, unitCellDimensions=None, periodicBoxVectors=None):
    positions = positions[:len(list(self._topology.atoms()))]
    writeModel(self,
               positions,
               unitCellDimensions=unitCellDimensions,
               periodicBoxVectors=periodicBoxVectors)
DCDFile.writeModel = writeModelPatched

# Load the topology and coordinate files.
prmtop = AmberPrmtopFile('test.prm7')
inpcrd = AmberInpcrdFile('test.rst7')

# Initialise the molecular system.
system = prmtop.createSystem(nonbondedMethod=PME,
                             nonbondedCutoff=1*nanometer,
                             constraints=HBonds)

# Define the integrator.
integrator = LangevinMiddleIntegrator(300.0*kelvin,
                                1/picosecond,
                                0.002*picoseconds)

# Set the simulation platform.
platform = Platform.getPlatformByName('CPU')
properties = {}

# Initialise and configure the simulation object.
simulation = Simulation(prmtop.topology,
                        system,
                        integrator,
                        platform,
                        properties)
simulation.context.setPositions(inpcrd.positions)
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

# Setting initial system velocities.
simulation.context.setVelocitiesToTemperature(300.0)

# Add reporters.
simulation.reporters.append(DCDReporter('test.dcd', 100, append=False))
log_file = open('test.log', 'a')
simulation.reporters.append(StateDataReporter(log_file,
                                              100,
                                              step=True,
                                              time=True,
                                              potentialEnergy=True,
                                              kineticEnergy=True,
                                              totalEnergy=True,
                                              temperature=True,
                                              volume=True,
                                              totalSteps=True,
                                              speed=True,
                                              remainingTime=True,
                                              separator=' '))

# Run the simulation.
simulation.step(500)

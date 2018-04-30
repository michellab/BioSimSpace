import BioSimSpace as BSS

import matplotlib.pyplot as plt

# Glob the input files.
files = BSS.IO.glob("amber/ala/*")

# Load the molecular system.
print("\nLoading molecules...")
system = BSS.IO.readMolecules(files)

# Create a minimisation protocol.
protocol = BSS.Protocol.Minimisation(steps=1000)

# Initialise the AMBER process.
print("\nInitialising minimisation process...")
process = BSS.Process.Amber(system, protocol, name="minimise", work_dir="minimise")

# Get the list of auto-generated input files.
filenames = process.inputFiles()
print("\nCreated AMBER input files: %s" % filenames)

# Start the minimisation.
print("\nStarting minimisation...")
process.start()

# Get the minimised molecular structure.
minimised = process.getSystem()

# Write the minimised structure to file.
filenames = BSS.IO.saveMolecules("minimised", minimised, system.fileFormat())
print("\nWritten minimised molecular structure to: %s" % filenames)

# Print final energy and timing information.
print("\nMinimised energy is %.2f kcal/mol." % process.getTotalEnergy())
print("Minimisation took %.2f minutes." % process.runTime())

# Create a short heating equilibration protocol.
protocol = BSS.Protocol.Equilibration(runtime=0.01, temperature_start=0, temperature_end=300)

# Initialise the AMBER process.
print("\nInitialising equilibration process...")
process = BSS.Process.Amber(minimised, protocol, name="equilibrate", work_dir="equilibrate")

# Get the list of auto-generated input files.
filenames = process.inputFiles()
print("\nCreated AMBER input files: %s" % filenames)

# Start the equlibration.
print("\nStarting equlibration...")
process.start()

# Get the minimised molecular structure.
equilibrated = process.getSystem()

# Write the equilibrated structure to file.
filenames = BSS.IO.saveMolecules("equilibrated", equilibrated, system.fileFormat())
print("\nWritten equilibrated molecular structure to: %s" % filenames)

# Print final energy and timing information.
print("\nEquilibrated energy is %.2f kcal/mol." % process.getTotalEnergy())
print("Equilibration took %.2f minutes." % process.runTime())

# Create a production protocol.
protocol = BSS.Protocol.Production(runtime=0.01, restart=True)

# Initialise the AMBER process.
print("\nInitialising production process...")
process = BSS.Process.Amber(equilibrated, protocol, name="production", work_dir="production")

# Get the list of auto-generated input files.
filenames = process.inputFiles()
print("\nCreated AMBER input files: %s" % filenames)

# Start the production run.
print("\nStarting production run...")
process.start()

# Get the final molecular structure.
final = process.getSystem()

# Get an MDTraj trajectory object.
traj = process.getTrajectory()

# Get a selection of trajectory frames as a list of systems.
frames = traj.getFrames([0, 5, 10])

# Write the final structure to file.
filenames = BSS.IO.saveMolecules("final", final, system.fileFormat())
print("\nWritten final molecular structure to: %s" % filenames)

# Print final timing information.
print("\nFinal energy is %.2f kcal/mol." % process.getTotalEnergy())
print("Production run took %.2f minutes." % process.runTime())

# Get a list of the time records and the corresponding total energies.
time = process.getTime(time_series=True)
energy = process.getTotalEnergy(time_series=True)

print("Plotting total energy vs time.")

# Create a plot of the total energy vs time.
fig, ax = plt.subplots()
ax.plot(time, energy, '-bo')

ax.set(xlabel="Time (ns)", ylabel="Total Energy (kcal/mol)")
ax.grid()

fig.savefig("energy.pdf", bbox_inches="tight")
plt.show()

print("\nDone!")

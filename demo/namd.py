import BioSimSpace as BSS

from glob import glob
from sys import exit

try:
    from Sire import try_import
    matplotlib = try_import("matplotlib")
    import matplotlib.pyplot as plt
except ImportError:
    raise ImportError("Matplotlib is not installed. Please install matplotlib in order to use BioSimSpace.")

# Glob the input files.
files = glob("namd/ubiquitin/*")

# Load the molecular system.
print("\nLoading molecules...")
system = BSS.readMolecules(files)

# Create a minimisation protocol.
protocol = BSS.Protocol.Minimisation(steps=1000)

# Initialise the NAMD process.
print("\nInitialising minimisation process...")
proc = BSS.Process.Namd(system, protocol, name="minimise",
        work_dir="minimise", charmm_params=True)

# Get the list of auto-generated input files.
filenames = proc.inputFiles()
print("\nCreated NAMD input files: %s" % filenames)

# Start the minimisation.
print("\nStarting minimisation...")
proc.start()

# Wait for the process to end.
proc.wait()

# Check for errors.
if proc.isError():
    print("The process failed!")
    exit()

# Get the minimised molecular structure.
minimised = proc.getSystem()

# Write the minimised structure to file.
filenames = BSS.saveMolecules("minimised", minimised, system.fileFormat())
print("\nWritten minimised molecular structure to: %s" % filenames)

# Print final energy and timing information.
print("\nMinimised energy is %.2f kcal/mol." % proc.getTotalEnergy())
print("Minimisation took %.2f minutes." % proc.runTime())

# Create a short equilibration protocol.
protocol = BSS.Protocol.Equilibration(runtime=0.01)

# Initialise the NAMD process.
print("\nInitialising equilibration process...")
proc = BSS.Process.Namd(minimised, protocol, name="equilibrate",
        work_dir="equilibrate", charmm_params=True)

# Get the list of auto-generated input files.
filenames = proc.inputFiles()
print("\nCreated NAMD input files: %s" % filenames)

# Start the equlibration.
print("\nStarting equlibration...")
proc.start()

# Wait for the process to end.
proc.wait()

# Check for errors.
if proc.isError():
    print("The process failed!")
    exit()

# Get the minimised molecular structure.
equilibrated = proc.getSystem()

# Write the equilibrated structure to file.
filenames = BSS.saveMolecules("equilibrated", equilibrated, system.fileFormat())
print("\nWritten equilibrated molecular structure to: %s" % filenames)

# Print final energy and timing information.
print("\nEquilibrated energy is %.2f kcal/mol." % proc.getTotalEnergy())
print("Equilibration took %.2f minutes." % proc.runTime())

# Create a production protocol.
protocol = BSS.Protocol.Production(runtime=0.01)

# Initialise the NAMD process.
print("\nInitialising production process...")
proc = BSS.Process.Namd(equilibrated, protocol, name="production",
        work_dir="production", charmm_params=True)

# Get the list of auto-generated input files.
filenames = proc.inputFiles()
print("\nCreated NAMD input files: %s" % filenames)

# Start the production run.
print("\nStarting production run...")
proc.start()

# Wait for the process to end.
proc.wait()

# Check for errors.
if proc.isError():
    print("The process failed!")
    exit()

# Get the final molecular structure.
final = proc.getSystem()

# Write the final structure to file.
filenames = BSS.saveMolecules("final", final, system.fileFormat())
print("\nWritten final molecular structure to: %s" % filenames)

# Print final timing information.
print("\nProduction run took %.2f minutes." % proc.runTime())

# Get a list of the time-steps and the corresponding total energies.
time = proc.getTime(time_series=True)
energy = proc.getTotalEnergy(time_series=True)

print("Plotting total energy vs time.")

# Create a plot of the total energy vs time step.
fig, ax = plt.subplots()
ax.plot(time, energy, '-bo')

ax.set(xlabel="Time (ns)", ylabel="Total Energy (kcal/mol)")
ax.grid()

fig.savefig("energy.pdf", bbox_inches="tight")
plt.show()

print("\nDone!")

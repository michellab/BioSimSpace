import BioSimSpace as BSS

from sys import exit

# Load the molecular system.
print("\nLoading molecules...")
system = BSS.readMolecules( ["namd/alanin/alanin.psf",
                             "namd/alanin/alanin.pdb",
                             "namd/alanin/alanin.params"] )

# Create a minimisation protocol.
protocol = BSS.Protocol.Minimisation(steps=1000)

# Initialise the NAMD process.
print("\nInitialising minimisation process...")
proc = BSS.Sample.NamdProcess(system, protocol, name="minimise",
        work_dir="minimise", charmm_params=False)

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

# Create a default equilibration protocol.
protocol = BSS.Protocol.Equilibration()

# Initialise the NAMD process.
print("\nInitialising equilibration process...")
proc = BSS.Sample.NamdProcess(minimised, protocol, name="equilibrate",
        work_dir="equilibrate", charmm_params=False)

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

# Create a default production protocol.
protocol = BSS.Protocol.Production()

# Initialise the NAMD process.
print("\nInitialising production process...")
proc = BSS.Sample.NamdProcess(equilibrated, protocol, name="production",
        work_dir="production", charmm_params=False)

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

# Print final energy and timing information.
print("\nFinal energy is %.2f kcal/mol." % proc.getTotalEnergy())
print("Production run took %.2f minutes." % proc.runTime())

print("\nDone!")

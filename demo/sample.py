import BioSimSpace as BSS

# Load the molecular system.
print("\nLoading molecules...")
system = BSS.readMolecules( ["namd/tiny/tiny.psf",
                             "namd/tiny/tiny.pdb",
                             "namd/tiny/par_all22_prot.inp"] )

# Create a default minimisation protocol.
protocol = BSS.Protocol.Minimisation()

# Initialise the NAMD process.
print("\nInitialising NAMD process...")
proc = BSS.Sample.NamdProcess(system, protocol, name="test")

# Get the list of auto-generated input files.
filenames = proc.input_files()
print("\nCreated NAMD input files: %s" % filenames)

# Start the NAMD simulation.
print("\nStarting simulation...")
proc.start()

# Wait for the process to end.
proc.wait()

# Check for errors.
if proc.isError():
    print("The process failed!")

else:
    # Get the minimised molecular structure.
    min_structure = proc.getSystem()

    # Write the minimised structure to file.
    filenames = BSS.saveMolecules("minimised", min_structure, system.fileFormat())
    print("\nWritten minimised structure to: %s" % filenames)

print("\nDone!")

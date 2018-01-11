import BioSimSpace as BSS

# Load the molecular system.
print("\nLoading molecules...")
system = BSS.readMolecules( ["namd/alanin/alanin.psf",
                             "namd/alanin/alanin.pdb",
                             "namd/alanin/alanin.params"] )

# Create a default minimisation protocol.
protocol = BSS.Protocol.Equilibration()

# Initialise the NAMD process.
print("\nInitialising NAMD process...")
proc = BSS.Sample.NamdProcess(system, protocol, name="test", work_dir="tmp", charmm_params=False)

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
    # Get the final molecular structure.
    new_system = proc.getSystem()

    # Write the minimised structure to file.
    filenames = BSS.saveMolecules("final", new_system, system.fileFormat())
    print("\nWritten final molecular system to: %s" % filenames)
    print("\nFinal energy is %.2f kcal/mol." % proc.getTotalEnergy())

print("\nDone!")

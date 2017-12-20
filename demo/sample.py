import BioSimSpace as BSS

# Load the molecular system.
print("\nLoading molecules...")
system = BSS.readMolecules( ["namd/tiny/tiny.psf", \
                             "namd/tiny/tiny.pdb", \
                             "namd/tiny/par_all22_prot.inp"] )
# Initialise the NAMD process.
print("\nInitialising NAMD process...")
proc = BSS.Sample.NamdProcess(system, BSS.Protocol.NamdProtocol, name="test")

# Get the list of auto-generated input files.
filenames = proc.input_files()
print("\nCreated NAMD input files: %s" % filenames)

# Start the NAMD simulation.
print("\nStarting simulation...")
proc.start()

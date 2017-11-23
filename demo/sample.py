import BioSimSpace as BSS

# Load the molecular system.
print("\nLoading molecules...")
system = BSS.readMolecules( ["namd/tiny/tiny.psf", \
                             "namd/tiny/tiny.pdb", \
                             "namd/tiny/par_all22_prot.inp"] )

# Initialise the NAMD process.
print("\nInitialising NAMD process...")
proc = BSS.Sample.NamdProcess(system, BSS.Protocol.NamdProtocol, "test")

# Setup the simulation workspace.
print("\nSetting up simulation...")
filenames = proc.setup()
print("\nCreated NAMD input files: %s" % filenames)

"""
@package biosimspace
@author Lester Hedges
@brief A class for running simulations using NAMD.
"""

import Sire.Base
import Sire.IO

import tempfile

class NamdProcess(Sire.Base.Process):
    """ A class for running simulations using NAMD. """

    def __init__(self, system, protocol, name = "namd"):
        """ Constructor. """

        # Copy the passed system, protocol, and process name.
        self.system = system
        self.protocol = protocol
        self.name = name

        # Create a temporary working directory and store the directory name.
        self.tmp_dir = tempfile.TemporaryDirectory()

        # The names of the input files.
        self.namd_file = "%s/%s.namd" % (self.tmp_dir.name, name)
        self.psf_file = "%s/%s.psf" % (self.tmp_dir.name, name)
        self.pdb_file = "%s/%s.pdb" % (self.tmp_dir.name, name)
        self.param_file = "%s/%s.params" % (self.tmp_dir.name, name)

        # Create the list of input files.
        self.input_files = [self.namd_file, self.psf_file, self.pdb_file, self.param_file]

    def setup(self):
        """ Setup the input files and working directory ready for simulation. """

        # Create the input files...

        # PSF and parameter files.
        psf = Sire.IO.CharmmPSF(self.system)
        psf.writeToFile(self.psf_file)

        # PDB file.
        pdb = Sire.IO.CharmmPSF(self.system)
        pdb.writeToFile(self.pdb_file)

        # NAMD requires donor and acceptor record entries in the PSF file.
        # We check that these are present and append blank records if not.
        has_donors = False
        has_acceptors = False

        # Open the PSF file for reading.
        with open(self.psf_file) as f:

            # Read the next line.
            line = f.readline()

            # There are donor records.
            if "!NDON" in line:
                has_donors = True

            # There are acceptor records
            elif "NACC" in line:
                has_acceptors = True

        # Append empty donor record.
        if not has_donors:
            f = open(self.psf_file, "a")
            f.write("\n%8d !NDON: donors\n" % 0)
            f.close()

        # Append empty acceptor record.
        if not has_acceptors:
            f = open(self.psf_file, "a")
            f.write("\n%8d !NACC: acceptors\n" % 0)
            f.close()

        # Return the list of input files.
        return self.input_files

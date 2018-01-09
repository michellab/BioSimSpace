"""
@package biosimspace
@author  Lester Hedges
@brief   A class for running simulations using NAMD.
"""

import Sire.Base
import Sire.IO

from . import process
from ..Protocol.protocol_type import ProtocolType

from contextlib import redirect_stdout
import __main__ as main
import os

class NamdProcess(process.Process):
    """A class for running simulations using NAMD."""

    def __init__(self, system, protocol, exe=None, name="namd", work_dir=None):
        """Constructor.

           Keyword arguments:

           system   -- The molecular system.
           protocol -- The protocol for the NAMD process.
           exe      -- The full path to the NAMD executable.
           name     -- The name of the process.
           work_dir -- The working directory for the process.
        """

        # Call the base class constructor.
        super().__init__(system, protocol, name, work_dir)

        # If the path to the executable wasn't specified, then search
        # for it in $PATH.
        if exe is None:
            self._exe = Sire.Base.findExe("namd2").absoluteFilePath()
        else:
            self._exe = exe

        # The names of the input files.
        self._namd_file = "%s/%s.namd" % (self._work_dir, name)
        self._psf_file = "%s/%s.psf" % (self._work_dir, name)
        self._pdb_file = "%s/%s.pdb" % (self._work_dir, name)
        self._param_file = "%s/%s.params" % (self._work_dir, name)
        self._velocity_file = None

        # Create the list of input files.
        self._input_files = [self._namd_file, self._psf_file, self._pdb_file, self._param_file]

        # Determine the box size and origin from the atomic coordinates.
        self._box_size, self._box_origin = process._compute_box_size(self._system)

        # Now set up the working directory for the process.
        self._setup()

    def _setup(self):
        """Setup the input files and working directory ready for simulation."""

        # Create the input files...

        # PSF and parameter files.
        psf = Sire.IO.CharmmPSF(self._system)
        psf.writeToFile(self._psf_file)

        # PDB file.
        pdb = Sire.IO.PDB2(self._system)
        pdb.writeToFile(self._pdb_file)

        # Try to write a PDB "velocity" restart file.
        # The file will only be generated if all atoms in self._system have
        # a "velocity" property.

        # First generate a name for the velocity file.
        velocity_file = os.path.splitext(self._pdb_file)[0] + ".vel"

        # Write the velocity file.
        has_velocities = pdb.writeVelocityFile(velocity_file)

        # If a file was written, store the name of the file and update the
        # list of input files.
        if has_velocities:
            self._velocity_file = velocity_file
            self._input_files.append(self._velocity_file)

        # NAMD requires donor, acceptor, and non-bonded exclusion
        # record entries in the PSF file. We check that these
        # are present and append blank records if not.
        has_donors = False
        has_acceptors = False
        has_non_bonded = False

        # Open the PSF file for reading.
        with open(self._psf_file) as f:

            # Read the next line.
            line = f.readline()

            # There are donor records.
            if "!NDON" in line:
                has_donors = True

            # There are acceptor records.
            elif "NACC" in line:
                has_acceptors = True

            # There are non-bonded exclusion records.
            elif "NNB" in line:
                has_non_bonded = True

        # Append empty donor record.
        if not has_donors:
            f = open(self._psf_file, "a")
            f.write("\n%8d !NDON: donors\n" % 0)
            f.close()

        # Append empty acceptor record.
        if not has_acceptors:
            f = open(self._psf_file, "a")
            f.write("\n%8d !NACC: acceptors\n" % 0)
            f.close()

        # Append empty non-bonded exclusion record.
        if not has_acceptors:
            f = open(self._psf_file, "a")
            f.write("\n%8d !NNB: excluded\n" % 0)
            f.close()

        # Generate the NAMD configuration file.
        self._generate_config_file()

        # Return the list of input files.
        return self._input_files

    def _generate_config_file(self):
        """Generate a NAMD configuration file."""

        # Open the configuration file for writing.
        f = open(self._namd_file, "w")

        # Write generic configuration variables.

        # Topology.
        f.write("structure           %s.psf\n" % self._name)
        f.write("coordinates         %s.pdb\n\n" % self._name)

        # Velocities.
        if not self._velocity_file is None:
            f.write("velocities          %s.vel\n\n" % self._name)

        # Parameters.
        f.write("paraTypeCharmm      on\n")
        f.write("parameters          %s.params\n\n" % self._name)

        # Non-bonded potential parameters.
        f.write("exclude             scaled1-4\n")
        f.write("switching           on\n")
        f.write("switchdist          10.\n")
        f.write("cutoff              12.\n\n")

        # Output file parameters.
        f.write("outputName          %s_out\n" % self._name)
        f.write("binaryOutput        no\n\n")

        # Printing frequency.
        f.write("outputEnergies      100\n")
        f.write("outputTiming        1000\n\n")

        # Add configuration variables for a minimisation simulation.
        if self._protocol.type() == ProtocolType.MINIMISATION:
            f.write("temperature         %s\n\n" % self._protocol.temperature)
            f.write("minimize            %s\n" % self._protocol.steps)

        # Close the configuration file.
        f.close()

    def input_files(self):
        """Return the list of input files."""
        return self._input_files

    def start(self):
        """Start the NAMD simulation."""

        # Store the current working directory.
        dir = os.getcwd()

        # Change to the working directory for the process.
        # This avoid problems with relative paths.
        os.chdir(self._work_dir)

        # Start the simulation.
        self._process = Sire.Base.Process.run(self._exe, "%s.namd" % self._name,
                                                         "%s.out"  % self._name,
                                                         "%s.err"  % self._name)

        # Change back to the original working directory.
        os.chdir(dir)

    def getSystem(self):
        """Get the latest molecular configuration as a Sire system."""

        # Read the PDB coordinate file and construct a parameterised molecular
        # system using the original PSF and param files.

        # List of files.
        files = [ "%s/%s_out.coor" % (self._work_dir, self._name),
                  self._psf_file, self._pdb_file, self._param_file ]

        # Create and return the molecular system.
        return Sire.IO.MoleculeParser.read(files)

"""
@package biosimspace
@author  Lester Hedges
@brief   A class for running simulations using AMBER.
"""

from Sire.Base import findExe, Process
from Sire.IO import AmberPrm, AmberRst7, MoleculeParser

from . import process
from ..Protocol.protocol import Protocol, ProtocolType

from math import ceil, floor
from os import path
from timeit import default_timer as timer
from warnings import warn

import __main__ as main
import os

try:
    from Sire import try_import
    pygtail = try_import("pygtail")
except ImportError:
    raise ImportError("Pygtail is not installed. Please install pygtail in order to use BioSimSpace.")

class AmberProcess(process.Process):
    """A class for running simulations using AMBER."""

    def __init__(self, system, protocol, exe=None, name="amber",
            work_dir=None, seed=None):
        """Constructor.

           Keyword arguments:

           system        -- The molecular system.
           protocol      -- The protocol for the AMBER process.
           exe           -- The full path to the AMBER executable.
           name          -- The name of the process.
           work_dir      -- The working directory for the process.
           seed          -- A random number seed.
        """

        # Call the base class constructor.
        super().__init__(system, protocol, name, work_dir)

        # If the path to the executable wasn't specified, then search
        # for it in $PATH.
        if exe is None:
            self._exe = findExe("sander").absoluteFilePath()

        else:
            # Make sure executable exists.
            if path.isfile(exe):
                self._exe = protocol
            else:
                raise IOError(('AMBER executable doesn\'t exist: "{x}"').format(x=exe))

        # Initialise the energy dictionary and header.
        self._stdout_dict = process._MDict()
        self._stdout_title = None

        # The names of the input files.
        self._crd_file = "%s/%s.crd" % (self._work_dir, name)
        self._top_file = "%s/%s.top" % (self._work_dir, name)

        # Set the path for the AMBER configuration file.
        # The 'protocol' argument may contain the path to a custom file.

        # Set the config file name.
        if not self._is_custom:
            self._config_file = "%s/%s.amber" % (self._work_dir, name)

        # The user has supplied a custom config file.
        else:
            # Make sure the file exists.
            if path.isfile(protocol):
                self._config_file = protocol
            else:
                raise IOError(('AMBER configuration file doesn\'t exist: "{x}"').format(x=config_file))

        # Create the list of input files.
        self._input_files = [self._config_file, self._crd_file, self._top_file]

        # Now set up the working directory for the process.
        self._setup()

    def _setup(self):
        """Setup the input files and working directory ready for simulation."""

        # Create the input files...

        # CRD file (coordinates).
        crd = AmberRst7(self._system)
        crd.writeToFile(self._crd_file)

        # TOP file (topology).
        top = AmberPrm(self._system)
        top.writeToFile(self._top_file)

        # Return the list of input files.
        return self._input_files

    def _generate_config_file(self):
        """Generate an AMBER configuration file."""

    def start(self):
        """Start the AMBER simulation."""

    def getSystem(self):
        """Get the latest molecular configuration as a Sire system."""

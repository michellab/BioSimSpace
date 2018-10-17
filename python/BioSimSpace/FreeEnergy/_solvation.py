######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2018
#
# Authors: Lester Hedges <lester.hedges@gmail.com>
#
# BioSimSpace is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# BioSimSpace is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with BioSimSpace. If not, see <http://www.gnu.org/licenses/>.
#####################################################################

"""
Functionality for running solvation free energy calculations.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

from ..Gateway import ResourceManager as _ResourceManager
from ..Process import ProcessRunner as _ProcessRunner
from ..Process import Somd as _Somd
from ..Protocol import FreeEnergy as _FreeEnergy
from .._SireWrappers import System as _System

import os as _os
import tempfile as _tempfile

class Solvation():
    """A class for configuring and running solvation free energy simulations."""

    def __init__(self, system, protocol=None, work_dir=None):
        """Constructor.

           Positional arguments
           --------------------

           system : BioSimSpace._SireWrappers.System
               The molecular system.


           Keyword arguments
           -----------------

           protocol : BioSimSpace.Protocol.FreeEnergy
               The simulation protocol.

           work_dir : str
               The working directory for the simulation.
        """

        # Validate the input.

        if type(system) is not _System:
            raise TypeError("'system' must be of type 'BioSimSpace._SireWrappers.System'")
        else:
            # Store a copy of the solvated system.
            self._system_solvated = _System(system._getSireSystem().__deepcopy__())

            # The system must have a single perturbable molecule.
            if system.nPerturbableMolecules() != 1:
                raise ValueError("The system must contain a single perturbable molecule!")

            # The system must be solvated.
            if system.nWaterMolecules() == 0:
                raise ValueError("The system must be solvated!")

            # Create the vacuum system.
            self._system_vacuum = _System(system._getSireSystem().__deepcopy__())
            self._system_vacuum.removeWaterMolecules()

        if protocol is not None:
            if type(protocol) is not _FreeEnergy:
                raise TypeError("'protocol' must be of type 'BioSimSpace.Protocol.FreeEnergy'")
            else:
                self._protocol = protocol
        else:
            # Use a default protocol.
            self._protocol = _FreeEnergy()

        # Create a temporary working directory and store the directory name.
        if work_dir is None:
            self._tmp_dir = _tempfile.TemporaryDirectory()
            self._work_dir = self._tmp_dir.name

        # User specified working directory.
        else:
            self._work_dir = work_dir

            # Create the directory if it doesn't already exist.
            if not _os.path.isdir(work_dir):
                _os.makedirs(work_dir, exist_ok=True)

        # Initalise lists to hold all of the process objects for the free and
        # vacuum legs.
        free = []
        vacuum = []

        # Set the directories for the free and vacuum legs.
        free_dir = "%s/free" % self._work_dir
        vacuum_dir = "%s/vacuum" % self._work_dir

        # Get the lambda values from the protocol.
        lam_vals = self._protocol.getLambdaValues()

        # Loop over all of the lambda values.
        for lam in lam_vals:
            # Update the protocol lambda values.
            self._protocol.setLambdaValues(lam=lam, lam_vals=lam_vals)

            # Create and append the required processes for each leg.
            # Nest the working directories inside self._work_dir.

            free.append(_Somd(self._system_solvated, self._protocol,
                platform="CPU", work_dir="%s/lambda_%s" % (free_dir, lam)))

            vacuum.append(_Somd(self._system_vacuum, self._protocol,
                platform="CPU", work_dir="%s/lambda_%s" % (vacuum_dir, lam)))

        # Initialise the process runner. All processes have already been nested
        # inside the working directory of the Solvation object.
        self._runner = _ProcessRunner(free + vacuum, work_dir=self._work_dir, nest_dirs=False)

    def run(self):
        """Run the solvation free energy simulation."""
        self._runner.startAll()

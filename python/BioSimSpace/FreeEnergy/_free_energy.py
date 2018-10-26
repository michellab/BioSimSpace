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
Base class for free energy simulations.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

from Sire.Base import getBinDir as _getBinDir

from .._Exceptions import MissingSoftwareError as _MissingSoftwareError
from ..Gateway import ResourceManager as _ResourceManager
from ..Process import ProcessRunner as _ProcessRunner
from ..Process import Somd as _Somd
from ..Protocol import FreeEnergy as _FreeEnergy
from .._SireWrappers import System as _System

import BioSimSpace.Units as _Units

import math as _math
import os as _os
import subprocess as _subprocess
import tempfile as _tempfile

class FreeEnergy():
    """Base class for configuring and running free energy simulations."""

    # Check that the analyse_freenrg script exists.
    _analyse_freenrg = _getBinDir() + "/analyse_freenrg"

    if not _os.path.isfile(_analyse_freenrg):
        raise _MissingSoftwareError("'Cannot find free energy analysis script in expected location: '%s'" % _analyse_freenrg)

    def __init__(self, protocol=None, work_dir=None):
        """Constructor.

           Keyword arguments
           -----------------

           protocol : BioSimSpace.Protocol.FreeEnergy
               The simulation protocol.

           work_dir : str
               The working directory for the simulation.
        """

	# Don't allow user to create an instance of this base class.
        if type(self) is FreeEnergy:
            raise Exception("<FreeEnergy> must be subclassed.")

        # Validate the input.

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

    def run(self):
        """Run the simulation."""
        self._runner.startAll()

    def analyse(self):
        """Analyse the solvation free energy data.

           Returns
           -------

           pmf0 : [ ( float, BioSimSpace.Types.Energy, BioSimSpace.Types.Energy ) ]
               The potential of mean force (PMF) for the first leg of the
               simulation. The data is a list of tuples, where each tuple
               contains the lambda value, the PMF, and the standard error.

           pmf1 : [ ( float, BioSimSpace.Types.Energy, BioSimSpace.Types.Energy ) ]
               The potential of mean force (PMF) for the second leg of the
               simulation. The data is a list of tuples, where each tuple
               contains the lambda value, the PMF, and the standard error.

           free_energy : ( BioSimSpace.Types.Energy, BioSimSpace.Types.Energy )
               The free energy difference and its associated error.
        """

        # Create the commands for the two legs.
        command0 = "%s mbar -i %s/lambda_*/simfile.dat -o %s/mbar_leg0.txt" % (self._analyse_freenrg, self._dir0, self._work_dir)
        command1 = "%s mbar -i %s/lambda_*/simfile.dat -o %s/mbar_leg1.txt" % (self._analyse_freenrg, self._dir1, self._work_dir)

        # Run the first command.
        proc = _subprocess.run(command0, shell=True, stdout=_subprocess.PIPE, stderr=_subprocess.PIPE)
        if proc.returncode != 0:
            return None

        # Run the second command.
        proc = _subprocess.run(command1, shell=True, stdout=_subprocess.PIPE, stderr=_subprocess.PIPE)
        if proc.returncode != 0:
            return None

        # Initialise lists to hold the data from each leg.
        leg0 = []
        leg1 = []

        # Extract the data from the output files.

        # First leg.
        with open("%s/mbar_leg0.txt" % self._work_dir) as file:

            # Read all of the lines into a list.
            lines = []
            for line in file:
                lines.append(line.rstrip())

            # Find the MBAR data.
            for x, line in enumerate(lines):
                if "PMF from MBAR" in line:
                    # Increment the line index.
                    x += 1

                    # Loop until we hit the next comment.
                    while lines[x][0] != "#":
                        # Split the line.
                        data = lines[x].split()

                        # Append the data.
                        leg0.append((float(data[0]),
                                     float(data[1]) * _Units.Energy.kcal_per_mol,
                                     float(data[2]) * _Units.Energy.kcal_per_mol))

                        # Increment the line index.
                        x += 1

                    break

        # Second leg.
        with open("%s/mbar_leg1.txt" % self._work_dir) as file:

            # Read all of the lines into a list.
            lines = []
            for line in file:
                lines.append(line.rstrip())

            # Find the MBAR data.
            for x, line in enumerate(lines):
                if "PMF from MBAR" in line:
                    # Increment the line index.
                    x += 1

                    # Loop until we hit the next comment.
                    while lines[x][0] != "#":
                        # Split the line.
                        data = lines[x].split()

                        # Append the data.
                        leg1.append((float(data[0]),
                                     float(data[1]) * _Units.Energy.kcal_per_mol,
                                     float(data[2]) * _Units.Energy.kcal_per_mol))

                        # Increment the line index.
                        x += 1

                    break

        # Work out the difference in free energy.
        free_energy = (leg0[-1][1] - leg0[0][1]) - (leg1[-1][1] - leg1[0][1])

        # Propagate the errors. (These add in quadrature.)

        # First leg.
        error0 = _math.sqrt((leg0[-1][2].magnitude() * leg0[-1][2].magnitude()) +
                            (leg0[0][2].magnitude()*leg0[0][2].magnitude()))

        # Second leg.
        error1 = _math.sqrt((leg1[-1][2].magnitude() * leg1[-1][2].magnitude()) +
                            (leg1[0][2].magnitude() * leg1[0][2].magnitude()))

        # Free energy difference.
        error = _math.sqrt((error0 * error0) + (error1 * error1)) * _Units.Energy.kcal_per_mol

        # Bundle the free energy and its associated error.
        free_energy = (free_energy, error)

        return (leg0, leg1, free_energy)

    def _initialise_runner(self, system0, system1):
        """Internral helper function to initialise the process runner.


           Positional arguments
           --------------------

           system0 : BioSimSpace._SireWrappers.System
               The system for the first free energy leg.

           system1 : BioSimSpace._SireWrappers.System
               The system for the second free energy leg.
        """

        if type(system0) is not _System:
            raise TypeError("'system0' must be of type 'BioSimSpace._SireWrappers.System'")

        if type(system1) is not _System:
            raise TypeError("'system1' must be of type 'BioSimSpace._SireWrappers.System'")

        # Initialise lists to store the processes for each leg.
        leg0 = []
        leg1 = []

        # Get the simulation type.
        sim_type = self.__class__.__name__

        # Store the working directories for the legs.

        if sim_type == "Solvation":
            self._dir0 = "%s/free" % self._work_dir
            self._dir1 = "%s/vacuum" % self._work_dir
        elif sim_type == "Binding":
            self._dir0 = "%s/bound" % self._work_dir
            self._dir1 = "%s/free" % self._work_dir
        else:
            raise TypeError("Unsupported FreeEnergy simulation: '%s'" % sim_type)

        # Get the lambda values from the protocol.
        lam_vals = self._protocol.getLambdaValues()

        # Loop over all of the lambda values.
        for lam in lam_vals:
            # Update the protocol lambda values.
            self._protocol.setLambdaValues(lam=lam, lam_vals=lam_vals)

            # Create and append the required processes for each leg.
            # Nest the working directories inside self._work_dir.

            leg0.append(_Somd(system0, self._protocol,
                platform="CUDA", work_dir="%s/lambda_%s" % (self._dir0, lam)))

            leg1.append(_Somd(system1, self._protocol,
                platform="CUDA", work_dir="%s/lambda_%s" % (self._dir1, lam)))

        # Initialise the process runner. All processes have already been nested
        # inside the working directory of the Solvation object.
        self._runner = _ProcessRunner(leg0 + leg1, work_dir=self._work_dir, nest_dirs=False)

######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2019
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
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["FreeEnergy"]

from collections import OrderedDict as _OrderedDict

import math as _math
import sys as _sys
import os as _os
import subprocess as _subprocess
import tempfile as _tempfile
import warnings as _warnings

from Sire.Base import getBinDir as _getBinDir
from Sire.Base import getShareDir as _getShareDir

from Sire import IO as _SireIO
from Sire import Mol as _SireMol

from BioSimSpace import _gmx_exe
from BioSimSpace._Exceptions import MissingSoftwareError as _MissingSoftwareError
from BioSimSpace._SireWrappers import Molecules as _Molecules
from BioSimSpace._SireWrappers import System as _System
from BioSimSpace import Process as _Process
from BioSimSpace import Protocol as _Protocol
from BioSimSpace import Units as _Units

class FreeEnergy():
    """Base class for configuring and running free energy simulations."""

    # Check that the analyse_freenrg script exists.
    if _sys.platform != "win32":
        _analyse_freenrg = _os.path.join(_getBinDir(), "analyse_freenrg")
    else:
        _analyse_freenrg = _os.path.join(_os.path.normpath(_getShareDir()), "scripts", "analyse_freenrg.py")
    if not _os.path.isfile(_analyse_freenrg):
        raise _MissingSoftwareError("Cannot find free energy analysis script in expected location: '%s'" % _analyse_freenrg)
    if _sys.platform == "win32":
        _analyse_freenrg = "%s %s" % (_os.path.join(_os.path.normpath(_getBinDir()), "sire_python.exe"), _analyse_freenrg)

    # Create a list of supported molecular dynamics engines.
    _engines = ["GROMACS", "SOMD"]

    def __init__(self, protocol=None, work_dir=None, engine=None):
        """Constructor.

           Parameters
           ----------

           protocol : :class:`Protocol.FreeEnergy <BioSimSpace.Protocol.FreeEnergy>`
               The simulation protocol.

           work_dir : str
               The working directory for the simulation.

           engine: str
               The molecular dynamics engine used to run the simulation. Available
               options are "GROMACS", or "SOMD". If this argument is omitted then
               BioSimSpace will choose an appropriate engine for you.
        """

	# Don't allow user to create an instance of this base class.
        if type(self) is FreeEnergy:
            raise Exception("<FreeEnergy> must be subclassed.")

        # Validate the input.

        if protocol is not None:
            if type(protocol) is not _Protocol.FreeEnergy:
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

        # Validate the user specified molecular dynamics engine.
        if engine is not None:
            if type(engine) is not str:
                raise Types("'engine' must be of type 'str'.")

            # Strip whitespace from engine and convert to upper case.
            engine = engine.replace(" ", "").upper()

            # Check that the engine is supported.
            if engine not in self._engines:
                raise ValueError("Unsupported molecular dynamics engine '%s'. "
                                 "Supported engines are: %r." % ", ".join(self._engines))

            # Make sure GROMACS is installed if GROMACS engine is selected.
            if engine == "GROMACS" and _gmx_exe is None:
                raise _MissingSoftwareError("Cannot use GROMACS engine as GROMACS is not installed!")
        else:
            # Use SOMD as a default.
            engine = "SOMD"

        # Set the engine.
        self._engine = engine

    def run(self):
        """Run the simulation."""
        self._runner.startAll()

    def _analyse_gromacs(self):
        """Analyse the GROMACS free energy data.

           Returns
           -------

           pmf0 : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
               The potential of mean force (PMF) for the first leg of the
               simulation. The data is a list of tuples, where each tuple
               contains the lambda value, the PMF, and the standard error.

           pmf0 : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
               The potential of mean force (PMF) for the second leg of the
               simulation. The data is a list of tuples, where each tuple
               contains the lambda value, the PMF, and the standard error.

           free_energy : (:class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)
               The free energy difference and its associated error.
        """

        # Create the commands for the two legs.
        command0 = "%s bar -f %s/lambda_*/*.xvg -o %s/bar_leg0.xvg" % (_gmx_exe, self._dir0, self._work_dir)
        command1 = "%s bar -f %s/lambda_*/*.xvg -o %s/bar_leg1.xvg" % (_gmx_exe, self._dir1, self._work_dir)

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
        with open("%s/bar_leg0.xvg" % self._work_dir) as file:

            # Read all of the lines into a list.
            lines = []
            for line in file:
                # Ignore comments and xmgrace directives.
                if line[0] != "#" and line[0] != "@":
                    lines.append(line.rstrip())

            # Store the initial free energy reading.
            leg0.append((0.0,
                         0.0 * _Units.Energy.kcal_per_mol,
                         0.0 * _Units.Energy.kcal_per_mol))

            # Zero the accumulated error.
            total_error = 0

            # Zero the accumulated free energy difference.
            total_freenrg = 0

            # Process the BAR data.
            for x, line in enumerate(lines):
                # Extract the data from the line.
                data = line.split()

                # Update the total free energy difference.
                total_freenrg += float(data[1])

                # Extract the error.
                error = float(data[2])

                # Update the accumulated error.
                total_error = _math.sqrt(total_error*total_error + error*error)

                # Append the data.
                leg0.append(((x + 1) / (len(lines)),
                             (total_freenrg * _Units.Energy.kt).kcal_per_mol(),
                             (total_error * _Units.Energy.kt).kcal_per_mol()))

        # Second leg.
        with open("%s/bar_leg1.xvg" % self._work_dir) as file:

            # Read all of the lines into a list.
            lines = []
            for line in file:
                # Ignore comments and xmgrace directives.
                if line[0] != "#" and line[0] != "@":
                    lines.append(line.rstrip())

            # Store the initial free energy reading.
            leg1.append((0.0,
                         0.0 * _Units.Energy.kcal_per_mol,
                         0.0 * _Units.Energy.kcal_per_mol))

            # Zero the accumulated error.
            total_error = 0

            # Zero the accumulated free energy difference.
            total_freenrg = 0

            # Process the BAR data.
            for x, line in enumerate(lines):
                # Extract the data from the line.
                data = line.split()

                # Update the total free energy difference.
                total_freenrg += float(data[1])

                # Extract the error.
                error = float(data[2])

                # Update the accumulated error.
                total_error = _math.sqrt(total_error*total_error + error*error)

                # Append the data.
                leg1.append(((x + 1) / (len(lines)),
                             (total_freenrg * _Units.Energy.kt).kcal_per_mol(),
                             (total_error * _Units.Energy.kt).kcal_per_mol()))

        # Work out the difference in free energy.
        free_energy = (leg0[-1][1] - leg0[0][1]) - (leg1[-1][1] - leg1[0][1])

        # Propagate the errors. (These add in quadrature.)

        # First leg.
        error0 = _math.sqrt((leg0[-1][2].magnitude() * leg0[-1][2].magnitude()) +
                            (leg0[0][2].magnitude()  * leg0[0][2].magnitude()))

        # Second leg.
        error1 = _math.sqrt((leg1[-1][2].magnitude() * leg1[-1][2].magnitude()) +
                            (leg1[0][2].magnitude()  * leg1[0][2].magnitude()))

        # Free energy difference.
        error = _math.sqrt((error0 * error0) + (error1 * error1)) * _Units.Energy.kcal_per_mol

        # Bundle the free energy and its associated error.
        free_energy = (free_energy, error)

        return (leg0, leg1, free_energy)

    def _analyse_somd(self):
        """Analyse the SOMD free energy data.

           Returns
           -------

           pmf0 : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
               The potential of mean force (PMF) for the first leg of the
               simulation. The data is a list of tuples, where each tuple
               contains the lambda value, the PMF, and the standard error.

           pmf0 : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
               The potential of mean force (PMF) for the second leg of the
               simulation. The data is a list of tuples, where each tuple
               contains the lambda value, the PMF, and the standard error.

           free_energy : (:class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)
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
        """Internal helper function to initialise the process runner.

           Parameters
           ----------

           system0 : :class:`System <BioSimSpace._SireWrappers.System>`
               The system for the first free energy leg.

           system1 : :class:`System <BioSimSpace._SireWrappers.System>`
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

        # Convert to an appropriate AMBER topology. (Required by SOMD.)
        if self._engine == "SOMD":
            # Try to get the water model used to solvate the system.
            try:
                water_model = system0._sire_object.property("water_model").toString()
                waters0 = _SireIO.setAmberWater(system0._sire_object.search("water"), water_model)
                waters1 = _SireIO.setAmberWater(system1._sire_object.search("water"), water_model)

            # If the system wasn't solvated by BioSimSpace, e.g. read from file, then try
            # to guess the water model from the topology.
            except:
                num_point = system0.getWaterMolecules()[0].nAtoms()

                if num_point == 3:
                    # TODO: Assume TIP3P. Not sure how to detect SPC/E.
                    waters0 = _SireIO.setAmberWater(system0._sire_object.search("water"), "TIP3P")
                    waters1 = _SireIO.setAmberWater(system1._sire_object.search("water"), "TIP3P")
                    water_model = "tip3p"
                elif num_point == 4:
                    waters0 = _SireIO.setAmberWater(system0._sire_object.search("water"), "TIP4P")
                    waters1 = _SireIO.setAmberWater(system1._sire_object.search("water"), "TIP4P")
                    water_model = "tip4p"
                elif num_point == 5:
                    waters0 = _SireIO.setAmberWater(system0._sire_object.search("water"), "TIP5P")
                    waters1 = _SireIO.setAmberWater(system1._sire_object.search("water"), "TIP5P")
                    water_model = "tip5p"
                else:
                    raise RuntimeError("Unsupported %d-point water model!" % num_point)

                # Warn the user that we've guessed the water topology.
                _warnings.warn("Guessed water topology: %r" % water_model)

            # Remove the existing water molecules from the systems.
            system0.removeWaterMolecules()
            system1.removeWaterMolecules()

            # Convert the waters to BioSimSpace molecule containers.
            waters0 = _Molecules(waters0.toMolecules())
            waters1 = _Molecules(waters1.toMolecules())

            # Add the updated water topology back into the systems.
            system0.addMolecules(waters0)
            system1.addMolecules(waters1)

        # Get the lambda values from the protocol.
        lam_vals = self._protocol.getLambdaValues()

        # Loop over all of the lambda values.
        for lam in lam_vals:
            # Update the protocol lambda values.
            self._protocol.setLambdaValues(lam=lam, lam_vals=lam_vals)

            # Create and append the required processes for each leg.
            # Nest the working directories inside self._work_dir.

            # SOMD.
            if self._engine == "SOMD":
                # Check for GPU support.
                if "CUDA_VISIBLE_DEVICES" in _os.environ:
                    platform = "CUDA"
                else:
                    platform = "CPU"

                leg0.append(_Process.Somd(system0, self._protocol,
                    platform=platform, work_dir="%s/lambda_%5.4f" % (self._dir0, lam)))

                leg1.append(_Process.Somd(system1, self._protocol,
                    platform=platform, work_dir="%s/lambda_%5.4f" % (self._dir1, lam)))

            # GROMACS.
            elif self._engine == "GROMACS":
                leg0.append(_Process.Gromacs(system0, self._protocol,
                    work_dir="%s/lambda_%5.4f" % (self._dir0, lam)))

                leg1.append(_Process.Gromacs(system1, self._protocol,
                    work_dir="%s/lambda_%5.4f" % (self._dir1, lam)))

        # Initialise the process runner. All processes have already been nested
        # inside the working directory so no need to re-nest.
        self._runner = _Process.ProcessRunner(leg0 + leg1, work_dir=self._work_dir, nest_dirs=False)

    def _update_run_args(self, args):
        """Internal function to update run arguments for all subprocesses.

           Parameters
           ----------

           args : dict, collections.OrderedDict
               A dictionary which contains the new command-line arguments
               for the process executable.
        """

        if type(args) is not dict and type(args) is not _OrderedDict:
            raise TypeError("'args' must be of type 'dict', or 'collections.OrderedDict'")

        for process in self._runner.processes():
            process.setArgs(args)

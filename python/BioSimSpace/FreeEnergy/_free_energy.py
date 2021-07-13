######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2021
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
__email__ = "lester.hedges@gmail.com"

__all__ = ["FreeEnergy", "analyse", "getData"]

from collections import OrderedDict as _OrderedDict
from glob import glob as _glob

import copy as _copy
import math as _math
import sys as _sys
import os as _os
import shutil as _shutil
import subprocess as _subprocess
import tempfile as _tempfile
import warnings as _warnings
import zipfile as _zipfile

from Sire.Base import getBinDir as _getBinDir
from Sire.Base import getShareDir as _getShareDir

from Sire import IO as _SireIO
from Sire import Mol as _SireMol

from BioSimSpace import _gmx_exe
from BioSimSpace import _is_notebook
from BioSimSpace._Exceptions import AnalysisError as _AnalysisError
from BioSimSpace._Exceptions import MissingSoftwareError as _MissingSoftwareError
from BioSimSpace._SireWrappers import Molecules as _Molecules
from BioSimSpace._SireWrappers import System as _System
from BioSimSpace._Utils import cd as _cd
from BioSimSpace import Process as _Process
from BioSimSpace import Protocol as _Protocol
from BioSimSpace import Units as _Units

if _is_notebook:
    from IPython.display import FileLink as _FileLink

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

    def __init__(self, protocol=None, work_dir=None, engine=None,
            setup_only=False, ignore_warnings=False, show_errors=True):
        """Constructor.

           Parameters
           ----------

           protocol : :class:`Protocol.FreeEnergy <BioSimSpace.Protocol.FreeEnergy>`, \
                     [:class:`Protocol.FreeEnergy <BioSimSpace.Protocol.FreeEnergy>`,
                      :class:`Protocol.FreeEnergy <BioSimSpace.Protocol.FreeEnergy>`]
               The simulation protocol. If one obect is passed, then the same
               protocol will be used for both legs of the simulation. Passing
               two objects enables a different protocol for each leg, e.g. a
               different lambda schedule, or run time.

           work_dir : str
               The working directory for the simulation.

           engine: str
               The molecular dynamics engine used to run the simulation. Available
               options are "GROMACS", or "SOMD". If this argument is omitted then
               BioSimSpace will choose an appropriate engine for you.

           setup_only: bool
               Whether to only support simulation setup. If True, then no
               simulation processes objects will be created, only the directory
               hierarchy and input files to run a simulation externally. This
               can be useful when you don't intend to use BioSimSpace to run
               the simulation. Note that a 'work_dir' must also be specified.

           ignore_warnings : bool
               Whether to ignore warnings when generating the binary run file.
               This option is specific to GROMACS and will be ignored when a
               different molecular dynamics engine is chosen.

           show_errors : bool
               Whether to show warning/error messages when generating the binary
               run file. This option is specific to GROMACS and will be ignored
               when a different molecular dynamics engine is chosen.
        """

        # Don't allow user to create an instance of this base class.
        if type(self) is FreeEnergy:
            raise Exception("<FreeEnergy> must be subclassed.")

        # Flag that this is a dual leg simulation (default).
        self._is_dual = True

        # Validate the input.

        if protocol is not None:
            # Convert tuple to list.
            if type(protocol) is tuple:
                protocol = list(protocol)

            if type(protocol) is _Protocol.FreeEnergy:
                self._protocol0 = protocol
                self._protocol1 = protocol
            elif type(protocol) is list:
                if not(all(isinstance(x, _Protocol.FreeEnergy) for x in protocol)):
                    raise TypeError("'protocol' must be a list of 'BioSimSpace.Protocol.FreeEnergy' types")
                # Same protocol for both legs.
                if len(protocol) == 1:
                    self._protocol0 = protocol[0]
                    self._protocol1 = protocol[0]
                elif len(protocol) == 2:
                    self._protocol0 = protocol[0]
                    self._protocol1 = protocol[1]
                else:
                    raise TypeError("'protocol' can only contain a maxmimum of two items, one for each leg.")
            else:
                raise TypeError("'protocol' must be of type 'BioSimSpace.Protocol.FreeEnergy' "
                                "or a list of 'BioSimSpace.Protocol.FreeEnergy' types." )
        else:
            # Use a default protocol.
            self._protocol0 = _FreeEnergy()
            self._protocol1 = _FreeEnergy()

        if type(setup_only) is not bool:
            raise TypeError("'setup_only' must be of type 'bool'.")
        else:
            self._setup_only = setup_only

        # Create a temporary working directory and store the directory name.
        if work_dir is None:
            if setup_only:
                raise ValueError("A 'work_dir' must be specified when 'setup_only' is True!")
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
            if engine == "GROMACS":
              if _gmx_exe is None:
                raise _MissingSoftwareError("Cannot use GROMACS engine as GROMACS is not installed!")

              if self._protocol0.getPerturbationType() != "full" or \
                 self._protocol1.getPerturbationType() != "full":
                    raise NotImplementedError("GROMACS currently only supports the 'full' perturbation "
                                              "type. Please use engine='SOMD' when running multistep "
                                              "perturbation types.")
        else:
            # Use SOMD as a default.
            engine = "SOMD"

        # Set the engine.
        self._engine = engine

        if type(ignore_warnings) is not bool:
            raise ValueError("'ignore_warnings' must be of type 'bool.")
        self._ignore_warnings = ignore_warnings

        if type(show_errors) is not bool:
            raise ValueError("'show_errors' must be of type 'bool.")
        self._show_errors = show_errors

    def run(self, serial=True):
        """Run the simulation.

           Parameters
           ----------

           serial : bool
               Whether to run the individual processes for the lambda windows
               in serial.
        """
        if type(serial) is not bool:
            raise TypeError("'serial' must be of type 'bool'.")

        if self._setup_only:
            _warnings.warn("No proceses exist! Object created in 'setup_only' mode.")
        else:
            self._runner.startAll(serial=serial)

    def wait(self):
        """Wait for the simulation to finish."""
        if self._setup_only:
            _warnings.warn("No proceses exist! Object created in 'setup_only' mode.")
        else:
            self._runner.wait()

    def workDir(self):
        """Return the working directory.

           Returns
           -------

           work_dir : str
               The path of the working directory.
        """
        return self._work_dir

    @classmethod
    def getData(cls, name="data", file_link=False, work_dir=None):
        """Return a link to a zip file containing the data files required for
           post-simulation analysis.

           Parameters
           ----------

           name : str
               The name of the zip file.

           file_link : bool
               Whether to return a FileLink when working in Jupyter.

           work_dir : str
               The working directory for the simulation.

           Returns
           -------

           ouput : str, IPython.display.FileLink
               A path, or file link, to an archive of the process input.
        """

        if work_dir is None and cls._work_dir is None:
            raise ValueError("'work_dir' must be set!")
        elif work_dir is None:
            work_dir = cls._work_dir
        else:
            if type(work_dir) is not str:
                raise TypeError("'work_dir' must be of type 'str'.")
            if not _os.path.isdir(work_dir):
                raise ValueError("'work_dir' doesn't exist!")

        if type(name) is not str:
            raise TypeError("'name' must be of type 'str'")

        # Generate the zip file name.
        zipname = "%s.zip" % name

        # Get the current working directory.
        cwd = _os.getcwd()

        # Change into the working directory.
        with _cd(work_dir):
            # Glob all of the analysis files.

            # First try SOMD data.
            files = _glob("*/*/gradients.dat")

            if len(files) == 0:
                files = _glob("*/*/gromacs.xvg")

                if len(files) == 0:
                    raise ValueError(f"Couldn't find any analysis files in '{work_dir}'")

            # Write to the zip file.
            with _zipfile.ZipFile(cwd + f"/{zipname}", "w") as zip:
                for file in files:
                    zip.write(file)

        # Return a link to the archive.
        if _is_notebook:
            if file_link:
                # Create a FileLink to the archive.
                f_link = _FileLink(zipname)

                # Set the download attribute so that JupyterLab doesn't try to open the file.
                f_link.html_link_str = f"<a href='%s' target='_blank' download='{zipname}'>%s</a>"

                # Return a link to the archive.
                return f_link
            else:
                return zipname
        # Return the path to the archive.
        else:
            return zipname

    @staticmethod
    def analyse(work_dir, simulation_type=None):
        """Analyse existing free-energy data from a simulation working directory.

           For the return values, leg0 refers to the simulation leg with the
           larger number of molecules, i.e. the "bound" leg for a "binding"
           free-energy simulation, or the "free" leg for a "solvation"
           free-energy simulation.

           Parameters
           ----------

           work_dir : str
               The working directory for the simulation.

           simulation_type : str
               The type of free-energy perturbation simulation_type. Options are:
               "solvation", or "binding". This allows use of a single "work_dir"
               for combined "solvation" and "binding" free-energy simulations.
               If None, then BioSimSpace will attempt to figure out the simulation
               type from the directory structure within "work_dir".

           Returns
           -------

           pmf0 : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
               The potential of mean force (PMF) for the first leg of the
               simulation. The data is a list of tuples, where each tuple
               contains the lambda value, the PMF, and the standard error.

           pmf1 : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
               The potential of mean force (PMF) for the second leg of the
               simulation. The data is a list of tuples, where each tuple
               contains the lambda value, the PMF, and the standard error.

           free_energy : (:class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)
               The free energy difference and its associated error.

           overlap0 : [ [ float, float, ... ] ]
               The overlap matrix. This gives the overlap between each window
               of the first leg. This parameter is only computed for the SOMD
               engine and will be None when GROMACS is used.

           overlap1 : [ [ float, float, ... ] ]
               The overlap matrix. This gives the overlap between each window
               of the second leg. This parameter is only computed for the SOMD
               engine and will be None when GROMACS is used.
        """

        if type(work_dir) is not str:
            raise TypeError("'work_dir' must be of type 'str'.")
        if not _os.path.isdir(work_dir):
            raise ValueError("'work_dir' doesn't exist!")

        if simulation_type is not None:
            if type(simulation_type) is not str:
                raise TypeError("'simulation_type' must be of type 'str'.")
            # Strip whitespace and convert to lower case.
            simulation_type = simulation_type.lower().replace(" ", "")
            if simulation_type not in ["solvation", "binding"]:
                raise ValueError("'simulation_type' must be either 'solvation' or 'binding'.")

        # Whether this is a dual-leg simulation.
        is_dual = False

        # First work out whether this is a binding or solvation simulation.

        if simulation_type is None:
            # Binding.
            if _os.path.isdir(work_dir + "/bound"):
                dir0 = work_dir + "/bound"
                dir1 = work_dir + "/free"
                if _os.path.isdir(dir1):
                    is_dual = True

            # Solvation.
            elif _os.path.isdir(work_dir + "/free"):
                dir0 = work_dir + "/free"
                dir1 = work_dir + "/vacuum"
                if _os.path.isdir(dir1):
                    is_dual = True

            # Invalid directory structure.
            else:
                msg = (f"Could not find '{work_dir}/bound' or "
                       f"'{work_dir}/free'?")
                raise ValueError(msg)

        else:
            if simulation_type == "binding":
                if _os.path.isdir(work_dir + "/bound"):
                    dir0 = work_dir + "/bound"
                    dir1 = work_dir + "/free"
                    if _os.path.isdir(dir1):
                        is_dual = True

                # Invalid directory structure.
                else:
                    msg = (f"Could not find '{work_dir}/bound'")
                    raise ValueError(msg)

            elif simulation_type == "solvation":
                if _os.path.isdir(work_dir + "/free"):
                    dir0 = work_dir + "/free"
                    dir1 = work_dir + "/vacuum"
                    if _os.path.isdir(dir1):
                        is_dual = True

                # Invalid directory structure.
                else:
                    msg = (f"Could not find '{work_dir}/free'")
                    raise ValueError(msg)

        # First test for SOMD files.
        data = _glob(dir0 + "/lambda_*/gradients.dat")

        # SOMD.
        if len(data) > 0:
            return FreeEnergy._analyse_somd(work_dir, dir0, dir1, is_dual)

        # Now check for GROMACS output.
        else:
            data = _glob(dir0 + "/lambda_*/gromacs.xvg")
            if len(data) == 0:
                raise ValueError("Couldn't find any SOMD or GROMACS free-energy output?")
            return FreeEnergy._analyse_gromacs(work_dir, dir0, dir1, is_dual)

    @staticmethod
    def _analyse_gromacs(work_dir=None, dir0=None, dir1=None, is_dual=True):
        """Analyse the GROMACS free energy data.

           Parameters
           ----------

           work_dir : str
               The path to the working directory.

           dir0 : str
               The path to the directory of the first leg.

           dir1 : str
               The path to the directory of the second leg.

           is_dual : bool
               Whether this is a dual-leg simulation.

           Returns
           -------

           pmf0 : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
               The potential of mean force (PMF) for the first leg of the
               simulation. The data is a list of tuples, where each tuple
               contains the lambda value, the PMF, and the standard error.

           pmf1 : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
               The potential of mean force (PMF) for the second leg of the
               simulation. The data is a list of tuples, where each tuple
               contains the lambda value, the PMF, and the standard error.

           free_energy : (:class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)
               The free energy difference and its associated error.
        """

        if type(work_dir) is not str:
            raise TypeError("'work_dir' must be of type 'str'.")
        if not _os.path.isdir(work_dir):
            raise ValueError("'work_dir' doesn't exist!")

        if type(is_dual) is not bool:
            raise TypeError("'is_dual' must be of type 'bool'.")

        if type(dir0) is not str:
            raise TypeError("'dir0' must be of type 'str'.")
        if not _os.path.isdir(dir0):
            raise ValueError("'dir0' doesn't exist!")

        if is_dual:
            if type(dir1) is not str:
                raise TypeError("'dir1' must be of type 'str'.")
            if not _os.path.isdir(dir1):
                raise ValueError("'dir1' doesn't exist!")

        # Create the commands for the two legs.
        command0 = "%s bar -f %s/lambda_*/*.xvg -o %s/bar_leg0.xvg" % (_gmx_exe, dir0, work_dir)
        command1 = "%s bar -f %s/lambda_*/*.xvg -o %s/bar_leg1.xvg" % (_gmx_exe, dir1, work_dir)

        # Run the first command.
        proc = _subprocess.run(command0, shell=True, stdout=_subprocess.PIPE, stderr=_subprocess.PIPE)
        if proc.returncode != 0:
            raise _AnalysisError("GROMACS free-energy analysis failed!")

        # Run the second command.
        if is_dual:
            proc = _subprocess.run(command1, shell=True, stdout=_subprocess.PIPE, stderr=_subprocess.PIPE)
            if proc.returncode != 0:
                raise _AnalysisError("GROMACS free-energy analysis failed!")

        # Initialise lists to hold the data from each leg.
        leg0 = []
        leg1 = []

        # Extract the data from the output files.

        # First leg.
        with open("%s/bar_leg0.xvg" % work_dir) as file:

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
        if is_dual:
            with open("%s/bar_leg1.xvg" % work_dir) as file:

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
        if is_dual:
            free_energy = (leg0[-1][1] - leg0[0][1]) - (leg1[-1][1] - leg1[0][1])
        else:
            free_energy = leg0[-1][1] - leg0[0][1]

        # Propagate the errors. (These add in quadrature.)

        # First leg.
        error0 = _math.sqrt((leg0[-1][2].magnitude() * leg0[-1][2].magnitude()) +
                            (leg0[0][2].magnitude()  * leg0[0][2].magnitude()))

        # Second leg.
        if is_dual:
            error1 = _math.sqrt((leg1[-1][2].magnitude() * leg1[-1][2].magnitude()) +
                                (leg1[0][2].magnitude()  * leg1[0][2].magnitude()))
        else:
            error1 = 0

        # Free energy difference.
        error = _math.sqrt((error0 * error0) + (error1 * error1)) * _Units.Energy.kcal_per_mol

        # Bundle the free energy and its associated error.
        free_energy = (free_energy, error)

        return (leg0, leg1, free_energy, None, None)

    @classmethod
    def _analyse_somd(cls, work_dir=None, dir0=None, dir1=None, is_dual=True):
        """Analyse the SOMD free energy data.

           Parameters
           ----------

           work_dir : str
               The path to the working directory.

           dir0 : str
               The path to the directory of the first leg.

           dir1 : str
               The path to the directory of the second leg.

           is_dual : bool
               Whether this is a dual-leg simulation.

           Returns
           -------

           pmf0 : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
               The potential of mean force (PMF) for the first leg of the
               simulation. The data is a list of tuples, where each tuple
               contains the lambda value, the PMF, and the standard error.

           pmf1 : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
               The potential of mean force (PMF) for the second leg of the
               simulation. The data is a list of tuples, where each tuple
               contains the lambda value, the PMF, and the standard error.

           free_energy : (:class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)
               The free energy difference and its associated error.

           overlap0 : [ [ float, float, ... ] ]
               The overlap matrix. This gives the overlap between each window
               of the first leg. This parameter is only computed for the SOMD
               engine and will be None when GROMACS is used.

           overlap1 : [ [ float, float, ... ] ]
               The overlap matrix. This gives the overlap between each window
               of the second leg. This parameter is only computed for the SOMD
               engine and will be None when GROMACS is used.
        """

        if type(work_dir) is not str:
            raise TypeError("'work_dir' must be of type 'str'.")
        if not _os.path.isdir(work_dir):
            raise ValueError("'work_dir' doesn't exist!")

        if type(is_dual) is not bool:
            raise TypeError("'is_dual' must be of type 'bool'.")

        if type(dir0) is not str:
            raise TypeError("'dir0' must be of type 'str'.")
        if not _os.path.isdir(dir0):
            raise ValueError("'dir0' doesn't exist!")

        if is_dual:
            if type(dir1) is not str:
                raise TypeError("'dir1' must be of type 'str'.")
            if not _os.path.isdir(dir1):
                raise ValueError("'dir1' doesn't exist!")

        # Create the commands for the two legs.
        command0 = "%s mbar -i %s/lambda_*/simfile.dat -o %s/mbar_leg0.txt --overlap --subsampling" % (cls._analyse_freenrg, dir0, work_dir)
        command1 = "%s mbar -i %s/lambda_*/simfile.dat -o %s/mbar_leg1.txt --overlap --subsampling" % (cls._analyse_freenrg, dir1, work_dir)

        # Run the first command.
        proc = _subprocess.run(command0, shell=True, stdout=_subprocess.PIPE, stderr=_subprocess.PIPE)
        if proc.returncode != 0:
            raise _AnalysisError("SOMD free-energy analysis failed!")

        # Run the second command.
        if is_dual:
            proc = _subprocess.run(command1, shell=True, text=True, stdout=_subprocess.PIPE, stderr=_subprocess.PIPE)
            if proc.returncode != 0:
                raise _AnalysisError("SOMD free-energy analysis failed!")

        # Initialise lists to hold the data from each leg.
        leg0 = []
        leg1 = []

        # Initialise lists to hold the overlap matrix for each leg.
        overlap0 = []
        overlap1 = []

        # Extract the data from the output files.

        # First leg.
        with open("%s/mbar_leg0.txt" % work_dir) as file:

            # Process the MBAR data.
            for line in file:
                # Process the overlap matrix.
                if "#Overlap matrix" in line:

                    # Get the next row.
                    row = next(file)

                    # Loop until we hit the next section.
                    while not row.startswith("#DG"):
                        # Extract the data for this row.
                        data = [float(x) for x in row.split()]

                        # Append to the overlap matrix.
                        overlap0.append(data)

                        # Get the next line.
                        row = next(file)

                # Process the PMF.
                elif "PMF from MBAR" in line:
                    # Get the next row.
                    row = next(file)

                    # Loop until we hit the next section.
                    while not row.startswith("#TI"):
                        # Split the line.
                        data = row.split()

                        # Append the data.
                        leg0.append((float(data[0]),
                                     float(data[1]) * _Units.Energy.kcal_per_mol,
                                     float(data[2]) * _Units.Energy.kcal_per_mol))

                        # Get the next line.
                        row = next(file)


        # Second leg.
        if is_dual:
            with open("%s/mbar_leg1.txt" % work_dir) as file:
                # Process the MBAR data.
                for line in file:
                    # Process the overlap matrix.
                    if "#Overlap matrix" in line:

                        # Get the next row.
                        row = next(file)

                        # Loop until we hit the next section.
                        while not row.startswith("#DG"):
                            # Extract the data for this row.
                            data = [float(x) for x in row.split()]

                            # Append to the overlap matrix.
                            overlap1.append(data)

                            # Get the next line.
                            row = next(file)

                    # Process the PMF.
                    elif "PMF from MBAR" in line:
                        # Get the next row.
                        row = next(file)

                        # Loop until we hit the next section.
                        while not row.startswith("#TI"):
                            # Split the line.
                            data = row.split()

                            # Append the data.
                            leg1.append((float(data[0]),
                                        float(data[1]) * _Units.Energy.kcal_per_mol,
                                        float(data[2]) * _Units.Energy.kcal_per_mol))

                            # Get the next line.
                            row = next(file)

        # Work out the difference in free energy.
        if is_dual:
            free_energy = (leg0[-1][1] - leg0[0][1]) - (leg1[-1][1] - leg1[0][1])
        else:
            free_energy = leg0[-1][1] - leg0[0][1]

        # Propagate the errors. (These add in quadrature.)

        # First leg.
        error0 = _math.sqrt((leg0[-1][2].magnitude() * leg0[-1][2].magnitude()) +
                            (leg0[0][2].magnitude()*leg0[0][2].magnitude()))

        # Second leg.
        if is_dual:
            error1 = _math.sqrt((leg1[-1][2].magnitude() * leg1[-1][2].magnitude()) +
                                (leg1[0][2].magnitude() * leg1[0][2].magnitude()))
        else:
            error1 = 0

        # Free energy difference.
        error = _math.sqrt((error0 * error0) + (error1 * error1)) * _Units.Energy.kcal_per_mol

        # Bundle the free energy and its associated error.
        free_energy = (free_energy, error)

        return (leg0, leg1, free_energy, overlap0, overlap1)

    def _initialise_runner(self, system0, system1=None):
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

        if system1 is None:
            is_dual = False
        else:
            is_dual = True
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
            if is_dual:
                self._dir1 = "%s/vacuum" % self._work_dir
        elif sim_type == "Binding":
            self._dir0 = "%s/bound" % self._work_dir
            if is_dual:
                self._dir1 = "%s/free" % self._work_dir
        else:
            raise TypeError("Unsupported FreeEnergy simulation: '%s'" % sim_type)

        # Convert to an appropriate AMBER topology. (Required by SOMD for its
        # FEP setup.)
        if self._engine == "SOMD":
            system0._set_water_topology("AMBER")
            if is_dual:
                system1._set_water_topology("AMBER")

        # Setup all of the simulation processes for each leg.

        # Get the lambda values from the protocol for the first leg.
        lam_vals = self._protocol0.getLambdaValues()

        # Create a process for the first lambda value.
        lam = lam_vals[0]

        # Update the protocol lambda values.
        self._protocol0.setLambdaValues(lam=lam, lam_vals=lam_vals)

        # Create and append the required processes for each leg.
        # Nest the working directories inside self._work_dir.

        # Name the first directory.
        first_dir = "%s/lambda_%5.4f" % (self._dir0, lam)

        # SOMD.
        if self._engine == "SOMD":
            # Check for GPU support.
            if "CUDA_VISIBLE_DEVICES" in _os.environ:
                platform = "CUDA"
            else:
                platform = "CPU"

            first_process = _Process.Somd(system0, self._protocol0,
                platform=platform, work_dir=first_dir)
            if self._setup_only:
                del(first_process)
            else:
                leg0.append(first_process)

        # GROMACS.
        elif self._engine == "GROMACS":
            first_process = _Process.Gromacs(system0, self._protocol0,
                work_dir=first_dir, ignore_warnings=self._ignore_warnings,
                show_errors=self._show_errors)
            if self._setup_only:
                del(first_process)
            else:
                leg0.append(first_process)

        # Loop over the rest of the lambda values.
        for x, lam in enumerate(lam_vals[1:]):
            # Name the directory.
            new_dir = "%s/lambda_%5.4f" % (self._dir0, lam)

            # Use the full path.
            if new_dir[0] != "/":
                new_dir = _os.getcwd() + "/" + new_dir

            # Delete any existing directories.
            if _os.path.isdir(new_dir):
                _shutil.rmtree(new_dir, ignore_errors=True)

            # Copy the first directory to that of the current lambda value.
            _shutil.copytree(first_dir, new_dir)

            # Update the protocol lambda values.
            self._protocol0.setLambdaValues(lam=lam, lam_vals=lam_vals)

            # Now update the lambda values in the config files.

            # SOMD.
            if self._engine == "SOMD":
                new_config = []
                with open(new_dir + "/somd.cfg", "r") as f:
                    for line in f:
                        if "lambda_val" in line:
                            new_config.append("lambda_val = %s\n" % lam)
                        else:
                            new_config.append(line)
                with open(new_dir + "/somd.cfg", "w") as f:
                    for line in new_config:
                        f.write(line)

                # Create a copy of the process and update the working
                # directory.
                if not self._setup_only:
                    process                 = _copy.copy(first_process)
                    process._system         = first_process._system.copy()
                    process._protocol       = self._protocol0
                    process._work_dir       = new_dir
                    process._std_out_file   = new_dir + "/somd.out"
                    process._std_err_file   = new_dir + "/somd.err"
                    process._rst_file       = new_dir + "/somd.rst7"
                    process._top_file       = new_dir + "/somd.prm7"
                    process._traj_file      = new_dir + "/traj000000001.dcd"
                    process._restart_file   = new_dir + "/latest.rst"
                    process._config_file    = new_dir + "/somd.cfg"
                    process._pert_file      = new_dir + "/somd.pert"
                    process._gradients_file = new_dir + "/gradients.dat"
                    process._input_files    = [process._config_file,
                                                process._rst_file,
                                                process._top_file,
                                                process._pert_file]
                    leg0.append(process)

            # GROMACS.
            elif self._engine == "GROMACS":
                # Delete the existing binary run file.
                _os.remove(new_dir + "/gromacs.tpr")
                new_config = []
                with open(new_dir + "/gromacs.mdp", "r") as f:
                    for line in f:
                        if "init-lambda-state" in line:
                            new_config.append("init-lambda-state = %d\n" % (x+1))
                        else:
                            new_config.append(line)
                with open(new_dir + "/gromacs.mdp", "w") as f:
                    for line in new_config:
                        f.write(line)

                mdp     = new_dir + "/gromacs.mdp"
                mdp_out = new_dir + "/gromacs.out.mdp"
                gro     = new_dir + "/gromacs.gro"
                top     = new_dir + "/gromacs.top"
                tpr     = new_dir + "/gromacs.tpr"

                # Use grompp to generate the portable binary run input file.
                command = "%s grompp -f %s -po %s -c %s -p %s -r %s -o %s" \
                    % (_gmx_exe, mdp, mdp_out, gro, top, gro, tpr)

                # Run the command. If this worked for the first lambda value,
                # then it should work for all others.
                proc = _subprocess.run(command, shell=True, text=True,
                    stdout=_subprocess.PIPE, stderr=_subprocess.PIPE)

                # Create a copy of the process and update the working
                # directory.
                if not self._setup_only:
                    process                 = _copy.copy(first_process)
                    process._system         = first_process._system.copy()
                    process._protocol       = self._protocol0
                    process._work_dir       = new_dir
                    process._std_out_file   = new_dir + "/gromacs.out"
                    process._std_err_file   = new_dir + "/gromacs.err"
                    process._gro_file       = new_dir + "/gromacs.gro"
                    process._top_file       = new_dir + "/gromacs.top"
                    process._traj_file      = new_dir + "/gromacs.trr"
                    process._config_file    = new_dir + "/gromacs.mdp"
                    process._tpr_file       = new_dir + "/gromacs.tpr"
                    process._input_files    = [process._config_file,
                                                process._gro_file,
                                                process._top_file,
                                                process._tpr_file]
                    leg0.append(process)

        if is_dual:
            # Get the lambda values from the protocol for the first leg.
            lam_vals = self._protocol1.getLambdaValues()

            # Create a process for the first lambda value.
            lam = lam_vals[0]

            # Update the protocol lambda values.
            self._protocol1.setLambdaValues(lam=lam, lam_vals=lam_vals)

            # Create and append the required processes for each leg.
            # Nest the working directories inside self._work_dir.

            # Name the first directory.
            first_dir = "%s/lambda_%5.4f" % (self._dir1, lam)

            # SOMD.
            if self._engine == "SOMD":
                # Check for GPU support.
                if "CUDA_VISIBLE_DEVICES" in _os.environ:
                    platform = "CUDA"
                else:
                    platform = "CPU"

                first_process = _Process.Somd(system1, self._protocol1,
                    platform=platform, work_dir=first_dir)
                if self._setup_only:
                    del(first_process)
                else:
                    leg1.append(first_process)

            # GROMACS.
            elif self._engine == "GROMACS":
                first_process = _Process.Gromacs(system1, self._protocol1,
                    work_dir=first_dir, ignore_warnings=self._ignore_warnings,
                    show_errors=self._show_errors)
                if self._setup_only:
                    del(first_process)
                else:
                    leg1.append(first_process)

            # Loop over the rest of the lambda values.
            for x, lam in enumerate(lam_vals[1:]):
                # Name the directory.
                new_dir = "%s/lambda_%5.4f" % (self._dir1, lam)

                # Use the full path.
                if new_dir[0] != "/":
                    new_dir = _os.getcwd() + "/" + new_dir

                # Delete any existing directories.
                if _os.path.isdir(new_dir):
                    _shutil.rmtree(new_dir, ignore_errors=True)

                # Copy the first directory to that of the current lambda value.
                _shutil.copytree(first_dir, new_dir)

                # Update the protocol lambda values.
                self._protocol1.setLambdaValues(lam=lam, lam_vals=lam_vals)

                # Now update the lambda values in the config files.

                # SOMD.
                if self._engine == "SOMD":
                    new_config = []
                    with open(new_dir + "/somd.cfg", "r") as f:
                        for line in f:
                            if "lambda_val" in line:
                                new_config.append("lambda_val = %s\n" % lam)
                            else:
                                new_config.append(line)
                    with open(new_dir + "/somd.cfg", "w") as f:
                        for line in new_config:
                            f.write(line)

                    # Create a copy of the process and update the working
                    # directory.
                    if not self._setup_only:
                        process                 = _copy.copy(first_process)
                        process._system         = first_process._system.copy()
                        process._protocol       = self._protocol1
                        process._work_dir       = new_dir
                        process._std_out_file   = new_dir + "/somd.out"
                        process._std_err_file   = new_dir + "/somd.err"
                        process._rst_file       = new_dir + "/somd.rst7"
                        process._top_file       = new_dir + "/somd.prm7"
                        process._traj_file      = new_dir + "/traj000000001.dcd"
                        process._restart_file   = new_dir + "/latest.rst"
                        process._config_file    = new_dir + "/somd.cfg"
                        process._pert_file      = new_dir + "/somd.pert"
                        process._gradients_file = new_dir + "/gradients.dat"
                        process._input_files    = [process._config_file,
                                                    process._rst_file,
                                                    process._top_file,
                                                    process._pert_file]
                        leg1.append(process)

                # GROMACS.
                elif self._engine == "GROMACS":
                    # Delete the existing binary run file.
                    _os.remove(new_dir + "/gromacs.tpr")
                    new_config = []
                    with open(new_dir + "/gromacs.mdp", "r") as f:
                        for line in f:
                            if "init-lambda-state" in line:
                                new_config.append("init-lambda-state = %d\n" % (x+1))
                            else:
                                new_config.append(line)
                    with open(new_dir + "/gromacs.mdp", "w") as f:
                        for line in new_config:
                            f.write(line)

                    mdp     = new_dir + "/gromacs.mdp"
                    mdp_out = new_dir + "/gromacs.out.mdp"
                    gro     = new_dir + "/gromacs.gro"
                    top     = new_dir + "/gromacs.top"
                    tpr     = new_dir + "/gromacs.tpr"

                    # Use grompp to generate the portable binary run input file.
                    command = "%s grompp -f %s -po %s -c %s -p %s -r %s -o %s" \
                        % (_gmx_exe, mdp, mdp_out, gro, top, gro, tpr)

                    # Run the command. If this worked for the first lambda value,
                    # then it should work for all others.
                    proc = _subprocess.run(command, shell=True, text=True,
                        stdout=_subprocess.PIPE, stderr=_subprocess.PIPE)

                    # Create a copy of the process and update the working
                    # directory.
                    if not self._setup_only:
                        process                 = _copy.copy(first_process)
                        process._system         = first_process._system.copy()
                        process._protocol       = self._protocol1
                        process._work_dir       = new_dir
                        process._std_out_file   = new_dir + "/gromacs.out"
                        process._std_err_file   = new_dir + "/gromacs.err"
                        process._gro_file       = new_dir + "/gromacs.gro"
                        process._top_file       = new_dir + "/gromacs.top"
                        process._traj_file      = new_dir + "/gromacs.trr"
                        process._config_file    = new_dir + "/gromacs.mdp"
                        process._tpr_file       = new_dir + "/gromacs.tpr"
                        process._input_files    = [process._config_file,
                                                process._gro_file,
                                                process._top_file,
                                                process._tpr_file]
                        leg1.append(process)

        if not self._setup_only:
            # Initialise the process runner. All processes have already been nested
            # inside the working directory so no need to re-nest.
            self._runner = _Process.ProcessRunner(leg0 + leg1)

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

def analyse(work_dir, simulation_type=None):
    """Analyse existing free-energy data from a simulation working directory.

       For the return values, leg0 refers to the simulation leg with the
       larger number of molecules, i.e. the "bound" leg for a "binding"
       free-energy simulation, or the "free" leg for a "solvation"
       free-energy simulation.

       Parameters
       ----------

       work_dir : str
           The working directory for the simulation.

       simulation_type : str
           The type of free-energy perturbation simulation_type. Options are:
           "solvation", or "binding". This allows use of a single "work_dir"
           for combined "solvation" and "binding" free-energy simulations.
           If None, then BioSimSpace will attempt to figure out the simulation
           type from the directory structure within "work_dir".

       Returns
       -------

       pmf0 : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
           The potential of mean force (PMF) for the first leg of the
           simulation. The data is a list of tuples, where each tuple
           contains the lambda value, the PMF, and the standard error.

       pmf1 : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
           The potential of mean force (PMF) for the second leg of the
           simulation. The data is a list of tuples, where each tuple
           contains the lambda value, the PMF, and the standard error.

       free_energy : (:class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)
           The free energy difference and its associated error.

       overlap0 : [ [ float, float, ... ] ]
           The overlap matrix. This gives the overlap between each window
           of the first leg. This parameter is only computed for the SOMD
           engine and will be None when GROMACS is used.

       overlap1 : [ [ float, float, ... ] ]
           The overlap matrix. This gives the overlap between each window
           of the second leg. This parameter is only computed for the SOMD
           engine and will be None when GROMACS is used.
    """

    return FreeEnergy.analyse(work_dir, simulation_type)

def getData(name="data", file_link=False, work_dir=None):
    """Return a link to a zip file containing the data files required for
       post-simulation analysis.

       Parameters
       ----------

       name : str
           The name of the zip file.

       file_link : bool
           Whether to return a FileLink when working in Jupyter.

       work_dir : str
           The working directory for the simulation.

       Returns
       -------

       ouput : str, IPython.display.FileLink
           A path, or file link, to an archive of the process input.
    """
    return FreeEnergy.getData(name=name, file_link=file_link, work_dir=work_dir)

######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2023
#
# Authors: Lester Hedges <lester.hedges@gmail.com>
#
# BioSimSpace is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
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

"""Functionality for relative free-energy simulations."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Relative", "getData"]

from collections import OrderedDict as _OrderedDict
from glob import glob as _glob

import copy as _copy
import math as _math
import shlex as _shlex
import sys as _sys
import os as _os
import shutil as _shutil
import subprocess as _subprocess
import warnings as _warnings
import zipfile as _zipfile

from sire.legacy.Base import getBinDir as _getBinDir
from sire.legacy.Base import getShareDir as _getShareDir

from sire.legacy import IO as _SireIO
from sire.legacy import Mol as _SireMol

from .. import _gmx_exe
from .. import _is_notebook
from .._Exceptions import AnalysisError as _AnalysisError
from .._Exceptions import MissingSoftwareError as _MissingSoftwareError
from .._SireWrappers import Molecules as _Molecules
from .._SireWrappers import System as _System
from .._Utils import cd as _cd
from .. import Process as _Process
from .. import Protocol as _Protocol
from .. import Types as _Types
from .. import Units as _Units
from .. import _Utils

if _is_notebook:
    from IPython.display import FileLink as _FileLink

# Check that the analyse_freenrg script exists.
if _sys.platform != "win32":
    _analyse_freenrg = _os.path.join(_getBinDir(), "analyse_freenrg")
else:
    _analyse_freenrg = _os.path.join(
        _os.path.normpath(_getShareDir()), "scripts", "analyse_freenrg.py"
    )
if not _os.path.isfile(_analyse_freenrg):
    raise _MissingSoftwareError(
        "Cannot find free energy analysis script in expected location: '%s'"
        % _analyse_freenrg
    )
if _sys.platform == "win32":
    _analyse_freenrg = "%s %s" % (
        _os.path.join(_os.path.normpath(_getBinDir()), "sire_python.exe"),
        _analyse_freenrg,
    )


class Relative:
    """Class for configuring and running relative free-energy perturbation simulations."""

    # Create a list of supported molecular dynamics engines.
    _engines = ["GROMACS", "SOMD"]

    def __init__(
        self,
        system,
        protocol=None,
        work_dir=None,
        engine=None,
        setup_only=False,
        ignore_warnings=False,
        show_errors=True,
        property_map={},
    ):
        """
        Constructor.

        Parameters
        ----------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system for the perturbation leg. This must contain
            a single perturbable molecule and is assumed to have already
            been equilibrated.

        protocol : :class:`Protocol.FreeEnergyMinimisation <BioSimSpace.Protocol.FreeEnergyMinimisation>`, \
                   :class:`Protocol.FreeEnergyEquilibration <BioSimSpace.Protocol.FreeEnergyEquilibration>`, \
                   :class:`Protocol.FreeEnergyProduction <BioSimSpace.Protocol.FreeEnergyProduction>`
            The simulation protocol.

        work_dir : str
            The working directory for the free-energy perturbation
            simulation.

        engine : str
            The molecular dynamics engine used to run the simulation. Available
            options are "GROMACS", or "SOMD". If this argument is omitted then
            BioSimSpace will choose an appropriate engine for you.

        setup_only : bool
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

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Validate the input.

        if not isinstance(system, _System):
            raise TypeError(
                "'system' must be of type 'BioSimSpace._SireWrappers.System'"
            )
        else:
            # Store a copy of solvated system.
            self._system = system.copy()

        # Validate the user specified molecular dynamics engine.
        if engine is not None:
            if not isinstance(engine, str):
                raise TypeError("'engine' must be of type 'str'.")

            # Strip whitespace from engine and convert to upper case.
            engine = engine.replace(" ", "").upper()

            # Check that the engine is supported.
            if engine not in self._engines:
                raise ValueError(
                    f"Unsupported molecular dynamics engine {engine}. "
                    "Supported engines are: %r." % ", ".join(self._engines)
                )

            # Make sure GROMACS is installed if GROMACS engine is selected.
            if engine == "GROMACS":
                if _gmx_exe is None:
                    raise _MissingSoftwareError(
                        "Cannot use GROMACS engine as GROMACS is not installed!"
                    )

                # The system must have a perturbable molecule.
                if system.nPerturbableMolecules() == 0:
                    raise ValueError(
                        "The system must contain a perturbable molecule! "
                        "Use the 'BioSimSpace.Align' package to map and merge molecules."
                    )

        else:
            # Use SOMD as a default.
            engine = "SOMD"

            # The system must have a single perturbable molecule.
            if system.nPerturbableMolecules() != 1:
                raise ValueError(
                    "The system must contain a single perturbable molecule! "
                    "Use the 'BioSimSpace.Align' package to map and merge molecules."
                )

        # Set the engine.
        self._engine = engine

        # Validate the protocol.
        if protocol is not None:
            from ..Protocol._free_energy_mixin import _FreeEnergyMixin

            if isinstance(protocol, _FreeEnergyMixin):
                if engine == "SOMD" and not isinstance(
                    protocol, _Protocol.FreeEnergyProduction
                ):
                    raise ValueError(
                        "Currently SOMD only supports protocols of type 'BioSimSpace.Protocol.FreeEnergyProduction'"
                    )

                self._protocol = protocol
            else:
                raise TypeError(
                    "'protocol' must be of type 'BioSimSpace.Protocol.FreeEnergy'"
                )
        else:
            # Use a default protocol.
            self._protocol = _Protocol.FreeEnergy()

        # Check that multi-step perturbation isn't specified if GROMACS is the chosen engine.
        if engine == "GROMACS":
            if self._protocol.getPerturbationType() != "full":
                raise NotImplementedError(
                    "GROMACS currently only supports the 'full' perturbation "
                    "type. Please use engine='SOMD' when running multistep "
                    "perturbation types."
                )

        if not isinstance(setup_only, bool):
            raise TypeError("'setup_only' must be of type 'bool'.")
        else:
            self._setup_only = setup_only

        if work_dir is None and setup_only:
            raise ValueError(
                "A 'work_dir' must be specified when 'setup_only' is True!"
            )

        # Create the working directory.
        self._work_dir = _Utils.WorkDir(work_dir)

        if not isinstance(ignore_warnings, bool):
            raise ValueError("'ignore_warnings' must be of type 'bool.")
        self._ignore_warnings = ignore_warnings

        if not isinstance(show_errors, bool):
            raise ValueError("'show_errors' must be of type 'bool.")
        self._show_errors = show_errors

        # Check that the map is valid.
        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'")
        self._property_map = property_map

        # Create fake instance methods for 'analyse' and 'difference'. These
        # pass instance data through to the staticmethod versions.
        self.analyse = self._analyse
        self.difference = self._difference

        # Initialise the process runner.
        self._initialise_runner(self._system)

    def run(self, serial=True):
        """
        Run the simulation.

        Parameters
        ----------

        serial : bool
            Whether to run the individual processes for the lambda windows
            in serial.
        """
        if not isinstance(serial, bool):
            raise TypeError("'serial' must be of type 'bool'.")

        if self._setup_only:
            _warnings.warn("No processes exist! Object created in 'setup_only' mode.")
        else:
            self._runner.startAll(serial=serial)

    def wait(self):
        """Wait for the simulation to finish."""
        if self._setup_only:
            _warnings.warn("No processes exist! Object created in 'setup_only' mode.")
        else:
            self._runner.wait()

    def kill(self, index):
        """
        Kill a process for a specific lambda window.

        Parameters
        ----------

        index : int
            The index of the lambda window.
        """
        self._runner.kill(index)

    def killAll(self):
        """Kill any running processes for all lambda windows."""

        self._runner.killAll()

    def workDir(self):
        """
        Return the working directory.

        Returns
        -------

        work_dir : str
            The path of the working directory.
        """
        return str(self._work_dir)

    def getData(self, name="data", file_link=False, work_dir=None):
        """
        Return a link to a zip file containing the data files required for
        post-simulation analysis.

        Parameters
        ----------

        name : str
            The name of the zip file.

        file_link : bool
            Whether to return a FileLink when working in Jupyter.

        work_dir : str
            The working directory for the free-energy perturbation
            simulation.

        Returns
        -------

        output : str, IPython.display.FileLink
            A path, or file link, to an archive of the process input.
        """

        if self._work_dir is None:
            raise ValueError("'work_dir' must be set!")
        else:
            if not isinstance(work_dir, str):
                raise TypeError("'work_dir' must be of type 'str'.")
            if not _os.path.isdir(work_dir):
                raise ValueError("'work_dir' doesn't exist!")

        if not isinstance(name, str):
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
                    raise ValueError(
                        f"Couldn't find any analysis files in '{work_dir}'"
                    )

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
                f_link.html_link_str = (
                    f"<a href='%s' target='_blank' download='{zipname}'>%s</a>"
                )

                # Return a link to the archive.
                return f_link
            else:
                return zipname
        # Return the path to the archive.
        else:
            return zipname

    @staticmethod
    def analyse(work_dir):
        """
        Analyse existing free-energy data from a simulation working directory.

        Parameters
        ----------

        work_dir : str
            The working directory for the simulation.

        Returns
        -------

        pmf : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
            The potential of mean force (PMF). The data is a list of tuples,
            where each tuple contains the lambda value, the PMF, and the
            standard error.

        overlap : [ [ float, float, ... ] ]
            The overlap matrix. This gives the overlap between each lambda
            window.  This parameter is only computed for the SOMD engine and
            will be None when GROMACS is used.
        """

        if not isinstance(work_dir, str):
            raise TypeError("'work_dir' must be of type 'str'.")
        if not _os.path.isdir(work_dir):
            raise ValueError("'work_dir' doesn't exist!")

        # First test for SOMD files.
        data = _glob(work_dir + "/lambda_*/gradients.dat")

        # SOMD.
        if len(data) > 0:
            return Relative._analyse_somd(work_dir)

        # Now check for GROMACS output.
        else:
            data = _glob(work_dir + "/lambda_*/gromacs.xvg")
            if len(data) == 0:
                raise ValueError(
                    "Couldn't find any SOMD or GROMACS free-energy output?"
                )
            return Relative._analyse_gromacs(work_dir)

    def _analyse(self):
        """
        Analyse free-energy data for this object.

        Returns
        -------

        pmf : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
            The potential of mean force (PMF). The data is a list of tuples,
            where each tuple contains the lambda value, the PMF, and the
            standard error.

        overlap : [ [ float, float, ... ] ]
            The overlap matrix. This gives the overlap between each lambda
            window.  This parameter is only computed for the SOMD engine and
            will be None when GROMACS is used.
        """

        # Return the result of calling the staticmethod, passing in the working
        # directory of this object.
        return Relative.analyse(str(self._work_dir))

    @staticmethod
    def _analyse_gromacs(work_dir=None):
        """
        Analyse the GROMACS free energy data.

        Parameters
        ----------

        work_dir : str
            The path to the working directory.

        Returns
        -------

        pmf : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
            The potential of mean force (PMF). The data is a list of tuples,
            where each tuple contains the lambda value, the PMF, and the
            standard error.
        """

        if not isinstance(work_dir, str):
            raise TypeError("'work_dir' must be of type 'str'.")
        if not _os.path.isdir(work_dir):
            raise ValueError("'work_dir' doesn't exist!")

        # Create the command.
        xvg_files = _glob(f"{work_dir}/lambda_*/*.xvg")
        command = "%s bar -f %s -o %s/bar.xvg" % (
            _gmx_exe,
            " ".join(xvg_files),
            work_dir,
        )

        # Run the first command.
        proc = _subprocess.run(
            _Utils.command_split(command),
            shell=False,
            stdout=_subprocess.PIPE,
            stderr=_subprocess.PIPE,
        )
        if proc.returncode != 0:
            raise _AnalysisError("GROMACS free-energy analysis failed!")

        # Initialise list to hold the data.
        data = []

        # Extract the data from the output files.

        # First leg.
        with open("%s/bar.xvg" % work_dir) as file:
            # Read all of the lines into a list.
            lines = []
            for line in file:
                # Ignore comments and xmgrace directives.
                if line[0] != "#" and line[0] != "@":
                    lines.append(line.rstrip())

            # Store the initial free energy reading.
            data.append(
                (
                    0.0,
                    0.0 * _Units.Energy.kcal_per_mol,
                    0.0 * _Units.Energy.kcal_per_mol,
                )
            )

            # Zero the accumulated error.
            total_error = 0

            # Zero the accumulated free energy difference.
            total_freenrg = 0

            # Process the BAR data.
            for x, line in enumerate(lines):
                # Extract the data from the line.
                records = line.split()

                # Update the total free energy difference.
                total_freenrg += float(records[1])

                # Extract the error.
                error = float(records[2])

                # Update the accumulated error.
                total_error = _math.sqrt(total_error * total_error + error * error)

                # Append the data.
                data.append(
                    (
                        (x + 1) / (len(lines)),
                        (total_freenrg * _Units.Energy.kt).kcal_per_mol(),
                        (total_error * _Units.Energy.kt).kcal_per_mol(),
                    )
                )

        return (data, None)

    @staticmethod
    def _analyse_somd(work_dir=None):
        """
        Analyse the SOMD free energy data.

        Parameters
        ----------

        work_dir : str
            The path to the working directory.

        Returns
        -------

        pmf : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
            The potential of mean force (PMF). The data is a list of tuples,
            where each tuple contains the lambda value, the PMF, and the
            standard error.

        overlap : [ [ float, float, ... ] ]
            The overlap matrix. This gives the overlap between each lambda
            window.  This parameter is only computed for the SOMD engine and
            will be None when GROMACS is used.
        """

        if not isinstance(work_dir, str):
            raise TypeError("'work_dir' must be of type 'str'.")
        if not _os.path.isdir(work_dir):
            raise ValueError("'work_dir' doesn't exist!")

        # Create the command.
        command = (
            "%s mbar -i %s/lambda_*/simfile.dat -o %s/mbar.txt --overlap --subsampling"
            % (_analyse_freenrg, work_dir, work_dir)
        )

        # Run the first command.
        proc = _subprocess.run(
            _Utils.command_split(command),
            shell=False,
            stdout=_subprocess.PIPE,
            stderr=_subprocess.PIPE,
        )
        if proc.returncode != 0:
            raise _AnalysisError("SOMD free-energy analysis failed!")

        # Re-run without subsampling if the subsampling has resulted in less than 50 samples.
        with open("%s/mbar.txt" % work_dir) as file:
            for line in file:
                if (
                    "#WARNING SUBSAMPLING ENERGIES RESULTED IN LESS THAN 50 SAMPLES"
                    in line
                ):
                    _warnings.warn(
                        "Subsampling resulted in less than 50 samples, "
                        f"re-running without subsampling for '{work_dir}'"
                    )
                    command = (
                        "%s mbar -i %s/lambda_*/simfile.dat -o %s/mbar.txt --overlap"
                        % (_analyse_freenrg, work_dir, work_dir)
                    )
                    proc = _subprocess.run(
                        _Utils.command_split(command),
                        shell=False,
                        stdout=_subprocess.PIPE,
                        stderr=_subprocess.PIPE,
                    )
                    if proc.returncode != 0:
                        raise _AnalysisError("SOMD free-energy analysis failed!")
                    break

        # Initialise list to hold the data.
        data = []

        # Initialise list to hold the overlap matrix.
        overlap = []

        # Extract the data from the output files.

        # First leg.
        with open("%s/mbar.txt" % work_dir) as file:
            # Process the MBAR data.
            for line in file:
                # Process the overlap matrix.
                if "#Overlap matrix" in line:
                    # Get the next row.
                    row = next(file)

                    # Loop until we hit the next section.
                    while not row.startswith("#DG"):
                        # Extract the data for this row.
                        records = [float(x) for x in row.split()]

                        # Append to the overlap matrix.
                        overlap.append(records)

                        # Get the next line.
                        row = next(file)

                # Process the PMF.
                elif "PMF from MBAR" in line:
                    # Get the next row.
                    row = next(file)

                    # Loop until we hit the next section.
                    while not row.startswith("#TI"):
                        # Split the line.
                        records = row.split()

                        # Append the data.
                        data.append(
                            (
                                float(records[0]),
                                float(records[1]) * _Units.Energy.kcal_per_mol,
                                float(records[2]) * _Units.Energy.kcal_per_mol,
                            )
                        )

                        # Get the next line.
                        row = next(file)

        return (data, overlap)

    @staticmethod
    def difference(pmf, pmf_ref):
        """
        Compute the relative free-energy difference between two perturbation
        legs.

        Parameters
        ----------

        pmf : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
            The potential of mean force (PMF) of interest. The data is a list
            of tuples, where each tuple contains the lambda value, the PMF,
            and the standard error.

        pmf_ref : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
            The reference potential of mean force (PMF). The data is a list
            of tuples, where each tuple contains the lambda value, the PMF,
            and the standard error.

        Returns
        -------

        free_energy : (:class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)
            The relative free-energy difference and its associated error.
        """

        if not isinstance(pmf, list):
            raise TypeError("'pmf' must be of type 'list'.")
        if not isinstance(pmf_ref, list):
            raise TypeError("'pmf_ref' must be of type 'list'.")

        for rec in pmf:
            if not isinstance(rec, tuple):
                raise TypeError(
                    "'pmf1' must contain tuples containing lambda "
                    "values and the associated free-energy and error."
                )
            else:
                if len(rec) != 3:
                    raise ValueError(
                        "Each tuple in 'pmf1' must contain three items: "
                        "a lambda value and the associated free energy "
                        "and error."
                    )
                for val in rec[1:]:
                    if not isinstance(val, _Types.Energy):
                        raise TypeError(
                            "'pmf' must contain 'BioSimSpace.Types.Energy' types."
                        )

        for rec in pmf_ref:
            if not isinstance(rec, tuple):
                raise TypeError(
                    "'pmf_ref' must contain tuples containing lambda "
                    "values and the associated free-energy and error."
                )
            else:
                if len(rec) != 3:
                    raise ValueError(
                        "Each tuple in 'pmf_ref' must contain three items: "
                        "a lambda value and the associated free energy "
                        "and error."
                    )
                for val in rec[1:]:
                    if not isinstance(val, _Types.Energy):
                        raise TypeError(
                            "'pmf_ref' must contain 'BioSimSpace.Types.Energy' types."
                        )

        # Work out the difference in free energy.
        free_energy = (pmf[-1][1] - pmf[0][1]) - (pmf_ref[-1][1] - pmf_ref[0][1])

        # Propagate the errors. (These add in quadrature.)

        # Measure.
        error0 = _math.sqrt(
            (pmf[-1][2].value() * pmf[-1][2].value())
            + (pmf[0][2].value() * pmf[0][2].value())
        )

        # Reference.
        error1 = _math.sqrt(
            (pmf_ref[-1][2].value() * pmf_ref[-1][2].value())
            + (pmf_ref[0][2].value() * pmf_ref[0][2].value())
        )

        # Error for free-energy difference.
        error = (
            _math.sqrt((error0 * error0) + (error1 * error1))
            * _Units.Energy.kcal_per_mol
        )

        return (free_energy, error)

    def _difference(self, pmf_ref):
        """
        Compute the relative free-energy difference between two perturbation
        legs.

        Parameters
        ----------

        pmf_ref : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
            The reference potential of mean force (PMF). The data is a list
            of tuples, where each tuple contains the lambda value, the PMF,
            and the standard error.

        Returns
        -------

        free_energy : (:class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)
            The relative free-energy difference and its associated error.
        """

        # Calculate the PMF for this object.
        pmf, _ = self.analyse()

        # Now call the staticmethod passing in both PMFs.
        return Relative.difference(pmf, pmf_ref)

    def _initialise_runner(self, system):
        """
        Internal helper function to initialise the process runner.

        Parameters
        ----------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system.
        """

        # Initialise list to store the processe
        processes = []

        # Convert to an appropriate water topology.
        if self._engine == "SOMD":
            system._set_water_topology("AMBER", property_map=self._property_map)
        elif self._engine == "GROMACS":
            system._set_water_topology("GROMACS", property_map=self._property_map)

        # Setup all of the simulation processes for each leg.

        # Get the lambda values from the protocol for the first leg.
        lam_vals = self._protocol.getLambdaValues()

        # Create a process for the first lambda value.
        lam = lam_vals[0]

        # Update the protocol lambda values.
        self._protocol.setLambdaValues(lam=lam, lam_vals=lam_vals)

        # Create and append the required processes for each leg.
        # Nest the working directories inside self._work_dir.

        # Name the first directory.
        first_dir = "%s/lambda_%5.4f" % (self._work_dir, lam)

        # SOMD.
        if self._engine == "SOMD":
            # Check for GPU support.
            if "CUDA_VISIBLE_DEVICES" in _os.environ:
                platform = "CUDA"
            else:
                platform = "CPU"

            first_process = _Process.Somd(
                system,
                self._protocol,
                platform=platform,
                work_dir=first_dir,
                property_map=self._property_map,
            )
            if self._setup_only:
                del first_process
            else:
                processes.append(first_process)

        # GROMACS.
        elif self._engine == "GROMACS":
            first_process = _Process.Gromacs(
                system,
                self._protocol,
                work_dir=first_dir,
                ignore_warnings=self._ignore_warnings,
                show_errors=self._show_errors,
            )
            if self._setup_only:
                del first_process
            else:
                processes.append(first_process)

        # Loop over the rest of the lambda values.
        for x, lam in enumerate(lam_vals[1:]):
            # Name the directory.
            new_dir = "%s/lambda_%5.4f" % (self._work_dir, lam)

            # Use absolute path.
            if not _os.path.isabs(new_dir):
                new_dir = _os.path.abspath(new_dir)

            # Delete any existing directories.
            if _os.path.isdir(new_dir):
                _shutil.rmtree(new_dir, ignore_errors=True)

            # Copy the first directory to that of the current lambda value.
            _shutil.copytree(first_dir, new_dir)

            # Update the protocol lambda values.
            self._protocol.setLambdaValues(lam=lam, lam_vals=lam_vals)

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
                    process = _copy.copy(first_process)
                    process._system = first_process._system.copy()
                    process._protocol = self._protocol
                    process._work_dir = new_dir
                    process._std_out_file = new_dir + "/somd.out"
                    process._std_err_file = new_dir + "/somd.err"
                    process._rst_file = new_dir + "/somd.rst7"
                    process._top_file = new_dir + "/somd.prm7"
                    process._traj_file = new_dir + "/traj000000001.dcd"
                    process._restart_file = new_dir + "/latest.rst"
                    process._config_file = new_dir + "/somd.cfg"
                    process._pert_file = new_dir + "/somd.pert"
                    process._gradients_file = new_dir + "/gradients.dat"
                    process._input_files = [
                        process._config_file,
                        process._rst_file,
                        process._top_file,
                        process._pert_file,
                    ]
                    processes.append(process)

            # GROMACS.
            elif self._engine == "GROMACS":
                # Delete the existing binary run file.
                _os.remove(new_dir + "/gromacs.tpr")
                new_config = []
                with open(new_dir + "/gromacs.mdp", "r") as f:
                    for line in f:
                        if "init-lambda-state" in line:
                            new_config.append("init-lambda-state = %d\n" % (x + 1))
                        else:
                            new_config.append(line)
                with open(new_dir + "/gromacs.mdp", "w") as f:
                    for line in new_config:
                        f.write(line)

                mdp = new_dir + "/gromacs.mdp"
                mdp_out = new_dir + "/gromacs.out.mdp"
                gro = new_dir + "/gromacs.gro"
                top = new_dir + "/gromacs.top"
                tpr = new_dir + "/gromacs.tpr"

                # Use grompp to generate the portable binary run input file.
                command = "%s grompp -f %s -po %s -c %s -p %s -r %s -o %s" % (
                    _gmx_exe,
                    mdp,
                    mdp_out,
                    gro,
                    top,
                    gro,
                    tpr,
                )

                # Run the command. If this worked for the first lambda value,
                # then it should work for all others.
                proc = _subprocess.run(
                    _Utils.command_split(command),
                    shell=False,
                    text=True,
                    stdout=_subprocess.PIPE,
                    stderr=_subprocess.PIPE,
                )

                # Create a copy of the process and update the working
                # directory.
                if not self._setup_only:
                    process = _copy.copy(first_process)
                    process._system = first_process._system.copy()
                    process._protocol = self._protocol
                    process._work_dir = new_dir
                    process._std_out_file = new_dir + "/gromacs.out"
                    process._std_err_file = new_dir + "/gromacs.err"
                    process._gro_file = new_dir + "/gromacs.gro"
                    process._top_file = new_dir + "/gromacs.top"
                    process._traj_file = new_dir + "/gromacs.trr"
                    process._config_file = new_dir + "/gromacs.mdp"
                    process._tpr_file = new_dir + "/gromacs.tpr"
                    process._input_files = [
                        process._config_file,
                        process._gro_file,
                        process._top_file,
                        process._tpr_file,
                    ]
                    processes.append(process)

        if not self._setup_only:
            # Initialise the process runner. All processes have already been nested
            # inside the working directory so no need to re-nest.
            self._runner = _Process.ProcessRunner(processes)

    def _update_run_args(self, args):
        """
        Internal function to update run arguments for all subprocesses.

        Parameters
        ----------

        args : dict, collections.OrderedDict
            A dictionary which contains the new command-line arguments
            for the process executable.
        """

        if not isinstance(args, dict):
            raise TypeError("'args' must be of type 'dict'")

        for process in self._runner.processes():
            process.setArgs(args)


def getData(name="data", file_link=False, work_dir=None):
    """
    Return a link to a zip file containing the data files required for
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

    output : str, IPython.display.FileLink
        A path, or file link, to an archive of the process input.
    """
    return Relative.getData(name=name, file_link=file_link, work_dir=work_dir)
